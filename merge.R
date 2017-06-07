# merge.R
#
# Read in raw data for household count & size (`h`), real household
# disposable income (`rhdi`), general elections (`ge`), UK military
# fatalities (`mf`), and UK population (`pop`), and mould them into usable
# data frames.

h <- read.table("uk_h.dat", header=TRUE)
rhdi <- read.table("uk_rhdi.dat", header=TRUE)
ge <- read.table("uk_ge.dat", header=TRUE)
mf <- read.table("uk_mf.dat", header=TRUE)
pop <- read.table("uk_pop.dat", header=TRUE)

# Compute the incumbent's share of the two-party popular vote in each
# general election for which `ge` has data. Then compute each election's
# year & quarter.
ge$PVI2 <- round(100 * ge$PVI / (ge$PVI + ge$PVC), 2)
ge$DAY <- as.Date(ge$DAY)
ge$YR <- as.integer(format(ge$DAY, "%Y"))
ge$Q <- floor((as.numeric(format(ge$DAY, "%m")) + 2) / 3)  # quarter
ge$Y <- ge$YR + (as.integer(format(ge$DAY, "%j")) / 365.25)

# Extrapolate the last few years of UK population estimates as far into the
# future as the military-fatality counts go.
pop <- pop[order(pop$YR),]  # force recent data to the bottom of `pop`
if (max(pop$YR) < max(mf$YR)) {
	pop_extrap <- data.frame(YR=seq(1+max(pop$YR), max(mf$YR)))
	pop_extrap$POP <- predict(lm(POP ~ YR + I(YR^2), tail(pop, 9)),
	                          pop_extrap)
	pop <- rbind(pop, pop_extrap)
}

# Define `mfs` as the relevant subset of military fatalities: those in
# conflicts where UK forces were in unprovoked, foreign, hostile deployments.
mfs <- mf[(mf$U * mf$F * mf$H) == 1, c("YR", "NAME", "NMF")]

# Merge the military-fatality counts and total UK population estimates.
# Do so by spline interpolating the annual population estimates into
# quarterly estimates; filling in gaps in the relevant fatality estimates
# with zeroes; and merging in those fatality estimates, dividing them by 4
# and linearly interpolating them to make them quarterly too.
fp <- data.frame(YR=rep(pop$YR, 4),
                 Q=rep(c(2,1,3,4), each=length(pop$YR)),
                 POP=c(pop$POP, rep(NA, 3 * length(pop$YR))))
fp <- fp[order(fp$YR, fp$Q),]
fp$Y <- fp$YR + ((fp$Q - 0.5) / 4)
fp$POP[is.na(fp$POP)] <- splinefun(fp$Y, fp$POP)(fp$Y[is.na(fp$POP)])
mfsa <- aggregate(NMF ~ YR, mfs, sum)  # aggregate relevant fatalities by year
empty_years <- (min(pop$YR):max(pop$YR))
empty_years <- empty_years[!(empty_years %in% mfsa$YR)]
if (length(empty_years) > 0) {
	mfsa <- rbind(mfsa, data.frame(YR=empty_years, NMF=0))
}
mfsa <- mfsa[order(mfsa$YR),]
fp$F[fp$Q == 2] <- mfsa$NMF / 4
fp$F[is.na(fp$F)] <- approxfun(fp$Y, fp$F, rule=2)(fp$Y[is.na(fp$F)])

# Calculate the ratio of (relevant) military fatalities to the UK population
# in millions, by year and quarter.
fp$FTP <- fp$F / fp$POP

# The first plot's coming up, so set nice parameters.
par(las=1, mar=c(4.5, 4.1, 0.2, 0.2), mfrow=c(2,1))

# Plot the fatalities-to-population time series, throwing in a label for
# every included conflict, and (just for visual comparison) adding an
# analogous fatalities-to-population time series for all conflicts in the
# original data set (even those not counted by Hibbs's model).
plot(fp$Y, fp$FTP, type="l",
     xlab="year", ylab="relevant military fatalities per million popul.")
fp_all_conflicts <- aggregate(NMF ~ YR, mf, sum)
fp_all_conflicts$NMF <- 0.25 * fp_all_conflicts$NMF / pop$POP
fp_all_conflicts$YR <- fp_all_conflicts$YR + 0.4
lines(fp_all_conflicts, col="#00000030", lty="dashed")
grid()
for (conflict_name in unique(mfs$NAME)) {
	mfss <- mfs[mfs$NAME == conflict_name,]
	conflict_biggest_yr <- mfss[(mfss$NAME == conflict_name)
	                            & (mfss$NMF == max(mfss$NMF)),]
	conflict_biggest_yr <- conflict_biggest_yr[1,]
	text(0.5 + conflict_biggest_yr$YR,
	     0.33 * conflict_biggest_yr$NMF / mean(fp$POP),
	     conflict_name, cex=max(0.4, sqrt(conflict_biggest_yr$NMF) / 6),
	     srt=-8)
}

# Let's try working out the cumulative quarterly fatalities-to-population
# ratio, resetting the cumulative ratio to zero when power changes hands.

ftp_carryover <- 0

for (i in 1:length(fp$FTP)) {

	fp$CFTP[i] <- ftp_carryover + fp$FTP[i]
	ftp_carryover <- fp$CFTP[i]

	# Was there an election this quarter? If so, did the incumbent party
	# lose? If it did, reset the rolling FTP total for next quarter.
	elec_this_q <- ge[(ge$YR == fp$YR[i]) & (ge$Q == fp$Q[i]),]
	if (nrow(elec_this_q)) {
		if (elec_this_q$VIC != elec_this_q$INC) {
			ftp_carryover <- 0
		}
	}

}

# Plot the cumulative fatalities time series against time.

plot(fp$Y, fp$CFTP, type="l",
     xlab="year", ylab="cumulative military-fatalities-to-popul. ratio")
grid()
abline(v=ge$Y, lty="dotted")

# Fill in gaps in the household numbers from my favoured spreadsheet by
# extrapolating from the England & Wales (and sometimes Scotland) numbers.
# Do that in a function to nicely modularize and minimize namespace pollution.

impute_household_counts <- function(h)
{
	impute1 <- lm(NNI ~ NE + NW + NS, h)
	impute2 <- lm(NS ~ NE + NW, h)
	impute3 <- lm(NNI ~ NE + NW, h)

	h$NNI[is.na(h$NNI)] <- predict(impute1, h[is.na(h$NNI),])
	h$NS[is.na(h$NS)] <- predict(impute2, h[is.na(h$NS),])
	h$NNI[is.na(h$NNI)] <- predict(impute3, h[is.na(h$NNI),])
	h$NUK1[is.na(h$NUK1)] <-
		rowSums(h[c("NE", "NW", "NS", "NNI")])[is.na(h$NUK1)]
	h[c("NS", "NNI", "NUK1")] <- round(h[c("NS", "NNI", "NUK1")])

	return(h)
}

h <- impute_household_counts(h)

# Where `h` lacks a UK household size average, but that average can be
# computed from known averages for the four constituent countries, calculate
# the UK-level average.
for (i in 1:nrow(h)) {
	hi <- h[i,]
	if (is.na(hi$SUK)) {
		h$SUK[i] <- weighted.mean(hi[c("SE", "SW", "SS", "SNI")],
		                          hi[c("NE", "NW", "NS", "NNI")])
		h$SUK[i] <- round(h$SUK[i], 2)
	}
}

# Define a function to plot the average-household-size time series.

plo_s <- function(y_lab="average household size")
{
	co <- "#00000070"
	plot(h$YR, h$SE, type="b", ylim=c(2.1, 4.0), col=co,
	     xlab="year", ylab=y_lab)
	points(h$YR[!is.na(h$SW)], h$SW[!is.na(h$SW)], type="b",
	       pch=2, col=co)
	points(h$YR[!is.na(h$SS)], h$SS[!is.na(h$SS)], type="b",
	       pch=3, col=co)
	points(h$YR[!is.na(h$SNI)], h$SNI[!is.na(h$SNI)], type="b", pch=4)
	points(h$YR[!is.na(h$SUK)], h$SUK[!is.na(h$SUK)], type="b", cex=2)
	grid()
	for (col_name in c("SE", "SW", "SS", "SNI", "SUK")) {
		x <- min(h$YR[!is.na(h[[col_name]])])
		if (!(col_name %in% c("SUK", "SW"))) {
			text(x, 0.08 + (h[[col_name]])[h$YR == x], col_name)
		} else if (col_name == "SW") {
			text(1991, 0.1 + (h[["SW"]])[h$YR == 1991], "SW")
		} else {
			text(2001, 0.11 + (h[["SUK"]])[h$YR == 2001], "SUK")
		}
	}
}

# Plot average household sizes before imputation.
plo_s()

# Time to adjust & impute like mad.
# 1. There's a really dodgy-looking jump in the English average for 1972
#    through 1980. Push those years' averages down to bring them in line with
#    the surrounding data.
# 2. Make the hopeful assumption that the gap between Welsh & English average
#    household sizes was constant before 1991, imputing Welsh AHS backwards
#    decadally from that.
# 3. Wedge the years 1952 through 1960 and 1962 through 1970 into `h` by
#    copying rows and blanking them, then linearly interpolate an annual time
#    series for each country in the UK from 1951 onwards.
# 4. Now extrapolate averages for NI, Scotland, & Wales by extending the
#    straight best-fit line to the 2001-2016 data through that whole period.
# 5. Turn to household counts, which are a prerequisite for recomputing
#    average household size. Extrapolate missing household counts for NI,
#    Scotland, England, & Wales by straightforward linear interpolation.
# 6. Fill in remaining gaps in the UK average household size series by
#    taking a household-count-weighted average of the individual countries'
#    average household sizes.

h$SE[(h$YR > 1971) & (h$YR < 1981)] <-
	h$SE[(h$YR > 1971) & (h$YR < 1981)] - 0.03

h$SW[(h$YR < 1991) & ((h$YR %% 10) == 1)] <-
	h$SE[(h$YR < 1991) & ((h$YR %% 10) == 1)] +
	(h$SW[h$YR == 1991] - h$SE[h$YR == 1991])

h <- rbind(h[1,], head(h, 9), h[2,], head(h, 9), h[3:nrow(h),])
h[c(2:10, 12:20),] <- NA
h$YR[c(2:10, 12:20)] <- c(1952:1960, 1962:1970)
h <- h[order(h$YR),]
for (col_name in c("SE", "SW", "SS", "SNI")) {
	h[[col_name]] <- approxfun(h$YR, h[[col_name]])(h$YR)
}

for (col_name in c("SE", "SW", "SS", "SNI")) {
	h_sub <- h[h$YR > 2000,]
	temp_df <- data.frame(Y=h_sub[[col_name]], X=h_sub$YR)
	subse <- (h$YR > 2000) & is.na(h[[col_name]])
	df <- data.frame(X=h$YR[subse])
	h[[col_name]][subse] <- predict(lm(Y ~ X, temp_df), df)
}

for (col_name in c("NE", "NW", "NS", "NNI")) {
	h[[col_name]] <- round(approxfun(h$YR, h[[col_name]])(h$YR))
}

for (i in 1:nrow(h)) {
	hi <- h[i,]
	if (is.na(hi$SUK)) {
		h$SUK[i] <- weighted.mean(hi[c("SE", "SW", "SS", "SNI")],
		                          hi[c("NE", "NW", "NS", "NNI")])
		h$SUK[i] <- round(h$SUK[i], 2)
	}
}

# Plot household-size averages after imputation.
plo_s("imputed average household size")

# Finally `SUK` is a complete time series of (mostly imputed!) average
# UK household size. Now return to the UK household number time series
# and fill in the 1962-1970 gap.

h$NUK1[is.na(h$NUK1)] <- rowSums(h[c("NE", "NW", "NS", "NNI")])[is.na(h$NUK1)]

# Multiply the average UK household size series by the UK household number
# series to produce a time series of the number of people in households.
# (Divide by a thousand to make the final number have units of millions of
# people.)
h$PIH <- h$NUK1 * h$SUK / 1e3

# Now to merge this with the real household disposable income data. Merge it
# into the data frame `rpdi`, then use spline interpolation to make a smooth
# people-in-households time series. Notice that I take the annual PIH
# estimates to correspond to quarter 2. That's because census results are
# most often presented as representative of Q2 in the year they were taken.
rpdi <- rhdi[rhdi$Y >= min(h$YR),]
rpdi$Y <- rpdi$YR + ((rpdi$Q - 0.5) / 4)
rpdi$PIH <- NA
rpdi$PIH[rpdi$Q == 2] <- h$PIH[h$YR %in% unique(rpdi$YR)]
rpdi$PIH <- splinefun(rpdi$Y, rpdi$PIH)(rpdi$Y)
plot(rpdi$Y, rpdi$PIH, type="l",
     xlab="year", ylab="people in households (millions)")
lines(pop$Y, pop$POP, type="l", lty="dotted")
grid(col="darkgrey")

# Finally I can just divide the `RHDI` column by the `PIH` column in the
# `rpdi` data frame to get the time series I originally hoped to find online:
# quarterly real disposable income per person.
rpdi$RPDI <- rpdi$RHDI / rpdi$PIH
plot(rpdi$Y, rpdi$RPDI, type="l",
     xlab="year", ylab="real, personal dispos. income (£ per quarter)")
grid(col="darkgrey")

# But Hibbs's model is of course dependent on changes in the logarithm of
# the RDI per person, multiplied by 400 (to make it "the quarter-to-quarter
# log-percentage change expressed at annual rates"). So work that out.
# Let's call it delta ln R.
rpdi$DLR <- c(NA, 400 * diff(log(rpdi$RPDI)))

# Round the relevant columns in `rpdi` and write it to a file.
rpdi[c("PIH", "RPDI", "DLR")] <- signif(rpdi[c("PIH", "RPDI", "DLR")], 4)
write.table(rpdi, "uk_rpdi.dat", FALSE, FALSE, "\t", row.names=FALSE)
