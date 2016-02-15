# Read in raw data for household count & size, real household disposable
# income, and general elections, and mould them into usable data frames.

h <- read.table("uk_h.dat", header=TRUE)
rhdi <- read.table("uk_rhdi.dat", header=TRUE)
el <- read.table("uk_ge.dat", header=TRUE)
#mf <- read.table("uk_mf.dat", header=TRUE)
#pop <- read.table("uk_pop.dat", header=TRUE)

# Compute the incumbent's share of the two-party popular vote in each
# general election for which `el` has data. Then compute each election's
# year & quarter.
el$PVI2 <- round(100 * el$PVI / (el$PVI + el$PVC), 2)
el$DAY <- as.Date(el$DAY)
el$YR <- as.integer(format(el$DAY, "%Y"))
el$Q <- floor((as.numeric(format(el$DAY, "%m")) + 2) / 3)  # quarter
el$Y <- el$YR + (as.integer(format(el$DAY, "%j")) / 365.25)

# Fill in gaps in the household numbers from my favoured spreadsheet by
# extrapolating from the England & Wales (and sometimes Scotland) numbers.

impute1 <- lm(NNI ~ NE + NW + NS, h)
impute2 <- lm(NS ~ NE + NW, h)
impute3 <- lm(NNI ~ NE + NW, h)

h$NNI[is.na(h$NNI)] <- predict(impute1, h[is.na(h$NNI),])
h$NS[is.na(h$NS)] <- predict(impute2, h[is.na(h$NS),])
h$NNI[is.na(h$NNI)] <- predict(impute3, h[is.na(h$NNI),])
h[c("NS", "NNI")] <- round(h[c("NS", "NNI")])
h$NUK1[is.na(h$NUK1)] <- rowSums(h[c("NE", "NW", "NS", "NNI")])[is.na(h$NUK1)]

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

# Show the average-household-size time series.

plo_s <- function(y_lab="average household size")
{
	co <- "#00000070"
	plot(h$YR, h$SE, type="b", ylim=c(2.2, 3.8), col=co,
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
			text(x, 0.07 + (h[[col_name]])[h$YR == x], col_name)
		} else if (col_name == "SW") {
			text(1991, 0.08 + (h[["SW"]])[h$YR == 1991], "SW")
		} else {
			text(2001, 0.08 + (h[["SUK"]])[h$YR == 2001], "SUK")
		}
	}
}

par(las=1, mar=c(4.5, 4.1, 0.2, 0.2), mfrow=c(2,1))

# Do a pre-imputation plot of average household sizes.
plo_s()

# Time to adjust & impute like mad.
# 1. There's a really dodgy-looking jump in the English average for 1972
#    through 1980. Push those years' averages down to bring them in line with
#    the surrounding data.
# 2. Make the hopeful assumption that the gap between Welsh & English average
#    household sizes was constant before 1991, imputing Welsh AHS backwards
#    decadally from that.
# 3. Wedge the years 1962 through 1970 into `h` by copying rows and blanking
#    them, then linearly interpolate an annual time series for each country
#    in the UK for 1961 onwards.
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

h <- rbind(head(h, 9), h)
h[1:9,] <- NA
h$YR[1:9] <- 1962:1970
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

# Do a post-imputation plot of household sizes.
plo_s("imputed average household size")

# Finally the `SUK` is a complete time series of (mostly imputed!) average
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
grid(col="darkgrey")

# Finally I can just divide the RHDI column by the PIH column in the `rpdi`
# data frame to obtain the time series I originally hoped to find online:
# quarterly real disposable income per person.
rpdi$RPDI <- rpdi$RHDI / rpdi$PIH
plot(rpdi$Y, rpdi$RPDI, type="l",
     xlab="year", ylab="real, personal dispos. income (£ per quarter)")
grid(col="darkgrey")

# But Hibbs's model is of course dependent on changes in the logarithm of
# the RDI per person, multiplied by 400 (to make it "the quarter-to-quarter
# log-percentage change expressed at annual rates". So work that out.
# Let's call it delta ln R.
rpdi$DLR <- c(NA, 400 * diff(log(rpdi$RPDI)))

# Round the relevant columns in `rpdi` and write it to a file.
rpdi[c("PIH", "RPDI", "DLR")] <- signif(rpdi[c("PIH", "RPDI", "DLR")], 4)
write.table(rpdi, "uk_rpdi.dat", FALSE, FALSE, "\t", row.names=FALSE)
