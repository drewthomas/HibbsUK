uk_rpdi.dat: uk_rhdi.dat uk_h.dat merge.R
	R -q --vanilla < merge.R

uk_rhdi.dat: UKEA_CSDB_DS.csdb.csv Makefile
	/bin/echo -e "YR\tQ\tRHDI" > uk_rhdi.dat
	grep " Q" UKEA_CSDB_DS.csdb.csv | tr -d 'Q"' | tr ', ' ' \t' \
		>> uk_rhdi.dat
