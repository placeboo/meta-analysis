#!/bin/sh
#
#
outfile="crt.n.$1.seed$2"

R CMD BATCH --vanilla "--args $1 $2" crt.R $outfile

echo "Outfile is `pwd`/$outfile" | mail -v -s "crt finished" jiaqiyin@uw.edu
