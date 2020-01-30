#!/bin/sh
#
#
outfile="meta-analysis.n1.$1.n2.$2.seed$3"

R CMD BATCH --vanilla "--args $1 $2 $3" meta-analysis.R $outfile

echo "Outfile is `pwd`/$outfile" | mail -v -s "meta finished" jiaqiyin@uw.edu
