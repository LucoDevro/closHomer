#! /bin/bash

genomes=$1
outdir=$2

mkdir -p $outdir
dir -1 $genomes | grep -E "\.fna$" | xargs -I % echo $genomes/% > genomes_list
fastANI --ql genomes_list --rl genomes_list -t 12 -o $outdir/out --minFraction 0.05
rm genomes_list
