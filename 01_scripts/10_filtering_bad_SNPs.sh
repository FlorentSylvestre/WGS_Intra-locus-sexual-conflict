#!/bin/bash

invcf=$1
badsnps="/97_SNPs_lists/Sex_duplications/BAD_Snps"

bcftools view -m2 -M2 -T ^$badsnps -Oz -o 04_datasets/Full_genome/CLEAN_$(basename $i)) $i
