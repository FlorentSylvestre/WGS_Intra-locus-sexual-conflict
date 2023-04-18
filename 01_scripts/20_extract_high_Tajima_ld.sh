#!/bin/bash

##User input
Global_vcf=$1

#Fixed variable
List_win="11_Res_permuts/windows_balancing_selection"
LD_path="14_ld"
Output_vcf="04_datasets/High_Tajima"


while read line
 do
        echo $line
        region=$( echo $line | sed "s/ /_/" | sed "s/ /-/")
        echo $region
        bcftools view -r ${region/_/:} -Oz -o $Output_vcf/$region.vcf.gz $Global_vcf
        plink --double-id\
              --ld-window 999999\
              --ld-window-kb 1000000\
              --ld-window-r2 0\
              --out $LD_path/$region\
              --r2\
              --vcf $Output_vcf/$region.vcf.gz\


 done<<(grep Shared $List_win)
