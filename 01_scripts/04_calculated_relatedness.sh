/bin/bash

module load vcftools #To run if running on a slurm server

INPUT_VCF="04_datasets/Full_genome/NO_INV_COV_MISS_QUAL_MAF010_filtered_99_autosomes_biallelic_SOR_sorted_renamed_NoUn.vcf.gz"
OUTPUT="06_relatedness/NO_INV_COV_MISS_QUAL_MAF010_filtered_99_autosomes_biallelic_SOR_sorted_renamed_NoUn

vcftools --gzvcf $INPUT_VCF --out $OUTPUT --relatedness2

