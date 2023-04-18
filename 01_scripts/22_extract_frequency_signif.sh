#!/bon/bash

vcfin=$1

vcftools --gzvcf $vcfin \
  --keep 00_maps/00_maps/sex/M\
  --freq2\
  --out 15_exploring_SNP-by-SNP/Shared_two_method_male\
  --positions 98_SNP_lists/Shared_signif_snps

vcftools --gzvcf $vcfin \
  --keep 00_maps/00_maps/sex/F\
  --freq2\
  --out 15_exploring_SNP-by-SNP/Shared_two_method_female\
  --positions 98_SNP_lists/Shared_signif_snp
  
  vcftools --gzvcf $vcfin \
  --keep 00_maps/00_maps/sex/F\
  --freq2\
  --out 15_exploring_SNP-by-SNP/All_two_method_female\
  --positions 98_SNP_lists/All_signif_snps
  
  vcftools --gzvcf $vcfin \
  --keep 00_maps/00_maps/sex/M\
  --freq2\
  --out 15_exploring_SNP-by-SNP/All_two_method_male\
  --positions 98_SNP_lists/All_signif_snps
