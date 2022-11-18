gatk=/media/software/gatk-4.2.6.1/gatk

#$gatk SelectVariants \
#--variant ApplyVSQR.recalibrated.vcf \
#--output Woman2_filtered_ApplyVSQR.recalibrated.vcf \
#--exclude-filtered TRUE \
#--restrict-alleles-to BIALLELIC \
#-sn Woman2

#$gatk SelectVariants \
#--variant ApplyVSQR.recalibrated.vcf \
#--output Woman1_filtered_ApplyVSQR.recalibrated.vcf \
#--exclude-filtered TRUE \
#--restrict-alleles-to BIALLELIC \
#-sn Woman1

#$gatk SelectVariants \
#--variant ApplyVSQR.recalibrated.vcf \
#--output Woman3_filtered_ApplyVSQR.recalibrated.vcf \
#--exclude-filtered TRUE \
#--restrict-alleles-to BIALLELIC \
#-sn Woman3

#$gatk SelectVariants \
#--variant ApplyVSQRControl.recalibrated.vcf \
#--output Control_filtered_ApplyVSQR.recalibrated.vcf \
#--exclude-filtered TRUE \
#--restrict-alleles-to BIALLELIC

# applies a filter to the samples and control so only the SNPs and Indels that are expressed on one 1 X-chr
#bcftools filter -i 'GT="het"' Control_filtered_ApplyVSQR.recalibrated.vcf > het_Control.vcf

#bcftools filter -i 'GT="het"' Woman1_filtered_ApplyVSQR.recalibrated.vcf > het_Woman1.vcf

#bcftools filter -i 'GT="het"' Woman2_filtered_ApplyVSQR.recalibrated.vcf > het_Woman2.vcf

#bcftools filter -i 'GT="het"' Woman3_filtered_ApplyVSQR.recalibrated.vcf > het_Woman3.vcf

#Switch to bed format for easier handling in R
#bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD\t%SAMPLE]\n' het_Woman1.vcf > het_Woman1.bed
#bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD\t%SAMPLE]\n' het_Woman2.vcf > het_Woman2.bed
#bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD\t%SAMPLE]\n' het_Woman3vcf > het_Woman3.bed

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD\t%SAMPLE]\n' het_Control.vcf > het_Control.bed

#$gatk SelectVariants \
#--variant ApplyVSQRControl.recalibrated.vcf \
#--output filtered_ApplyVSQR.Control.recalibrated.vcf \
#--exclude-filtered TRUE
