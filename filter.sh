gatk=/media/god/data1/software/gatk-4.2.6.1/gatk

#$gatk SelectVariants \
#--variant ApplyVSQR.recalibrated.vcf \
#--output 13PLJ_filtered_ApplyVSQR.recalibrated.vcf \
#--exclude-filtered TRUE \
#--restrict-alleles-to BIALLELIC \
#-sn GTEX-13PLJ-0003-SM-6WSSCN

#$gatk SelectVariants \
#--variant ApplyVSQR.recalibrated.vcf \
#--output UPIC_filtered_ApplyVSQR.recalibrated.vcf \
#--exclude-filtered TRUE \
#--restrict-alleles-to BIALLELIC \
#-sn GTEX-UPIC-0004-SM-5SOEF

#$gatk SelectVariants \
#--variant ApplyVSQRControl.recalibrated.vcf \
#--output Control_filtered_ApplyVSQR.recalibrated.vcf \
#--exclude-filtered TRUE \
#--restrict-alleles-to BIALLELIC


#bcftools filter -i 'GT="het"' Control_filtered_ApplyVSQR.recalibrated.vcf > het_Control.vcf

#bcftools filter -i 'GT="het"' UPIC_filtered_ApplyVSQR.recalibrated.vcf > het_UPIC.vcf

#bcftools filter -i 'GT="het"' ZZPU_filtered_ApplyVSQR.recalibrated.vcf > het_ZZPU.vcf


#bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD\t%SAMPLE]\n' het_13PLJ.vcf > het_13PLJ.bed
#bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD\t%SAMPLE]\n' het_UPIC.vcf > het_UPIC.bed
#bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD\t%SAMPLE]\n' het_ZZPU.vcf > het_ZZPU.bed

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD\t%SAMPLE]\n' het_Control.vcf > het_Control.bed

#$gatk SelectVariants \
#--variant ApplyVSQRControl.recalibrated.vcf \
#--output filtered_ApplyVSQR.Control.recalibrated.vcf \
#--exclude-filtered TRUE
