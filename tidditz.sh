#### References and variants #####
ref_GATK=/media/user/redacted/indexes_annotations/bwa/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna



cat /media/user/redacted/wgs/runfile.txt | while read participant    BAM_file
do
tiddit --sv --threads 100 -o /media/user/redacted/wgs/tiddit_out/$participant. --bam /media/user/redacted/wgs/$BAM_file --ref $ref_GATK
done



cat /media/user/redacted/wgs/runfile.txt | while read participant    BAM_file
do
tiddit --cov -o /media/user/redacted/wgs/tiddit_out/cov/$participant'_'cov. --bam /media/user/redacted/wgs/$BAM_file
done
