tiddit=/media/god/programms/tiddit
sampel=/media/god/wgs/samps.txt
ref_bwa_mem=/media/god/data1/genome_indexes/bwa/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna


$tiddit --cov -o $sample -z 5000 $ref_bwa_mem $sample 
