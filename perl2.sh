#!/usr/bin/perl

# The raw .g.vcf file is then run through GenotypeGVCFs which performs joint genotyping on all samples together. step 2
####################################################################################################################

use warnings;

use strict;

my $gatkcommand = "/media/god/data1/software/gatk-4.2.6.1/gatk --java-options '-Xmx64g -Djava.io.tmpdir=/tmp' GenotypeGVCFs -R /media/god/data1/genome_indexes/bwa/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -G StandardAnnotation -G AS_StandardAnnotation";

my $contigfile = 'contigs.txt';

open my $contig, $contigfile or die "Could not open $contigfile: $!";

while( my $line = <$contig>)  {

    chomp($line);

        my $gatkfinalcommand = $gatkcommand." --variant gendb://".$line." -O ./fin/All_".$line.".vcf.gz";

        system($gatkfinalcommand);

}
