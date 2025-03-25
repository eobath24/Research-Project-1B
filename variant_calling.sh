#!/bin/bash


# Variant calling commands
# Arg 1 is assembly file
# Arg 2 & 3 are illumina reads
# Arg 4 is name of strain


# Mapping and calling varaints
minimap2 -a -x sr $1 $2 $3 | samtools view -h -F 0x900 - | samtools sort -O bam > $4_mapped_reads.bam
samtools flagstat $4_mapped_reads.bam
bcftools mpileup -Ou -f $1 $4_mapped_reads.bam| bcftools call -vc -Ov > $4_variants.vcf

# Consensus sequence generation of 10kb region of the tgr locus
bcftools view -I $4_variants.vcf -Oz -o $4_variants.vcf.gz
bcftools index $4_variants.vcf.gz
samtools faidx $1 contig_1224:5000-14999 | bcftools consensus $4_variants.vcf.gz > $4_consensus_1224_10kb.fasta



# Depths and coverage of each $4
samtools depth -a $4_mapped_reads.bam > $4_depth.txt
python get_depths.py -i $4_depth.txt -o $4_coverage.bed



# View the mapped reads and the specific SNPs
samtools index $4_mapped_reads.bam
samtools tview $4_mapped_reads.bam $1

