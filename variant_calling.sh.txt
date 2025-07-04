#################################################################################################################
Variant calling script
#################################################################################################################


# Mapping and calling variants
minimap2 -a -x sr $assembly $illuminaR1 $illuminaR2 | samtools view -h -F 0x900 - | samtools sort -O bam > $filename_mapped_reads.bam
samtools flagstat $filename_mapped_reads.bam
bcftools mpileup -Ou -f $assembly $filename_mapped_reads.bam| bcftools call -vc -Ov > $filename_variants.vcf



# Consensus sequence generation of 10kb region of the tgr locus
bcftools view -I $filename_variants.vcf -Oz -o $filename_variants.vcf.gz
bcftools index $filename_variants.vcf.gz
samtools faidx $assembly contig_1224:5000-14999 | bcftools consensus $filename_variants.vcf.gz > $filename_consensus_1224_10kb.fasta



# Depths and coverage of each strain
samtools depth -a $filename_mapped_reads.bam > $filename_depth.txt
python get_depths.py -i $filename_depth.txt -o $filename_coverage.bed



# View the mapped reads and the specific SNPs
samtools index $filename_mapped_reads.bam
samtools tview $filename_mapped_reads.bam $1





