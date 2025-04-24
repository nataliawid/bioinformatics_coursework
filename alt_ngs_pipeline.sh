#!/bin/bash

R1="NGS0001.R1.fastq.qz"
R2="NGS0001.R2.fastq.qz"


TRIMMED_R1="NGS0001.trimmed.R1.fastq"
TRIMMED_R2="NGS0001.trimmed.R2.fastq"

REF="hg19.fa"       
BED="annotation.bed"
SAMPLE="NGS0001"
RGID="NGS0001"


FASTP_REPORT_DIR="fastp_reports"
FASTQC_REPORT_DIR="fastqc_reports"
OUTDIR="alignment_output"

mkdir -p $FASTP_REPORT_DIR $FASTQC_REPORT_DIR $OUTDIR

echo "Running fastp..."
fastp \
  -i $R1 -I $R2 \
  -o $TRIMMED_R1 -O $TRIMMED_R2 \
  --html $FASTP_REPORT_DIR/fastp_report.html \
  --json $FASTP_REPORT_DIR/fastp_report.json \
  --thread 4


echo "Running FastQC..."
fastqc $TRIMMED_R1 $TRIMMED_R2 -o $FASTQC_REPORT_DIR

echo "Running BWA MEM..."
bwa mem -R "@RG\tID:$RGID\tSM:$SAMPLE\tPL:ILLUMINA" $REF $TRIMMED_R1 $TRIMMED_R2 | \
  samtools sort -o $OUTDIR/aligned.bam


echo "Marking duplicates..."
picard MarkDuplicates \
  I=$OUTDIR/aligned.bam \
  O=$OUTDIR/dedup.bam \
  M=$OUTDIR/dedup.metrics.txt

echo "Filtering BAM by MAPQ >= 20..."
samtools view -b -q 20 $OUTDIR/dedup.bam > $OUTDIR/dedup.qf.bam

echo "Generating alignment stats..."
samtools flagstat $OUTDIR/dedup.qf.bam > $OUTDIR/flagstats.txt
samtools idxstats $OUTDIR/dedup.qf.bam > $OUTDIR/idxstats.txt
samtools depth -a $OUTDIR/dedup.qf.bam > $OUTDIR/depth.txt
picard CollectInsertSizeMetrics \
  I=$OUTDIR/dedup.qf.bam \
  O=$OUTDIR/insert_size_metrics.txt \
  H=$OUTDIR/insert_size_histogram.pdf \
  M=0.5


echo "Calling variants with bcftools..."
bcftools mpileup -f $REF -R $BED -Ou $OUTDIR/dedup.qf.bam | \
  bcftools call -mv -Ob -o $OUTDIR/variants.bcf

bcftools view $OUTDIR/variants.bcf > $OUTDIR/variants.vcf

echo "Filtering variants..."
bcftools filter -e 'QUAL<20 || DP<10' $OUTDIR/variants.vcf > $OUTDIR/variants.filtered.vcf

echo "Annotating variants..."
table_annovar.pl $OUTDIR/variants.filtered.vcf humandb/ \
  -buildver hg19 \
  -out $OUTDIR/annovar_output \
  -remove \
  -protocol refGene \
  -operation g \
  -nastring . \
  -vcfinput

snpEff hg19 $OUTDIR/variants.filtered.vcf > $OUTDIR/snpeff_output.vcf

echo "Prioritizing variants..."
grep -v "dbSNP" $OUTDIR/annovar_output.hg19_multianno.txt | grep "exonic" > $OUTDIR/prioritized_variants.txt

echo "Pipeline completed successfully."