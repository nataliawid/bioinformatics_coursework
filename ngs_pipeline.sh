#!/bin/bash

# Files
R1="NGS0001.R1.fastq.qz"
R2="NGS0001.R2.fastq.qz"

# Output
TRIMMED_R1="NGS0001.trimmed.R1.fastq"
TRIMMED_R2="NGS0001.trimmed.R2.fastq"

# Output directories
FASTP_REPORT_DIR="fastp_reports"
FASTQC_REPORT_DIR="fastqc_reports"

mkdir -p $FASTP_REPORT_DIR $FASTQC_REPORT_DIR

# Perform quality assessment and trimming using fastp
echo "Running fastp..."
fastp \
  -i $R1 -I $R2 \
  -o $TRIMMED_R1 -O $TRIMMED_R2 \
  --html $FASTP_REPORT_DIR/fastp_report.html \
  --json $FASTP_REPORT_DIR/fastp_report.json \
  --threads 4

# Perform basic quality assessment of paired trimmed sequencing data using FastQC
echo "Running FastQC..."
fastqc $TRIMMED_R1 $TRIMMED_R2 -o $FASTQC_REPORT_DIR

# Variables
R1="NGS0001.trimmed.R1.fastq"
R2="NGS0001.trimmed.R2.fastq"
REF="hg19.fa"  # Make sure the BWA index is built
SAMPLE="NGS0001"
RGID="NGS0001"  # Read group ID
OUTDIR="alignment_output"
BED="annotation.bed"

mkdir -p $OUTDIR

# Align with BWA MEM including read group information
echo "Aligning with BWA..."
bwa mem -M -R "@RG\tID:$RGID\tSM:$SAMPLE\tPL:ILLUMINA" $REF $R1 $R2 \
  | samtools view -Sb - > $OUTDIR/${SAMPLE}.raw.bam

# Sort and mark duplicates
echo "Sorting and marking duplicates..."
samtools sort -o $OUTDIR/${SAMPLE}.sorted.bam $OUTDIR/${SAMPLE}.raw.bam
picard MarkDuplicates \
  I=$OUTDIR/${SAMPLE}.sorted.bam \
  O=$OUTDIR/${SAMPLE}.dedup.bam \
  M=$OUTDIR/${SAMPLE}.dup_metrics.txt \
  CREATE_INDEX=true

# Quality filtering (MAPQ >= 30)
echo "Filtering BAM file (MAPQ >= 30)..."
samtools view -b -q 30 $OUTDIR/${SAMPLE}.dedup.bam > $OUTDIR/${SAMPLE}.filtered.bam
samtools index $OUTDIR/${SAMPLE}.filtered.bam

# Generate alignment statistics
echo "Generating statistics..."

# Flagstats
samtools flagstat $OUTDIR/${SAMPLE}.filtered.bam > $OUTDIR/${SAMPLE}.flagstat.txt

# Idxstats
samtools idxstats $OUTDIR/${SAMPLE}.filtered.bam > $OUTDIR/${SAMPLE}.idxstats.txt

# Depth of coverage
if [[ -f "$BED" ]]; then
  bedtools coverage -a $BED -b $OUTDIR/${SAMPLE}.filtered.bam > $OUTDIR/${SAMPLE}.coverage.txt
else
  echo "No annotation.bed found, skipping BED coverage."
fi

# Insert size metrics
picard CollectInsertSizeMetrics \
  I=$OUTDIR/${SAMPLE}.filtered.bam \
  O=$OUTDIR/${SAMPLE}.insertsize.txt \
  H=$OUTDIR/${SAMPLE}.insertsize_histogram.pdf \
  M=0.5

# Variant Calling using Freebayes (restricted to BED regions)
echo "Calling variants with Freebayes..."
freebayes -f $REF -t $BED $OUTDIR/${SAMPLE}.filtered.bam > $OUTDIR/${SAMPLE}.raw.vcf

# Quality Filtering of Variants
echo "Filtering variants (QUAL >= 30, DP >= 10)..."
bcftools filter -i 'QUAL>=30 && DP>=10' $OUTDIR/${SAMPLE}.raw.vcf -o $OUTDIR/${SAMPLE}.filtered.vcf

# Variant Annotation with ANNOVAR
echo "Annotating with ANNOVAR..."
# Convert VCF to ANNOVAR input
convert2annovar.pl -format vcf4 $OUTDIR/${SAMPLE}.filtered.vcf -outfile $OUTDIR/${SAMPLE}.avinput

# Annotate using ANNOVAR (you'll need to download the proper hg19 databases)
table_annovar.pl $OUTDIR/${SAMPLE}.avinput humandb/ -buildver hg19 \
  -out $OUTDIR/${SAMPLE}.annovar \
  -remove -protocol refGene,dbnsfp42a,avsnp150 \
  -operation g,f,f -nastring . -vcfinput

# Variant Annotation with SnpEff
echo "Annotating with SnpEff..."
snpEff hg19 $OUTDIR/${SAMPLE}.filtered.vcf > $OUTDIR/${SAMPLE}.snpEff.vcf

# Variant Prioritization – Exonic, not in dbSNP
echo "Filtering for exonic variants not in dbSNP..."
# Assumes avsnp150 column contains dbSNP ID. We'll pull exonic and no-snp entries from ANNOVAR output.
awk '($2 == "exonic" || $2 ~ /^exonic/) && $6 == "."' $OUTDIR/${SAMPLE}.annovar.hg19_multianno.txt > $OUTDIR/${SAMPLE}.prioritized.txt

echo "All alignment steps completed."