#!/bin/bash
THREADS=${SLURM_CPUS_PER_TASK:-1}
set -e # stop on errors
set -x # # print each commmand before executing 

# Step 1: Prepare the reference genome
cd /hb/groups/kay_lab/popoolation2/test
mkdir -p ref

# if the reference fasta isn't already in ref/, move it there
if [ ! -f ref/refgenome.bwt ]; then
  mv refgenome ref/
fi

# if no BWA index exists, create one
if [ ! -f ref/refgenome.bwt ]; then
  bwa index ref/refgenome
fi

# step 2: map the reads to the reference genome
mkdir -p map

# map pop2p
if [ ! -f map/pop2p.sam ]; then
  bwa mem -t $THREADS ref/refgenome pop2p.fastq > map/pop2p.sam
fi

# map pop2w
if [ ! -f map/pop2w.sam ]; then
  bwa mem -t $THREADS ref/refgenome pop2w.fastq > map/pop2w.sam
fi

# step 3: remove ambiguously mapped reads
if [ ! -f map/pop2p.sorted.bam ]; then
  samtools view -@ $THREADS -q 20 -b map/pop2p.sam | samtools sort -@ $THREADS -o map/pop2p.sorted.bam -T map/pop2p.temp
fi

if [ ! -f map/pop2w.sorted.bam ]; then
  samtools view -@ $THREADS -q 20 -b map/pop2w.sam | samtools sort -@ $THREADS -o map/pop2w.sorted.bam -T map/pop2w.temp
fi

# Generate mpileup file from both populations  
if [ ! -f p1_p2.mpileup ]; then
  samtools mpileup -f ref/refgenome map/pop2p.sorted.bam map/pop2w.sorted.bam > p1_p2.mpileup
fi

# Step 4: Create a syncronized file
if [ ! -f p1_p2.sync ]; then
  java -ea -Xmx7g -jar /hb/groups/kay_lab/popoolation2/pool2/popoolation2_1201/mpileup2sync.jar \
  --input p1_p2.mpileup \
  --output p1_p2.sync \
  --fastq-type sanger \
  --min-qual 20 \
  --threads 8
fi

# Step 5: Calculate allele frequency differences
if [ ! -f p1_p2.pwc ]; then
  perl /hb/groups/kay_lab/popoolation2/pool2/popoolation2_1201/snp-frequency-diff.pl \
  --input p1_p2.sync \
  --output-prefix p1_p2 \
  --min-count 6 \
  --min-coverage 50 \
  --max-coverage 200
fi

# This step results in two output files, _rc, Read Counts which contains the major and minor alleles per SNP, counts, and other data (chr, position) and _pwc (Pariwise Comparison, contains allele frequency differences between all pairs of populations).

# Step 6: Fst-values: measure differentiation between populations, Calculate Fst for every SNP
if [ ! -f p1_p2.fst ]; then
  perl /hb/groups/kay_lab/popoolation2/pool2/popoolation2_1201/fst-sliding.pl \
  --input p1_p2.sync \
  --output p1_p2.fst \
  --suppress-noninformative \
  --min-count 6 \
  --min-coverage 50 \
  --max-coverage 200 \
  --min-covered-fraction 1 \
  --window-size 1 \
  --step-size 1 \
  --pool-size 500
fi

if [ ! -f p1_p2_w500.fst ]; then
  perl /hb/groups/kay_lab/popoolation2/pool2/popoolation2_1201/fst-sliding.pl \
  --input p1_p2.sync \
  --output p1_p2_w500.fst \
  --min-count 6 \
  --min-coverage 50 \
  --max-coverage 200 \
  --min-covered-fraction 1 \
  --window-size 500 \
  --step-size 500 \
  --pool-size 500 
fi

# Step 7: index BAMs and prep for IGV 
if [ ! -f map/pop2p.sorted.bam.bai ]; then
samtools index map/pop2p.sorted.bam
fi
if [ ! -f map/pop2w.sorted.bam.bai ]; then
samtools index map/pop2w.sorted.bam
fi

if [ ! -f p1_p2.igv ]; then
  perl /hb/groups/kay_lab/popoolation2/pool2/popoolation2_1201/export/pwc2igv.pl \
  --input p1_p2.fst \
  --output p1_p2.igv
  fi

# Step 8: Estimate the significance of allele frequency differences using Fisher's Exact Test 
if [ ! -f p1_p2.fet ]; then
  perl /hb/groups/kay_lab/popoolation2/pool2/popoolation2_1201/fisher-test.pl \
  --input p1_p2.sync \
  --output p1_p2.fet \
  --min-count 6 \
  --min-coverage 50 \
  --max-coverage 200 \
  --suppress-noninformative
fi

# Load the above results into IGV file
if [ ! -f p1_p2_fet.igv ]; then
  perl /hb/groups/kay_lab/popoolation2/pool2/popoolation2_1201/export/pwc2igv.pl \
  --input p1_p2.fet \
  --output p1_p2_fet.igv
fi
