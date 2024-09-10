# SNP Calling Pipeline using HPCC System at TTU

This pipeline details the SNP calling process for S.nigra, using the TTU HPCC system. 

---

## 1. MD5 Check

```bash
find -name "*.gz" -exec md5sum '{}' \; > md5sum.txt
```

---

## 2. Unzip All Files

```bash
gunzip -d $(find ./ -type f -name '*.gz')
```

---

## 3. Copy Unzipped Files to a New Folder

```bash
cp -d $(find ./ -type f -name '*.gz') /lustre/scratch/all_in_one
```

---

## 4. Create a List of File Names and Remove Duplicates

```bash
ls *.fastq | sed 's/.fastq//g' | sed 's/_R1//g' | sed 's/_R2//g' | sed 's/reads\///g' > allnames_1.txt
sort allnames_1.txt | uniq > allnames.txt
wc -l allnames.txt
```

---

## 5. Convert FASTQ to Unmapped BAM (uBAM)

```bash
#!/bin/bash
#SBATCH -J Pmex_snp
#SBATCH -N 1
#SBATCH --ntasks-per-node=64
#SBATCH -o log/01/array-%A_%a.out
#SBATCH -e log/01/array-%A_%a.err
#SBATCH -p nocona
#SBATCH -a 1-20:1

prefixNum="$SLURM_ARRAY_TASK_ID"
prefixNumP="$prefixNum"p
prefix=$(sed -n "$prefixNumP" allnames.txt)

gatk FastqToSam --java-options "-Xmx8G" \
    --FASTQ ./"$prefix"_R1.fastq \
    --FASTQ2 ./"$prefix"_R2.fastq \
    --OUTPUT uBAM/"$prefix".unmapped.bam \
    --SAMPLE_NAME "$prefix" \
    --LIBRARY_NAME NEBNext-"$prefix" \
    --PLATFORM illumina \
    --SEQUENCING_CENTER OMRF
```

---

## 6. Align FASTQ to Reference Genome (Aligned BAM)

```bash
#!/bin/bash
#SBATCH -J Pmex_snp
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH -o log/02/array-%A_%a.out
#SBATCH -e log/02/array-%A_%a.err
#SBATCH -p nocona
#SBATCH -a 1-20:1

prefixNum="$SLURM_ARRAY_TASK_ID"
prefixNumP="$prefixNum"p
prefix=$(sed -n "$prefixNumP" allnames.txt)

gatk --java-options "-Dsamjdk.compression_level=5 -Xms3000m" SamToFastq \
    --INPUT uBAM/"$prefix".unmapped.bam \
    --FASTQ /dev/stdout --INTERLEAVE true --NON_PF true | \
    bwa mem -K 100000000 -p -v 3 -t 16 -Y /lustre/work/minguo/genome/purpurea_v5/Salix_purpurea_var_94006.mainGenome-noZ.fasta \
    /dev/stdin - | \
    samtools view -1 - > aBAM/"$prefix".aligned.bam
```

---

## 7. Merge BAM Files

```bash
#!/bin/bash
#SBATCH -J Pmex_snp
#SBATCH -N 1
#SBATCH --ntasks-per-node=8
#SBATCH -o log/03/array-%A_%a.out
#SBATCH -e log/03/array-%A_%a.err
#SBATCH -p nocona
#SBATCH -a 1-20:1

prefixNum="$SLURM_ARRAY_TASK_ID"
prefixNumP="$prefixNum"p
prefix=$(sed -n "$prefixNumP" allnames.txt)

gatk --java-options "-Dsamjdk.compression_level=5 -Xms3000m" MergeBamAlignment \
    --VALIDATION_STRINGENCY SILENT --EXPECTED_ORIENTATIONS FR --ATTRIBUTES_TO_RETAIN X0 \
    --ALIGNED_BAM aBAM/"$prefix".aligned.bam \
    --UNMAPPED_BAM uBAM/"$prefix".unmapped.bam \
    --OUTPUT mBAM/"$prefix".merged.bam \
    --REFERENCE_SEQUENCE /lustre/work/minguo/genome/purpurea_v5/Salix_purpurea_var_94006.mainGenome-noZ.fasta
```

---

## 8. Sort and Fix BAM

```bash
#!/bin/bash
#SBATCH -J Pmex_snp
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH -o log/04/array-%A_%a.out
#SBATCH -e log/04/array-%A_%a.err
#SBATCH -p nocona
#SBATCH -a 1-20:1

prefixNum="$SLURM_ARRAY_TASK_ID"
prefixNumP="$prefixNum"p
prefix=$(sed -n "$prefixNumP" allnames.txt)

gatk --java-options "-Dsamjdk.compression_level=5 -Xms4000m" SortSam \
    --INPUT mBAM/"$prefix".merged.bam --OUTPUT /dev/stdout --SORT_ORDER 'coordinate' | \
    gatk --java-options "-Dsamjdk.compression_level=5 -Xms500m" SetNmMdAndUqTags \
    --INPUT /dev/stdin --OUTPUT sBAM/"$prefix".sorted.bam \
    --CREATE_INDEX true --CREATE_MD5_FILE true
```

---


## 9. Mark Duplicates

```bash
#!/bin/bash
#SBATCH -J Pmex_snp
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH -o log/05/array-%A_%a.out
#SBATCH -e log/05/array-%A_%a.err
#SBATCH -p nocona
#SBATCH -a 1-20:1

prefixNum="$SLURM_ARRAY_TASK_ID"
prefixNumP="$prefixNum"p
prefix=$(sed -n "$prefixNumP" allnames.txt)

gatk --java-options "-Dsamjdk.compression_level=5 -Xms4000m" MarkDuplicates \
    --INPUT sBAM/"$prefix".sorted.bam --OUTPUT dBAM/"$prefix".dedup.bam \
    --METRICS_FILE dBAM/"$prefix".dedup.metrics \
    --CREATE_INDEX true --CREATE_MD5_FILE true
```

---

## 10. Create Directories for Each Individual

```bash
while read -r i; do mkdir /lustre/scratch/minguo/pmex/raw/"$i"; done < allnames.txt
while read -r i; do mkdir /lustre/scratch/minguo/pmex/raw/log/"$i"; done < allnames.txt
```

---

## 11. Call Variants using HaplotypeCaller

```bash
#!/bin/bash
#SBATCH -J Pmex_snp
#SBATCH -N 1
#SBATCH --ntasks-per-node=32
#SBATCH -o log/try/array-%A_%a.out
#SBATCH -e log/try/array-%A_%a.err
#SBATCH -p nocona
#SBATCH -a 1-20:1

prefixNum="$SLURM_ARRAY_TASK_ID"
prefixNumP="$prefixNum"p
prefix=$(sed -n "$prefixNumP" allnames.txt)

gatk HaplotypeCaller -R /lustre/work/minguo/genome/purpurea_v5/Salix_purpurea_var_94006.mainGenome-noZ.fasta \
    -I dBAM/"$prefix".dedup.bam -ERC GVCF --heterozygosity 0.015 --indel-heterozygosity 0.01 \
    -O /lustre/scratch/minguo/pmex/raw/"$prefix".g.vcf
```

---

## 12. Combine GVCFs into a Genomics Database

```bash
while read -r i; do printf "$i\t/lustre/scratch/minguo/pmex/raw/$i.g.vcf\n"; done < allnames.txt > all.sample_map
```

```bash
#!/bin/bash
#SBATCH -J Pmex_snp
#SBATCH -N 1
#SBATCH --ntasks-per-node=32
#SBATCH -o log/06/array-%A_%a.out
#SBATCH -e log/06/array-%A_%a.err
#SBATCH -p nocona

gatk --java-options "-Xmx4g -Xms4g" GenomicsDBImport \
    --sample-name-map all.sample_map --genomicsdb-workspace-path /lustre/scratch/minguo/genodb \
    -L salixIntervals.list
```

---

## 13. Genotype the GVCFs

```bash
#!/bin/bash
#SBATCH -J Pmex_snp
#SBATCH -N 1
#SBATCH --ntasks-per-node=32
#SBATCH -o log/07/array-%A_%a.out
#SBATCH -e log/07/array-%A_%a.err
#SBATCH -p nocona

gatk --java-options "-Xmx5g -Xms5g" GenotypeGVCFs \
    -R /lustre/work/minguo/genome/purpurea_v5/Salix_purpurea_var_94006.mainGenome-noZ.fasta \
    -O /lustre/scratch/minguo/pmex/raw/raw_vcf/all.genotypes.raw.vcf \
    -G StandardAnnotation \
    --use-new-qual-calculator \
    -V gendb:///lustre/scratch/minguo/genodb \
    --heterozygosity 0.015 \
    --indel-heterozygosity 0.01 \
    --stand-call-conf 30.0
```
