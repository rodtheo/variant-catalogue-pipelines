# Benchmarking Variant Catalogue Pipelines

This repo contains code to compare pipelines designed to generate variant catalogues, a list of variants and their frequencies in a population, from whole genome sequences.

The following pipelines were compared:

1. Sarek (nextflow)
2. Variant Catalogue (nextflow)
3. CAFFEE (nextflow)

## Dataset

The test data used in this experiment was downloaded from public repositories. In order to accelerate the download rate we priorize the download of files stored in aws cloud buckets. The downloads were done using aws-cli tool.

1. Human Reference Genome (GRCh38) from GATK Bundle.
```
$ aws s3 --no-sign-request cp s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta* .
```
2. Databases for variant annotation: VEP version 110 and SNPeff version 105. 
```
https://annotation-cache.github.io/vep_cache/
```
3. Samples to call variants: as it's a test dataset, we extracted on-the-fly the reads mapped to locations 10,000,000 to 20,000,000 of chromossomes 12, X and Y from BAM files of GIAB NA24694, NA12878 and NA24385. After the BAM subset we convert it to paired-end fastq. This was done with:

First, create the following BED named `include_regions.bed`:

```
20	10000000	20000000
X	10000000	20000000
Y	10000000	20000000
```

- NA24694
```
$ samtools view --threads 4 -h -L include_regions_chr.bed --use-index s3://giab/data/ChineseTrio/HG006_NA24694-huCA017E_father/NA24694_Father_HiSeq100x/NHGRI_Illumina100X_Chinesetrio_novoalign_bams/HG006.GRCh38_full_plus_hs38d1_analysis_set_minus_alts.100x.bam | samtools view --threads 4 -Shb - | samtools sort --threads 4 -n -O bam -o test_sorted_NA24694.bam -

$ bedtools bamtofastq -i test_sorted_NA24694.bam -fq NA24694.R1.fastq -fq2 NA24694.R2.fastq
```


- NA12878
```
$ samtools view --threads 4 -h -L include_regions_chr.bed --use-index s3://giab/data/NA12878/NIST_NA12878_HG001_HiSeq_300x/RMNISTHS_30xdownsample.bam | samtools view --threads 4 -Shb - | samtools sort --threads 4 -n -O bam -o test_sorted_NA12878.bam -

$ bedtools bamtofastq -i test_sorted_NA12878.bam -fq NA12878.R1.fastq -fq2 NA12878.R2.fastq
```


- NA24385
```
$ samtools view --threads 4 -h -L include_regions_chr.bed --use-index s3://giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.GRCh38.60x.1.bam | samtools view --threads 4 -Shb - | samtools sort --threads 4 -n -O bam -o test_sorted_NA24385.bam -

$ bedtools bamtofastq -i test_sorted_NA24385.bam -fq NA24385.R1.fastq -fq2 NA24385.R2.fastq
```

## Pipelines Execution

### Sarek

```
git clone

cd 

nextflow run . -params-file params.yaml -profile docker --max_cpus 6 --max_memory 12GB
```

Where custom parameters can be found at `sarek/params.yaml`.

### Variant Catalogue


```
git clone

cd

nextflow run . --input assets/local_samplesheet.csv -c conf/test_local.conf -profile docker --outdir test_out_local -resume --cache_path /home/rodtheo/Bioinfo/project-data/nextflow/sarek/data/ref/vep_cache --cache_version 110
```
