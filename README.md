# Benchmarking Variant Catalogue Pipelines

This repo contains code to compare pipelines designed to generate variant catalogues, a list of variants and their frequencies in a population, from whole genome sequences.

The following pipelines were compared:

1. Sarek (nextflow)
2. Variant Catalogue (nextflow)
3. CAFFEE (nextflow)

## Dataset

The test data used in this experiment was downloaded from public repositories. In order to accelerate the download rate we priorize the download of files stored in aws cloud buckets. The downloads were done using aws-cli tool.

ascat_alleles           = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh38/Annotation/ASCAT/G1000_alleles_hg38.zip"
ascat_genome            = 'hg38'
ascat_loci              = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh38/Annotation/ASCAT/G1000_loci_hg38.zip"
ascat_loci_gc           = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh38/Annotation/ASCAT/GC_G1000_hg38.zip"
ascat_loci_rt           = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh38/Annotation/ASCAT/RT_G1000_hg38.zip"
bwa                     = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh38/Sequence/BWAIndex/"
bwamem2                 = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh38/Sequence/BWAmem2Index/"
cf_chrom_len            = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh38/Sequence/Length/Homo_sapiens_assembly38.len"
chr_dir                 = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh38/Sequence/Chromosomes"
dbsnp                   = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/dbsnp_146.hg38.vcf.gz"
dbsnp_tbi               = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/dbsnp_146.hg38.vcf.gz.tbi"
dbsnp_vqsr              = '--resource:dbsnp,known=false,training=true,truth=false,prior=2.0 dbsnp_146.hg38.vcf.gz'
dict                    = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.dict"
dragmap                 = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh38/Sequence/dragmap/"
fasta                   = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta"
fasta_fai               = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFast/Homo_sapiens_assembly38.fasta.fai"
germline_resource       = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/af-only-gnomad.hg38.vcf.gz"
germline_resource_tbi   = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/af-only-gnomad.hg38.vcf.gz.tbi"
intervals               = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh38/Annotation/interval/wgs_calling_regions_noseconds.hg38.bed"
known_indels            = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/{Mills_and_1000G_gold_standard.indels.hg38,beta/Homo_sapiens_assembly38.known_indels}.vcf.gz"
known_indels_tbi        = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/{Mills_and_1000G_gold_standard.indels.hg38,beta/Homo_sapiens_assembly38.known_indels}.vcf.gz.tbi"
known_indels_vqsr       = '--resource:gatk,known=false,training=true,truth=true,prior=10.0 Homo_sapiens_assembly38.known_indels.vcf.gz --resource:mills,known=false,training=true,truth=true,prior=10.0 Mills_and_1000G_gold_standard.indels.hg38.vcf.gz'
known_snps              = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/1000G_omni2.5.hg38.vcf.gz"
known_snps_tbi          = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/1000G_omni2.5.hg38.vcf.gz.tbi"
known_snps_vqsr         = '--resource:1000G,known=false,training=true,truth=true,prior=10.0 1000G_omni2.5.hg38.vcf.gz'
mappability             = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh38/Annotation/Control-FREEC/out100m2_hg38.gem"
ngscheckmate_bed        = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh38/Annotation/NGSCheckMate/SNP_GRCh38_hg38_wChr.bed"
pon                     = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/1000g_pon.hg38.vcf.gz"
pon_tbi                 = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/1000g_pon.hg38.vcf.gz.tbi"
sentieon_dnascope_model = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh38/Annotation/Sentieon/SentieonDNAscopeModel1.1.model"
snpeff_db               = 105
snpeff_genome           = 'GRCh38'
vep_cache_version       = 110
vep_genome              = 'GRCh38'
vep_species             = 'homo_sapiens'


1. Human Reference Genome (GRCh38) from GATK Bundle.
```
$ aws s3 --no-sign-request cp s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta* .
```
2. Databases for variant annotation: VEP version 110 and SNPeff version 105. 
```
https://annotation-cache.github.io/vep_cache/
```
3. Samples to call variants: as it's a test dataset, we extracted on-the-fly the reads mapped to chromossome 12, X and Y from BAM files of GIAB NA24694, NA12878 and NA24385. After the BAM subset we convert it to paired-end fastq. This was done with:

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
$ samtools view --threads 4 -h -L include_regions_chr.bed --use-index s3://giab/data/ChineseTrio/HG006_NA24694-huCA017E_father/NA24694_Father_HiSeq100x/NHGRI_Illumina100X_Chinesetrio_novoalign_bams/HG006.GRCh38_full_plus_hs38d1_analysis_set_minus_alts.100x.bam | samtools view --threads 4 -Shb - | samtools sort --threads 4 -n -O bam -o test_sorted_NA24694.bam -

$ bedtools bamtofastq -i test_sorted_NA24694.bam -fq NA24694.R1.fastq -fq2 NA24694.R2.fastq
```
