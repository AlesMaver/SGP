# Oblika zbranih podatkov in specifikacije analize - Predlog (SI)
## Vhodni podatki
Vhodni podatek | Format | Imenovanje | Opis
--- | --- | --- | ---
Nenalegana branja | uBAM | SGP00001.WGS.ubam | Zaporedja v binarni obliki (BAM) brez podatka o položaju na genomu (nenalegana branja). Format uBAM ima pred FASTQ prednost, da vsebuje več informacij o izvoru in skupinah branj. 
* Vhodne podatke je mogoče po opravljenih primarnih analizah zbrisati, saj datoteka CRAM vsebuje izvornih informacij in je mogoče regenerirati BAM in FASTQ datoteke. 

## Rezultati analiz
Rezultat | Format | Imenovanje | Opis
--- | --- | --- | ---
Nalegana branja | CRAM | SGP00001.WGS.hg38.cram (.crai) | Binarni format naleganih branj, ki uporablja vhodno genomsko referenco za opis razlik med naleganim branjem in referenčnim zapisom genoma.
Genomske različice | gVCF (gzip) | SGP00001.WGS.hg38.gvcf.gz | Oblika hranjenja podatkov o različicah, ki vsebuje podatke o zanesljivosti detekcije različic tako na vseh mestih (tako na mestih označenih različic kakor tudi na vseh ostalih mestih genoma). Ta format podatkov omogoča tudi detekcijo genetskih različic v sledečih kohortnih analizah. 

## Uporabljeni delovni tokovi
Delovni tok | Povezava | Opis
--- | --- | ---
WholeGenomeGermlineSingleSample_v2.0 | [link](https://github.com/broadinstitute/warp/tree/develop/pipelines/broad/dna_seq/germline/single_sample/wgs) | Iz opisa Broadovega repozitorija: The Whole Genome Germline Single Sample pipeline implements data pre-processing and initial variant calling according to the GATK Best Practices (June 2016) for germline SNP and Indel discovery in human whole-genome sequencing data.
* Repozitorij je zrcaljen v repozitoriju SGP in bo spreminjan po potrebi projekta (vendar z ohranitvijo ključnih korakov priporočil dobre prakse)

## Jezik delovnih tokov
Format | Opis | Specifikacije | Program za izvedbo
--- | --- | --- | ---
WDL | Workflow definition language | [WDL 1.0](https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md) | [Cromwell](https://github.com/broadinstitute/cromwell/)

# Referenčni podatki - Predlog
Za namene Slovenskega genomskega projekta bomo uporabljali genomski sestav v različici hg38

Pri analizah bomo uporabljali naslednje referenčne podatke, ki so na voljo na naslovu : gs://gcp-public-data--broad-references/hg38/v0/

Opis vira | Datoteka | V uporabi
--- | --- | ---
GATK resource bundle hg38 | contamination-resources/1000g/1000g.phase3.100k.b38.vcf.gz.dat.UD | DA
GATK resource bundle hg38 | contamination-resources/1000g/1000g.phase3.100k.b38.vcf.gz.dat.bed | DA
GATK resource bundle hg38 | contamination-resources/1000g/1000g.phase3.100k.b38.vcf.gz.dat.mu | DA
GATK resource bundle hg38 | contamination-resources/ | DA
GATK resource bundle hg38 | wgs_calling_regions.hg38.interval_list | DA
GATK resource bundle hg38 | Homo_sapiens_assembly38.dict | DA
GATK resource bundle hg38 | Homo_sapiens_assembly38.fasta | DA
GATK resource bundle hg38 | Homo_sapiens_assembly38.fasta.fai | DA
GATK resource bundle hg38 | Homo_sapiens_assembly38.fasta.64.alt | DA
GATK resource bundle hg38 | Homo_sapiens_assembly38.fasta.64.sa | DA
GATK resource bundle hg38 | Homo_sapiens_assembly38.fasta.64.amb | DA
GATK resource bundle hg38 | Homo_sapiens_assembly38.fasta.64.bwt | DA
GATK resource bundle hg38 | Homo_sapiens_assembly38.fasta.64.ann | DA
GATK resource bundle hg38 | Homo_sapiens_assembly38.fasta.64.pac | DA
GATK resource bundle hg38 | Mills_and_1000G_gold_standard.indels.hg38.vcf.gz | DA
GATK resource bundle hg38 | Homo_sapiens_assembly38.known_indels.vcf.gz | DA
GATK resource bundle hg38 | Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi | DA
GATK resource bundle hg38 | Homo_sapiens_assembly38.known_indels.vcf.gz.tbi | DA
GATK resource bundle hg38 | Homo_sapiens_assembly38.dbsnp138.vcf | DA
GATK resource bundle hg38 | Homo_sapiens_assembly38.dbsnp138.vcf.idx | DA
GATK resource bundle hg38 | wgs_evaluation_regions.hg38.interval_list | DA
GATK resource bundle hg38 | Homo_sapiens_assembly38.haplotype_database.txt | DA
GATK resource bundle hg38 | wgs_coverage_regions.hg38.interval_list | DA
GATK resource bundle hg38 | Homo_sapiens_assembly38.haplotype_database.txt ./wgs_test_reference/ | DA
GATK resource bundle hg38 | wgs_coverage_regions.hg38.interval_list ./wgs_test_reference/ | DA

#
#
# Format and analysis specifications - PROPOSAL (EN)
## Input files
Input | Format | Naming | Description
--- | --- | --- | ---
Unmapped reads | uBAM | SGP00001.WGS.ubam | Binary unaligned read files in a single BAM file, that keeps the read group and other information. Preferrable to FASTQ for data storage. 
* The input files may be deleted after the complete project analysis is complete as the relevant information is all contained in the output files

## Output files
Output | Format | Naming | Description
--- | --- | --- | ---
Aligned reads | CRAM | SGP00001.WGS.hg38.cram (.crai) | Binary alignment format that uses a genomic reference to describe differences between the aligned sequence fragments and the reference sequence
Variant calls | gVCF (gzipped) | SGP00001.WGS.hg38.gvcf.gz | A variant storage format that contains the information on variant calls and variation likelihoods across the whole genome (including likelihoods for vairant and invariant sites). This format also serves as the input for variant calling in cohort mode. 

## Production pipelines
Pipeline | Location | Description
--- | --- | ---
WholeGenomeGermlineSingleSample_v2.0 | [link](https://github.com/broadinstitute/warp/tree/develop/pipelines/broad/dna_seq/germline/single_sample/wgs) | From the Broad's repo: The Whole Genome Germline Single Sample pipeline implements data pre-processing and initial variant calling according to the GATK Best Practices (June 2016) for germline SNP and Indel discovery in human whole-genome sequencing data.
* The pipeline is mirrored in our GitHub repo and will be modified in accordance to the needs of our project (while maintaining compliance with GATK best practice guidelines)

## Pipeline language definition
Format | Description | Specs | Execution engine
--- | --- | --- | ---
WDL | Workflow definition language | [WDL 1.0](https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md) | [Cromwell](https://github.com/broadinstitute/cromwell/)
