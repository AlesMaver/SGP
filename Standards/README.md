# Data formats specification - PROPOSAL
## Input files
Input | Format | Description
--- | --- | ---
Unmapped reads | uBAM | Binary unaligned read files in a single BAM file, that keeps the read group and other information. Preferrable to FASTQ for data storage. 
* The input files may be deleted after the complete project analysis is complete as the relevant information is all contained in the output files

## Output files
Output | Format | Description
--- | --- | ---
Aligned reads | CRAM | Binary alignment format that uses a genomic reference to describe differences between the aligned sequence fragments and the reference sequence
Variant calls | gVCF (gzipped) | A variant storage format that contains the information on variant calls and variation likelihoods across the whole genome (including likelihoods for vairant and invariant sites). This format also serves as the input for variant calling in cohort mode. 

## Production pipelines
Pipeline | Location | Description
--- | --- | ---
WholeGenomeGermlineSingleSample_v2.0 | [link](https://github.com/broadinstitute/warp/tree/develop/pipelines/broad/dna_seq/germline/single_sample/wgs) | From the Broad's repo: The Whole Genome Germline Single Sample pipeline implements data pre-processing and initial variant calling according to the GATK Best Practices (June 2016) for germline SNP and Indel discovery in human whole-genome sequencing data.
* The pipeline is mirrored in our GitHub repo and will be modified in accordance to the needs of our project (while maintaining compliance with GATK best practice guidelines)

## Pipeline language definition
Format | Description | Specs | Execution engine
--- | --- | --- | ---
WDL | Workflow definition language | [WDL 1.0](https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md) | [Cromwell] (https://github.com/broadinstitute/cromwell/)
