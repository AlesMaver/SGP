# Slovenian Genome Project pipelines
This folder contains the proposed pipelines for primary data processing of the whole genome sequencing data performed within the Slovenian genome project. 

This repository is written predominantly in the workflow definition language (WDL format), which contains both - the code and the software execution environment definitions (dockers) for data analysis. The workflows implemented are compatible with local execution or clusters. The task requirements are in process of optimization for SLURM execution environment, but the Google Cloud compute definitions are maintained as in the original Broad's pipelines. 

[Cromwell engine](https://github.com/broadinstitute/cromwell) is advised to run the WDL worklows. 

The current contents are as follows:
* The Broad's warp pipelines which are proposed to serve as the best practice pipeline for processing of sequence reads and generating the CRAM and GVCF files
  * MODIFICATION: Corrected missing input_bam_index entries in several calls, which was causing the pipeline to fail
  * PROPOSED MODIFICATION: Allow input of two FASTQ files and their conversion to uBAM as inputs
  * PROPOSED MODIFICATION: Always output CRAM files rather than BAM. The genome reference should be included in the CRAM file name, ie. with "_hg38_" label

