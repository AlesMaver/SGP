# Slovenian Genome Project pipelines
This folder will contain the proposed pipelines for primary data processing of the whole genome sequencing data performed within the Slovenian genome project. 

The current contents are as follows:
* The Broad's warp pipelines which are proposed to serve as the best practice pipeline for processing of sequence reads and generating the CRAM and GVCF files
  * MODIFICATION: Corrected missing input_bam_index entries in several calls, which was causing the pipeline to fail
  * PROPOSED MODIFICATION: Allow input of two FASTQ files and their conversion to uBAM as inputs
  * PROPOSED MODIFICATION: Always output CRAM files rather than BAM. The genome reference should be included in the CRAM file name, ie. with "_hg38_" label
