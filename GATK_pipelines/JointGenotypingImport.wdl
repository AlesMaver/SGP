version 1.0

import "https://raw.githubusercontent.com/broadinstitute/warp/develop/tasks/broad/JointGenotypingTasks.wdl" as Tasks


# Joint Genotyping for hg38 Whole Genomes and Exomes (has not been tested on hg19)
workflow JointGenotyping {

  String pipeline_version = "1.6.1"

  input {
    File unpadded_intervals_file

    String callset_name
    File sample_name_map
    Array[File]? genomicsdb_tars

    File ref_fasta
    File ref_fasta_index
    File ref_dict

    File dbsnp_vcf
    File dbsnp_vcf_index

    Int small_disk
    Int medium_disk
    Int large_disk
    Int huge_disk

    Array[String] snp_recalibration_tranche_values
    Array[String] snp_recalibration_annotation_values
    Array[String] indel_recalibration_tranche_values
    Array[String] indel_recalibration_annotation_values

    File haplotype_database

    File eval_interval_list
    File hapmap_resource_vcf
    File hapmap_resource_vcf_index
    File omni_resource_vcf
    File omni_resource_vcf_index
    File one_thousand_genomes_resource_vcf
    File one_thousand_genomes_resource_vcf_index
    File mills_resource_vcf
    File mills_resource_vcf_index
    File axiomPoly_resource_vcf
    File axiomPoly_resource_vcf_index
    File dbsnp_resource_vcf = dbsnp_vcf
    File dbsnp_resource_vcf_index = dbsnp_vcf_index

    # ExcessHet is a phred-scaled p-value. We want a cutoff of anything more extreme
    # than a z-score of -4.5 which is a p-value of 3.4e-06, which phred-scaled is 54.69
    Float excess_het_threshold = 54.69
    Float snp_filter_level
    Float indel_filter_level
    Int SNP_VQSR_downsampleFactor

    Int? top_level_scatter_count
    Boolean? gather_vcfs
    Int snps_variant_recalibration_threshold = 500000
    Boolean rename_gvcf_samples = true
    Float unbounded_scatter_count_scale_factor = 0.15
    Int gnarly_scatter_count = 10
    Boolean use_gnarly_genotyper = false
    Boolean use_allele_specific_annotations = true
    Boolean cross_check_fingerprints = true
    Boolean scatter_cross_check_fingerprints = false
  }

  Boolean allele_specific_annotations = !use_gnarly_genotyper && use_allele_specific_annotations

  Array[Array[String]] sample_name_map_lines = read_tsv(sample_name_map)
  Int num_gvcfs = length(sample_name_map_lines)

  # Make a 2.5:1 interval number to samples in callset ratio interval list.
  # We allow overriding the behavior by specifying the desired number of vcfs
  # to scatter over for testing / special requests.
  # Zamboni notes say "WGS runs get 30x more scattering than Exome" and
  # exome scatterCountPerSample is 0.05, min scatter 10, max 1000

  # For small callsets (fewer than 1000 samples) we can gather the VCF shards and collect metrics directly.
  # For anything larger, we need to keep the VCF sharded and gather metrics collected from them.
  # We allow overriding this default behavior for testing / special requests.
  Boolean is_small_callset = select_first([gather_vcfs, num_gvcfs <= 1000])

  Int unbounded_scatter_count = select_first([top_level_scatter_count, round(unbounded_scatter_count_scale_factor * num_gvcfs)])
  Int scatter_count = if unbounded_scatter_count > 2 then unbounded_scatter_count else 2 #I think weird things happen if scatterCount is 1 -- IntervalListTools is noop?

  call Tasks.CheckSamplesUnique {
    input:
      sample_name_map = sample_name_map
  }

  call Tasks.SplitIntervalList {
    input:
      interval_list = unpadded_intervals_file,
      scatter_count = scatter_count,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      disk_size = small_disk,
      sample_names_unique_done = CheckSamplesUnique.samples_unique
  }

  Array[File] unpadded_intervals = SplitIntervalList.output_intervals

  scatter (idx in range(length(unpadded_intervals))) {
    # The batch_size value was carefully chosen here as it
    # is the optimal value for the amount of memory allocated
    # within the task; please do not change it without consulting
    # the Hellbender (GATK engine) team!
    if ( !defined(genomicsdb_tars) ) {
        call ImportGVCFs_import as ImportGVCFs {
          input:
            sample_name_map = sample_name_map,
            interval = unpadded_intervals[idx],
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,
            workspace_dir_name = "genomicsdb" + "_" + idx,
            index = idx,
            disk_size = medium_disk,
            batch_size = 50
        }
    }

    if ( defined(genomicsdb_tars) ) {
        Array[File] genomicsdb_tars_provided = select_first([genomicsdb_tars])

        call ImportGVCFs_update {
          input:
            sample_name_map = sample_name_map,
            interval = unpadded_intervals[idx],
            workspace_tar = genomicsdb_tars_provided[idx],
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,
            workspace_dir_name = "genomicsdb" + "_" + idx,
            index = idx,
            disk_size = medium_disk,
            batch_size = 50
        }
    }

    File genomicsdb_tar = select_first([ImportGVCFs.output_genomicsdb, ImportGVCFs_update.output_genomicsdb])

  }

  output {
    # Output the interval list generated/used by this run workflow.
    Array[File] output_workspaces = genomicsdb_tar
  }
  meta {
    allowNestedInputs: true
  }
}


# Override
task ImportGVCFs_import {

  input {
    File sample_name_map
    File interval
    File ref_fasta
    File ref_fasta_index
    File ref_dict

    String workspace_dir_name
    String? index

    Int disk_size
    Int batch_size

    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.2.6.1"
  }

  String idx = select_first([index, ""])

  command <<<
    set -euo pipefail

    rm -rf ~{workspace_dir_name}

    # We've seen some GenomicsDB performance regressions related to intervals, so we're going to pretend we only have a single interval
    # using the --merge-input-intervals arg
    # There's no data in between since we didn't run HaplotypeCaller over those loci so we're not wasting any compute

    # The memory setting here is very important and must be several GiB lower
    # than the total memory allocated to the VM because this tool uses
    # a significant amount of non-heap memory for native libraries.
    # Also, testing has shown that the multithreaded reader initialization
    # does not scale well beyond 5 threads, so don't increase beyond that.
    gatk --java-options "-Xms8000m -Xmx25000m" \
      GenomicsDBImport \
      --genomicsdb-workspace-path ~{workspace_dir_name} \
      --batch-size ~{batch_size} \
      -L ~{interval} \
      --sample-name-map ~{sample_name_map} \
      --reader-threads 5 \
      --merge-input-intervals \
      --consolidate

    tar -cf ~{workspace_dir_name}.tar ~{workspace_dir_name}
  >>>

  runtime {
    memory: "26000 MiB"
    cpu: 4
    bootDiskSizeGb: 15
    disks: "local-disk " + disk_size + " HDD"
    docker: gatk_docker
    preemptible: 1
  }

  output {
    File output_genomicsdb = "~{workspace_dir_name}.tar"
  }
}

# Override
task ImportGVCFs_update {

  input {
    File sample_name_map
    File interval
    File ref_fasta
    File ref_fasta_index
    File ref_dict

    File workspace_tar

    String workspace_dir_name
    String? index

    Int disk_size
    Int batch_size

    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.2.6.1"
  }

  String idx = select_first([index, ""])

  command <<<
    set -euo pipefail

    tar -xf ~{workspace_tar}
    WORKSPACE=$(basename ~{workspace_tar} .tar)

    # We've seen some GenomicsDB performance regressions related to intervals, so we're going to pretend we only have a single interval
    # using the --merge-input-intervals arg
    # There's no data in between since we didn't run HaplotypeCaller over those loci so we're not wasting any compute

    # The memory setting here is very important and must be several GiB lower
    # than the total memory allocated to the VM because this tool uses
    # a significant amount of non-heap memory for native libraries.
    # Also, testing has shown that the multithreaded reader initialization
    # does not scale well beyond 5 threads, so don't increase beyond that.
    gatk --java-options "-Xms8000m -Xmx25000m" \
      GenomicsDBImport \
      --genomicsdb-update-workspace-path $WORKSPACE \
      --batch-size ~{batch_size} \
      --sample-name-map ~{sample_name_map} \
      --reader-threads 5 \
      --merge-input-intervals \
      --consolidate

    tar -cf ~{workspace_dir_name}.tar ~{workspace_dir_name}
  >>>

  runtime {
    memory: "26000 MiB"
    cpu: 4
    bootDiskSizeGb: 15
    disks: "local-disk " + disk_size + " HDD"
    docker: gatk_docker
    preemptible: 1
  }

  output {
    File output_genomicsdb = "~{workspace_dir_name}.tar"
  }
}
