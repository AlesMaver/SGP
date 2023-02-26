version 1.0

import "warp/tasks/broad/JointGenotypingTasks.wdl" as Tasks


# Joint Genotyping for hg38 Whole Genomes and Exomes (has not been tested on hg19)
workflow JointGenotyping {

  String pipeline_version = "1.6.1"

  input {
    File unpadded_intervals_file

    String callset_name
    File sample_name_map
    Array[File] genomicsdb_tars

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

    if ( defined(genomicsdb_tars) ) {
        Array[File] genomicsdb_tars_provided = select_first([genomicsdb_tars])

        call ImportGDB {
          input:
            interval = unpadded_intervals[idx],
            workspace_tar = genomicsdb_tars_provided[idx],
            workspace_dir_name = "genomicsdb" + "_" + idx,
            index = idx,
            disk_size = medium_disk,
            batch_size = 50
        }
    }

    File genomicsdb_tar = select_first([ImportGDB.output_genomicsdb])

    if (use_gnarly_genotyper) {

      call Tasks.SplitIntervalList as GnarlyIntervalScatterDude {
        input:
          interval_list = unpadded_intervals[idx],
          scatter_count = gnarly_scatter_count,
          ref_fasta = ref_fasta,
          ref_fasta_index = ref_fasta_index,
          ref_dict = ref_dict,
          disk_size = small_disk,
          sample_names_unique_done = CheckSamplesUnique.samples_unique
      }

      Array[File] gnarly_intervals = GnarlyIntervalScatterDude.output_intervals

      scatter (gnarly_idx in range(length(gnarly_intervals))) {
        call Tasks.GnarlyGenotyper {
          input:
            workspace_tar = genomicsdb_tar,
            interval = gnarly_intervals[gnarly_idx],
            output_vcf_filename = callset_name + "." + idx + "." + gnarly_idx + ".vcf.gz",
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,
            dbsnp_vcf = dbsnp_vcf,
        }
      }

      Array[File] gnarly_gvcfs = GnarlyGenotyper.output_vcf

      call Tasks.GatherVcfs as TotallyRadicalGatherVcfs {
        input:
          input_vcfs = gnarly_gvcfs,
          output_vcf_name = callset_name + "." + idx + ".gnarly.vcf.gz",
          disk_size = large_disk
      }
    }

    if (!use_gnarly_genotyper) {
      call GenotypeGVCFs {
        input:
          workspace_tar = genomicsdb_tar,
          interval = unpadded_intervals[idx],
          output_vcf_filename = callset_name + "." + idx + ".vcf.gz",
          ref_fasta = ref_fasta,
          ref_fasta_index = ref_fasta_index,
          ref_dict = ref_dict,
          dbsnp_vcf = dbsnp_vcf,
          disk_size = medium_disk
      }
    }

    File genotyped_vcf = select_first([TotallyRadicalGatherVcfs.output_vcf, GenotypeGVCFs.output_vcf])
    File genotyped_vcf_index = select_first([TotallyRadicalGatherVcfs.output_vcf_index, GenotypeGVCFs.output_vcf_index])

    call Tasks.HardFilterAndMakeSitesOnlyVcf {
      input:
        vcf = genotyped_vcf,
        vcf_index = genotyped_vcf_index,
        excess_het_threshold = excess_het_threshold,
        ## variant_filtered_vcf_filename = callset_name + "." + idx + ".variant_filtered.vcf.gz",
        variant_filtered_vcf_filename = callset_name + "." + idx + ".vf.vcf.gz",
        ## sites_only_vcf_filename = callset_name + "." + idx + ".sites_only.variant_filtered.vcf.gz",
        sites_only_vcf_filename = callset_name + "." + idx + ".so.vcf.gz",
        disk_size = medium_disk
    }
  }

  call Tasks.GatherVcfs as SitesOnlyGatherVcf {
    input:
      input_vcfs = HardFilterAndMakeSitesOnlyVcf.sites_only_vcf,
      output_vcf_name = callset_name + ".sites_only.vcf.gz",
      disk_size = medium_disk
  }

  call Tasks.IndelsVariantRecalibrator {
    input:
      sites_only_variant_filtered_vcf = SitesOnlyGatherVcf.output_vcf,
      sites_only_variant_filtered_vcf_index = SitesOnlyGatherVcf.output_vcf_index,
      recalibration_filename = callset_name + ".indels.recal",
      tranches_filename = callset_name + ".indels.tranches",
      recalibration_tranche_values = indel_recalibration_tranche_values,
      recalibration_annotation_values = indel_recalibration_annotation_values,
      mills_resource_vcf = mills_resource_vcf,
      mills_resource_vcf_index = mills_resource_vcf_index,
      axiomPoly_resource_vcf = axiomPoly_resource_vcf,
      axiomPoly_resource_vcf_index = axiomPoly_resource_vcf_index,
      dbsnp_resource_vcf = dbsnp_resource_vcf,
      dbsnp_resource_vcf_index = dbsnp_resource_vcf_index,
      use_allele_specific_annotations = allele_specific_annotations,
      disk_size = small_disk
  }

  if (num_gvcfs > snps_variant_recalibration_threshold) {
    call Tasks.SNPsVariantRecalibratorCreateModel {
      input:
        sites_only_variant_filtered_vcf = SitesOnlyGatherVcf.output_vcf,
        sites_only_variant_filtered_vcf_index = SitesOnlyGatherVcf.output_vcf_index,
        recalibration_filename = callset_name + ".snps.recal",
        tranches_filename = callset_name + ".snps.tranches",
        recalibration_tranche_values = snp_recalibration_tranche_values,
        recalibration_annotation_values = snp_recalibration_annotation_values,
        downsampleFactor = SNP_VQSR_downsampleFactor,
        model_report_filename = callset_name + ".snps.model.report",
        hapmap_resource_vcf = hapmap_resource_vcf,
        hapmap_resource_vcf_index = hapmap_resource_vcf_index,
        omni_resource_vcf = omni_resource_vcf,
        omni_resource_vcf_index = omni_resource_vcf_index,
        one_thousand_genomes_resource_vcf = one_thousand_genomes_resource_vcf,
        one_thousand_genomes_resource_vcf_index = one_thousand_genomes_resource_vcf_index,
        dbsnp_resource_vcf = dbsnp_resource_vcf,
        dbsnp_resource_vcf_index = dbsnp_resource_vcf_index,
        use_allele_specific_annotations = allele_specific_annotations,
        disk_size = small_disk
    }

    scatter (idx in range(length(HardFilterAndMakeSitesOnlyVcf.sites_only_vcf))) {
      call SNPsVariantRecalibrator as SNPsVariantRecalibratorScattered {
        input:
          sites_only_variant_filtered_vcf = HardFilterAndMakeSitesOnlyVcf.sites_only_vcf[idx],
          sites_only_variant_filtered_vcf_index = HardFilterAndMakeSitesOnlyVcf.sites_only_vcf_index[idx],
          recalibration_filename = callset_name + ".snps." + idx + ".recal",
          tranches_filename = callset_name + ".snps." + idx + ".tranches",
          recalibration_tranche_values = snp_recalibration_tranche_values,
          recalibration_annotation_values = snp_recalibration_annotation_values,
          model_report = SNPsVariantRecalibratorCreateModel.model_report,
          hapmap_resource_vcf = hapmap_resource_vcf,
          hapmap_resource_vcf_index = hapmap_resource_vcf_index,
          omni_resource_vcf = omni_resource_vcf,
          omni_resource_vcf_index = omni_resource_vcf_index,
          one_thousand_genomes_resource_vcf = one_thousand_genomes_resource_vcf,
          one_thousand_genomes_resource_vcf_index = one_thousand_genomes_resource_vcf_index,
          dbsnp_resource_vcf = dbsnp_resource_vcf,
          dbsnp_resource_vcf_index = dbsnp_resource_vcf_index,
          use_allele_specific_annotations = allele_specific_annotations,
          disk_size = small_disk
        }
    }

    call Tasks.GatherTranches as SNPGatherTranches {
      input:
        tranches = SNPsVariantRecalibratorScattered.tranches,
        output_filename = callset_name + ".snps.gathered.tranches",
        mode = "SNP",
        disk_size = small_disk
    }
  }

  if (num_gvcfs <= snps_variant_recalibration_threshold) {
    call SNPsVariantRecalibrator as SNPsVariantRecalibratorClassic {
      input:
        sites_only_variant_filtered_vcf = SitesOnlyGatherVcf.output_vcf,
        sites_only_variant_filtered_vcf_index = SitesOnlyGatherVcf.output_vcf_index,
        recalibration_filename = callset_name + ".snps.recal",
        tranches_filename = callset_name + ".snps.tranches",
        recalibration_tranche_values = snp_recalibration_tranche_values,
        recalibration_annotation_values = snp_recalibration_annotation_values,
        hapmap_resource_vcf = hapmap_resource_vcf,
        hapmap_resource_vcf_index = hapmap_resource_vcf_index,
        omni_resource_vcf = omni_resource_vcf,
        omni_resource_vcf_index = omni_resource_vcf_index,
        one_thousand_genomes_resource_vcf = one_thousand_genomes_resource_vcf,
        one_thousand_genomes_resource_vcf_index = one_thousand_genomes_resource_vcf_index,
        dbsnp_resource_vcf = dbsnp_resource_vcf,
        dbsnp_resource_vcf_index = dbsnp_resource_vcf_index,
        use_allele_specific_annotations = allele_specific_annotations,
        disk_size = small_disk
    }
  }

  scatter (idx in range(length(HardFilterAndMakeSitesOnlyVcf.variant_filtered_vcf))) {
    #for really large callsets we give to friends, just apply filters to the sites-only
    call Tasks.ApplyRecalibration {
      input:
        recalibrated_vcf_filename = callset_name + ".filtered." + idx + ".vcf.gz",
        input_vcf = HardFilterAndMakeSitesOnlyVcf.variant_filtered_vcf[idx],
        input_vcf_index = HardFilterAndMakeSitesOnlyVcf.variant_filtered_vcf_index[idx],
        indels_recalibration = IndelsVariantRecalibrator.recalibration,
        indels_recalibration_index = IndelsVariantRecalibrator.recalibration_index,
        indels_tranches = IndelsVariantRecalibrator.tranches,
        snps_recalibration = if defined(SNPsVariantRecalibratorScattered.recalibration) then select_first([SNPsVariantRecalibratorScattered.recalibration])[idx] else select_first([SNPsVariantRecalibratorClassic.recalibration]),
        snps_recalibration_index = if defined(SNPsVariantRecalibratorScattered.recalibration_index) then select_first([SNPsVariantRecalibratorScattered.recalibration_index])[idx] else select_first([SNPsVariantRecalibratorClassic.recalibration_index]),
        snps_tranches = select_first([SNPGatherTranches.tranches_file, SNPsVariantRecalibratorClassic.tranches]),
        indel_filter_level = indel_filter_level,
        snp_filter_level = snp_filter_level,
        use_allele_specific_annotations = allele_specific_annotations,
        disk_size = medium_disk
    }

    # For large callsets we need to collect metrics from the shards and gather them later.
    if (!is_small_callset) {
      call Tasks.CollectVariantCallingMetrics as CollectMetricsSharded {
        input:
          input_vcf = ApplyRecalibration.recalibrated_vcf,
          input_vcf_index = ApplyRecalibration.recalibrated_vcf_index,
          metrics_filename_prefix = callset_name + "." + idx,
          dbsnp_vcf = dbsnp_vcf,
          dbsnp_vcf_index = dbsnp_vcf_index,
          interval_list = eval_interval_list,
          ref_dict = ref_dict,
          disk_size = medium_disk
      }
    }
  }

  # For small and large callsets we gather the VCF shards.
  call Tasks.GatherVcfs as FinalGatherVcf {
    input:
      input_vcfs = ApplyRecalibration.recalibrated_vcf,
      output_vcf_name = callset_name + ".vcf.gz",
      disk_size = huge_disk
  }

  # For small callsets we collect metrics on a gathered VCF
  if (is_small_callset) {
    call Tasks.CollectVariantCallingMetrics as CollectMetricsOnFullVcf {
      input:
        input_vcf = FinalGatherVcf.output_vcf,
        input_vcf_index = FinalGatherVcf.output_vcf_index,
        metrics_filename_prefix = callset_name,
        dbsnp_vcf = dbsnp_vcf,
        dbsnp_vcf_index = dbsnp_vcf_index,
        interval_list = eval_interval_list,
        ref_dict = ref_dict,
        disk_size = large_disk
    }
  }

  # For large callsets we still need to gather the sharded metrics.
  if (!is_small_callset) {
    call Tasks.GatherVariantCallingMetrics {
      input:
        input_details = select_all(CollectMetricsSharded.detail_metrics_file),
        input_summaries = select_all(CollectMetricsSharded.summary_metrics_file),
        output_prefix = callset_name,
        disk_size = medium_disk
    }
  }

  # WE DO NOT NEED FINGERPRINT CHECKS - comment out
  # # CrossCheckFingerprints takes forever on large callsets.
  # # We scatter over the input GVCFs to make things faster.
  # if (scatter_cross_check_fingerprints) {
  #   call Tasks.GetFingerprintingIntervalIndices {
  #     input:
  #       unpadded_intervals = unpadded_intervals,
  #       haplotype_database = haplotype_database
  #   }

  #   Array[Int] fingerprinting_indices = GetFingerprintingIntervalIndices.indices_to_fingerprint

  #   scatter (idx in fingerprinting_indices) {
  #     File vcfs_to_fingerprint = HardFilterAndMakeSitesOnlyVcf.variant_filtered_vcf[idx]
  #   }

  #   call Tasks.GatherVcfs as GatherFingerprintingVcfs {
  #     input:
  #       input_vcfs = vcfs_to_fingerprint,
  #       output_vcf_name = callset_name + ".gathered.fingerprinting.vcf.gz",
  #       disk_size = medium_disk
  #   }

  #   call Tasks.SelectFingerprintSiteVariants {
  #     input:
  #       input_vcf = GatherFingerprintingVcfs.output_vcf,
  #       base_output_name = callset_name + ".fingerprinting",
  #       haplotype_database = haplotype_database,
  #       disk_size = medium_disk
  #   }

  #   call Tasks.PartitionSampleNameMap {
  #     input:
  #       sample_name_map = sample_name_map,
  #       line_limit = 1000
  #   }

  #   scatter (idx in range(length(PartitionSampleNameMap.partitions))) {

  #     Array[File] files_in_partition = read_lines(PartitionSampleNameMap.partitions[idx])

  #     call Tasks.CrossCheckFingerprint as CrossCheckFingerprintsScattered {
  #       input:
  #         gvcf_paths = files_in_partition,
  #         vcf_paths = vcfs_to_fingerprint,
  #         sample_name_map = sample_name_map,
  #         haplotype_database = haplotype_database,
  #         output_base_name = callset_name + "." + idx,
  #         scattered = true
  #     }
  #   }

  #   call Tasks.GatherPicardMetrics as GatherFingerprintingMetrics {
  #     input:
  #       metrics_files = CrossCheckFingerprintsScattered.crosscheck_metrics,
  #       output_file_name = callset_name + ".fingerprintcheck",
  #       disk_size = small_disk
  #   }
  # }

  # if (!scatter_cross_check_fingerprints) {

  #   scatter (line in sample_name_map_lines) {
  #     File gvcf_paths = line[1]
  #   }

  #   call Tasks.CrossCheckFingerprint as CrossCheckFingerprintSolo {
  #     input:
  #       gvcf_paths = gvcf_paths,
  #       vcf_paths = ApplyRecalibration.recalibrated_vcf,
  #       sample_name_map = sample_name_map,
  #       haplotype_database = haplotype_database,
  #       output_base_name = callset_name
  #   }
  # }

  # Get the metrics from either code path
  File output_detail_metrics_file = select_first([CollectMetricsOnFullVcf.detail_metrics_file, GatherVariantCallingMetrics.detail_metrics_file])
  File output_summary_metrics_file = select_first([CollectMetricsOnFullVcf.summary_metrics_file, GatherVariantCallingMetrics.summary_metrics_file])

  # Get the VCFs from either code path
  Array[File?] output_vcf_files = if defined(FinalGatherVcf.output_vcf) then [FinalGatherVcf.output_vcf] else ApplyRecalibration.recalibrated_vcf
  Array[File?] output_vcf_index_files = if defined(FinalGatherVcf.output_vcf_index) then [FinalGatherVcf.output_vcf_index] else ApplyRecalibration.recalibrated_vcf_index

  output {
    # Metrics from either the small or large callset
    File detail_metrics_file = output_detail_metrics_file
    File summary_metrics_file = output_summary_metrics_file

    # Outputs from the small callset path through the wdl.
    Array[File] output_vcfs = select_all(output_vcf_files)
    Array[File] output_vcf_indices = select_all(output_vcf_index_files)

    # Output the interval list generated/used by this run workflow.
    Array[File] output_intervals = SplitIntervalList.output_intervals

    # Disable fingerprint check outputs
    # # Output the metrics from crosschecking fingerprints.
    # File crosscheck_fingerprint_check = select_first([CrossCheckFingerprintSolo.crosscheck_metrics, GatherFingerprintingMetrics.gathered_metrics])
  }
  meta {
    allowNestedInputs: true
  }
}


task ImportGDB {

  input {
    File interval
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

    # tar -xf ~{workspace_tar}
    # tar -cf ~{workspace_dir_name}.tar ~{workspace_dir_name}
    cp -lf ~{workspace_tar} ~{workspace_dir_name}.tar
  >>>

  runtime {
    memory: "3750 MiB"
    cpu: 1
    bootDiskSizeGb: 15
    disks: "local-disk " + disk_size + " HDD"
    docker: gatk_docker
    preemptible: 1
  }

  output {
    File output_genomicsdb = "~{workspace_dir_name}.tar"
  }
}


# Override: cpu, before 2, used 4
#           mem 26->32 GB, java mem 31 GB, default 25 GB
# Added: --max-alternate-alleles 4, default 6
#        --genomicsdb-max-alternate-alleles 7, default max-alternate-alleles + 3
#        --max-genotype-count 256 (== 2^8); default 2^10=1024
# Note: --disable-bam-index-caching: Caching is automatically disabled if there are no intervals specified.
task GenotypeGVCFs {

  input {
    File workspace_tar
    File interval

    String output_vcf_filename

    File ref_fasta
    File ref_fasta_index
    File ref_dict

    String dbsnp_vcf

    Int disk_size
    # This is needed for gVCFs generated with GATK3 HaplotypeCaller
    Boolean allow_old_rms_mapping_quality_annotation_data = true
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.2.6.1"
  }

  parameter_meta {
    interval: {
      localization_optional: true
    }
  }

  command <<<
    set -euo pipefail

    tar -xf ~{workspace_tar}
    WORKSPACE=$(basename ~{workspace_tar} .tar)

    gatk --java-options "-Xms8000m -Xmx31000m" \
      GenotypeGVCFs \
      -R ~{ref_fasta} \
      -O ~{output_vcf_filename} \
      -D ~{dbsnp_vcf} \
      -G StandardAnnotation -G AS_StandardAnnotation \
      --only-output-calls-starting-in-intervals \
      -V gendb://$WORKSPACE \
      -L ~{interval} \
      ~{true='--allow-old-rms-mapping-quality-annotation-data' false='' allow_old_rms_mapping_quality_annotation_data} \
      --merge-input-intervals \
      --max-alternate-alleles 4 \
      --genomicsdb-max-alternate-alleles 7 \
      --max-genotype-count 256
  >>>

  runtime {
    memory: "32000 MiB" # Increased from "26000 MiB" to fit default memory 8GB/core @ Vega; additionally increased -Xmx21000m -> -Xmx31000m
    cpu: 4 # Increased 2->4 in order to facilitate obtain enough memory
    bootDiskSizeGb: 15
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 1
    docker: gatk_docker
  }

  output {
    File output_vcf = "~{output_vcf_filename}"
    File output_vcf_index = "~{output_vcf_filename}.tbi"
  }
}

# Override: cpu: 2 => 16 in order to obtain enough memory (Vega: 4GB/thread, 8GB/core)
task SNPsVariantRecalibrator {

  input {
    String recalibration_filename
    String tranches_filename
    File? model_report

    Array[String] recalibration_tranche_values
    Array[String] recalibration_annotation_values

    File sites_only_variant_filtered_vcf
    File sites_only_variant_filtered_vcf_index

    File hapmap_resource_vcf
    File omni_resource_vcf
    File one_thousand_genomes_resource_vcf
    File dbsnp_resource_vcf
    File hapmap_resource_vcf_index
    File omni_resource_vcf_index
    File one_thousand_genomes_resource_vcf_index
    File dbsnp_resource_vcf_index
    Boolean use_allele_specific_annotations
    Int max_gaussians = 6

    Int disk_size
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.2.6.1"
    Int? machine_mem_mb

  }

  Int auto_mem = ceil(3 * size([sites_only_variant_filtered_vcf,
                              hapmap_resource_vcf,
                              omni_resource_vcf,
                              one_thousand_genomes_resource_vcf,
                              dbsnp_resource_vcf],
                      "MiB"))
  Int machine_mem = select_first([machine_mem_mb, if auto_mem < 7000 then 7000 else auto_mem])
  Int java_mem = machine_mem - 1000
  Int max_heap = machine_mem - 500


  String model_report_arg = if defined(model_report) then "--input-model $MODEL_REPORT --output-tranches-for-scatter" else ""

  command <<<
    set -euo pipefail

    MODEL_REPORT=~{model_report}

    gatk --java-options "-Xms~{java_mem}m -Xmx~{max_heap}m" \
      VariantRecalibrator \
      -V ~{sites_only_variant_filtered_vcf} \
      -O ~{recalibration_filename} \
      --tranches-file ~{tranches_filename} \
      --trust-all-polymorphic \
      -tranche ~{sep=' -tranche ' recalibration_tranche_values} \
      -an ~{sep=' -an ' recalibration_annotation_values} \
      ~{true='--use-allele-specific-annotations' false='' use_allele_specific_annotations} \
      -mode SNP \
      ~{model_report_arg} \
      --max-gaussians ~{max_gaussians} \
      -resource:hapmap,known=false,training=true,truth=true,prior=15 ~{hapmap_resource_vcf} \
      -resource:omni,known=false,training=true,truth=true,prior=12 ~{omni_resource_vcf} \
      -resource:1000G,known=false,training=true,truth=false,prior=10 ~{one_thousand_genomes_resource_vcf} \
      -resource:dbsnp,known=true,training=false,truth=false,prior=7 ~{dbsnp_resource_vcf}
  >>>

  runtime {
    memory: "~{machine_mem} MiB"
    cpu: 16
    bootDiskSizeGb: 15
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 1
    docker: gatk_docker
  }

  output {
    File recalibration = "~{recalibration_filename}"
    File recalibration_index = "~{recalibration_filename}.idx"
    File tranches = "~{tranches_filename}"
  }
}
