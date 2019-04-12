cwlVersion: v1.0
class: Workflow

requirements:
  SubworkflowFeatureRequirement: {}
  ScatterFeatureRequirement: {}
  StepInputExpressionRequirement: {}
  InlineJavascriptRequirement: {}
  MultipleInputFeatureRequirement: {}

### INPUT PART:
##################################################
inputs:
  sample_id:
    type: string
  fastq1:
    type: File[]
  fastq2: 
    type:
      type: array
      items: File
  reference:
    type: File
    secondaryFiles:
      - .fai
      - ^.1.bt2
      - ^.2.bt2
      - ^.3.bt2
      - ^.4.bt2
      - ^.rev.1.bt2
      - ^.rev.2.bt2
  reference_info:
    type: File
    default: "test"
  macs2_genome_size:
    type: string
  adapter1: 
    type: 
      - "null"
      - type: array
        items: [string, "null"]
  adapter2:
    type: 
      - "null"
      - type: array
        items: string
    default: "test"
  max_mapping_insert_length:
    type: long
    default: 2500
  shift_bps_upstream:
    type: float
    default: 37.5
  extend_bps_downstream:
    type: int
    default: 73
  
steps:
  #########################################################################################################
  ## upstream processing
  trim_and_map:
    run: "../workflow_modules/trim_and_map.cwl"
    scatter: [fastq1, fastq2, adapter1, adapter2]
    scatterMethod: 'dotproduct'
    in:
      fastq1:
        source: fastq1
      fastq2: 
        source: fastq2
      reference:
        source: reference
      adapter1: 
        source: adapter1
      adapter2:
        source: adapter2
      is_paired_end:
        default: true
      max_mapping_insert_length:
        source: max_mapping_insert_length
      sample_id:
        source: sample_id
    out:
      - pre_trim_fastqc_zip
      - pre_trim_fastqc_html
      - fastq1_trimmed
      - fastq2_trimmed
      - trim_galore_log
      - post_trim_fastqc_html
      - post_trim_fastqc_zip
      - bam
      - bowtie2_log
 
  merge_duprem_filter:
    run: "../workflow_modules/merge_duprem_filter.cwl"
    in:
      sample_id:
        source: sample_id
      bams: 
        source: trim_and_map/bam
      is_paired_end:
        default: true
    out:
      - post_filter_fastqc_zip
      - post_filter_fastqc_html
      - picard_markdup_stdout
      - bam

  name_sorting_filtered_bam:
      doc: samtools sort - sorting of filtered bam file by read name
      run: "../tools/samtools_sort_name.cwl"
      in:
        bam_unsorted:
          source: merge_duprem_filter/bam
      out:
        - bam_sorted
       
  converting_bam_to_bedpe:
    doc: bedtools bamtobed
    run: "../tools/bedtools_bamtobed_pe.cwl"
    in:
      bam:
        source: name_sorting_filtered_bam/bam_sorted
    out:
      - bedpe

  ##################################################################################################
  ## Generating ATAC signal (coverage tracks and peaks)
  ## The standard procedure is to shift the read ends to represent the center of the transposition
  ## event. However, not only the center of the transposition event represents accessible:
  ## This approach consideres the bp length which is covered by Tn5 enzyme 
  ## upon binding. Moreover, the span between the binding region is classified with respect to the
  ## binding region size of a nucleosome.
  ## Thereby, four signal tracks are generated:
  ## (1) accessible/"naked" DNA: 
  ##      - includes the Tn5 binding regions (TBR) of all reads
  ## (2) nucleosome free:
  ##      - includes the TBR of all reads
  ##      - plus the region between two TBR (=insert region) if it is < 147 bp (size 
  ##        of nucleosome binding region)
  ## (3) potentially nucleosome bound:
  ##      - exclude the TBR of all reads
  ##      - include insert region if it is > 147 bp
  ## (4) open chromatin / fragment tracks excluding tn5 binding region:
  ##      - include the TBR of all reads
  ##      - include insert region of all read pairs
  ## This approach precisely reflects the accessiblity of the genome.
  ## Moreover, it might be especially intresting for nucleosome positioning.
  ## However, the advantages of this approach are rarely tested and 
  ## there are also good reasons for using the established read shifting method 
  ## (e.g. better comparability towards other studies).
  ## Therefore, this workflow produces standard signal tags 
  ## (where reads are shifted towards the center of the transposition event)
  ## as proposed by the original ATAC paper (https://www.nature.com/articles/nmeth.2688)
  ## and calls peaks on them as discribed in the encode pipeline.

  generating_atac_signal_tags:
    doc: 
    run: "../tools/generate_atac_signal_tags.cwl"
    in:
      bedpe_alignm: # paired alignments in bedpe format
        source: converting_bam_to_bedpe/bedpe
      output_basename:
        source: sample_id
    out:
      # all output files are already sorted by coordinate
      - bed_tn5_bind_region_signal # (1)
      - bed_tn5_center_1bp_signal # (as original ATAC paper, only 1bp of the center of the transp. event)
      - bed_nucl_free_signal # (2)
      - bed_nucl_bound_signal # (3)
      - bed_fragments_tn5_incl_signal # (4)
      - fragment_sizes_tsv
      - filtering_stats_tsv
      - frag_size_stats_tsv
      - irreg_mappings_bedpe

  ## generate signal tracks:

  generating_tn5_bind_region_signal_tracks:
    doc:
    run: "../workflow_modules/bed_to_coverage_track.cwl"
    in:
      bed:
        source: generating_atac_signal_tags/bed_tn5_bind_region_signal
      reference_info:
        source: reference_info
    out:
      - bigwig
      - bam

  generating_tn5_center_1bp_signal_tracks:
    doc:
    run: "../workflow_modules/bed_to_coverage_track.cwl"
    in:
      bed:
        source: generating_atac_signal_tags/bed_tn5_center_1bp_signal
      reference_info:
        source: reference_info
    out:
      - bigwig
      - bam

  generating_nucl_free_signal_tracks:
    doc:
    run: "../workflow_modules/bed_to_coverage_track.cwl"
    in:
      bed:
        source: generating_atac_signal_tags/bed_nucl_free_signal
      reference_info:
        source: reference_info
    out:
      - bigwig
      - bam

  generating_nucl_bound_signal_tracks:
    doc:
    run: "../workflow_modules/bed_to_coverage_track.cwl"
    in:
      bed:
        source: generating_atac_signal_tags/bed_nucl_bound_signal
      reference_info:
        source: reference_info
    out:
      - bigwig
      - bam

  generating_fragments_tn5_excl_signal_tracks:
    doc:
    run: "../workflow_modules/bed_to_coverage_track.cwl"
    in:
      bed:
        source: generating_atac_signal_tags/bed_fragments_tn5_incl_signal
      reference_info:
        source: reference_info
    out:
      - bigwig
      - bam

  ############################################################################################################################
  ## Peak calling:
  
  ## peak calling and downstream QC - Tn5 binding region:

  peak_calling_macs2_tn5_bind_region:
    doc: peak calling using macs2
    run: "../tools/macs2_callpeak_atac_tn5_binding_region.cwl"
    in:
      treatment_bed:
        source: generating_atac_signal_tags/bed_tn5_bind_region_signal
      output_basename:
        source: sample_id
        valueFrom: $(self + "_tn5_bind_region.macs2")
      genome_size:
        source: macs2_genome_size
    out: 
      - broad_peaks_bed
      - gapped_peaks_bed
      - peaks_xls


  ## peak calling and downstream QC - Tn5 center:
  ## (following the peak calling approach of the encode pipeline)
  peak_calling_macs2_tn5_center:
    doc: peak calling using macs2
    run: "../tools/macs2_callpeak_atac_tn5_center.cwl"
    in:
      treatment_bed:
        source: generating_atac_signal_tags/bed_tn5_center_1bp_signal
      output_basename:
        source: sample_id
        valueFrom: $(self + "_tn5_center_shift_ext.macs2")
      genome_size:
        source: macs2_genome_size
    out: 
      - broad_peaks_bed
      - gapped_peaks_bed
      - peaks_xls
      - treat_pileup_bdg # coverage tracks of shifted and extended reads according to encode

  ## convert shifted and extended coverage tracks to bigwig:

  clip_shift_ext_bedgraph:
    doc: clips features exceeding the chromosome boundaries
    run: "../tools/bedtools_slop_clip.cwl"
    in:
      bed:
        source: peak_calling_macs2_tn5_center/treat_pileup_bdg
      reference_info:
        source: reference_info
    out:
      - bed_clipped

  sorting_shift_ext_bedgraph:
    doc: LC_COLLATE=C sort -k1,1 -k2,2n 
    run: "../tools/bedgraph_sort.cwl"
    in:
      bedgraph:
        source: clip_shift_ext_bedgraph/bed_clipped
      output_name:
        source: sample_id
        valueFrom: $(self + "_tn5_center_shift_ext.bedgraph")
    out:
      - bedgraph_sorted

  converting_shift_ext_bedgraph_to_bigwig:
    doc: bedGraphToBigWig (kentUtils)
    run: "../tools/kentutils_bedGraphToBigWig.cwl"
    in:
      bedgraph_sorted:
        source: sorting_shift_ext_bedgraph/bedgraph_sorted
      reference_info:
        source: reference_info
    out:
      - bigwig

  #########################################################################################################
  ## Nucleosome position calling using NucleoATAC:
  
  nucl_position_calling:
    doc: NucleoATAC
    run: "../tools/nucleoatac.cwl"
    in:
      bam:
        source: merge_duprem_filter/bam
      bed:
        source: peak_calling_macs2_tn5_center/broad_peaks_bed
      fasta:
        source: reference
      output_basename:
        source: sample_id
        valueFrom: $(self + ".nucleoatac")
    out:
      - nucl_occ_tracks
      - nucl_occ_lower_bound_tracks
      - nucl_occ_upper_bound_tracks
      - nucl_dist_txt
      - nucl_dist_plot
      - fragsize_in_peaks_txt
      - nucl_occ_fit_txt
      - nucl_occ_fit_plot
      - nucl_occ_peaks_bed
      - nucl_vplot_data
      - nucl_pos_bed
      - nucl_pos_redundant_bed
      - nucl_norm_crosscor_tracks
      - nucl_norm_smooth_crosscor_tracks
      - combined_nucl_pos_bed
      - nfr_pos_bed
      - nucleoatac_stderr
      - nucleoatac_stdout

  ##########################################################################################################################
  ## Quality Controls:

  plot_fragment_size_distribution:
    run: "../tools/plot_frag_size_distr.cwl"
    in:
      fragment_sizes_tsv:
        source: generating_atac_signal_tags/fragment_sizes_tsv
      output_basename:
        source: sample_id
    out:
      - frag_size_distr_plot
      - frag_size_distr_tsv

  ## plot genome coverage:
  # (with respect to the complete fragment between a reads pair - i.e. the "open chrom" tracks)
  qc_plot_coverage_fragments_tn5_excl:
    doc: |
      deeptools plotCoverage - plots how many times a certain fraction of the 
      genome was covered (consideres the complete fragment between a reads pair).
    run: "../tools/deeptools_plotCoverage.cwl"
    in:
      bam:
        source: generating_fragments_tn5_excl_signal_tracks/bam
      sample_id:
        source: sample_id
    out:
      - qc_plot_coverage_plot  
      - qc_plot_coverage_tsv

  qc_plot_fingerprint:
    run: "../tools/deeptools_plotFingerprint.cwl"
    in:
      bam:
        source: generating_tn5_bind_region_signal_tracks/bam
      sample_id:
        source: sample_id
      is_paired_end:
        default: true
    out:
      - qc_plot_fingerprint_plot  
      - qc_plot_fingerprint_tsv
      - qc_plot_fingerprint_stderr

  qc_phantompeakqualtools:
    run: "../tools/phantompeakqualtools.cwl"
    in:
      bam:
        source: merge_duprem_filter/bam
    out:
      - qc_crosscorr_summary  
      - qc_crosscorr_plot
      - qc_phantompeakqualtools_stderr
      - qc_phantompeakqualtools_stdout

  create_summary_qc_report:
    doc: |
      multiqc summarizes the qc results from fastqc 
      and other tools
    run: "../tools/multiqc_hack.cwl"
    in:
      qc_files_array_of_array:
        source:
          - trim_and_map/pre_trim_fastqc_zip
          - trim_and_map/pre_trim_fastqc_html
          - trim_and_map/post_trim_fastqc_html
          - trim_and_map/post_trim_fastqc_zip
          - trim_and_map/trim_galore_log
        linkMerge: merge_flattened
      qc_files_array:
        source:
          - trim_and_map/bowtie2_log
          - merge_duprem_filter/post_filter_fastqc_zip
          - merge_duprem_filter/post_filter_fastqc_html
          - generating_atac_signal_tags/frag_size_stats_tsv
          - generating_atac_signal_tags/fragment_sizes_tsv
          - generating_atac_signal_tags/filtering_stats_tsv
          - plot_fragment_size_distribution/frag_size_distr_tsv
          - peak_calling_macs2_tn5_bind_region/peaks_xls
          - peak_calling_macs2_tn5_center/peaks_xls
          - qc_plot_coverage_fragments_tn5_excl/qc_plot_coverage_tsv
          - qc_plot_fingerprint/qc_plot_fingerprint_tsv
          - qc_phantompeakqualtools/qc_phantompeakqualtools_stdout
          - qc_phantompeakqualtools/qc_crosscorr_summary
          - merge_duprem_filter/picard_markdup_stdout
        linkMerge: merge_flattened
      report_name:
        source: sample_id
    out:
      - multiqc_zip
      - multiqc_html

  ##########################################################################################################################

outputs:
  pre_trim_fastqc_zip:
    type:
      type: array
      items: 
        type: array
        items: File
    outputSource: trim_and_map/pre_trim_fastqc_zip
  pre_trim_fastqc_html:
    type:
      type: array
      items: 
        type: array
        items: File
    outputSource: trim_and_map/pre_trim_fastqc_html
  trim_galore_log:
    type:
      type: array
      items: 
        type: array
        items: File
    outputSource: trim_and_map/pre_trim_fastqc_zip
  post_trim_fastqc_html:
    type:
      type: array
      items: 
        type: array
        items: File
    outputSource: trim_and_map/post_trim_fastqc_html
  post_trim_fastqc_zip:
    type:
      type: array
      items: 
        type: array
        items: File
    outputSource: trim_and_map/post_trim_fastqc_zip
  bowtie2_log:
    type:
      type: array
      items: File
    outputSource: trim_and_map/bowtie2_log

  post_filter_fastqc_zip:
    type:
      type: array
      items: File
    outputSource: merge_duprem_filter/post_filter_fastqc_zip
  post_filter_fastqc_html:
    type:
      type: array
      items: File
    outputSource: merge_duprem_filter/post_filter_fastqc_html
  bam:
    type: File
    secondaryFiles: .bai
    outputSource: merge_duprem_filter/bam
  picard_markdup_stdout:
    type: File
    outputSource: merge_duprem_filter/picard_markdup_stdout

  frag_size_stats_tsv:
    type: File
    outputSource: generating_atac_signal_tags/frag_size_stats_tsv
  filtering_stats_tsv:
    type: File
    outputSource: generating_atac_signal_tags/filtering_stats_tsv
  fragment_sizes_tsv:
    type: File
    outputSource: generating_atac_signal_tags/fragment_sizes_tsv
  irreg_mappings_bedpe:
    type: File
    outputSource: generating_atac_signal_tags/irreg_mappings_bedpe

  bigwig_tn5_bind_region_signal:
    type: File
    outputSource: generating_tn5_bind_region_signal_tracks/bigwig 
  bigwig_tn5_center_1bp_signal:
    type: File
    outputSource: generating_tn5_center_1bp_signal_tracks/bigwig 
  bigwig_nucl_free_signal:
    type: File
    outputSource: generating_nucl_free_signal_tracks/bigwig  
  bigwig_nucl_bound_signal:
    type: File
    outputSource: generating_nucl_bound_signal_tracks/bigwig
  bigwig_fragments_tn5_excl_signal:
    type: File
    outputSource: generating_fragments_tn5_excl_signal_tracks/bigwig

  bam_tn5_bind_region_signal:
    type: File
    outputSource: generating_tn5_bind_region_signal_tracks/bam 
  bam_tn5_center_1bp_signal:
    type: File
    outputSource: generating_tn5_center_1bp_signal_tracks/bam 
#  bam_nucl_free_signal:
#    type: File
#    outputSource: generating_nucl_free_signal_tracks/bam  
#  bam_nucl_bound_signal:
#    type: File
#    outputSource: generating_nucl_bound_signal_tracks/bam
#  bam_fragments_tn5_excl_signal:
#    type: File
#    outputSource: generating_fragments_tn5_excl_signal_tracks/bam

  broad_peaks_bed_tn5_bind_region_signal:
    type: File
    outputSource: peak_calling_macs2_tn5_bind_region/broad_peaks_bed
  gapped_peaks_bed_tn5_bind_region_signal:
    type: File
    outputSource: peak_calling_macs2_tn5_bind_region/gapped_peaks_bed
  peaks_xls_tn5_bind_region_signal:
    type: File
    outputSource: peak_calling_macs2_tn5_bind_region/peaks_xls
    
  broad_peaks_bed_tn5_center_shift_ext_signal:
    type: File
    outputSource: peak_calling_macs2_tn5_center/broad_peaks_bed
  gapped_peaks_bed_tn5_center_shift_ext_signal:
    type: File
    outputSource: peak_calling_macs2_tn5_center/gapped_peaks_bed
  peaks_xls_tn5_center_shift_ext_signal:
    type: File
    outputSource: peak_calling_macs2_tn5_center/peaks_xls
  bigwig_tn5_center_shift_ext_signal:
    type: File
    outputSource: converting_shift_ext_bedgraph_to_bigwig/bigwig

  nucl_occ_tracks:
    type: File?
    outputSource: nucl_position_calling/nucl_occ_tracks
  nucl_occ_lower_bound_tracks:
    type: File?
    outputSource: nucl_position_calling/nucl_occ_lower_bound_tracks
  nucl_occ_upper_bound_tracks:
    type: File?
    outputSource: nucl_position_calling/nucl_occ_upper_bound_tracks
  nucl_dist_txt:
    type: File?
    outputSource: nucl_position_calling/nucl_dist_txt
  nucl_dist_plot:
    type: File?
    outputSource: nucl_position_calling/nucl_dist_plot
  fragsize_in_peaks_txt:
    type: File?
    outputSource: nucl_position_calling/fragsize_in_peaks_txt
  nucl_occ_fit_txt:
    type: File?
    outputSource: nucl_position_calling/nucl_occ_fit_txt
  nucl_occ_fit_plot:
    type: File?
    outputSource: nucl_position_calling/nucl_occ_fit_plot
  nucl_occ_peaks_bed:
    type: File?
    outputSource: nucl_position_calling/nucl_occ_peaks_bed
  nucl_vplot_data:
    type: File?
    outputSource: nucl_position_calling/nucl_vplot_data
  nucl_pos_bed:
    type: File?
    outputSource: nucl_position_calling/nucl_pos_bed
  nucl_pos_redundant_bed:
    type: File?
    outputSource: nucl_position_calling/nucl_pos_redundant_bed
  nucl_norm_crosscor_tracks:
    type: File?
    outputSource: nucl_position_calling/nucl_norm_crosscor_tracks
  nucl_norm_smooth_crosscor_tracks:
    type: File?
    outputSource: nucl_position_calling/nucl_norm_smooth_crosscor_tracks
  combined_nucl_pos_bed:
    type: File?
    outputSource: nucl_position_calling/combined_nucl_pos_bed
  nfr_pos_bed:
    type: File?
    outputSource: nucl_position_calling/nfr_pos_bed
  nucleoatac_stderr:
    type: File
    outputSource: nucl_position_calling/nucleoatac_stderr
  nucleoatac_stdout:
    type: File
    outputSource: nucl_position_calling/nucleoatac_stdout

  frag_size_distr_plot:
    type: File
    outputSource: plot_fragment_size_distribution/frag_size_distr_plot
  frag_size_distr_tsv:
    type: File
    outputSource: plot_fragment_size_distribution/frag_size_distr_tsv
  coverage_plot_fragments_tn5_excl:
    type: File
    outputSource: qc_plot_coverage_fragments_tn5_excl/qc_plot_coverage_plot
  coverage_counts_fragments_tn5_excl:
    type: File
    outputSource: qc_plot_coverage_fragments_tn5_excl/qc_plot_coverage_tsv
  qc_plot_fingerprint_plot:
    type: File?
    outputSource: qc_plot_fingerprint/qc_plot_fingerprint_plot
  qc_plot_fingerprint_tsv:
    type: File?
    outputSource: qc_plot_fingerprint/qc_plot_fingerprint_tsv
  qc_plot_fingerprint_stderr:
    type: File
    outputSource: qc_plot_fingerprint/qc_plot_fingerprint_stderr
  qc_crosscorr_summary:
    type: File?
    outputSource: qc_phantompeakqualtools/qc_crosscorr_summary
  qc_crosscorr_plot:
    type: File?
    outputSource: qc_phantompeakqualtools/qc_crosscorr_plot
  qc_phantompeakqualtools_stderr:
    type: File?
    outputSource: qc_phantompeakqualtools/qc_phantompeakqualtools_stderr

  multiqc_zip:
    type: File
    outputSource: create_summary_qc_report/multiqc_zip
  multiqc_html:
    type: File
    outputSource: create_summary_qc_report/multiqc_html