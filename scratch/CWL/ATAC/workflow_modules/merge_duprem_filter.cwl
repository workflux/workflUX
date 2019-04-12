
cwlVersion: v1.0
class: Workflow

requirements:
  SubworkflowFeatureRequirement: {}
  ScatterFeatureRequirement: {}
  StepInputExpressionRequirement: {}
  InlineJavascriptRequirement: {}

### INPUT PART:
##################################################
inputs:
  sample_id:
    type: string
  bams:
    type:
      type: array
      items: File
  is_paired_end:
    type: boolean
        
### WORKFLOW STEPS:
##################################################
steps:
  lane_replicate_merging:
    doc: samtools merge - merging bam files of lane replicates
    run: "../tools/samtools_merge.cwl"
    in:
      bams:
        source: bams
      output_name:
        source: sample_id
        valueFrom: $(self + ".bam")
  
    out:
       - bam_merged

  sorting_merged_bam:
    doc: samtools sort - sorting of merged bam
    run: "../tools/samtools_sort.cwl"
    in:
      bam_unsorted:
        source: lane_replicate_merging/bam_merged
    out:
       - bam_sorted

  remove_duplicates:
    doc: picard markdup - emoves duplicates from a single sorted bam file.
    run: "../tools/picard_markdup.cwl"
    in:
      bam_sorted:
        source: sorting_merged_bam/bam_sorted
    out:
      - bam_duprem
      - picard_markdup_stdout

  filter_by_mapq:
    doc: samtools view
    run: "../tools/samtools_view_filter.cwl"
    in:
      bam:
        source: remove_duplicates/bam_duprem
      is_paired_end:
        source: is_paired_end
    out:
      - bam_filtered

  sorting_filtered_bam:
    doc: samtools sort - sorting of filtered bam
    run: "../tools/samtools_sort.cwl"
    in:
      bam_unsorted:
        source: filter_by_mapq/bam_filtered
    out:
       - bam_sorted

  indexing_filtered_bam:
    doc: |
      samtools index - indexes sorted bam
    run: "../tools/samtools_index_hack.cwl"
    in:
      bam_sorted:
        source: sorting_filtered_bam/bam_sorted
    out:
       - bam_sorted_indexed

  qc_post_filtering:
    doc: fastqc - quality control for reads directly after mapping
    run: "../tools/fastqc.cwl"
    in:
      bam:
        source: indexing_filtered_bam/bam_sorted_indexed
    out:
      - fastqc_zip
      - fastqc_html
      
### OUTPUTS:
##################################################
outputs:
  #bam_merged:
  #  type: File
  #  outputSource: lane_replicate_merging/bam_merged
  #bam_merged_sorted:
  #  type: File
  #  outputSource: sorting_merged_bam/bam_sorted
  #bam_merged_duprem:
  #  type: File
  #  outputSource: remove_duplicates/bam_duprem
  #bam_merged_duprem_filtered:
  #  type: File
  #  outputSource: filter_by_mapq/bam_filtered
  post_filter_fastqc_zip:
    type: 
      type: array
      items: File
    outputSource: qc_post_filtering/fastqc_zip
  post_filter_fastqc_html:
    type: 
      type: array
      items: File
    outputSource: qc_post_filtering/fastqc_html
  bam:
    type: File
    secondaryFiles: .bai
    outputSource: indexing_filtered_bam/bam_sorted_indexed
  picard_markdup_stdout:
    type: File
    outputSource: remove_duplicates/picard_markdup_stdout
    
  
    
    
    