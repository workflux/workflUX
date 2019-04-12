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
  fastq1:
    type: File
  fastq2: 
    type: File?
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
  adapter1: 
    type: [string, "null"]
  adapter2:
    type: [string, "null"]
  is_paired_end:
    type: boolean
  max_mapping_insert_length:
    type: long?
  sample_id:
    type: string
        
### WORKFLOW STEPS:
##################################################
steps:
  qc_pre_trim:
    doc: fastqc - quality control for trimmed fastq
    run: "../tools/fastqc.cwl"
    in:
      fastq1:
        source: fastq1
      fastq2:
        source: fastq2
    out:
      - fastqc_zip
      - fastqc_html

  adaptor_trimming_and_qc_post_trim:
    doc: trim galore - adapter trimming using trim_galore
    run: "../tools/trim_galore.cwl"
    in:
      fastq1:
        source: fastq1
      fastq2:
        source: fastq2
      adapter1:
        source: adapter1
      adapter2:
        source: adapter2   
    out:
      - fastq1_trimmed
      - fastq2_trimmed
      - fastq1_trimmed_unpaired
      - fastq2_trimmed_unpaired
      - trim_galore_log
      - post_trim_fastqc_html
      - post_trim_fastqc_zip

  mapping:
    doc: bowite2 - mapper, produces sam file
    run: "../tools/bowtie2.cwl"
    in:
      fastq1:
        source: adaptor_trimming_and_qc_post_trim/fastq1_trimmed
      fastq2:
        source: adaptor_trimming_and_qc_post_trim/fastq2_trimmed
      reference_index:
        source: reference
      is_paired_end:
        source: is_paired_end
      max_mapping_insert_length:
        source: max_mapping_insert_length
      output_basename:
        source: sample_id
    out:
      - sam
      - bowtie2_log
  sam2bam:
    doc: samtools view - convert sam to bam
    run: "../tools/samtools_view_sam2bam.cwl"
    in:
      sam:
        source: mapping/sam
    out:
      - bam_unsorted
  sort_bam: 
    doc: samtools sort - sorts unsorted bam file by coordinates.
    run: "../tools/samtools_sort.cwl"
    in:
      bam_unsorted:
        source: sam2bam/bam_unsorted
    out:
      - bam_sorted
      
### OUTPUTS:
##################################################
outputs:
  pre_trim_fastqc_zip:
    type:
      type: array
      items: File
    outputSource: qc_pre_trim/fastqc_zip
  pre_trim_fastqc_html:
    type:
      type: array
      items: File
    outputSource: qc_pre_trim/fastqc_html
  fastq1_trimmed:
    type: File
    outputSource: adaptor_trimming_and_qc_post_trim/fastq1_trimmed
  fastq2_trimmed:
    type: File
    outputSource: adaptor_trimming_and_qc_post_trim/fastq2_trimmed
  trim_galore_log:
    type:
      type: array
      items: File
    outputSource: adaptor_trimming_and_qc_post_trim/trim_galore_log
  post_trim_fastqc_html:
    type:
      type: array
      items: File
    outputSource: adaptor_trimming_and_qc_post_trim/post_trim_fastqc_html
  post_trim_fastqc_zip:
    type:
      type: array
      items: File
    outputSource: adaptor_trimming_and_qc_post_trim/post_trim_fastqc_zip
  bowtie2_log:
    type: File
    outputSource: mapping/bowtie2_log
  #sam:
  #  type: File
  #  outputSource: mapping/sam
  #bam_unsorted:
  #  type: File
  #  outputSource: sam2bam/bam_unsorted
  bam:
    type: File
    outputSource: sort_bam/bam_sorted
    