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
  fastq1_trimmed:
    type: 
      type: array
      items: File
  fastq2_trimmed: 
    type: 
      type: array
      items: [File, "null"]
  reference_spike_in:
    type: File
    secondaryFiles:
      - .fai
      - ^.1.bt2
      - ^.2.bt2
      - ^.3.bt2
      - ^.4.bt2
      - ^.rev.1.bt2
      - ^.rev.2.bt2
  is_paired_end:
    type: boolean
  sample_id:
    type: string
        
### WORKFLOW STEPS:
##################################################
steps:
  mapping:
    doc: bowite2 - mapper, produces sam file
    run: "../tools/bowtie2.cwl"
    scatter: [fastq1, fastq2]
    scatterMethod: 'dotproduct'
    in:
      fastq1:
        source: fastq1_trimmed
      fastq2:
        source: fastq2_trimmed
      reference_index:
        source: reference_spike_in
      is_paired_end:
        source: is_paired_end
    out:
      - sam
      - bowtie2_log

  sam2bam:
    doc: samtools view - convert sam to bam
    run: "../tools/samtools_view_sam2bam.cwl"
    scatter: [sam]
    scatterMethod: 'dotproduct'
    in:
      sam:
        source: mapping/sam
    out:
      - bam_unsorted  

  lane_replicate_merging:
    doc: samtools merge - merging bam files of lane replicates
    run: "../tools/samtools_merge.cwl"
    in:
      bams:
        source: sam2bam/bam_unsorted
      output_name:
        source: sample_id
        valueFrom: $(self + "_spike_in.bam")
  
    out:
       - bam_merged

  sort_bam: 
    doc: samtools sort - sorts unsorted bam file by coordinates.
    run: "../tools/samtools_sort.cwl"
    in:
      bam_unsorted:
        source: lane_replicate_merging/bam_merged
    out:
      - bam_sorted

  duplicate_removal:
    doc: remove duplicate reads
    run: "../tools/picard_markdup.cwl"
    in:
      bam_sorted:
        source: sort_bam/bam_sorted
    out:
      - bam_duprem

  count_aligned_reads:
    run: "../tools/samtools_view_count_alignments.cwl"
    in:
      bam:
        source: duplicate_removal/bam_duprem
      is_paired_end:
        source: is_paired_end
    out:
      - aln_read_count_file
      - aln_read_count

      
### OUTPUTS:
##################################################
outputs:
  bam:
    type: File
    outputSource: duplicate_removal/bam_duprem
  aln_read_count:
    type: long
    outputSource: count_aligned_reads/aln_read_count
  aln_read_count_file:
    type: File
    outputSource: count_aligned_reads/aln_read_count_file
  
    