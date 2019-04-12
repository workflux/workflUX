cwlVersion: v1.0
class: Workflow

### INPUT PART:
##################################################
inputs:
  bed:
    type: File
  reference_info:
    type: File

### WORKFLOW STEPS:
##################################################
steps:
  converting_bed_to_bam:
    doc: |
      bedtools bedtobam - converts bed to bam;
      as most tools can handle bam as input but not always
      the bed format(e.g. deeptools), moreover, bam is compressed;
      therefore, it will be used as final output (instead of the bed file)
    run: "../tools/bedtools_bedtobam.cwl"
    in:
      bed:
        source: bed
      reference_info:
        source: reference_info
    out:
      - bam

  sorting_bam:
    doc: samtools sort - sorting of merged bam
    run: "../tools/samtools_sort.cwl"
    in:
      bam_unsorted:
        source: converting_bed_to_bam/bam
    out:
       - bam_sorted
       
  indexing_bam:
    doc: |
      samtools index - indexes sorted bam
    run: "../tools/samtools_index_hack.cwl"
    in:
      bam_sorted:
        source: sorting_bam/bam_sorted
    out:
       - bam_sorted_indexed

  converting_bed_to_bedgraph:
    doc:  bedtools genomeCov -bg
    run: "../tools/bedtools_genomecov.cwl"
    in:
      bed:
        source: bed
      reference_info:
        source: reference_info
    out:
      - bedgraph

  sorting_bedgraph:
    doc: LC_COLLATE=C sort -k1,1 -k2,2n 
    run: "../tools/bedgraph_sort.cwl"
    in:
      bedgraph:
        source: converting_bed_to_bedgraph/bedgraph
    out:
      - bedgraph_sorted

  converting_bedgraph_to_bigwig:
    doc: bedGraphToBigWig (kentUtils)
    run: "../tools/kentutils_bedGraphToBigWig.cwl"
    in:
      bedgraph_sorted:
        source: sorting_bedgraph/bedgraph_sorted
      reference_info:
        source: reference_info
    out:
      - bigwig

### OUTPUTS:
##################################################
outputs:
  bigwig:
    type: File
    outputSource: converting_bedgraph_to_bigwig/bigwig
  bam:
    type: File
    outputSource: indexing_bam/bam_sorted_indexed
    
  
    
    
    