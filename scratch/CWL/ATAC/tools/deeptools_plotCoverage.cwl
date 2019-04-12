
cwlVersion: v1.0
class: CommandLineTool
requirements:
  InlineJavascriptRequirement: {}
  StepInputExpressionRequirement: {}
hints:
  ResourceRequirement:
    coresMin: 1
    ramMin: 15000
    #tmpdirMin: 10000
  DockerRequirement:
    dockerPull: kerstenbreuer/deeptools:3.1.1
  
baseCommand: ["plotCoverage"]
arguments:
  - valueFrom: |
      ${
        if ( inputs.is_paired_end ){
           return null;
        }
        else {
          return "--extendReads";
        }
      }
    position: 1
  - valueFrom: |
      ${
        if ( inputs.is_paired_end ){
           return null;
        }
        else {
          return inputs.fragment_size;
        }
      }
    position: 2
  - valueFrom: $(inputs.sample_id)
    prefix: --labels
    position: 10
  - valueFrom: $(inputs.sample_id + ".plot_cov.png")
    prefix: --plotFile
    position: 10
  - valueFrom: $(inputs.sample_id + ".plot_cov.tsv")
    prefix: --outRawCounts
    position: 10
  
inputs:
  bam:
    doc: must be indexed
    type: File
    secondaryFiles: .bai
    inputBinding:
        position: 100
        prefix: --bamfiles
  fragment_size:
    type: int?
  sample_id:
    type: string
  is_paired_end:
    doc: if paired end, reads are extended
    type: boolean
    default: true
    
outputs:
  qc_plot_coverage_plot:
    type: File
    outputBinding:
      glob: $(inputs.sample_id + ".plot_cov.png")
  qc_plot_coverage_tsv:
    type: File
    outputBinding:
      glob: $(inputs.sample_id + ".plot_cov.tsv")
    