cwlVersion: v1.0
class: CommandLineTool
requirements:
  - class: InlineJavascriptRequirement
  - class: StepInputExpressionRequirement
hints:
  ResourceRequirement:
    coresMin: 1
    ramMin: 20000
    #tmpdirMin: 10000
  DockerRequirement:
    dockerPull: kerstenbreuer/deeptools:3.1.1

baseCommand: ["bamCoverage"]
arguments:  
  - valueFrom: --extendReads
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
  - valueFrom: $(inputs.bam.nameroot + ".bigwig")
    prefix: --outFileName
    position: 10
  - valueFrom: "bigwig"
    prefix: --outFileFormat
    position: 10
  - valueFrom: "RPGC"
    prefix: --normalizeUsing
    position: 10
  
inputs:
  bam:
    doc: bam file as input; needs bai index file in the same directory
    type: File
    secondaryFiles: .bai
    inputBinding:
        position: 100
        prefix: --bam
  is_paired_end:
    doc: if false, reads are extended by fragment_size
    type: boolean
  fragment_size:
    doc: mean library fragment size; used to extend the reads
    type: long?
  effective_genome_size:
    doc: |
      the effectively mappable genome size, 
      see: https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html
    type: long
    inputBinding:
        position: 10
        prefix: --effectiveGenomeSize
  bin_size:
    type: int
    default: 10
    inputBinding:
      prefix: --binSize
      position: 10
  ignoreForNormalization:
    type:
      type: array
      items: string
    default: ["chrX", "chrY", "chrM"]
    inputBinding:
      prefix: --ignoreForNormalization
      position: 10
  spike_in_count:
    doc: number of reads aligned to the spike in reference, optional
    type: long?
    inputBinding:
      position: 10
      prefix: --scaleFactor
      valueFrom: |
        ${ 
          if( self == null ){
            return null
          }
          else{
            (1.0 / parseFloat(self)) 
          }
        }

outputs:
  bigwig:
    type: File
    outputBinding:
      glob: $(inputs.bam.nameroot + ".bigwig")
    