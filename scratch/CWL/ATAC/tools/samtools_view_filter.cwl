# for single end data

cwlVersion: v1.0
class: CommandLineTool
hints:
  ResourceRequirement:
    coresMin: 1
    ramMin: 10000
  DockerRequirement:
    dockerPull: kerstenbreuer/samtools:1.7

baseCommand: ["samtools", "view"]

inputs:
  bam:
    doc: aligned reads to be checked in bam format
    type: File
    inputBinding:
      position: 10
  is_paired_end:
    doc: if paired end, only properly paired reads pass
    type: boolean
    default: true

arguments:
  - valueFrom: -h
    position: 1
    # include the headers
  - valueFrom: -b
    position: 1
    # output in bam format
  - valueFrom: "4"
    prefix: -F
    position: 1
  - valueFrom: "20"
    # only include reads with mapping quality >= 20
    prefix: -q
    position: 1
  # when paired end, "-f 3" is added to the command line
  - valueFrom: |
      ${
        if ( inputs.is_paired_end ){
           return "-f";
        }
        else {
          return null;
        }
      }
    position: 2
  - valueFrom: |
      ${
        if ( inputs.is_paired_end ){
           return "3";
        }
        else {
          return null;
        }
      }
    position: 3
stdout: $(inputs.bam.nameroot + "_filt.bam")

outputs:
  bam_filtered:
    type: stdout
  
  
