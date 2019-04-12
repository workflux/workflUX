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
  sam:
    doc: aligned reads to be checked in sam or bam format
    type: File
    inputBinding:
      position: 2

arguments:
  - valueFrom: -h
    position: 1
    # include the headers
  - valueFrom: -b
    position: 1
    # output in bam format

stdout: $(inputs.sam.nameroot).bam

outputs:
  bam_unsorted:
    type: stdout
  
  
