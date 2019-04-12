# clips features exceeding the chromosome boundaries

cwlVersion: v1.0
class: CommandLineTool
requirements:
  InlineJavascriptRequirement: {}
  StepInputExpressionRequirement: {}
hints:
  ResourceRequirement:
    coresMin: 1
    ramMin: 15000
    tmpdirMin: 10000
  DockerRequirement:
    dockerPull: biocontainers/bedtools:2.25.0

### BASE COMMAND AND ARGUMENTS:
##################################################
baseCommand: ["bedtools", "slop"]
arguments:
  - valueFrom: "0"
    prefix: -b
    position: 1
stdout: $(inputs.bed.nameroot + "_clipped.bed")
  

### INPUT PART:
##################################################
inputs:
  bed:
    type: File
    inputBinding:
      prefix: "-i"
      position: 2
  reference_info:
    type: File
    inputBinding:
      prefix: "-g"
      position: 3
 
### OUTPUT PART:
##################################################
outputs:
  bed_clipped:
    type: stdout
    