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
baseCommand: ["bedtools", "genomecov"]
arguments:
  - valueFrom: "-bg"
    position: 1
stdout: $(inputs.bed.nameroot + ".bedgraph")
  

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
  bedgraph:
    type: stdout
    