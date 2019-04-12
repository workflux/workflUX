cwlVersion: v1.0
class: CommandLineTool
requirements:
  InlineJavascriptRequirement: {}
hints:
  ResourceRequirement:
    coresMin: 1
    ramMin: 15000
    #tmpdirMin: 10000
  DockerRequirement:
    dockerPull: biocontainers/bedtools:2.25.0
  
### BASE COMMAND AND ARGUMENTS:
##################################################
baseCommand: ["bedtools", "bamtobed"]
arguments:
  - valueFrom: -bedpe
    position: 1
stdout: $(inputs.bam.nameroot + ".bedpe")

### INPUT PART:
##################################################
inputs:
  bam:
    type: File
    inputBinding:
      prefix: -i
      position: 10
 
### OUTPUT PART:
##################################################
outputs:
  bedpe:
    type: stdout