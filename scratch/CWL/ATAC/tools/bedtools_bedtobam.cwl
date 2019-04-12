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
baseCommand: ["bedtools", "bedtobam"]
stdout: $(inputs.bed.nameroot + ".bam")      

### INPUT PART:
##################################################
inputs:
  bed:
    type: File
    inputBinding:
      position: 1
      prefix: "-i"
  reference_info:
    type: File
    inputBinding:
      position: 2
      prefix: "-g"
    
 
### OUTPUT PART:
##################################################
outputs:
  bam:
    type: stdout