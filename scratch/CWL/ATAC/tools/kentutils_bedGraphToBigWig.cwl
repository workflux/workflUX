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
    dockerPull: kerstenbreuer/kent-bedgraphtobigwig:latest
  
### BASE COMMAND AND ARGUMENTS:
##################################################
baseCommand: ["bedGraphToBigWig"]
arguments:
  - valueFrom: $(inputs.bedgraph_sorted.nameroot + ".bigwig")
    position: 3
  

### INPUT PART:
##################################################
inputs:
  bedgraph_sorted:
    type: File
    inputBinding:
      position: 1
  reference_info:
    type: File
    inputBinding:
      position: 2

    
 
### OUTPUT PART:
##################################################
outputs:
  bigwig:
    type: File
    outputBinding:
      glob: $(inputs.bedgraph_sorted.nameroot + ".bigwig")
    