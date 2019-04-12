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
    dockerPull: kerstenbreuer/samtools:1.7
  
### BASE COMMAND AND ARGUMENTS:
##################################################
baseCommand: ["bash", "-c"]
arguments:
  - valueFrom: $("LC_COLLATE=C sort -k1,1 -k2,2n " + inputs.bedgraph.path)
stdout: |
  ${
    if( inputs.output_name == null ){
      return inputs.bedgraph.basename;
    }
    else{
      return inputs.output_name;
    }
  }

### INPUT PART:
##################################################
inputs:
  bedgraph:
    type: File
  output_name:
    doc: optional, if not specified, output file will be named as input file
    type: string?

 
### OUTPUT PART:
##################################################
outputs:
  bedgraph_sorted:
    type: stdout
    