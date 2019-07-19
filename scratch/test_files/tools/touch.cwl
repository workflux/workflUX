cwlVersion: v1.0
class: CommandLineTool

hints:
  ResourceRequirement:
    coresMin: 1
    ramMin: 100
    #tmpdirMin: 10000
  
baseCommand: ["touch"]

inputs:
  filename:
    type: string
    inputBinding:
      position: 1
    
outputs:
  file:
    type: File
    outputBinding:
      glob: $(inputs.filename)