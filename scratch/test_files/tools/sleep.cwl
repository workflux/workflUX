cwlVersion: v1.0
class: CommandLineTool
requirements:
  InlineJavascriptRequirement: {}

hints:
  ResourceRequirement:
    coresMin: 1
    ramMin: 100
    #tmpdirMin: 10000
  
baseCommand: ["sleep"]
arguments:
  - valueFrom: ${return(Math.random()*10)}
    position: 1

inputs:
  file:
    type: File
    
outputs: []