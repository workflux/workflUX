cwlVersion: v1.0
class: CommandLineTool
hints:
  DockerRequirement:
    dockerPull: ubuntu:latest
  ResourceRequirement:
    coresMin: 1
    ramMin: 5000
    tmpdirMin: 1000

baseCommand: ["sleep"]

inputs:
  sleep_time:
    type: int
    inputBinding:
      position: 1

outputs: []
  
  
