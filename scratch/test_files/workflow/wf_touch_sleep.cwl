cwlVersion: v1.0
class: Workflow

inputs:
  filename:
    type: string
    
steps:
  touch:
    run: "../tools/touch.cwl"
    in:
      filename:
        source: filename
    out:
      - file

  sleep:
    run: "../tools/sleep.cwl"
    in:
      file:
        source: touch/file
    out: []
       
outputs:
  file:
    type: File
    outputSource: touch/file
    
  
    
    
    