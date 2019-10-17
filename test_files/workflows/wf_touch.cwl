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
      - test_file

outputs:
  test_file:
    type: File
    outputSource: touch/test_file
