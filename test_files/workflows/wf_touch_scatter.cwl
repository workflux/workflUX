cwlVersion: v1.0
class: Workflow

requirements:
  ScatterFeatureRequirement: {}

inputs:
  filename:
    type: 
      type: array
      items: string
    
steps:
  touch:
    scatter: [filename]
    scatterMethod: 'dotproduct'
    run: "../tools/touch.cwl"
    in:
      filename:
        source: filename
    out:
      - test_file

outputs:
  test_file:
    type: 
      type: array
      items: File
    outputSource: touch/test_file
