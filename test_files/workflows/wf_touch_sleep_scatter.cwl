cwlVersion: v1.0
class: Workflow

requirements:
  ScatterFeatureRequirement: {}

inputs:
  sleep_time:
    type: 
      type: array
      items: int
  filename:
    type: 
      type: array
      items: string
    
steps:
  sleep:
    scatter: [sleep_time]
    scatterMethod: 'dotproduct'
    run: "../tools/sleep.cwl"
    in:
      sleep_time:
        source: sleep_time
    out: []
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
