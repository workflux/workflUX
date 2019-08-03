cwlVersion: v1.0
class: Workflow

inputs:
  sleep_time:
    type: int
    
steps:
  sleep:
    run: "../tools/sleep.cwl"
    in:
      sleep_time:
        source: sleep_time
    out: []

outputs: []
