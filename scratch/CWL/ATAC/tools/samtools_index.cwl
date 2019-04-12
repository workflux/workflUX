cwlVersion: v1.0
class: CommandLineTool
requirements:
  InitialWorkDirRequirement:
    listing: 
      - $(inputs.bam_sorted)
hints:
  ResourceRequirement:
    coresMin: 1
    ramMin: 20000
    #tmpdirMin: 10000
  DockerRequirement:
    dockerPull: kerstenbreuer/samtools:1.7
  
### BASE COMMAND AND ARGUMENTS:
##################################################
baseCommand: ["samtools", "index"]
arguments:
  - valueFrom: -b  # specifies that index is created in bai format
    position: 1

### INPUT PART:
##################################################
inputs:
  bam_sorted:
    doc: sorted bam input file
    type: File
    inputBinding:
      position: 2
 
### OUTPUT PART:
##################################################
outputs:
  bam_sorted_indexed:
    type: File
    secondaryFiles: .bai
    outputBinding:
      glob: $(inputs.bam_sorted.basename)
      
    