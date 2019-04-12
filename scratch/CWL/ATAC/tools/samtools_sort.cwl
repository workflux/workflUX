doc: Sort a bam file by read names.
cwlVersion: v1.0
class: CommandLineTool
hints:
  ResourceRequirement:
    coresMin: 4
    ramMin: 15000
    #ramMin: 200 for testing on small device
  DockerRequirement:
    dockerPull: kerstenbreuer/samtools:1.7

baseCommand: ["samtools", "sort"]
arguments:
  - valueFrom: $(runtime.cores)
    prefix: -@

inputs:
  bam_unsorted:
    doc: aligned reads to be checked in sam or bam format
    type: File
    inputBinding:
      position: 2

#arguments:
  #- prefix: -m
    #valueFrom: ${ return (parseInt(15000/hints.ResourceRequirement.ramMin-100)) + "M" }
    #position: 1
    # specifies the allowed maximal memory usage per thread before
    # samtools start to outsource memory to temporary files

stdout: $(inputs.bam_unsorted.basename)

outputs:
  bam_sorted:
    type: stdout
  
  
