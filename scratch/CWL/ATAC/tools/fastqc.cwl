cwlVersion: v1.0
class: CommandLineTool
requirements:
  InlineJavascriptRequirement: {}
  StepInputExpressionRequirement: {}
hints:
  ResourceRequirement:
    coresMin: 1
    ramMin: 5000
  DockerRequirement:
    dockerPull: kerstenbreuer/trim_galore:0.4.4_1.14_0.11.7
  
baseCommand: "fastqc"
arguments: 
  - valueFrom: $(runtime.outdir)
    prefix: "-o"
    # specifies output directory
  - valueFrom: "--noextract"
    # reported file will be zipped

#stdout: $(inputs.log_file_name)

inputs:
  fastq1:
    type: File?
    inputBinding:
      position: 1
  fastq2:
    type: File?
    inputBinding:
      position: 2
  bam:
    type: File?
    inputBinding:
      position: 1

  #log_file_name:
  #  type: string
  #  doc: |
  #    only used for the log file; naming of the rest 
  #    refers to the fastq basename
 
outputs:
  fastqc_zip:
    doc: all data e.g. figures
    type:
      type: array
      items: File
    outputBinding:
      glob: "*_fastqc.zip"
  fastqc_html:
    doc: html report showing results from zip
    type:
      type: array
      items: File
    outputBinding:
      glob: "*_fastqc.html"
  #fastqc_log:
  #  doc: stdout log
  #  type: File
  #  outputBinding:
  #    glob: "*_fastqc.log"
    
