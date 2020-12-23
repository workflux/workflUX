cwlVersion: v1.0
class: Workflow
doc: |
  bla

requirements:
  SubworkflowFeatureRequirement: {}
  ScatterFeatureRequirement: {}
  StepInputExpressionRequirement: {}
  InlineJavascriptRequirement: {}

inputs:
  fastq1:
    type: File
    doc: test
  fastq2: 
    type: File
        
steps:
  fastqc:
    doc: fastqc - quality control for trimmed fastq
    run: "../tools/fastqc.cwl"
    in:
      fastq1:
        source: fastq1
      fastq2:
        source: fastq2
    out:
      - fastqc_zip
      - fastqc_html

outputs:
  fastqc_zip:
    type:
      type: array
      items: File
    outputSource: fastqc/fastqc_zip
  fastqc_html:
    type:
      type: array
      items: File
    outputSource: fastqc/fastqc_html