cwlVersion: v1.0
class: CommandLineTool
requirements:
  InlineJavascriptRequirement: {}
  StepInputExpressionRequirement: {}
hints:
  ResourceRequirement:
    coresMin: 4
    ramMin: 30000
    #tmpdirMin: 20000
  DockerRequirement:
    dockerPull: kerstenbreuer/bowtie2:2.2.6-2

baseCommand: ["bowtie2"]
arguments:
  - valueFrom: --very-sensitive
    position: 1
  - valueFrom: $(runtime.cores) # set the number of threads
    prefix: "-p"
    position: 1
  - position: 10 # prefix for fastq1, differs for paired/single end
    valueFrom: |
      ${
        if ( inputs.is_paired_end ){
           return "-1";
        }
        else {
          return "-U";
        }
      }
  - valueFrom: $(inputs.output_basename + ".sam") # set the number of threads
    prefix: "-S"
    position: 6
stderr: $( inputs.output_basename + ".bowtie2_stderr") # log file
  

inputs:
  reference_index:
    doc: path to the FM-index files for the chosen reference genome
    type: File
    secondaryFiles:
      - .fai
      - ^.1.bt2
      - ^.2.bt2
      - ^.3.bt2
      - ^.4.bt2
      - ^.rev.1.bt2
      - ^.rev.2.bt2
    inputBinding:
      position: 2
      prefix: "-x"
      valueFrom: $(self.path.replace(/\.fa/i,""))
  fastq1:
    type: File
    inputBinding:
      position: 11
  is_paired_end:
    type: boolean
  fastq2:
    type: File?
    inputBinding:
      valueFrom: |
        ${
            if ( inputs.is_paired_end ){
                return self;
            }
            else {
              return null;
            }
        }  
      position: 12
      prefix: "-2"
  max_mapping_insert_length:
    doc: usefull for very long fragments, as expected for ATAC
    type: long?
    inputBinding:
      prefix: --maxins
      position: 1
  output_basename:
    type: string

      
outputs:
  sam:
    type: File
    outputBinding:
      glob: "*.sam"
  bowtie2_log:
    type: stderr
    
