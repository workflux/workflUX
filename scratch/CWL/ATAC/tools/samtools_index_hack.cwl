# Please note: this version is a workaround for following problem: https://github.com/common-workflow-language/cwltool/issues/775
# If this will be fixed in the future, please use the version named without "_hack".

cwlVersion: v1.0
class: CommandLineTool
hints:
  ResourceRequirement:
    coresMin: 1
    ramMin: 20000
    #tmpdirMin: 10000
  DockerRequirement:
    dockerPull: kerstenbreuer/samtools:1.7
  
### BASE COMMAND AND ARGUMENTS:
##################################################
baseCommand: ["bash", "-c"]
arguments:
  - valueFrom: $("cp " + inputs.bam_sorted.path + " . && samtools index -b " + inputs.bam_sorted.basename )

### INPUT PART:
##################################################
inputs:
  bam_sorted:
    doc: sorted bam input file
    type: File
 
### OUTPUT PART:
##################################################
outputs:
  bam_sorted_indexed:
    type: File
    secondaryFiles: .bai
    outputBinding:
      glob: $(inputs.bam_sorted.basename)
      
    