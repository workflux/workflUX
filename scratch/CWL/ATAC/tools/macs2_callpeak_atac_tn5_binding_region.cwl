cwlVersion: v1.0
class: CommandLineTool
requirements:
  InlineJavascriptRequirement: {}
  StepInputExpressionRequirement: {}
hints:
  ResourceRequirement:
    coresMin: 1
    ramMin: 10000
  DockerRequirement:
    dockerPull: genomicpariscentre/macs2:2.1.0.20140616
 
### BASE COMMAND AND ARGUMENTS:
##################################################
baseCommand: ["macs2", "callpeak"]
arguments:  
  - valueFrom: "BED"
    prefix: "--format"
    position: 1
  
  - valueFrom: "--broad"
    position: 2
  - valueFrom: "--nomodel"
    position: 2
  - valueFrom: "all"
    prefix: "--keep-dup"
    position: 2

  # following settings remain the tn5 binding region bed
  # as it is, as all entries are already 29 bp long.
  # (However, these parameters are expected by MACS2, 
  # if --nomodel is set)
  - valueFrom: "0"
    prefix: "--shift"
    position: 2
  - valueFrom: "29"
    prefix: "--extsize"
    position: 2

  - valueFrom: "--bdg"
    position: 2
#stdout: $(inputs.output_basename + "_macs2.log")
#stderr: $(inputs.output_basename + "_macs2.logErr")

### INPUT PART:
##################################################
inputs:
  treatment_bed:
    type: File
    inputBinding:
        position: 100
        prefix: "--treatment"

  output_basename:
    doc: gives the base name for all output files
    type: string
    inputBinding:
        position: 101
        prefix: "--name"
  
  genome_size:
    doc: can be "mm", "hs", "ce", "dm", or the total number of genomic bp 
    type: string
    inputBinding:
        position: 3
        prefix: "--gsize"

 
### OUTPUT PART:
##################################################
outputs:
  broad_peaks_bed:    
    type: File
    outputBinding:
      glob: "*.broadPeak"
  gapped_peaks_bed:    
    type: File
    outputBinding:
      glob: "*.gappedPeak"
  #summits_bed:
  #  type: File
  #  outputBinding:
  #    glob: "*_summits.bed"
  peaks_xls:
    type: File
    outputBinding:
      glob: "*_peaks.xls"
  treat_pileup_bdg:
    type: File
    outputBinding:
      glob: "*_treat_pileup.bdg"
  #macs2_log:
  #  type: stdout
  #macs2_log_err:
  #  type: stderr
      
    
  
    