# Please note: this version is a workaround for following problem: https://github.com/common-workflow-language/cwltool/issues/775
# If this will be fixed in the future, please use the version named without "_hack".

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
    dockerPull: kerstenbreuer/multiqc:1.7
  

baseCommand: ["bash", "-c"]
arguments:
  - valueFrom: |
      ${
          var qc_files_array = inputs.qc_files_array;
          var qc_files_array_of_array = inputs.qc_files_array_of_array;
          var cmdline = "echo 'copying input file ...'";

          if ( qc_files_array != null ){
            for (var i=0; i<qc_files_array.length; i++){
              if( qc_files_array[i] != null ){
                cmdline += "; cp " + qc_files_array[i].path + " .";
              }
            }
          }

          if ( qc_files_array_of_array != null ){
            for (var i=0; i<qc_files_array_of_array.length; i++){ 
              for (var ii=0; ii<qc_files_array_of_array[i].length; ii++){
                if( qc_files_array_of_array[i][ii] != null ){
                  cmdline += "; cp " + qc_files_array_of_array[i][ii].path + " .";
                }
              }
            }
          }
          
          cmdline += "; echo \'copying done\'" +
              "; multiqc --zip-data-dir --cl_config \'log_filesize_limit: 100000000\' " +
              "--outdir " + runtime.outdir +
              " --filename " + inputs.report_name + "_report .";

          return cmdline
        }
  
inputs:
  qc_files_array:
    doc: |
      qc files which shall be part of the multiqc summary;
      optional, only one of qc_files_array or qc_files_array_of_array 
      must be provided
    type:
      - "null"
      - type: array
        items: [File, "null"]
  qc_files_array_of_array:
    doc: |
      qc files which shall be part of the multiqc summary;
      optional, only one of qc_files_array or qc_files_array_of_array 
      must be provided
    type:
      - "null"
      - type: array
        items: 
          type: array
          items: [File, "null"]
  report_name:
    doc: name used for the html report and the corresponding zip file
    type: string
    default: multiqc
      
outputs:
  multiqc_zip:
    type: File
    outputBinding:
      glob: $(inputs.report_name + "_report_data.zip")
  multiqc_html:
    type: File
    outputBinding:
      glob: $(inputs.report_name + "_report.html")
  