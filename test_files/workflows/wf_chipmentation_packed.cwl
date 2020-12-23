{
    "cwlVersion": "v1.0", 
    "$graph": [
        {
            "inputs": [
                {
                    "inputBinding": {
                        "position": 4
                    }, 
                    "type": "File", 
                    "id": "#bowtie2.cwl/fastq1"
                }, 
                {
                    "inputBinding": {
                        "position": 5, 
                        "valueFrom": "${\n    if ( inputs.is_paired_end ){\n        return self;\n    }\n    else {\n      return null;\n    }\n}  \n", 
                        "prefix": "-2"
                    }, 
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "id": "#bowtie2.cwl/fastq2"
                }, 
                {
                    "type": "boolean", 
                    "id": "#bowtie2.cwl/is_paired_end"
                }, 
                {
                    "doc": "path to the FM-index files for the chosen reference genome", 
                    "inputBinding": {
                        "position": 2, 
                        "prefix": "-x", 
                        "valueFrom": "$(self.path.replace(/\\.fa/i,\"\"))"
                    }, 
                    "type": "File", 
                    "id": "#bowtie2.cwl/reference_index", 
                    "secondaryFiles": [
                        ".fai", 
                        "^.1.bt2", 
                        "^.2.bt2", 
                        "^.3.bt2", 
                        "^.4.bt2", 
                        "^.rev.1.bt2", 
                        "^.rev.2.bt2"
                    ]
                }
            ], 
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }, 
                {
                    "class": "StepInputExpressionRequirement"
                }
            ], 
            "outputs": [
                {
                    "type": "stderr", 
                    "id": "#bowtie2.cwl/bowtie2_log"
                }, 
                {
                    "outputBinding": {
                        "glob": "*.sam"
                    }, 
                    "type": "File", 
                    "id": "#bowtie2.cwl/sam"
                }
            ], 
            "baseCommand": [
                "bowtie2"
            ], 
            "id": "#bowtie2.cwl", 
            "arguments": [
                {
                    "position": 1, 
                    "valueFrom": "--very-sensitive"
                }, 
                {
                    "position": 1, 
                    "valueFrom": "$(runtime.cores)", 
                    "prefix": "-p"
                }, 
                {
                    "position": 3, 
                    "valueFrom": "${\n  if ( inputs.is_paired_end ){\n     return \"-1\";\n  }\n  else {\n    return \"-U\";\n  }\n}\n"
                }, 
                {
                    "position": 6, 
                    "valueFrom": "$(inputs.fastq1.nameroot + \".sam\")", 
                    "prefix": "-S"
                }
            ], 
            "stderr": "$( inputs.fastq1.nameroot + \"_bowtie2.log\")", 
            "class": "CommandLineTool", 
            "hints": [
                {
                    "dockerPull": "kerstenbreuer/bowtie2:2.2.6-2", 
                    "class": "DockerRequirement"
                }, 
                {
                    "coresMin": 4, 
                    "ramMin": 30000, 
                    "class": "ResourceRequirement"
                }
            ]
        }, 
        {
            "inputs": [
                {
                    "doc": "bam file as input; needs bai index file in the same directory", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "--bam"
                    }, 
                    "type": "File", 
                    "id": "#deeptools_bamCoverage.cwl/bam", 
                    "secondaryFiles": ".bai"
                }, 
                {
                    "doc": "the effectively mappable genome size, \nsee: https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html\n", 
                    "inputBinding": {
                        "position": 4, 
                        "prefix": "--effectiveGenomeSize"
                    }, 
                    "type": "long", 
                    "id": "#deeptools_bamCoverage.cwl/effective_genome_size"
                }, 
                {
                    "doc": "mean library fragment size; used to extend the reads", 
                    "inputBinding": {
                        "position": 4, 
                        "prefix": "--extendReads"
                    }, 
                    "type": "int", 
                    "id": "#deeptools_bamCoverage.cwl/fragment_size"
                }
            ], 
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }, 
                {
                    "class": "StepInputExpressionRequirement"
                }
            ], 
            "outputs": [
                {
                    "outputBinding": {
                        "glob": "$(inputs.bam.nameroot + \".bigwig\")"
                    }, 
                    "type": "File", 
                    "id": "#deeptools_bamCoverage.cwl/bigwig"
                }
            ], 
            "baseCommand": [
                "bamCoverage"
            ], 
            "id": "#deeptools_bamCoverage.cwl", 
            "arguments": [
                {
                    "position": 2, 
                    "valueFrom": "$(inputs.bam.nameroot + \".bigwig\")", 
                    "prefix": "--outFileName"
                }, 
                {
                    "position": 3, 
                    "valueFrom": "bigwig", 
                    "prefix": "--outFileFormat"
                }, 
                {
                    "position": 5, 
                    "valueFrom": "RPGC", 
                    "prefix": "--normalizeUsing"
                }, 
                {
                    "position": 5, 
                    "valueFrom": "chrX", 
                    "prefix": "--ignoreForNormalization"
                }, 
                {
                    "position": 5, 
                    "valueFrom": "10", 
                    "prefix": "--binSize"
                }
            ], 
            "class": "CommandLineTool", 
            "hints": [
                {
                    "dockerPull": "kerstenbreuer/deeptools:3.1.1", 
                    "class": "DockerRequirement"
                }, 
                {
                    "coresMin": 1, 
                    "ramMin": 20000, 
                    "class": "ResourceRequirement"
                }
            ]
        }, 
        {
            "inputs": [
                {
                    "inputBinding": {
                        "position": 1
                    }, 
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "id": "#fastqc.cwl/bam"
                }, 
                {
                    "inputBinding": {
                        "position": 1
                    }, 
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "id": "#fastqc.cwl/fastq1"
                }, 
                {
                    "inputBinding": {
                        "position": 2
                    }, 
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "id": "#fastqc.cwl/fastq2"
                }
            ], 
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }, 
                {
                    "class": "StepInputExpressionRequirement"
                }
            ], 
            "outputs": [
                {
                    "doc": "html report showing results from zip", 
                    "outputBinding": {
                        "glob": "*_fastqc.html"
                    }, 
                    "type": {
                        "items": "File", 
                        "type": "array"
                    }, 
                    "id": "#fastqc.cwl/fastqc_html"
                }, 
                {
                    "doc": "all data e.g. figures", 
                    "outputBinding": {
                        "glob": "*_fastqc.zip"
                    }, 
                    "type": {
                        "items": "File", 
                        "type": "array"
                    }, 
                    "id": "#fastqc.cwl/fastqc_zip"
                }
            ], 
            "baseCommand": "fastqc", 
            "id": "#fastqc.cwl", 
            "arguments": [
                {
                    "valueFrom": "$(runtime.outdir)", 
                    "prefix": "-o"
                }, 
                {
                    "valueFrom": "--noextract"
                }
            ], 
            "class": "CommandLineTool", 
            "hints": [
                {
                    "dockerPull": "kerstenbreuer/trim_galore:0.4.4_1.14_0.11.7", 
                    "class": "DockerRequirement"
                }, 
                {
                    "coresMin": 1, 
                    "ramMin": 5000, 
                    "class": "ResourceRequirement"
                }
            ]
        }, 
        {
            "inputs": [
                {
                    "doc": "qc files which shall be part of the multiqc summary;\noptional, only one of qc_files_array or qc_files_array_of_array \nmust be provided\n", 
                    "type": [
                        "null", 
                        {
                            "items": "File", 
                            "type": "array"
                        }
                    ], 
                    "id": "#multiqc_hack.cwl/qc_files_array"
                }, 
                {
                    "doc": "qc files which shall be part of the multiqc summary;\noptional, only one of qc_files_array or qc_files_array_of_array \nmust be provided\n", 
                    "type": [
                        "null", 
                        {
                            "items": {
                                "items": "File", 
                                "type": "array"
                            }, 
                            "type": "array"
                        }
                    ], 
                    "id": "#multiqc_hack.cwl/qc_files_array_of_array"
                }, 
                {
                    "default": "multiqc", 
                    "doc": "name used for the html report and the corresponding zip file", 
                    "type": "string", 
                    "id": "#multiqc_hack.cwl/report_name"
                }
            ], 
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }, 
                {
                    "class": "StepInputExpressionRequirement"
                }
            ], 
            "outputs": [
                {
                    "outputBinding": {
                        "glob": "$(inputs.report_name + \"_report.html\")"
                    }, 
                    "type": "File", 
                    "id": "#multiqc_hack.cwl/multiqc_html"
                }, 
                {
                    "outputBinding": {
                        "glob": "$(inputs.report_name + \"_report_data.zip\")"
                    }, 
                    "type": "File", 
                    "id": "#multiqc_hack.cwl/multiqc_zip"
                }
            ], 
            "baseCommand": [
                "bash", 
                "-c"
            ], 
            "id": "#multiqc_hack.cwl", 
            "arguments": [
                {
                    "valueFrom": "${\n    var qc_files_array = inputs.qc_files_array;\n    var qc_files_array_of_array = inputs.qc_files_array_of_array;\n    var cmdline = \"echo 'copying input file ...'\";\n\n    if ( qc_files_array != null ){\n      for (var i=0; i<qc_files_array.length; i++){\n        cmdline += \"; cp \" + qc_files_array[i].path + \" .\";\n      }\n    }\n\n    if ( qc_files_array_of_array != null ){\n      for (var i=0; i<qc_files_array_of_array.length; i++){ \n        for (var ii=0; ii<qc_files_array_of_array[i].length; ii++){\n          cmdline += \"; cp \" + qc_files_array_of_array[i][ii].path + \" .\";\n        }\n      }\n    }\n    \n    cmdline += \"; echo \\'copying done\\'\" +\n        \"; multiqc --zip-data-dir --cl_config \\'log_filesize_limit: 100000000\\' \" +\n        \"--outdir \" + runtime.outdir +\n        \" --filename \" + inputs.report_name + \"_report .\";\n\n    return cmdline\n  }\n"
                }
            ], 
            "class": "CommandLineTool", 
            "hints": [
                {
                    "dockerPull": "kerstenbreuer/multiqc:1.7", 
                    "class": "DockerRequirement"
                }, 
                {
                    "coresMin": 1, 
                    "ramMin": 10000, 
                    "class": "ResourceRequirement"
                }
            ]
        }, 
        {
            "inputs": [
                {
                    "doc": "sorted bam input file", 
                    "inputBinding": {
                        "position": 11, 
                        "prefix": "INPUT=", 
                        "separate": false
                    }, 
                    "type": "File", 
                    "id": "#picard_markdup.cwl/bam_sorted"
                }, 
                {
                    "default": "/bin/picard.jar", 
                    "inputBinding": {
                        "position": 1
                    }, 
                    "type": "string", 
                    "id": "#picard_markdup.cwl/path_to_picards"
                }
            ], 
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }, 
                {
                    "class": "StepInputExpressionRequirement"
                }
            ], 
            "outputs": [
                {
                    "outputBinding": {
                        "glob": "$(inputs.bam_sorted.nameroot + \"_duprem.bam\")"
                    }, 
                    "type": "File", 
                    "id": "#picard_markdup.cwl/bam_duprem"
                }
            ], 
            "baseCommand": [
                "java", 
                "-jar"
            ], 
            "id": "#picard_markdup.cwl", 
            "arguments": [
                {
                    "position": 2, 
                    "valueFrom": "MarkDuplicates"
                }, 
                {
                    "position": 13, 
                    "valueFrom": "$(inputs.bam_sorted.nameroot + \"_duprem.bam\")", 
                    "prefix": "OUTPUT=", 
                    "separate": false
                }, 
                {
                    "position": 13, 
                    "valueFrom": "$(inputs.bam_sorted.nameroot + \"_duprem.log\")", 
                    "prefix": "METRICS_FILE=", 
                    "separate": false
                }, 
                {
                    "position": 14, 
                    "valueFrom": "REMOVE_DUPLICATES=TRUE"
                }, 
                {
                    "position": 15, 
                    "valueFrom": "ASSUME_SORTED=TRUE"
                }, 
                {
                    "position": 16, 
                    "valueFrom": "VALIDATION_STRINGENCY=SILENT"
                }
            ], 
            "class": "CommandLineTool", 
            "hints": [
                {
                    "dockerPull": "kerstenbreuer/picard_tools:2.17.4", 
                    "class": "DockerRequirement"
                }, 
                {
                    "coresMin": 1, 
                    "ramMin": 20000, 
                    "class": "ResourceRequirement"
                }
            ]
        }, 
        {
            "inputs": [
                {
                    "doc": "sorted bam input file", 
                    "type": "File", 
                    "id": "#samtools_index_hack.cwl/bam_sorted"
                }
            ], 
            "outputs": [
                {
                    "secondaryFiles": ".bai", 
                    "outputBinding": {
                        "glob": "$(inputs.bam_sorted.basename)"
                    }, 
                    "type": "File", 
                    "id": "#samtools_index_hack.cwl/bam_sorted_indexed"
                }
            ], 
            "baseCommand": [
                "bash", 
                "-c"
            ], 
            "id": "#samtools_index_hack.cwl", 
            "arguments": [
                {
                    "valueFrom": "$(\"cp \" + inputs.bam_sorted.path + \" . && samtools index -b \" + inputs.bam_sorted.basename )"
                }
            ], 
            "class": "CommandLineTool", 
            "hints": [
                {
                    "dockerPull": "kerstenbreuer/samtools:1.7", 
                    "class": "DockerRequirement"
                }, 
                {
                    "coresMin": 1, 
                    "ramMin": 20000, 
                    "class": "ResourceRequirement"
                }
            ]
        }, 
        {
            "inputs": [
                {
                    "doc": "name of merged bam file", 
                    "inputBinding": {
                        "position": 1
                    }, 
                    "type": "string", 
                    "id": "#samtools_merge.cwl/output_name"
                }, 
                {
                    "doc": "bam files to be merged", 
                    "inputBinding": {
                        "position": 2
                    }, 
                    "type": {
                        "items": "File", 
                        "type": "array"
                    }, 
                    "id": "#samtools_merge.cwl/bams"
                }
            ], 
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }, 
                {
                    "class": "StepInputExpressionRequirement"
                }
            ], 
            "outputs": [
                {
                    "outputBinding": {
                        "glob": "$(inputs.output_name)"
                    }, 
                    "type": "File", 
                    "id": "#samtools_merge.cwl/bam_merged"
                }
            ], 
            "baseCommand": [
                "samtools", 
                "merge"
            ], 
            "id": "#samtools_merge.cwl", 
            "class": "CommandLineTool", 
            "hints": [
                {
                    "dockerPull": "kerstenbreuer/samtools:1.7", 
                    "class": "DockerRequirement"
                }, 
                {
                    "coresMin": 1, 
                    "ramMin": 20000, 
                    "class": "ResourceRequirement"
                }
            ]
        }, 
        {
            "inputs": [
                {
                    "doc": "aligned reads to be checked in sam or bam format", 
                    "inputBinding": {
                        "position": 2
                    }, 
                    "type": "File", 
                    "id": "#samtools_sort.cwl/bam_unsorted"
                }
            ], 
            "stdout": "$(inputs.bam_unsorted.basename)", 
            "doc": "Sort a bam file by read names.", 
            "baseCommand": [
                "samtools", 
                "sort"
            ], 
            "id": "#samtools_sort.cwl", 
            "arguments": [
                {
                    "valueFrom": "$(runtime.cores)", 
                    "prefix": "-@"
                }
            ], 
            "outputs": [
                {
                    "type": "stdout", 
                    "id": "#samtools_sort.cwl/bam_sorted"
                }
            ], 
            "class": "CommandLineTool", 
            "hints": [
                {
                    "dockerPull": "kerstenbreuer/samtools:1.7", 
                    "class": "DockerRequirement"
                }, 
                {
                    "coresMin": 4, 
                    "ramMin": 15000, 
                    "class": "ResourceRequirement"
                }
            ]
        }, 
        {
            "inputs": [
                {
                    "doc": "aligned reads to be checked in bam format", 
                    "inputBinding": {
                        "position": 2
                    }, 
                    "type": "File", 
                    "id": "#samtools_view_filter.cwl/bam"
                }
            ], 
            "stdout": "$(inputs.bam.nameroot + \"_filt.bam\")", 
            "outputs": [
                {
                    "type": "stdout", 
                    "id": "#samtools_view_filter.cwl/bam_filtered"
                }
            ], 
            "baseCommand": [
                "samtools", 
                "view"
            ], 
            "id": "#samtools_view_filter.cwl", 
            "arguments": [
                {
                    "position": 1, 
                    "valueFrom": "-h"
                }, 
                {
                    "position": 1, 
                    "valueFrom": "-b"
                }, 
                {
                    "position": 1, 
                    "valueFrom": "4", 
                    "prefix": "-F"
                }, 
                {
                    "position": 1, 
                    "valueFrom": "20", 
                    "prefix": "-q"
                }
            ], 
            "class": "CommandLineTool", 
            "hints": [
                {
                    "dockerPull": "kerstenbreuer/samtools:1.7", 
                    "class": "DockerRequirement"
                }, 
                {
                    "coresMin": 1, 
                    "ramMin": 10000, 
                    "class": "ResourceRequirement"
                }
            ]
        }, 
        {
            "inputs": [
                {
                    "doc": "aligned reads to be checked in sam or bam format", 
                    "inputBinding": {
                        "position": 2
                    }, 
                    "type": "File", 
                    "id": "#samtools_view_sam2bam.cwl/sam"
                }
            ], 
            "stdout": "$(inputs.sam.nameroot).bam", 
            "outputs": [
                {
                    "type": "stdout", 
                    "id": "#samtools_view_sam2bam.cwl/bam_unsorted"
                }
            ], 
            "baseCommand": [
                "samtools", 
                "view"
            ], 
            "id": "#samtools_view_sam2bam.cwl", 
            "arguments": [
                {
                    "position": 1, 
                    "valueFrom": "-h"
                }, 
                {
                    "position": 1, 
                    "valueFrom": "-b"
                }
            ], 
            "class": "CommandLineTool", 
            "hints": [
                {
                    "dockerPull": "kerstenbreuer/samtools:1.7", 
                    "class": "DockerRequirement"
                }, 
                {
                    "coresMin": 1, 
                    "ramMin": 10000, 
                    "class": "ResourceRequirement"
                }
            ]
        }, 
        {
            "inputs": [
                {
                    "type": "File", 
                    "id": "#tn5_overhang_correction.cwl/bam"
                }, 
                {
                    "type": "boolean", 
                    "id": "#tn5_overhang_correction.cwl/is_paired_end"
                }, 
                {
                    "default": "tn5correct", 
                    "type": "string", 
                    "id": "#tn5_overhang_correction.cwl/out_suffix"
                }
            ], 
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }, 
                {
                    "class": "StepInputExpressionRequirement"
                }
            ], 
            "outputs": [
                {
                    "outputBinding": {
                        "glob": "$(inputs.bam.nameroot + \"_\" + inputs.out_suffix + \".bam\")"
                    }, 
                    "type": "File", 
                    "id": "#tn5_overhang_correction.cwl/bam_tn5_corrected"
                }
            ], 
            "baseCommand": [
                "bash", 
                "-c"
            ], 
            "id": "#tn5_overhang_correction.cwl", 
            "arguments": [
                {
                    "valueFrom": "${\n  var cmd_line = \"\";\n  \n  if ( inputs.is_paired_end ){ // for paired end data\n                                // unpaired will be removed\n  \n    ////// shift + strand reads\n    cmd_line += \"samtools view -h -f 3 -F 16 \" + inputs.bam.path; // only properly paired reads with the\n                                          // first read on the + strand are output;\n                                          // the header is included\n    \n    cmd_line += \" | awk  \\'BEGIN {OFS=\\\"\\\\t\\\"} \" +\n                \"{ \" + \n                \"if ( $1 ~ /^@/) { print }\" + //header lines are printed unmodified\n                \"else if ($9>=38) { \" + \n                \"$4=$4+4; $8=$8-5; $9=$9-9; $11=\\\"*\\\"; if ($8>0){ print } else {$8=1; print}}\" + // read start positions are shifted\n                \"else if ($9<=-38) { \" + \n                \"$4=$4+4; $8=$8-5; $9=$9+9; $11=\\\"*\\\"; if ($8>0){ print } else {$8=1; print}}\" + // read start positions are shifted\n                \"}\\' > correcting.sam\";\n    ///// shift - strand reads      \n    cmd_line += \" ; samtools view -f 19 \" + inputs.bam.path; // only properly paired reads with the\n                                          // first read on the - strand are output;\n                                          // the header is excluded\n    \n    cmd_line += \" | awk  \\'BEGIN {OFS=\\\"\\\\t\\\"} \" +\n                \"{ \" + \n                \"if ($9>=38) { \" + \n                \"$4=$4-5; $8=$8+4; $9=$9-9; $11=\\\"*\\\";  if ($4>0){ print } else {$4=1; print}}\" + // read start positions are shifted\n                \"else if ($9<=-38) { \" + \n                \"$4=$4-5; $8=$8+4; $9=$9+9; $11=\\\"*\\\"; if ($4>0){ print } else {$4=1; print}}\" + // read start positions are shifted\n                \"}\\' >> correcting.sam\";\n                \n    \n  }\n  else { // for single end data\n  \n    ////// shift + strand reads\n    cmd_line += \"samtools view -h -F 16 \" + inputs.bam.path; // paired end as well as\n                                          // - strand reads are excluded\n                                          // the header is included\n    \n    cmd_line += \" | awk  \\'BEGIN {OFS=\\\"\\\\t\\\"} \" +\n                \"{ \" + \n                \"if ( $1 ~ /^@/ ) { print }\" + //header lines are printed unmodified\n                \"else { $4=$4+4; $7=\\\"*\\\"; $8=0; $9=0; $11=\\\"*\\\"; print}\" + // read start positions are shifted\n                \"}\\' > correcting.sam\";\n    ///// shift - strand reads      \n    cmd_line += \" ; samtools view -f 16 \" + inputs.bam.path; // paired end as well as\n                                          // + strand reads are excluded\n                                          // the header is included\n    \n    cmd_line += \" | awk  \\'BEGIN {OFS=\\\"\\\\t\\\"} \" +\n                \"{ \" + \n                \"$4=$4-5; $7=\\\"*\\\"; $8=0; $9=0; $11=\\\"*\\\"; print \" + // read start positions are shifted\n                \"}\\' >> correcting.sam\";\n  }\n  \n  cmd_line += \" ; samtools sort -@ \" + runtime.cores + \" -O bam -T sorting.bam -o \" + inputs.bam.nameroot + \"_\" + inputs.out_suffix + \".bam correcting.sam\";\n  \n  return cmd_line;\n}\n  \n"
                }
            ], 
            "class": "CommandLineTool", 
            "hints": [
                {
                    "dockerPull": "kerstenbreuer/samtools:1.7", 
                    "class": "DockerRequirement"
                }, 
                {
                    "coresMin": 4, 
                    "ramMin": 20000, 
                    "class": "ResourceRequirement"
                }
            ]
        }, 
        {
            "inputs": [
                {
                    "doc": "Adapter sequence for first reads.\nif not specified, trim_galore will try to autodetect whether ...\n- Illumina universal adapter (AGATCGGAAGAGC)\n- Nextera adapter (CTGTCTCTTATA)\n- Illumina Small RNA 3' Adapter (TGGAATTCTCGG)\n... was used.\nYou can directly choose one of the above configurations\nby setting the string to \"illumina\", \"nextera\", or \"small_rna\".\n", 
                    "type": [
                        "null", 
                        "string"
                    ], 
                    "id": "#trim_galore.cwl/adapter1"
                }, 
                {
                    "doc": "Adapter sequence for second reads - only for paired end data.\nif not specified, trim_galore will try to autodetect whether ...\n- Illumina universal adapter (AGATCGGAAGAGC)\n- Nextera adapter (CTGTCTCTTATA)\n- Illumina Small RNA 3' Adapter (TGGAATTCTCGG)\n... was used.\nYou can directly choose one of the above configurations\nby setting the adapter1 string to \"illumina\", \"nextera\", or \"small_rna\".\n", 
                    "type": [
                        "null", 
                        "string"
                    ], 
                    "id": "#trim_galore.cwl/adapter2"
                }, 
                {
                    "doc": "raw reads in fastq format; can be gzipped;\nif paired end, the file contains the first reads;\nif single end, the file contains all reads\n", 
                    "inputBinding": {
                        "position": 10
                    }, 
                    "type": "File", 
                    "id": "#trim_galore.cwl/fastq1"
                }, 
                {
                    "doc": "(optional) raw reads in fastq format; can be gzipped;\nif paired end, the file contains the second reads;\nif single end, the file does not exist\n", 
                    "inputBinding": {
                        "position": 11
                    }, 
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "id": "#trim_galore.cwl/fastq2"
                }, 
                {
                    "default": 1, 
                    "doc": "minimum overlap with adapter seq in bp needed to trim", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "--stringency"
                    }, 
                    "type": "int", 
                    "id": "#trim_galore.cwl/min_adapter_overlap"
                }, 
                {
                    "default": 20, 
                    "doc": "discard reads that get shorter than this value", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "--length"
                    }, 
                    "type": "int", 
                    "id": "#trim_galore.cwl/min_read_length"
                }, 
                {
                    "default": 35, 
                    "doc": "if only one read of a pair passes the qc and adapter trimming,\nit needs at least this length to be rescued\n", 
                    "type": "int", 
                    "id": "#trim_galore.cwl/min_unpaired_read_rescue_length"
                }, 
                {
                    "default": 20, 
                    "doc": "trim all base with a phred score lower than this valueFrom", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "--quality"
                    }, 
                    "type": "int", 
                    "id": "#trim_galore.cwl/qual_trim_cutoff"
                }
            ], 
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ], 
            "outputs": [
                {
                    "outputBinding": {
                        "glob": "${\n    if ( inputs.fastq2 == null  ){ return \"*trimmed.fq*\" }\n    else { return \"*val_1.fq*\" }\n}\n"
                    }, 
                    "type": "File", 
                    "id": "#trim_galore.cwl/fastq1_trimmed"
                }, 
                {
                    "outputBinding": {
                        "glob": "*unpaired_1.fq*"
                    }, 
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "id": "#trim_galore.cwl/fastq1_trimmed_unpaired"
                }, 
                {
                    "outputBinding": {
                        "glob": "*val_2.fq*"
                    }, 
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "id": "#trim_galore.cwl/fastq2_trimmed"
                }, 
                {
                    "outputBinding": {
                        "glob": "*unpaired_2.fq*"
                    }, 
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "id": "#trim_galore.cwl/fastq2_trimmed_unpaired"
                }, 
                {
                    "doc": "html report of post-trimming fastqc", 
                    "outputBinding": {
                        "glob": "*fastqc.html"
                    }, 
                    "type": {
                        "items": "File", 
                        "type": "array"
                    }, 
                    "id": "#trim_galore.cwl/post_trim_fastqc_html"
                }, 
                {
                    "doc": "all data of post-trimming fastqc e.g. figures", 
                    "outputBinding": {
                        "glob": "*fastqc.zip"
                    }, 
                    "type": {
                        "items": "File", 
                        "type": "array"
                    }, 
                    "id": "#trim_galore.cwl/post_trim_fastqc_zip"
                }, 
                {
                    "outputBinding": {
                        "glob": "*trimming_report.txt"
                    }, 
                    "type": {
                        "items": "File", 
                        "type": "array"
                    }, 
                    "id": "#trim_galore.cwl/trim_galore_log"
                }
            ], 
            "baseCommand": "trim_galore", 
            "id": "#trim_galore.cwl", 
            "arguments": [
                {
                    "position": 1, 
                    "prefix": "--fastqc_args", 
                    "valueFrom": "\"--noextract\""
                }, 
                {
                    "position": 1, 
                    "prefix": "--gzip"
                }, 
                {
                    "position": 1, 
                    "valueFrom": "${\n  if ( inputs.adapter1 == \"illumina\" ){ return \"--illumina\" }\n  else if ( inputs.adapter1 == \"nextera\" ){ return \"--nextera\" }\n  else if ( inputs.adapter1 == \"small_rna\" ){ return \"--small_rna\" }\n  else { return null }\n}\n"
                }, 
                {
                    "position": 1, 
                    "prefix": "--adapter", 
                    "valueFrom": "${\n  if ( inputs.apdater1 != null && inputs.adapter1 != \"illumina\" && inputs.adapter1 != \"nextera\" && inputs.adapter1 != \"small_rna\" ){\n    return inputs.adapter1\n  } else {\n    return null\n  }\n}\n"
                }, 
                {
                    "position": 1, 
                    "prefix": "--adapter2", 
                    "valueFrom": "${\n  if ( inputs.fastq2 != null && inputs.apdater2 != null && inputs.adapter1 != \"illumina\" && inputs.adapter1 != \"nextera\" && inputs.adapter1 != \"small_rna\" ){\n    return inputs.adapter2\n  } else {\n    return null\n  }\n}\n"
                }, 
                {
                    "position": 1, 
                    "valueFrom": "${\n  if ( inputs.fastq2 == null ){ return null }\n  else { return \"--paired\" }\n}\n"
                }, 
                {
                    "position": 1, 
                    "valueFrom": "${\n  if ( inputs.fastq2 == null ){ return null }\n  else { return \"--retain_unpaired\" }\n}\n"
                }, 
                {
                    "position": 1, 
                    "prefix": "--length_1", 
                    "valueFrom": "${\n  if ( inputs.fastq2 == null ){ return null }\n  else { return inputs.min_unpaired_read_rescue_length }\n}\n"
                }, 
                {
                    "position": 1, 
                    "prefix": "--length_2", 
                    "valueFrom": "${\n  if ( inputs.fastq2 == null ){ return null }\n  else { return inputs.min_unpaired_read_rescue_length }\n}\n"
                }
            ], 
            "class": "CommandLineTool", 
            "hints": [
                {
                    "dockerPull": "kerstenbreuer/trim_galore:0.4.4_1.14_0.11.7", 
                    "class": "DockerRequirement"
                }, 
                {
                    "coresMin": 1, 
                    "ramMin": 7000, 
                    "class": "ResourceRequirement"
                }
            ]
        }, 
        {
            "inputs": [
                {
                    "type": {
                        "items": "File", 
                        "type": "array"
                    }, 
                    "id": "#merge_duprem_filter.cwl/bams"
                }, 
                {
                    "type": "string", 
                    "id": "#merge_duprem_filter.cwl/sample_id"
                }
            ], 
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }, 
                {
                    "class": "ScatterFeatureRequirement"
                }, 
                {
                    "class": "StepInputExpressionRequirement"
                }, 
                {
                    "class": "SubworkflowFeatureRequirement"
                }
            ], 
            "outputs": [
                {
                    "secondaryFiles": ".bai", 
                    "type": "File", 
                    "outputSource": "#merge_duprem_filter.cwl/indexing_filtered_bam/bam_sorted_indexed", 
                    "id": "#merge_duprem_filter.cwl/bam"
                }, 
                {
                    "outputSource": "#merge_duprem_filter.cwl/qc_post_filtering/fastqc_html", 
                    "type": {
                        "items": "File", 
                        "type": "array"
                    }, 
                    "id": "#merge_duprem_filter.cwl/post_filter_fastqc_html"
                }, 
                {
                    "outputSource": "#merge_duprem_filter.cwl/qc_post_filtering/fastqc_zip", 
                    "type": {
                        "items": "File", 
                        "type": "array"
                    }, 
                    "id": "#merge_duprem_filter.cwl/post_filter_fastqc_zip"
                }
            ], 
            "class": "Workflow", 
            "steps": [
                {
                    "doc": "samtools view", 
                    "out": [
                        "#merge_duprem_filter.cwl/filter_by_mapq/bam_filtered"
                    ], 
                    "run": "#samtools_view_filter.cwl", 
                    "id": "#merge_duprem_filter.cwl/filter_by_mapq", 
                    "in": [
                        {
                            "source": "#merge_duprem_filter.cwl/remove_duplicates/bam_duprem", 
                            "id": "#merge_duprem_filter.cwl/filter_by_mapq/bam"
                        }
                    ]
                }, 
                {
                    "doc": "samtools index - indexes sorted bam\n", 
                    "out": [
                        "#merge_duprem_filter.cwl/indexing_filtered_bam/bam_sorted_indexed"
                    ], 
                    "run": "#samtools_index_hack.cwl", 
                    "id": "#merge_duprem_filter.cwl/indexing_filtered_bam", 
                    "in": [
                        {
                            "source": "#merge_duprem_filter.cwl/sorting_filtered_bam/bam_sorted", 
                            "id": "#merge_duprem_filter.cwl/indexing_filtered_bam/bam_sorted"
                        }
                    ]
                }, 
                {
                    "doc": "samtools merge - merging bam files of lane replicates", 
                    "out": [
                        "#merge_duprem_filter.cwl/lane_replicate_merging/bam_merged"
                    ], 
                    "run": "#samtools_merge.cwl", 
                    "id": "#merge_duprem_filter.cwl/lane_replicate_merging", 
                    "in": [
                        {
                            "source": "#merge_duprem_filter.cwl/bams", 
                            "id": "#merge_duprem_filter.cwl/lane_replicate_merging/bams"
                        }, 
                        {
                            "source": "#merge_duprem_filter.cwl/sample_id", 
                            "valueFrom": "$(self + \".bam\")", 
                            "id": "#merge_duprem_filter.cwl/lane_replicate_merging/output_name"
                        }
                    ]
                }, 
                {
                    "doc": "fastqc - quality control for reads directly after mapping", 
                    "out": [
                        "#merge_duprem_filter.cwl/qc_post_filtering/fastqc_zip", 
                        "#merge_duprem_filter.cwl/qc_post_filtering/fastqc_html"
                    ], 
                    "run": "#fastqc.cwl", 
                    "id": "#merge_duprem_filter.cwl/qc_post_filtering", 
                    "in": [
                        {
                            "source": "#merge_duprem_filter.cwl/indexing_filtered_bam/bam_sorted_indexed", 
                            "id": "#merge_duprem_filter.cwl/qc_post_filtering/bam"
                        }
                    ]
                }, 
                {
                    "doc": "picard markdup - emoves duplicates from a single sorted bam file.", 
                    "out": [
                        "#merge_duprem_filter.cwl/remove_duplicates/bam_duprem"
                    ], 
                    "run": "#picard_markdup.cwl", 
                    "id": "#merge_duprem_filter.cwl/remove_duplicates", 
                    "in": [
                        {
                            "source": "#merge_duprem_filter.cwl/sorting_merged_bam/bam_sorted", 
                            "id": "#merge_duprem_filter.cwl/remove_duplicates/bam_sorted"
                        }
                    ]
                }, 
                {
                    "doc": "samtools sort - sorting of filtered bam", 
                    "out": [
                        "#merge_duprem_filter.cwl/sorting_filtered_bam/bam_sorted"
                    ], 
                    "run": "#samtools_sort.cwl", 
                    "id": "#merge_duprem_filter.cwl/sorting_filtered_bam", 
                    "in": [
                        {
                            "source": "#merge_duprem_filter.cwl/filter_by_mapq/bam_filtered", 
                            "id": "#merge_duprem_filter.cwl/sorting_filtered_bam/bam_unsorted"
                        }
                    ]
                }, 
                {
                    "doc": "samtools sort - sorting of merged bam", 
                    "out": [
                        "#merge_duprem_filter.cwl/sorting_merged_bam/bam_sorted"
                    ], 
                    "run": "#samtools_sort.cwl", 
                    "id": "#merge_duprem_filter.cwl/sorting_merged_bam", 
                    "in": [
                        {
                            "source": "#merge_duprem_filter.cwl/lane_replicate_merging/bam_merged", 
                            "id": "#merge_duprem_filter.cwl/sorting_merged_bam/bam_unsorted"
                        }
                    ]
                }
            ], 
            "id": "#merge_duprem_filter.cwl"
        }, 
        {
            "inputs": [
                {
                    "type": [
                        "string", 
                        "null"
                    ], 
                    "id": "#trim_and_map.cwl/adapter1"
                }, 
                {
                    "type": [
                        "string", 
                        "null"
                    ], 
                    "id": "#trim_and_map.cwl/adapter2"
                }, 
                {
                    "type": "File", 
                    "id": "#trim_and_map.cwl/fastq1"
                }, 
                {
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "id": "#trim_and_map.cwl/fastq2"
                }, 
                {
                    "type": "boolean", 
                    "id": "#trim_and_map.cwl/is_paired_end"
                }, 
                {
                    "secondaryFiles": [
                        ".fai", 
                        "^.1.bt2", 
                        "^.2.bt2", 
                        "^.3.bt2", 
                        "^.4.bt2", 
                        "^.rev.1.bt2", 
                        "^.rev.2.bt2"
                    ], 
                    "type": "File", 
                    "id": "#trim_and_map.cwl/reference"
                }
            ], 
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }, 
                {
                    "class": "ScatterFeatureRequirement"
                }, 
                {
                    "class": "StepInputExpressionRequirement"
                }, 
                {
                    "class": "SubworkflowFeatureRequirement"
                }
            ], 
            "outputs": [
                {
                    "outputSource": "#trim_and_map.cwl/sort_bam/bam_sorted", 
                    "type": "File", 
                    "id": "#trim_and_map.cwl/bam"
                }, 
                {
                    "outputSource": "#trim_and_map.cwl/mapping/bowtie2_log", 
                    "type": "File", 
                    "id": "#trim_and_map.cwl/bowtie2_log"
                }, 
                {
                    "outputSource": "#trim_and_map.cwl/adaptor_trimming_and_qc_post_trim/fastq1_trimmed", 
                    "type": "File", 
                    "id": "#trim_and_map.cwl/fastq1_trimmed"
                }, 
                {
                    "outputSource": "#trim_and_map.cwl/adaptor_trimming_and_qc_post_trim/fastq2_trimmed", 
                    "type": "File", 
                    "id": "#trim_and_map.cwl/fastq2_trimmed"
                }, 
                {
                    "outputSource": "#trim_and_map.cwl/adaptor_trimming_and_qc_post_trim/post_trim_fastqc_html", 
                    "type": {
                        "items": "File", 
                        "type": "array"
                    }, 
                    "id": "#trim_and_map.cwl/post_trim_fastqc_html"
                }, 
                {
                    "outputSource": "#trim_and_map.cwl/adaptor_trimming_and_qc_post_trim/post_trim_fastqc_zip", 
                    "type": {
                        "items": "File", 
                        "type": "array"
                    }, 
                    "id": "#trim_and_map.cwl/post_trim_fastqc_zip"
                }, 
                {
                    "outputSource": "#trim_and_map.cwl/qc_pre_trim/fastqc_html", 
                    "type": {
                        "items": "File", 
                        "type": "array"
                    }, 
                    "id": "#trim_and_map.cwl/pre_trim_fastqc_html"
                }, 
                {
                    "outputSource": "#trim_and_map.cwl/qc_pre_trim/fastqc_zip", 
                    "type": {
                        "items": "File", 
                        "type": "array"
                    }, 
                    "id": "#trim_and_map.cwl/pre_trim_fastqc_zip"
                }, 
                {
                    "outputSource": "#trim_and_map.cwl/adaptor_trimming_and_qc_post_trim/trim_galore_log", 
                    "type": {
                        "items": "File", 
                        "type": "array"
                    }, 
                    "id": "#trim_and_map.cwl/trim_galore_log"
                }
            ], 
            "class": "Workflow", 
            "steps": [
                {
                    "doc": "trim galore - adapter trimming using trim_galore", 
                    "out": [
                        "#trim_and_map.cwl/adaptor_trimming_and_qc_post_trim/fastq1_trimmed", 
                        "#trim_and_map.cwl/adaptor_trimming_and_qc_post_trim/fastq2_trimmed", 
                        "#trim_and_map.cwl/adaptor_trimming_and_qc_post_trim/fastq1_trimmed_unpaired", 
                        "#trim_and_map.cwl/adaptor_trimming_and_qc_post_trim/fastq2_trimmed_unpaired", 
                        "#trim_and_map.cwl/adaptor_trimming_and_qc_post_trim/trim_galore_log", 
                        "#trim_and_map.cwl/adaptor_trimming_and_qc_post_trim/post_trim_fastqc_html", 
                        "#trim_and_map.cwl/adaptor_trimming_and_qc_post_trim/post_trim_fastqc_zip"
                    ], 
                    "run": "#trim_galore.cwl", 
                    "id": "#trim_and_map.cwl/adaptor_trimming_and_qc_post_trim", 
                    "in": [
                        {
                            "source": "#trim_and_map.cwl/adapter1", 
                            "id": "#trim_and_map.cwl/adaptor_trimming_and_qc_post_trim/adapter1"
                        }, 
                        {
                            "source": "#trim_and_map.cwl/adapter2", 
                            "id": "#trim_and_map.cwl/adaptor_trimming_and_qc_post_trim/adapter2"
                        }, 
                        {
                            "source": "#trim_and_map.cwl/fastq1", 
                            "id": "#trim_and_map.cwl/adaptor_trimming_and_qc_post_trim/fastq1"
                        }, 
                        {
                            "source": "#trim_and_map.cwl/fastq2", 
                            "id": "#trim_and_map.cwl/adaptor_trimming_and_qc_post_trim/fastq2"
                        }
                    ]
                }, 
                {
                    "doc": "bowite2 - mapper, produces sam file", 
                    "out": [
                        "#trim_and_map.cwl/mapping/sam", 
                        "#trim_and_map.cwl/mapping/bowtie2_log"
                    ], 
                    "run": "#bowtie2.cwl", 
                    "id": "#trim_and_map.cwl/mapping", 
                    "in": [
                        {
                            "source": "#trim_and_map.cwl/adaptor_trimming_and_qc_post_trim/fastq1_trimmed", 
                            "id": "#trim_and_map.cwl/mapping/fastq1"
                        }, 
                        {
                            "source": "#trim_and_map.cwl/adaptor_trimming_and_qc_post_trim/fastq2_trimmed", 
                            "id": "#trim_and_map.cwl/mapping/fastq2"
                        }, 
                        {
                            "source": "#trim_and_map.cwl/is_paired_end", 
                            "id": "#trim_and_map.cwl/mapping/is_paired_end"
                        }, 
                        {
                            "source": "#trim_and_map.cwl/reference", 
                            "id": "#trim_and_map.cwl/mapping/reference_index"
                        }
                    ]
                }, 
                {
                    "doc": "fastqc - quality control for trimmed fastq", 
                    "out": [
                        "#trim_and_map.cwl/qc_pre_trim/fastqc_zip", 
                        "#trim_and_map.cwl/qc_pre_trim/fastqc_html"
                    ], 
                    "run": "#fastqc.cwl", 
                    "id": "#trim_and_map.cwl/qc_pre_trim", 
                    "in": [
                        {
                            "source": "#trim_and_map.cwl/fastq1", 
                            "id": "#trim_and_map.cwl/qc_pre_trim/fastq1"
                        }, 
                        {
                            "source": "#trim_and_map.cwl/fastq2", 
                            "id": "#trim_and_map.cwl/qc_pre_trim/fastq2"
                        }
                    ]
                }, 
                {
                    "doc": "samtools view - convert sam to bam", 
                    "out": [
                        "#trim_and_map.cwl/sam2bam/bam_unsorted"
                    ], 
                    "run": "#samtools_view_sam2bam.cwl", 
                    "id": "#trim_and_map.cwl/sam2bam", 
                    "in": [
                        {
                            "source": "#trim_and_map.cwl/mapping/sam", 
                            "id": "#trim_and_map.cwl/sam2bam/sam"
                        }
                    ]
                }, 
                {
                    "doc": "samtools sort - sorts unsorted bam file by coordinates.", 
                    "out": [
                        "#trim_and_map.cwl/sort_bam/bam_sorted"
                    ], 
                    "run": "#samtools_sort.cwl", 
                    "id": "#trim_and_map.cwl/sort_bam", 
                    "in": [
                        {
                            "source": "#trim_and_map.cwl/sam2bam/bam_unsorted", 
                            "id": "#trim_and_map.cwl/sort_bam/bam_unsorted"
                        }
                    ]
                }
            ], 
            "id": "#trim_and_map.cwl"
        }, 
        {
            "inputs": [
                {
                    "type": {
                        "items": [
                            "string", 
                            "null"
                        ], 
                        "type": "array"
                    }, 
                    "id": "#main/adapter1"
                }, 
                {
                    "type": {
                        "items": [
                            "string", 
                            "null"
                        ], 
                        "type": "array"
                    }, 
                    "id": "#main/adapter2"
                }, 
                {
                    "type": "long", 
                    "id": "#main/effective_genome_size"
                }, 
                {
                    "type": {
                        "items": [
                            "File", 
                            "null"
                        ], 
                        "type": "array"
                    }, 
                    "id": "#main/fastq1"
                }, 
                {
                    "type": {
                        "items": [
                            "File", 
                            "null"
                        ], 
                        "type": "array"
                    }, 
                    "id": "#main/fastq2"
                }, 
                {
                    "type": "int", 
                    "id": "#main/fragment_size"
                }, 
                {
                    "type": "boolean", 
                    "id": "#main/is_paired_end"
                }, 
                {
                    "secondaryFiles": [
                        ".fai", 
                        "^.1.bt2", 
                        "^.2.bt2", 
                        "^.3.bt2", 
                        "^.4.bt2", 
                        "^.rev.1.bt2", 
                        "^.rev.2.bt2"
                    ], 
                    "type": "File", 
                    "id": "#main/reference"
                }, 
                {
                    "type": "string", 
                    "id": "#main/sample_id"
                }
            ], 
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }, 
                {
                    "class": "MultipleInputFeatureRequirement"
                }, 
                {
                    "class": "ScatterFeatureRequirement"
                }, 
                {
                    "class": "StepInputExpressionRequirement"
                }, 
                {
                    "class": "SubworkflowFeatureRequirement"
                }
            ], 
            "outputs": [
                {
                    "outputSource": "#main/merge_duprem_filter/bam", 
                    "type": "File", 
                    "id": "#main/bam"
                }, 
                {
                    "outputSource": "#main/indexing_shifted_bam/bam_sorted_indexed", 
                    "type": "File", 
                    "id": "#main/bam_tn5_corrected"
                }, 
                {
                    "outputSource": "#main/generate_coverage_tracks/bigwig", 
                    "type": "File", 
                    "id": "#main/bigwig"
                }, 
                {
                    "outputSource": "#main/trim_and_map/bowtie2_log", 
                    "type": {
                        "items": "File", 
                        "type": "array"
                    }, 
                    "id": "#main/bowtie2_log"
                }, 
                {
                    "outputSource": "#main/create_summary_qc_report/multiqc_html", 
                    "type": "File", 
                    "id": "#main/multiqc_html"
                }, 
                {
                    "outputSource": "#main/create_summary_qc_report/multiqc_zip", 
                    "type": "File", 
                    "id": "#main/multiqc_zip"
                }, 
                {
                    "outputSource": "#main/merge_duprem_filter/post_filter_fastqc_html", 
                    "type": {
                        "items": "File", 
                        "type": "array"
                    }, 
                    "id": "#main/post_filter_fastqc_html"
                }, 
                {
                    "outputSource": "#main/merge_duprem_filter/post_filter_fastqc_zip", 
                    "type": {
                        "items": "File", 
                        "type": "array"
                    }, 
                    "id": "#main/post_filter_fastqc_zip"
                }, 
                {
                    "outputSource": "#main/trim_and_map/post_trim_fastqc_html", 
                    "type": {
                        "items": {
                            "items": "File", 
                            "type": "array"
                        }, 
                        "type": "array"
                    }, 
                    "id": "#main/post_trim_fastqc_html"
                }, 
                {
                    "outputSource": "#main/trim_and_map/post_trim_fastqc_zip", 
                    "type": {
                        "items": {
                            "items": "File", 
                            "type": "array"
                        }, 
                        "type": "array"
                    }, 
                    "id": "#main/post_trim_fastqc_zip"
                }, 
                {
                    "outputSource": "#main/trim_and_map/pre_trim_fastqc_html", 
                    "type": {
                        "items": {
                            "items": "File", 
                            "type": "array"
                        }, 
                        "type": "array"
                    }, 
                    "id": "#main/pre_trim_fastqc_html"
                }, 
                {
                    "outputSource": "#main/trim_and_map/pre_trim_fastqc_zip", 
                    "type": {
                        "items": {
                            "items": "File", 
                            "type": "array"
                        }, 
                        "type": "array"
                    }, 
                    "id": "#main/pre_trim_fastqc_zip"
                }, 
                {
                    "outputSource": "#main/trim_and_map/pre_trim_fastqc_zip", 
                    "type": {
                        "items": {
                            "items": "File", 
                            "type": "array"
                        }, 
                        "type": "array"
                    }, 
                    "id": "#main/trim_galore_log"
                }
            ], 
            "class": "Workflow", 
            "steps": [
                {
                    "doc": "multiqc summarizes the qc results from fastqc \nand other tools\n", 
                    "out": [
                        "#main/create_summary_qc_report/multiqc_zip", 
                        "#main/create_summary_qc_report/multiqc_html"
                    ], 
                    "run": "#multiqc_hack.cwl", 
                    "id": "#main/create_summary_qc_report", 
                    "in": [
                        {
                            "source": [
                                "#main/trim_and_map/bowtie2_log", 
                                "#main/merge_duprem_filter/post_filter_fastqc_zip", 
                                "#main/merge_duprem_filter/post_filter_fastqc_html"
                            ], 
                            "linkMerge": "merge_flattened", 
                            "id": "#main/create_summary_qc_report/qc_files_array"
                        }, 
                        {
                            "source": [
                                "#main/trim_and_map/pre_trim_fastqc_zip", 
                                "#main/trim_and_map/pre_trim_fastqc_html", 
                                "#main/trim_and_map/post_trim_fastqc_html", 
                                "#main/trim_and_map/post_trim_fastqc_zip", 
                                "#main/trim_and_map/trim_galore_log"
                            ], 
                            "linkMerge": "merge_flattened", 
                            "id": "#main/create_summary_qc_report/qc_files_array_of_array"
                        }, 
                        {
                            "source": "#main/sample_id", 
                            "id": "#main/create_summary_qc_report/report_name"
                        }
                    ]
                }, 
                {
                    "out": [
                        "#main/generate_coverage_tracks/bigwig"
                    ], 
                    "run": "#deeptools_bamCoverage.cwl", 
                    "id": "#main/generate_coverage_tracks", 
                    "in": [
                        {
                            "source": "#main/indexing_shifted_bam/bam_sorted_indexed", 
                            "id": "#main/generate_coverage_tracks/bam"
                        }, 
                        {
                            "source": "#main/effective_genome_size", 
                            "id": "#main/generate_coverage_tracks/effective_genome_size"
                        }, 
                        {
                            "source": "#main/fragment_size", 
                            "id": "#main/generate_coverage_tracks/fragment_size"
                        }
                    ]
                }, 
                {
                    "out": [
                        "#main/indexing_shifted_bam/bam_sorted_indexed"
                    ], 
                    "run": "#samtools_index_hack.cwl", 
                    "id": "#main/indexing_shifted_bam", 
                    "in": [
                        {
                            "source": "#main/tn5_overhang_correction/bam_tn5_corrected", 
                            "id": "#main/indexing_shifted_bam/bam_sorted"
                        }
                    ]
                }, 
                {
                    "out": [
                        "#main/merge_duprem_filter/post_filter_fastqc_zip", 
                        "#main/merge_duprem_filter/post_filter_fastqc_html", 
                        "#main/merge_duprem_filter/bam"
                    ], 
                    "run": "#merge_duprem_filter.cwl", 
                    "id": "#main/merge_duprem_filter", 
                    "in": [
                        {
                            "source": "#main/trim_and_map/bam", 
                            "id": "#main/merge_duprem_filter/bams"
                        }, 
                        {
                            "source": "#main/sample_id", 
                            "id": "#main/merge_duprem_filter/sample_id"
                        }
                    ]
                }, 
                {
                    "out": [
                        "#main/tn5_overhang_correction/bam_tn5_corrected"
                    ], 
                    "run": "#tn5_overhang_correction.cwl", 
                    "id": "#main/tn5_overhang_correction", 
                    "in": [
                        {
                            "source": "#main/merge_duprem_filter/bam", 
                            "id": "#main/tn5_overhang_correction/bam"
                        }, 
                        {
                            "source": "#main/is_paired_end", 
                            "id": "#main/tn5_overhang_correction/is_paired_end"
                        }
                    ]
                }, 
                {
                    "run": "#trim_and_map.cwl", 
                    "scatter": [
                        "#main/trim_and_map/fastq1", 
                        "#main/trim_and_map/fastq2", 
                        "#main/trim_and_map/adapter1", 
                        "#main/trim_and_map/adapter2"
                    ], 
                    "in": [
                        {
                            "source": "#main/adapter1", 
                            "id": "#main/trim_and_map/adapter1"
                        }, 
                        {
                            "source": "#main/adapter2", 
                            "id": "#main/trim_and_map/adapter2"
                        }, 
                        {
                            "source": "#main/fastq1", 
                            "id": "#main/trim_and_map/fastq1"
                        }, 
                        {
                            "source": "#main/fastq2", 
                            "id": "#main/trim_and_map/fastq2"
                        }, 
                        {
                            "source": "#main/is_paired_end", 
                            "id": "#main/trim_and_map/is_paired_end"
                        }, 
                        {
                            "source": "#main/reference", 
                            "id": "#main/trim_and_map/reference"
                        }
                    ], 
                    "scatterMethod": "dotproduct", 
                    "id": "#main/trim_and_map", 
                    "out": [
                        "#main/trim_and_map/pre_trim_fastqc_zip", 
                        "#main/trim_and_map/pre_trim_fastqc_html", 
                        "#main/trim_and_map/fastq1_trimmed", 
                        "#main/trim_and_map/fastq2_trimmed", 
                        "#main/trim_and_map/trim_galore_log", 
                        "#main/trim_and_map/post_trim_fastqc_html", 
                        "#main/trim_and_map/post_trim_fastqc_zip", 
                        "#main/trim_and_map/bam", 
                        "#main/trim_and_map/bowtie2_log"
                    ]
                }
            ], 
            "id": "#main"
        }
    ]
}