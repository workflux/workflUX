{
    "cwlVersion": "v1.0", 
    "$graph": [
        {
            "inputs": [
                {
                    "type": "File", 
                    "id": "#bedgraph_sort.cwl/bedgraph"
                }, 
                {
                    "doc": "optional, if not specified, output file will be named as input file", 
                    "type": [
                        "null", 
                        "string"
                    ], 
                    "id": "#bedgraph_sort.cwl/output_name"
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
            "stdout": "${\n  if( inputs.output_name == null ){\n    return inputs.bedgraph.basename;\n  }\n  else{\n    return inputs.output_name;\n  }\n}\n", 
            "outputs": [
                {
                    "type": "stdout", 
                    "id": "#bedgraph_sort.cwl/bedgraph_sorted"
                }
            ], 
            "baseCommand": [
                "bash", 
                "-c"
            ], 
            "id": "#bedgraph_sort.cwl", 
            "arguments": [
                {
                    "valueFrom": "$(\"LC_COLLATE=C sort -k1,1 -k2,2n \" + inputs.bedgraph.path)"
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
                    "ramMin": 15000, 
                    "class": "ResourceRequirement"
                }
            ]
        }, 
        {
            "inputs": [
                {
                    "inputBinding": {
                        "position": 10, 
                        "prefix": "-i"
                    }, 
                    "type": "File", 
                    "id": "#bedtools_bamtobed_pe.cwl/bam"
                }
            ], 
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ], 
            "stdout": "$(inputs.bam.nameroot + \".bedpe\")", 
            "outputs": [
                {
                    "type": "stdout", 
                    "id": "#bedtools_bamtobed_pe.cwl/bedpe"
                }
            ], 
            "baseCommand": [
                "bedtools", 
                "bamtobed"
            ], 
            "id": "#bedtools_bamtobed_pe.cwl", 
            "arguments": [
                {
                    "position": 1, 
                    "valueFrom": "-bedpe"
                }
            ], 
            "class": "CommandLineTool", 
            "hints": [
                {
                    "dockerPull": "biocontainers/bedtools:2.25.0", 
                    "class": "DockerRequirement"
                }, 
                {
                    "coresMin": 1, 
                    "ramMin": 15000, 
                    "class": "ResourceRequirement"
                }
            ]
        }, 
        {
            "inputs": [
                {
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "-i"
                    }, 
                    "type": "File", 
                    "id": "#bedtools_bedtobam.cwl/bed"
                }, 
                {
                    "inputBinding": {
                        "position": 2, 
                        "prefix": "-g"
                    }, 
                    "type": "File", 
                    "id": "#bedtools_bedtobam.cwl/reference_info"
                }
            ], 
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ], 
            "stdout": "$(inputs.bed.nameroot + \".bam\")", 
            "outputs": [
                {
                    "type": "stdout", 
                    "id": "#bedtools_bedtobam.cwl/bam"
                }
            ], 
            "baseCommand": [
                "bedtools", 
                "bedtobam"
            ], 
            "id": "#bedtools_bedtobam.cwl", 
            "class": "CommandLineTool", 
            "hints": [
                {
                    "dockerPull": "biocontainers/bedtools:2.25.0", 
                    "class": "DockerRequirement"
                }, 
                {
                    "coresMin": 1, 
                    "ramMin": 15000, 
                    "class": "ResourceRequirement"
                }
            ]
        }, 
        {
            "inputs": [
                {
                    "inputBinding": {
                        "position": 2, 
                        "prefix": "-i"
                    }, 
                    "type": "File", 
                    "id": "#bedtools_genomecov.cwl/bed"
                }, 
                {
                    "inputBinding": {
                        "position": 3, 
                        "prefix": "-g"
                    }, 
                    "type": "File", 
                    "id": "#bedtools_genomecov.cwl/reference_info"
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
            "stdout": "$(inputs.bed.nameroot + \".bedgraph\")", 
            "outputs": [
                {
                    "type": "stdout", 
                    "id": "#bedtools_genomecov.cwl/bedgraph"
                }
            ], 
            "baseCommand": [
                "bedtools", 
                "genomecov"
            ], 
            "id": "#bedtools_genomecov.cwl", 
            "arguments": [
                {
                    "position": 1, 
                    "valueFrom": "-bg"
                }
            ], 
            "class": "CommandLineTool", 
            "hints": [
                {
                    "dockerPull": "biocontainers/bedtools:2.25.0", 
                    "class": "DockerRequirement"
                }, 
                {
                    "coresMin": 1, 
                    "ramMin": 15000, 
                    "tmpdirMin": 10000, 
                    "class": "ResourceRequirement"
                }
            ]
        }, 
        {
            "inputs": [
                {
                    "inputBinding": {
                        "position": 2, 
                        "prefix": "-i"
                    }, 
                    "type": "File", 
                    "id": "#bedtools_slop_clip.cwl/bed"
                }, 
                {
                    "inputBinding": {
                        "position": 3, 
                        "prefix": "-g"
                    }, 
                    "type": "File", 
                    "id": "#bedtools_slop_clip.cwl/reference_info"
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
            "stdout": "$(inputs.bed.nameroot + \"_clipped.bed\")", 
            "outputs": [
                {
                    "type": "stdout", 
                    "id": "#bedtools_slop_clip.cwl/bed_clipped"
                }
            ], 
            "baseCommand": [
                "bedtools", 
                "slop"
            ], 
            "id": "#bedtools_slop_clip.cwl", 
            "arguments": [
                {
                    "position": 1, 
                    "valueFrom": "0", 
                    "prefix": "-b"
                }
            ], 
            "class": "CommandLineTool", 
            "hints": [
                {
                    "dockerPull": "biocontainers/bedtools:2.25.0", 
                    "class": "DockerRequirement"
                }, 
                {
                    "coresMin": 1, 
                    "ramMin": 15000, 
                    "tmpdirMin": 10000, 
                    "class": "ResourceRequirement"
                }
            ]
        }, 
        {
            "inputs": [
                {
                    "inputBinding": {
                        "position": 11
                    }, 
                    "type": "File", 
                    "id": "#bowtie2.cwl/fastq1"
                }, 
                {
                    "inputBinding": {
                        "position": 12, 
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
                    "doc": "usefull for very long fragments, as expected for ATAC", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "--maxins"
                    }, 
                    "type": [
                        "null", 
                        "long"
                    ], 
                    "id": "#bowtie2.cwl/max_mapping_insert_length"
                }, 
                {
                    "type": "string", 
                    "id": "#bowtie2.cwl/output_basename"
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
                    "position": 10, 
                    "valueFrom": "${\n  if ( inputs.is_paired_end ){\n     return \"-1\";\n  }\n  else {\n    return \"-U\";\n  }\n}\n"
                }, 
                {
                    "position": 6, 
                    "valueFrom": "$(inputs.output_basename + \".sam\")", 
                    "prefix": "-S"
                }
            ], 
            "stderr": "$( inputs.output_basename + \".bowtie2_stderr\")", 
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
                    "doc": "must be indexed", 
                    "inputBinding": {
                        "position": 100, 
                        "prefix": "--bamfiles"
                    }, 
                    "type": "File", 
                    "id": "#deeptools_plotCoverage.cwl/bam", 
                    "secondaryFiles": ".bai"
                }, 
                {
                    "type": [
                        "null", 
                        "int"
                    ], 
                    "id": "#deeptools_plotCoverage.cwl/fragment_size"
                }, 
                {
                    "default": true, 
                    "doc": "if paired end, reads are extended", 
                    "type": "boolean", 
                    "id": "#deeptools_plotCoverage.cwl/is_paired_end"
                }, 
                {
                    "type": "string", 
                    "id": "#deeptools_plotCoverage.cwl/sample_id"
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
                        "glob": "$(inputs.sample_id + \".plot_cov.png\")"
                    }, 
                    "type": "File", 
                    "id": "#deeptools_plotCoverage.cwl/qc_plot_coverage_plot"
                }, 
                {
                    "outputBinding": {
                        "glob": "$(inputs.sample_id + \".plot_cov.tsv\")"
                    }, 
                    "type": "File", 
                    "id": "#deeptools_plotCoverage.cwl/qc_plot_coverage_tsv"
                }
            ], 
            "baseCommand": [
                "plotCoverage"
            ], 
            "id": "#deeptools_plotCoverage.cwl", 
            "arguments": [
                {
                    "position": 1, 
                    "valueFrom": "${\n  if ( inputs.is_paired_end ){\n     return null;\n  }\n  else {\n    return \"--extendReads\";\n  }\n}\n"
                }, 
                {
                    "position": 2, 
                    "valueFrom": "${\n  if ( inputs.is_paired_end ){\n     return null;\n  }\n  else {\n    return inputs.fragment_size;\n  }\n}\n"
                }, 
                {
                    "position": 10, 
                    "valueFrom": "$(inputs.sample_id)", 
                    "prefix": "--labels"
                }, 
                {
                    "position": 10, 
                    "valueFrom": "$(inputs.sample_id + \".plot_cov.png\")", 
                    "prefix": "--plotFile"
                }, 
                {
                    "position": 10, 
                    "valueFrom": "$(inputs.sample_id + \".plot_cov.tsv\")", 
                    "prefix": "--outRawCounts"
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
                    "ramMin": 15000, 
                    "class": "ResourceRequirement"
                }
            ]
        }, 
        {
            "inputs": [
                {
                    "doc": "must be indexed", 
                    "inputBinding": {
                        "position": 100, 
                        "prefix": "--bamfiles"
                    }, 
                    "type": "File", 
                    "id": "#deeptools_plotFingerprint.cwl/bam", 
                    "secondaryFiles": ".bai"
                }, 
                {
                    "type": [
                        "null", 
                        "int"
                    ], 
                    "id": "#deeptools_plotFingerprint.cwl/fragment_size"
                }, 
                {
                    "default": true, 
                    "doc": "if paired end, reads are extended", 
                    "type": "boolean", 
                    "id": "#deeptools_plotFingerprint.cwl/is_paired_end"
                }, 
                {
                    "type": "string", 
                    "id": "#deeptools_plotFingerprint.cwl/sample_id"
                }
            ], 
            "permanentFailCodes": [], 
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }, 
                {
                    "class": "StepInputExpressionRequirement"
                }
            ], 
            "successCodes": [
                0, 
                1, 
                2
            ], 
            "outputs": [
                {
                    "outputBinding": {
                        "glob": "$(inputs.sample_id + \".plot_fingerp.png\")"
                    }, 
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "id": "#deeptools_plotFingerprint.cwl/qc_plot_fingerprint_plot"
                }, 
                {
                    "type": "stderr", 
                    "id": "#deeptools_plotFingerprint.cwl/qc_plot_fingerprint_stderr"
                }, 
                {
                    "outputBinding": {
                        "glob": "$(inputs.sample_id + \".plot_fingerp.tsv\")"
                    }, 
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "id": "#deeptools_plotFingerprint.cwl/qc_plot_fingerprint_tsv"
                }
            ], 
            "baseCommand": [
                "plotFingerprint"
            ], 
            "id": "#deeptools_plotFingerprint.cwl", 
            "arguments": [
                {
                    "position": 1, 
                    "valueFrom": "${\n  if ( inputs.is_paired_end ){\n     return null;\n  }\n  else {\n    return \"--extendReads\";\n  }\n}\n"
                }, 
                {
                    "position": 2, 
                    "valueFrom": "${\n  if ( inputs.is_paired_end ){\n     return null;\n  }\n  else {\n    return inputs.fragment_size;\n  }\n}\n"
                }, 
                {
                    "position": 10, 
                    "valueFrom": "$(inputs.sample_id)", 
                    "prefix": "--labels"
                }, 
                {
                    "position": 10, 
                    "valueFrom": "$(inputs.sample_id + \".plot_fingerp.png\")", 
                    "prefix": "--plotFile"
                }, 
                {
                    "position": 10, 
                    "valueFrom": "$(inputs.sample_id + \".plot_fingerp.tsv\")", 
                    "prefix": "--outRawCounts"
                }
            ], 
            "stderr": "$( inputs.sample_id + \".plot_fingerp.stderr\")", 
            "temporaryFailCodes": [], 
            "class": "CommandLineTool", 
            "hints": [
                {
                    "dockerPull": "kerstenbreuer/deeptools:3.1.1", 
                    "class": "DockerRequirement"
                }, 
                {
                    "coresMin": 1, 
                    "ramMin": 15000, 
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
                    "type": "File", 
                    "id": "#generate_atac_signal_tags.cwl/bedpe_alignm"
                }, 
                {
                    "type": "string", 
                    "id": "#generate_atac_signal_tags.cwl/output_basename"
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
                        "glob": "$(inputs.output_basename + \"_fragments_tn5_incl_tags.bed\")"
                    }, 
                    "type": "File", 
                    "id": "#generate_atac_signal_tags.cwl/bed_fragments_tn5_incl_signal"
                }, 
                {
                    "outputBinding": {
                        "glob": "$(inputs.output_basename + \"_pot_nucl_bound_tags.bed\")"
                    }, 
                    "type": "File", 
                    "id": "#generate_atac_signal_tags.cwl/bed_nucl_bound_signal"
                }, 
                {
                    "outputBinding": {
                        "glob": "$(inputs.output_basename + \"_nucl_free_tags.bed\")"
                    }, 
                    "type": "File", 
                    "id": "#generate_atac_signal_tags.cwl/bed_nucl_free_signal"
                }, 
                {
                    "outputBinding": {
                        "glob": "$(inputs.output_basename + \"_tn5_bind_region.bed\")"
                    }, 
                    "type": "File", 
                    "id": "#generate_atac_signal_tags.cwl/bed_tn5_bind_region_signal"
                }, 
                {
                    "outputBinding": {
                        "glob": "$(inputs.output_basename + \"_tn5_center_1bp.bed\")"
                    }, 
                    "type": "File", 
                    "id": "#generate_atac_signal_tags.cwl/bed_tn5_center_1bp_signal"
                }, 
                {
                    "outputBinding": {
                        "glob": "$(inputs.output_basename + \"_filtering_stats_mqc.tsv\")"
                    }, 
                    "type": "File", 
                    "id": "#generate_atac_signal_tags.cwl/filtering_stats_tsv"
                }, 
                {
                    "outputBinding": {
                        "glob": "$(inputs.output_basename + \"_frag_size_classification_mqc.tsv\")"
                    }, 
                    "type": "File", 
                    "id": "#generate_atac_signal_tags.cwl/frag_size_stats_tsv"
                }, 
                {
                    "outputBinding": {
                        "glob": "$(inputs.output_basename + \"_fragment_sizes.txt\")"
                    }, 
                    "type": "File", 
                    "id": "#generate_atac_signal_tags.cwl/fragment_sizes_tsv"
                }, 
                {
                    "outputBinding": {
                        "glob": "$(inputs.output_basename + \"_irreg_mappings.bedpe\")"
                    }, 
                    "type": "File", 
                    "id": "#generate_atac_signal_tags.cwl/irreg_mappings_bedpe"
                }
            ], 
            "baseCommand": [
                "bash", 
                "-c"
            ], 
            "id": "#generate_atac_signal_tags.cwl", 
            "arguments": [
                {
                    "valueFrom": "${\n  var cmd_line = \"touch \\\"\" + inputs.output_basename + \"_irreg_mappings.bedpe\\\" \" +\n      \"\\\"\" + inputs.output_basename + \"_tn5_bind_region_unsorted.bed\\\" \" +\n      \"\\\"\" + inputs.output_basename + \"_fragments_tn5_incl_tags_unsorted.bed\\\" \" +\n      \"\\\"\" + inputs.output_basename + \"_tn5_center_1bp_unsorted.bed\\\" \" +\n      \"\\\"\" + inputs.output_basename + \"_pot_nucl_bound_tags_unsorted.bed\\\" \" +\n      \"\\\"\" + inputs.output_basename + \"_nucl_free_tags_unsorted.bed\\\" \" +\n      \"\\\"\" + inputs.output_basename + \"_irreg_mappings.bedpe\\\"; \" +\n      \" awk \\'\" +\n              \"BEGIN {\" +\n                  \"OFS=\\\"\\\\t\\\";\" +\n                  \"chrM_read_count=0;\" +\n                  \"interchrom_map_read_count=0;\" +\n                  \"regular_read_count=0;\" +\n                  \"irregular_read_count=0;\" +\n                  \"too_small_fragment_count=0;\" +\n                  \"nucl_free_fragment_count=0;\" +\n                  \"nucl_bound_fragment_count=0;\" +\n                  \"wrong_strand_orient_count=0;\" +\n              \"}\" +\n              \"{\" +\n              \"\tif ( $1==\\\"chrM\\\" || $4==\\\"chrM\\\") {\" +\n                      \"irregular_read_count += 2;\" +\n                      \"if ( $1==$4 ) {\" +\n                          \"chrM_read_count += 2;\" +\n                          \"print $0, \\\"chrM\\\" >\\\"\" + inputs.output_basename + \"_irreg_mappings.bedpe\\\";\" +\n                      \"}\" +\n                      \"else {\" +\n                          \"chrM_read_count += 1;\" +\n                          \"interchrom_map_read_count += 2;\" +\n                          \"print $0, \\\"interchrom_map\\\" >\\\"\" + inputs.output_basename + \"_irreg_mappings.bedpe\\\";\" +\n                      \"}\" +\n                  \"}\" +\n                  \"else if ( $1!=$4 ) {\" +\n                      \"irregular_read_count += 2;\" +\n                      \"interchrom_map_read_count += 2;\" +\n                      \"print $0, \\\"interchrom_map\\\" >\\\"\" + inputs.output_basename + \"_irreg_mappings.bedpe\\\";\" +\n                  \"}\" +\n                  \"else if ( $9==$10 ) {\" +\n                      \"irregular_read_count += 2;\" +\n                      \"wrong_strand_orient_count += 2;\" +\n                      \"print $0, \\\"wrong_read_pair_orientation\\\" >\\\"\" + inputs.output_basename + \"_irreg_mappings.bedpe\\\";\" +\n                  \"}\" +\n                  \"else {\" +\n                      \"right_orientation=0;\" +\n                      \"if ( $9==\\\"+\\\" && $2<=$5 && $3<=$6 ){\" +\n                          \"right_orientation=1;\" +\n                          \"fragment_size=$6-$2;\" +\n                          \"start=$2;\" +\n                          \"end=$6;\" +\n                          \"print $1, $2+4, $2+5, $7, $8 > \\\"\" + inputs.output_basename + \"_tn5_center_1bp_unsorted.bed\\\";\" +\n                          \"print $1, $6-5, $6-4, ($7 \\\"(mate)\\\"), $8 > \\\"\" + inputs.output_basename + \"_tn5_center_1bp_unsorted.bed\\\";\" +\n                      \"}\" +\n                      \"else if ( $9==\\\"-\\\" && $2>=$5 && $3>=$6 ){\" +\n                          \"right_orientation=1;\" +\n                          \"fragment_size=$3-$5;\" +\n                          \"start=$5;\" +\n                          \"end=$3;\" +\n                          \"print $1, $5+4, $5+5, $7, $8 > \\\"\" + inputs.output_basename + \"_tn5_center_1bp_unsorted.bed\\\";\" +\n                          \"print $1, $3-5, $3-4, ($7 \\\"(mate)\\\"), $8 > \\\"\" + inputs.output_basename + \"_tn5_center_1bp_unsorted.bed\\\";\" +\n                      \"}\" +\n                      \"else {\" +\n                          \"wrong_strand_orient_count += 2;\" +\n                          \"print $0, \\\"wrong_read_pair_orientation\\\" >\\\"\" + inputs.output_basename + \"_irreg_mappings.bedpe\\\";\" +\n                      \"}\" +\n                      \"print fragment_size > \\\"\" + inputs.output_basename + \"_fragment_sizes.txt\\\";\" +\n                      \"if ( fragment_size>=38 && right_orientation ) {\" +\n                          \"regular_read_count += 2;\" +\n                          \"tag_start = start-10;\" +\n                          \"tag_end = start+19;\" +\n                          \"if ( tag_start < 0) { tag_start = 0 }\" +\n                          \"print $1, tag_start, tag_end, $7, $8 > \\\"\" + inputs.output_basename + \"_tn5_bind_region_unsorted.bed\\\";\" +\n                          \"tag_start = end-19;\" +\n                          \"tag_end = end+10;\" +\n                          \"if ( tag_start < 0) { tag_start = 0 }\" +\n                          \"print $1, tag_start, tag_end, ($7 \\\"(mate)\\\"), $8 > \\\"\" + inputs.output_basename + \"_tn5_bind_region_unsorted.bed\\\";\" +\n                          \"tag_start = start-10;\" +\n                          \"tag_end = end+10;\" +\n                          \"if ( tag_start < 0) { tag_start = 0 }\" +\n                          \"print $1, tag_start, tag_end, $7, $8 > \\\"\" + inputs.output_basename + \"_fragments_tn5_incl_tags_unsorted.bed\\\";\" +\n                          \"if (fragment_size>=185) { \" +\n                              \"nucl_bound_fragment_count += 1;\" +\n                              \"tag_start = start+19;\" +\n                              \"tag_end = end-19;\" +\n                              \"if ( tag_end < 0) { tag_end = 0 }\" +\n                              \"print $1, tag_start, tag_end, $7, $8 > \\\"\" + inputs.output_basename + \"_pot_nucl_bound_tags_unsorted.bed\\\";\" +\n                              \"tag_start = start-10;\" +\n                              \"tag_end = start+19;\" +\n                              \"if ( tag_start < 0) { tag_start = 0 }\" +\n                              \"print $1, tag_start, tag_end, $7, $8 > \\\"\" + inputs.output_basename + \"_nucl_free_tags_unsorted.bed\\\";\" +\n                              \"tag_start = end-19;\" +\n                              \"tag_end = end+10;\" +\n                              \"if ( tag_start < 0) { tag_start = 0 }\" +\n                              \"print $1, tag_start, tag_end, ($7 \\\"(mate)\\\"), $8 > \\\"\" + inputs.output_basename + \"_nucl_free_tags_unsorted.bed\\\";\" +\n                          \"}\" +\n                          \"else {\" +\n                              \"nucl_free_fragment_count += 1;\" +\n                              \"tag_start = start-10;\" +\n                              \"tag_end = end+10;\" +\n                              \"if ( tag_start < 0) { tag_start = 0 }\" +\n                              \"print $1, tag_start, tag_end, $7, $8 > \\\"\" + inputs.output_basename + \"_nucl_free_tags_unsorted.bed\\\";\" +\n                          \"}\" +\n                      \"}\" +\n                      \"else {\" +\n                          \"if ( right_orientation ){\" +\n                              \"irregular_read_count += 2;\" +\n                              \"too_small_fragment_count += 1;\" +\n                              \"print $0, \\\"too_small_frag_size\\\" >\\\"\" + inputs.output_basename + \"_irreg_mappings.bedpe\\\";\" +\n                          \"}\" +\n                      \"}\" +\n              \"\t}\" +\n              \"}\" +\n              \"END {\" +\n                  \"print \\\"\\# id: \\\\\\\"Filtering Statistics\\\\\\\"\\\" > \\\"\" + inputs.output_basename + \"_filtering_stats_mqc.tsv\\\";\" +\n                  \"print \\\"\\# description: \\\\\\\"- This section shows statistics on read filtering: (1) reads pairs that are mapping to different chromosomes as well as \" + \n                      \"(2) read pairs located on ChrM are filtered out; (3) read pairs that have a wrong orientation towards each other \" +\n                      \"(e.g. both reads on same strand, or reads pointing to different direction) are removed, too. \" +\n                      \"(4) Only the remaining regular reads are used for the fragment size analysis and the generation of atac signal tracks. \" +\n                      \"In addition to the filtering shown here, reads were also selected for high mapping quality and to be mapped in a proper pair.\" +\n                      \"\\\\\\\"\\\" > \\\"\" + inputs.output_basename + \"_filtering_stats_mqc.tsv\\\";\" +\n                  \"print \\\"\\# plot_type: \\\\\\\"bargraph\\\\\\\"\\\" > \\\"\" + inputs.output_basename + \"_filtering_stats_mqc.tsv\\\";\" +\n                  \"print \\\"chrM_reads\\\", chrM_read_count > \\\"\" + inputs.output_basename + \"_filtering_stats_mqc.tsv\\\";\" +\n                  \"print \\\"interchrom_map_reads\\\", interchrom_map_read_count > \\\"\" + inputs.output_basename + \"_filtering_stats_mqc.tsv\\\";\" +\n                  \"print \\\"wrong_read_pair_orientation\\\", wrong_strand_orient_count > \\\"\" + inputs.output_basename + \"_filtering_stats_mqc.tsv\\\";\" +\n                  \"print \\\"regular_reads\\\", regular_read_count > \\\"\" + inputs.output_basename + \"_filtering_stats_mqc.tsv\\\";\" +\n                  \"print \\\"\\# id: \\\\\\\"Fragment Length Classification\\\\\\\"\\\" > \\\"\" + inputs.output_basename + \"_frag_size_classification_mqc.tsv\\\";\" +\n                  \"print \\\"\\# description: \\\\\\\"Fragments are classified by their size: (1) nucleosome free fragements are smaller than a typical a nucleosome binding region while \"+\n                      \"(2) potentially nucleosome bound fragments are larger;(3) fragments that are classified as too small are shorter than expected for a Tn5 digestion. \" +\n                      \"For these calculations, the DNA span that is covered by the Tn5 enzyme during the trasposition is taken into account.\" + \n                      \"\\\\\\\"\\\" > \\\"\" + inputs.output_basename + \"_frag_size_classification_mqc.tsv\\\";\" +\n                  \"print \\\"\\# plot_type: \\\\\\\"bargraph\\\\\\\"\\\" > \\\"\" + inputs.output_basename + \"_frag_size_classification_mqc.tsv\\\";\" +\n                  \"print \\\"too_small_fragments\\\", too_small_fragment_count > \\\"\" + inputs.output_basename + \"_frag_size_classification_mqc.tsv\\\";\" +\n                  \"print \\\"nucl_free_fragments\\\", nucl_free_fragment_count > \\\"\" + inputs.output_basename + \"_frag_size_classification_mqc.tsv\\\";\" +\n                  \"print \\\"pot_nucl_bound_fragments\\\", nucl_bound_fragment_count > \\\"\" + inputs.output_basename + \"_frag_size_classification_mqc.tsv\\\";\" +\n              \"}\" +\n              \"\\' \" + inputs.bedpe_alignm.path;\n  cmd_line += \";LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n \" + inputs.output_basename + \"_tn5_bind_region_unsorted.bed > \" + inputs.output_basename + \"_tn5_bind_region.bed \" +\n              \";LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n \" + inputs.output_basename + \"_tn5_center_1bp_unsorted.bed > \" + inputs.output_basename + \"_tn5_center_1bp.bed \" +\n              \";LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n \" + inputs.output_basename + \"_nucl_free_tags_unsorted.bed > \" + inputs.output_basename + \"_nucl_free_tags.bed \" +\n              \";LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n \" + inputs.output_basename + \"_pot_nucl_bound_tags_unsorted.bed > \" + inputs.output_basename + \"_pot_nucl_bound_tags.bed \" +\n              \";LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n \" + inputs.output_basename + \"_fragments_tn5_incl_tags_unsorted.bed > \" + inputs.output_basename + \"_fragments_tn5_incl_tags.bed \" ;\n\n  return cmd_line;\n}\n  \n"
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
                    "ramMin": 15000, 
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
                    "type": "File", 
                    "id": "#kentutils_bedGraphToBigWig.cwl/bedgraph_sorted"
                }, 
                {
                    "inputBinding": {
                        "position": 2
                    }, 
                    "type": "File", 
                    "id": "#kentutils_bedGraphToBigWig.cwl/reference_info"
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
                        "glob": "$(inputs.bedgraph_sorted.nameroot + \".bigwig\")"
                    }, 
                    "type": "File", 
                    "id": "#kentutils_bedGraphToBigWig.cwl/bigwig"
                }
            ], 
            "baseCommand": [
                "bedGraphToBigWig"
            ], 
            "id": "#kentutils_bedGraphToBigWig.cwl", 
            "arguments": [
                {
                    "position": 3, 
                    "valueFrom": "$(inputs.bedgraph_sorted.nameroot + \".bigwig\")"
                }
            ], 
            "class": "CommandLineTool", 
            "hints": [
                {
                    "dockerPull": "kerstenbreuer/kent-bedgraphtobigwig:latest", 
                    "class": "DockerRequirement"
                }, 
                {
                    "coresMin": 1, 
                    "ramMin": 15000, 
                    "class": "ResourceRequirement"
                }
            ]
        }, 
        {
            "inputs": [
                {
                    "doc": "can be \"mm\", \"hs\", \"ce\", \"dm\", or the total number of genomic bp", 
                    "inputBinding": {
                        "position": 3, 
                        "prefix": "--gsize"
                    }, 
                    "type": "string", 
                    "id": "#macs2_callpeak_atac_tn5_binding_region.cwl/genome_size"
                }, 
                {
                    "doc": "gives the base name for all output files", 
                    "inputBinding": {
                        "position": 101, 
                        "prefix": "--name"
                    }, 
                    "type": "string", 
                    "id": "#macs2_callpeak_atac_tn5_binding_region.cwl/output_basename"
                }, 
                {
                    "inputBinding": {
                        "position": 100, 
                        "prefix": "--treatment"
                    }, 
                    "type": "File", 
                    "id": "#macs2_callpeak_atac_tn5_binding_region.cwl/treatment_bed"
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
                        "glob": "*.broadPeak"
                    }, 
                    "type": "File", 
                    "id": "#macs2_callpeak_atac_tn5_binding_region.cwl/broad_peaks_bed"
                }, 
                {
                    "outputBinding": {
                        "glob": "*.gappedPeak"
                    }, 
                    "type": "File", 
                    "id": "#macs2_callpeak_atac_tn5_binding_region.cwl/gapped_peaks_bed"
                }, 
                {
                    "outputBinding": {
                        "glob": "*_peaks.xls"
                    }, 
                    "type": "File", 
                    "id": "#macs2_callpeak_atac_tn5_binding_region.cwl/peaks_xls"
                }, 
                {
                    "outputBinding": {
                        "glob": "*_treat_pileup.bdg"
                    }, 
                    "type": "File", 
                    "id": "#macs2_callpeak_atac_tn5_binding_region.cwl/treat_pileup_bdg"
                }
            ], 
            "baseCommand": [
                "macs2", 
                "callpeak"
            ], 
            "id": "#macs2_callpeak_atac_tn5_binding_region.cwl", 
            "arguments": [
                {
                    "position": 1, 
                    "valueFrom": "BED", 
                    "prefix": "--format"
                }, 
                {
                    "position": 2, 
                    "valueFrom": "--broad"
                }, 
                {
                    "position": 2, 
                    "valueFrom": "--nomodel"
                }, 
                {
                    "position": 2, 
                    "valueFrom": "all", 
                    "prefix": "--keep-dup"
                }, 
                {
                    "position": 2, 
                    "valueFrom": "0", 
                    "prefix": "--shift"
                }, 
                {
                    "position": 2, 
                    "valueFrom": "29", 
                    "prefix": "--extsize"
                }, 
                {
                    "position": 2, 
                    "valueFrom": "--bdg"
                }
            ], 
            "class": "CommandLineTool", 
            "hints": [
                {
                    "dockerPull": "genomicpariscentre/macs2:2.1.0.20140616", 
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
                    "doc": "can be \"mm\", \"hs\", \"ce\", \"dm\", or the total number of genomic bp", 
                    "inputBinding": {
                        "position": 3, 
                        "prefix": "--gsize"
                    }, 
                    "type": "string", 
                    "id": "#macs2_callpeak_atac_tn5_center.cwl/genome_size"
                }, 
                {
                    "doc": "gives the base name for all output files", 
                    "inputBinding": {
                        "position": 101, 
                        "prefix": "--name"
                    }, 
                    "type": "string", 
                    "id": "#macs2_callpeak_atac_tn5_center.cwl/output_basename"
                }, 
                {
                    "inputBinding": {
                        "position": 100, 
                        "prefix": "--treatment"
                    }, 
                    "type": "File", 
                    "id": "#macs2_callpeak_atac_tn5_center.cwl/treatment_bed"
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
                        "glob": "*.broadPeak"
                    }, 
                    "type": "File", 
                    "id": "#macs2_callpeak_atac_tn5_center.cwl/broad_peaks_bed"
                }, 
                {
                    "outputBinding": {
                        "glob": "*.gappedPeak"
                    }, 
                    "type": "File", 
                    "id": "#macs2_callpeak_atac_tn5_center.cwl/gapped_peaks_bed"
                }, 
                {
                    "outputBinding": {
                        "glob": "*_peaks.xls"
                    }, 
                    "type": "File", 
                    "id": "#macs2_callpeak_atac_tn5_center.cwl/peaks_xls"
                }, 
                {
                    "outputBinding": {
                        "glob": "*_treat_pileup.bdg"
                    }, 
                    "type": "File", 
                    "id": "#macs2_callpeak_atac_tn5_center.cwl/treat_pileup_bdg"
                }
            ], 
            "baseCommand": [
                "macs2", 
                "callpeak"
            ], 
            "id": "#macs2_callpeak_atac_tn5_center.cwl", 
            "arguments": [
                {
                    "position": 1, 
                    "valueFrom": "BED", 
                    "prefix": "--format"
                }, 
                {
                    "position": 2, 
                    "valueFrom": "-37", 
                    "prefix": "--shift"
                }, 
                {
                    "position": 2, 
                    "valueFrom": "73", 
                    "prefix": "--extsize"
                }, 
                {
                    "position": 2, 
                    "valueFrom": "--broad"
                }, 
                {
                    "position": 2, 
                    "valueFrom": "--nomodel"
                }, 
                {
                    "position": 2, 
                    "valueFrom": "all", 
                    "prefix": "--keep-dup"
                }, 
                {
                    "position": 2, 
                    "valueFrom": "--bdg"
                }
            ], 
            "class": "CommandLineTool", 
            "hints": [
                {
                    "dockerPull": "genomicpariscentre/macs2:2.1.0.20140616", 
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
                    "doc": "qc files which shall be part of the multiqc summary;\noptional, only one of qc_files_array or qc_files_array_of_array \nmust be provided\n", 
                    "type": [
                        "null", 
                        {
                            "items": [
                                "File", 
                                "null"
                            ], 
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
                                "items": [
                                    "File", 
                                    "null"
                                ], 
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
                    "valueFrom": "${\n    var qc_files_array = inputs.qc_files_array;\n    var qc_files_array_of_array = inputs.qc_files_array_of_array;\n    var cmdline = \"echo 'copying input file ...'\";\n\n    if ( qc_files_array != null ){\n      for (var i=0; i<qc_files_array.length; i++){\n        if( qc_files_array[i] != null ){\n          cmdline += \"; cp \" + qc_files_array[i].path + \" .\";\n        }\n      }\n    }\n\n    if ( qc_files_array_of_array != null ){\n      for (var i=0; i<qc_files_array_of_array.length; i++){ \n        for (var ii=0; ii<qc_files_array_of_array[i].length; ii++){\n          if( qc_files_array_of_array[i][ii] != null ){\n            cmdline += \"; cp \" + qc_files_array_of_array[i][ii].path + \" .\";\n          }\n        }\n      }\n    }\n    \n    cmdline += \"; echo \\'copying done\\'\" +\n        \"; multiqc --zip-data-dir --cl_config \\'log_filesize_limit: 100000000\\' \" +\n        \"--outdir \" + runtime.outdir +\n        \" --filename \" + inputs.report_name + \"_report .\";\n\n    return cmdline\n  }\n"
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
                    "doc": "aligned and filtered reads, not shifted", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "--bam"
                    }, 
                    "type": "File", 
                    "id": "#nucleoatac.cwl/bam"
                }, 
                {
                    "doc": "regions of open chromatin (ATAC peaks)", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "--bed"
                    }, 
                    "type": "File", 
                    "id": "#nucleoatac.cwl/bed"
                }, 
                {
                    "doc": "reference genome in fasta format + samtools faidx index", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "--fasta"
                    }, 
                    "type": "File", 
                    "id": "#nucleoatac.cwl/fasta", 
                    "secondaryFiles": [
                        ".fai"
                    ]
                }, 
                {
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "--out"
                    }, 
                    "type": "string", 
                    "id": "#nucleoatac.cwl/output_basename"
                }
            ], 
            "permanentFailCodes": [], 
            "successCodes": [
                0, 
                1, 
                2
            ], 
            "stdout": "$( inputs.output_basename + \".stout\" )", 
            "outputs": [
                {
                    "outputBinding": {
                        "glob": "*.nucmap_combined.bed.gz"
                    }, 
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "id": "#nucleoatac.cwl/combined_nucl_pos_bed"
                }, 
                {
                    "outputBinding": {
                        "glob": "*.fragmentsizes.txt"
                    }, 
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "id": "#nucleoatac.cwl/fragsize_in_peaks_txt"
                }, 
                {
                    "outputBinding": {
                        "glob": "*.nfrpos.bed.gz"
                    }, 
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "id": "#nucleoatac.cwl/nfr_pos_bed"
                }, 
                {
                    "outputBinding": {
                        "glob": "*.nuc_dist.eps"
                    }, 
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "id": "#nucleoatac.cwl/nucl_dist_plot"
                }, 
                {
                    "outputBinding": {
                        "glob": "*.nuc_dist.txt"
                    }, 
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "id": "#nucleoatac.cwl/nucl_dist_txt"
                }, 
                {
                    "outputBinding": {
                        "glob": "*.nucleoatac_signal.bedgraph.gz"
                    }, 
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "id": "#nucleoatac.cwl/nucl_norm_crosscor_tracks"
                }, 
                {
                    "outputBinding": {
                        "glob": "*.nucleoatac_signal.smooth.bedgraph.gz"
                    }, 
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "id": "#nucleoatac.cwl/nucl_norm_smooth_crosscor_tracks"
                }, 
                {
                    "outputBinding": {
                        "glob": "*.occ_fit.eps"
                    }, 
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "id": "#nucleoatac.cwl/nucl_occ_fit_plot"
                }, 
                {
                    "outputBinding": {
                        "glob": "*.occ_fit.txt"
                    }, 
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "id": "#nucleoatac.cwl/nucl_occ_fit_txt"
                }, 
                {
                    "outputBinding": {
                        "glob": "*.occ.lower_bound.bedgraph.gz"
                    }, 
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "id": "#nucleoatac.cwl/nucl_occ_lower_bound_tracks"
                }, 
                {
                    "outputBinding": {
                        "glob": "*.occpeaks.bed.gz"
                    }, 
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "id": "#nucleoatac.cwl/nucl_occ_peaks_bed"
                }, 
                {
                    "outputBinding": {
                        "glob": "*.occ.bedgraph.gz"
                    }, 
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "id": "#nucleoatac.cwl/nucl_occ_tracks"
                }, 
                {
                    "outputBinding": {
                        "glob": "*.occ.upper_bound.bedgraph.gz"
                    }, 
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "id": "#nucleoatac.cwl/nucl_occ_upper_bound_tracks"
                }, 
                {
                    "outputBinding": {
                        "glob": "*.nucpos.bed.gz"
                    }, 
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "id": "#nucleoatac.cwl/nucl_pos_bed"
                }, 
                {
                    "outputBinding": {
                        "glob": "*.nucpos.redundant.bed.gz"
                    }, 
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "id": "#nucleoatac.cwl/nucl_pos_redundant_bed"
                }, 
                {
                    "outputBinding": {
                        "glob": "*.VMat"
                    }, 
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "id": "#nucleoatac.cwl/nucl_vplot_data"
                }, 
                {
                    "type": "stderr", 
                    "id": "#nucleoatac.cwl/nucleoatac_stderr"
                }, 
                {
                    "type": "stdout", 
                    "id": "#nucleoatac.cwl/nucleoatac_stdout"
                }
            ], 
            "baseCommand": [
                "nucleoatac", 
                "run"
            ], 
            "id": "#nucleoatac.cwl", 
            "stderr": "$( inputs.output_basename + \".stderr\" )", 
            "temporaryFailCodes": [], 
            "class": "CommandLineTool", 
            "hints": [
                {
                    "dockerPull": "kerstenbreuer/nucleoatac:0.3.4", 
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
                    "inputBinding": {
                        "position": 10, 
                        "prefix": "-c=", 
                        "separate": false
                    }, 
                    "type": "File", 
                    "id": "#phantompeakqualtools.cwl/bam"
                }
            ], 
            "permanentFailCodes": [], 
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }, 
                {
                    "class": "StepInputExpressionRequirement"
                }
            ], 
            "successCodes": [
                0, 
                1, 
                2
            ], 
            "outputs": [
                {
                    "outputBinding": {
                        "glob": "*.pdf"
                    }, 
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "id": "#phantompeakqualtools.cwl/qc_crosscorr_plot"
                }, 
                {
                    "outputBinding": {
                        "glob": "*.spp.out"
                    }, 
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "id": "#phantompeakqualtools.cwl/qc_crosscorr_summary"
                }, 
                {
                    "type": "stderr", 
                    "id": "#phantompeakqualtools.cwl/qc_phantompeakqualtools_stderr"
                }, 
                {
                    "type": "stdout", 
                    "id": "#phantompeakqualtools.cwl/qc_phantompeakqualtools_stdout"
                }
            ], 
            "baseCommand": [
                "Rscript", 
                "--verbose", 
                "--max-ppsize=500000", 
                "/usr/bin/phantompeakqualtools-1.2/run_spp.R"
            ], 
            "id": "#phantompeakqualtools.cwl", 
            "arguments": [
                {
                    "position": 10, 
                    "valueFrom": "$(runtime.tmpdir)", 
                    "prefix": "-tmpdir=", 
                    "separate": false
                }, 
                {
                    "position": 10, 
                    "valueFrom": "$(runtime.outdir)", 
                    "prefix": "-odir=", 
                    "separate": false
                }, 
                {
                    "position": 100, 
                    "valueFrom": "$(inputs.bam.nameroot + \".crosscor.pdf\")", 
                    "prefix": "-savp=", 
                    "separate": false
                }, 
                {
                    "position": 100, 
                    "valueFrom": "$(inputs.bam.nameroot + \".spp.out\")", 
                    "prefix": "-out=", 
                    "separate": false
                }
            ], 
            "stderr": "$(inputs.bam.nameroot + \".phantompeakqualtools_stdout\")", 
            "temporaryFailCodes": [], 
            "class": "CommandLineTool", 
            "hints": [
                {
                    "dockerPull": "kerstenbreuer/phantompeakqualtools:1.2", 
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
            "stdout": "$(inputs.bam_sorted.nameroot + \".picard_markdup.stdout\")", 
            "outputs": [
                {
                    "outputBinding": {
                        "glob": "$(inputs.bam_sorted.nameroot + \"_duprem.bam\")"
                    }, 
                    "type": "File", 
                    "id": "#picard_markdup.cwl/bam_duprem"
                }, 
                {
                    "type": "stdout", 
                    "id": "#picard_markdup.cwl/picard_markdup_stdout"
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
                    "type": "File", 
                    "id": "#plot_frag_size_distr.cwl/fragment_sizes_tsv"
                }, 
                {
                    "type": "string", 
                    "id": "#plot_frag_size_distr.cwl/output_basename"
                }
            ], 
            "outputs": [
                {
                    "outputBinding": {
                        "glob": "$(inputs.output_basename + \"_frag_size_distr.png\")"
                    }, 
                    "type": "File", 
                    "id": "#plot_frag_size_distr.cwl/frag_size_distr_plot"
                }, 
                {
                    "outputBinding": {
                        "glob": "$(inputs.output_basename + \"_fragm_sizes_mqc.tsv\")"
                    }, 
                    "type": "File", 
                    "id": "#plot_frag_size_distr.cwl/frag_size_distr_tsv"
                }
            ], 
            "baseCommand": [
                "Rscript", 
                "-e"
            ], 
            "id": "#plot_frag_size_distr.cwl", 
            "arguments": [
                {
                    "valueFrom": "${\n  var r_cmd = \"library(tidyverse); \" +\n            \"library(gridExtra); \" +\n            \"size_freqs <- read.table('\" + \n                inputs.fragment_sizes_tsv.path +\n                  \"', header=F)[,1] %>% as.integer() %>% table(); \" +\n            \"size_freq_tb <- tibble(frag_size=1:max(as.integer(names(size_freqs)))) %>% \" +\n              \"full_join( tibble( frag_size=as.integer(names(size_freqs)), freq=size_freqs), by='frag_size') %>% \" +\n              \"mutate(freq=sapply(freq,function(f) ifelse(is.na(f), 0, f))); \" +\n            \"plot_ <- ggplot(size_freq_tb) + \" +\n              \"geom_line(mapping=aes(x=frag_size,y=freq)) + \" +\n              \"ggtitle('Fragment Size Distribution') + \" +\n              \"xlab('') + ylab('frequency'); \" +\n            \"plot_log <- ggplot(size_freq_tb) + \" +\n              \"geom_line(mapping=aes(x=frag_size,y=freq)) + \" +\n              \"scale_y_log10() + \" +\n              \"xlab('fragment size [bp]') + ylab('frequency'); \" +\n            \"png( '\" + inputs.output_basename + \"_frag_size_distr.png', \" +\n                \"width = 850, height = 850); \" +\n            \"grid.arrange(plot_, plot_log, nrow=2);\" +\n            \"dev.off();\"+\n            \"mqc_file_headers <- c(\"+\n              \"'# id: \\\\\\\"Fragment Size Distribution\\\\\\\"', \" +\n              \"'# description: \\\\\\\"Distribution of fragment sizes (<1000) which typically shows the DNA pitch and multimers of nucleosomes. \" + \n                              \"For a better visualization, which also includes higher fragment sizes, please see the files ending with _frag_size_distr.png\\\\\\\"', \" +\n              \"'# pconfig:', \" +\n              \"'#    xmax: 1000'\" +\n              \");\" +\n            \"mqc_file <- file('\" + inputs.output_basename + \"_fragm_sizes_mqc.tsv');\" +\n            \"writeLines(mqc_file_headers, mqc_file);\" +\n            \"close(mqc_file);\" +\n            \"write.table(size_freq_tb, file=\\'\" + inputs.output_basename + \"_fragm_sizes_mqc.tsv\\', \" +\n              \"sep='\\\\t', row.names=F, col.names=F, append=T);\" +\n            \"print('done')\"\n  return r_cmd\n}\n"
                }
            ], 
            "class": "CommandLineTool", 
            "hints": [
                {
                    "dockerPull": "kerstenbreuer/tidyverse:1.2.1", 
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
                    "doc": "aligned reads to be checked in sam or bam format", 
                    "inputBinding": {
                        "position": 2
                    }, 
                    "type": "File", 
                    "id": "#samtools_sort_name.cwl/bam_unsorted"
                }
            ], 
            "stdout": "$(inputs.bam_unsorted.basename)", 
            "doc": "Sort a bam file by read names.", 
            "baseCommand": [
                "samtools", 
                "sort"
            ], 
            "id": "#samtools_sort_name.cwl", 
            "arguments": [
                {
                    "position": 1, 
                    "valueFrom": "$(runtime.cores)", 
                    "prefix": "-@"
                }, 
                {
                    "position": 1, 
                    "valueFrom": "-n"
                }
            ], 
            "outputs": [
                {
                    "type": "stdout", 
                    "id": "#samtools_sort_name.cwl/bam_sorted"
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
                        "position": 10
                    }, 
                    "type": "File", 
                    "id": "#samtools_view_filter.cwl/bam"
                }, 
                {
                    "default": true, 
                    "doc": "if paired end, only properly paired reads pass", 
                    "type": "boolean", 
                    "id": "#samtools_view_filter.cwl/is_paired_end"
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
                }, 
                {
                    "position": 2, 
                    "valueFrom": "${\n  if ( inputs.is_paired_end ){\n     return \"-f\";\n  }\n  else {\n    return null;\n  }\n}\n"
                }, 
                {
                    "position": 3, 
                    "valueFrom": "${\n  if ( inputs.is_paired_end ){\n     return \"3\";\n  }\n  else {\n    return null;\n  }\n}\n"
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
            "id": "#bed_to_coverage_track.cwl", 
            "inputs": [
                {
                    "type": "File", 
                    "id": "#bed_to_coverage_track.cwl/bed"
                }, 
                {
                    "type": "File", 
                    "id": "#bed_to_coverage_track.cwl/reference_info"
                }
            ], 
            "steps": [
                {
                    "doc": "bedtools bedtobam - converts bed to bam;\nas most tools can handle bam as input but not always\nthe bed format(e.g. deeptools), moreover, bam is compressed;\ntherefore, it will be used as final output (instead of the bed file)\n", 
                    "out": [
                        "#bed_to_coverage_track.cwl/converting_bed_to_bam/bam"
                    ], 
                    "run": "#bedtools_bedtobam.cwl", 
                    "id": "#bed_to_coverage_track.cwl/converting_bed_to_bam", 
                    "in": [
                        {
                            "source": "#bed_to_coverage_track.cwl/bed", 
                            "id": "#bed_to_coverage_track.cwl/converting_bed_to_bam/bed"
                        }, 
                        {
                            "source": "#bed_to_coverage_track.cwl/reference_info", 
                            "id": "#bed_to_coverage_track.cwl/converting_bed_to_bam/reference_info"
                        }
                    ]
                }, 
                {
                    "doc": "bedtools genomeCov -bg", 
                    "out": [
                        "#bed_to_coverage_track.cwl/converting_bed_to_bedgraph/bedgraph"
                    ], 
                    "run": "#bedtools_genomecov.cwl", 
                    "id": "#bed_to_coverage_track.cwl/converting_bed_to_bedgraph", 
                    "in": [
                        {
                            "source": "#bed_to_coverage_track.cwl/bed", 
                            "id": "#bed_to_coverage_track.cwl/converting_bed_to_bedgraph/bed"
                        }, 
                        {
                            "source": "#bed_to_coverage_track.cwl/reference_info", 
                            "id": "#bed_to_coverage_track.cwl/converting_bed_to_bedgraph/reference_info"
                        }
                    ]
                }, 
                {
                    "doc": "bedGraphToBigWig (kentUtils)", 
                    "out": [
                        "#bed_to_coverage_track.cwl/converting_bedgraph_to_bigwig/bigwig"
                    ], 
                    "run": "#kentutils_bedGraphToBigWig.cwl", 
                    "id": "#bed_to_coverage_track.cwl/converting_bedgraph_to_bigwig", 
                    "in": [
                        {
                            "source": "#bed_to_coverage_track.cwl/sorting_bedgraph/bedgraph_sorted", 
                            "id": "#bed_to_coverage_track.cwl/converting_bedgraph_to_bigwig/bedgraph_sorted"
                        }, 
                        {
                            "source": "#bed_to_coverage_track.cwl/reference_info", 
                            "id": "#bed_to_coverage_track.cwl/converting_bedgraph_to_bigwig/reference_info"
                        }
                    ]
                }, 
                {
                    "doc": "samtools index - indexes sorted bam\n", 
                    "out": [
                        "#bed_to_coverage_track.cwl/indexing_bam/bam_sorted_indexed"
                    ], 
                    "run": "#samtools_index_hack.cwl", 
                    "id": "#bed_to_coverage_track.cwl/indexing_bam", 
                    "in": [
                        {
                            "source": "#bed_to_coverage_track.cwl/sorting_bam/bam_sorted", 
                            "id": "#bed_to_coverage_track.cwl/indexing_bam/bam_sorted"
                        }
                    ]
                }, 
                {
                    "doc": "samtools sort - sorting of merged bam", 
                    "out": [
                        "#bed_to_coverage_track.cwl/sorting_bam/bam_sorted"
                    ], 
                    "run": "#samtools_sort.cwl", 
                    "id": "#bed_to_coverage_track.cwl/sorting_bam", 
                    "in": [
                        {
                            "source": "#bed_to_coverage_track.cwl/converting_bed_to_bam/bam", 
                            "id": "#bed_to_coverage_track.cwl/sorting_bam/bam_unsorted"
                        }
                    ]
                }, 
                {
                    "doc": "LC_COLLATE=C sort -k1,1 -k2,2n", 
                    "out": [
                        "#bed_to_coverage_track.cwl/sorting_bedgraph/bedgraph_sorted"
                    ], 
                    "run": "#bedgraph_sort.cwl", 
                    "id": "#bed_to_coverage_track.cwl/sorting_bedgraph", 
                    "in": [
                        {
                            "source": "#bed_to_coverage_track.cwl/converting_bed_to_bedgraph/bedgraph", 
                            "id": "#bed_to_coverage_track.cwl/sorting_bedgraph/bedgraph"
                        }
                    ]
                }
            ], 
            "class": "Workflow", 
            "outputs": [
                {
                    "outputSource": "#bed_to_coverage_track.cwl/indexing_bam/bam_sorted_indexed", 
                    "type": "File", 
                    "id": "#bed_to_coverage_track.cwl/bam"
                }, 
                {
                    "outputSource": "#bed_to_coverage_track.cwl/converting_bedgraph_to_bigwig/bigwig", 
                    "type": "File", 
                    "id": "#bed_to_coverage_track.cwl/bigwig"
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
                    "type": "boolean", 
                    "id": "#merge_duprem_filter.cwl/is_paired_end"
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
                    "outputSource": "#merge_duprem_filter.cwl/remove_duplicates/picard_markdup_stdout", 
                    "type": "File", 
                    "id": "#merge_duprem_filter.cwl/picard_markdup_stdout"
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
                        }, 
                        {
                            "source": "#merge_duprem_filter.cwl/is_paired_end", 
                            "id": "#merge_duprem_filter.cwl/filter_by_mapq/is_paired_end"
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
                        "#merge_duprem_filter.cwl/remove_duplicates/bam_duprem", 
                        "#merge_duprem_filter.cwl/remove_duplicates/picard_markdup_stdout"
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
                    "type": [
                        "null", 
                        "long"
                    ], 
                    "id": "#trim_and_map.cwl/max_mapping_insert_length"
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
                }, 
                {
                    "type": "string", 
                    "id": "#trim_and_map.cwl/sample_id"
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
                            "source": "#trim_and_map.cwl/max_mapping_insert_length", 
                            "id": "#trim_and_map.cwl/mapping/max_mapping_insert_length"
                        }, 
                        {
                            "source": "#trim_and_map.cwl/sample_id", 
                            "id": "#trim_and_map.cwl/mapping/output_basename"
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
                    "default": 73, 
                    "type": "int", 
                    "id": "#main/extend_bps_downstream"
                }, 
                {
                    "type": {
                        "items": "File", 
                        "type": "array"
                    }, 
                    "id": "#main/fastq1"
                }, 
                {
                    "type": {
                        "items": "File", 
                        "type": "array"
                    }, 
                    "id": "#main/fastq2"
                }, 
                {
                    "type": "string", 
                    "id": "#main/macs2_genome_size"
                }, 
                {
                    "default": 2500, 
                    "type": "long", 
                    "id": "#main/max_mapping_insert_length"
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
                    "type": "File", 
                    "id": "#main/reference_info"
                }, 
                {
                    "type": "string", 
                    "id": "#main/sample_id"
                }, 
                {
                    "default": 37, 
                    "type": "int", 
                    "id": "#main/shift_bps_upstream"
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
                    "secondaryFiles": ".bai", 
                    "type": "File", 
                    "outputSource": "#main/merge_duprem_filter/bam", 
                    "id": "#main/bam"
                }, 
                {
                    "outputSource": "#main/generating_tn5_bind_region_signal_tracks/bam", 
                    "type": "File", 
                    "id": "#main/bam_tn5_bind_region_signal"
                }, 
                {
                    "outputSource": "#main/generating_tn5_center_1bp_signal_tracks/bam", 
                    "type": "File", 
                    "id": "#main/bam_tn5_center_1bp_signal"
                }, 
                {
                    "outputSource": "#main/generating_fragments_tn5_excl_signal_tracks/bigwig", 
                    "type": "File", 
                    "id": "#main/bigwig_fragments_tn5_excl_signal"
                }, 
                {
                    "outputSource": "#main/generating_nucl_bound_signal_tracks/bigwig", 
                    "type": "File", 
                    "id": "#main/bigwig_nucl_bound_signal"
                }, 
                {
                    "outputSource": "#main/generating_nucl_free_signal_tracks/bigwig", 
                    "type": "File", 
                    "id": "#main/bigwig_nucl_free_signal"
                }, 
                {
                    "outputSource": "#main/generating_tn5_bind_region_signal_tracks/bigwig", 
                    "type": "File", 
                    "id": "#main/bigwig_tn5_bind_region_signal"
                }, 
                {
                    "outputSource": "#main/generating_tn5_center_1bp_signal_tracks/bigwig", 
                    "type": "File", 
                    "id": "#main/bigwig_tn5_center_1bp_signal"
                }, 
                {
                    "outputSource": "#main/converting_shift_ext_bedgraph_to_bigwig/bigwig", 
                    "type": "File", 
                    "id": "#main/bigwig_tn5_center_shift_ext_signal"
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
                    "outputSource": "#main/peak_calling_macs2_tn5_bind_region/broad_peaks_bed", 
                    "type": "File", 
                    "id": "#main/broad_peaks_bed_tn5_bind_region_signal"
                }, 
                {
                    "outputSource": "#main/peak_calling_macs2_tn5_center/broad_peaks_bed", 
                    "type": "File", 
                    "id": "#main/broad_peaks_bed_tn5_center_shift_ext_signal"
                }, 
                {
                    "outputSource": "#main/nucl_position_calling/combined_nucl_pos_bed", 
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "id": "#main/combined_nucl_pos_bed"
                }, 
                {
                    "outputSource": "#main/qc_plot_coverage_fragments_tn5_excl/qc_plot_coverage_tsv", 
                    "type": "File", 
                    "id": "#main/coverage_counts_fragments_tn5_excl"
                }, 
                {
                    "outputSource": "#main/qc_plot_coverage_fragments_tn5_excl/qc_plot_coverage_plot", 
                    "type": "File", 
                    "id": "#main/coverage_plot_fragments_tn5_excl"
                }, 
                {
                    "outputSource": "#main/generating_atac_signal_tags/filtering_stats_tsv", 
                    "type": "File", 
                    "id": "#main/filtering_stats_tsv"
                }, 
                {
                    "outputSource": "#main/plot_fragment_size_distribution/frag_size_distr_plot", 
                    "type": "File", 
                    "id": "#main/frag_size_distr_plot"
                }, 
                {
                    "outputSource": "#main/plot_fragment_size_distribution/frag_size_distr_tsv", 
                    "type": "File", 
                    "id": "#main/frag_size_distr_tsv"
                }, 
                {
                    "outputSource": "#main/generating_atac_signal_tags/frag_size_stats_tsv", 
                    "type": "File", 
                    "id": "#main/frag_size_stats_tsv"
                }, 
                {
                    "outputSource": "#main/generating_atac_signal_tags/fragment_sizes_tsv", 
                    "type": "File", 
                    "id": "#main/fragment_sizes_tsv"
                }, 
                {
                    "outputSource": "#main/nucl_position_calling/fragsize_in_peaks_txt", 
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "id": "#main/fragsize_in_peaks_txt"
                }, 
                {
                    "outputSource": "#main/peak_calling_macs2_tn5_bind_region/gapped_peaks_bed", 
                    "type": "File", 
                    "id": "#main/gapped_peaks_bed_tn5_bind_region_signal"
                }, 
                {
                    "outputSource": "#main/peak_calling_macs2_tn5_center/gapped_peaks_bed", 
                    "type": "File", 
                    "id": "#main/gapped_peaks_bed_tn5_center_shift_ext_signal"
                }, 
                {
                    "outputSource": "#main/generating_atac_signal_tags/irreg_mappings_bedpe", 
                    "type": "File", 
                    "id": "#main/irreg_mappings_bedpe"
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
                    "outputSource": "#main/nucl_position_calling/nfr_pos_bed", 
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "id": "#main/nfr_pos_bed"
                }, 
                {
                    "outputSource": "#main/nucl_position_calling/nucl_dist_plot", 
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "id": "#main/nucl_dist_plot"
                }, 
                {
                    "outputSource": "#main/nucl_position_calling/nucl_dist_txt", 
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "id": "#main/nucl_dist_txt"
                }, 
                {
                    "outputSource": "#main/nucl_position_calling/nucl_norm_crosscor_tracks", 
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "id": "#main/nucl_norm_crosscor_tracks"
                }, 
                {
                    "outputSource": "#main/nucl_position_calling/nucl_norm_smooth_crosscor_tracks", 
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "id": "#main/nucl_norm_smooth_crosscor_tracks"
                }, 
                {
                    "outputSource": "#main/nucl_position_calling/nucl_occ_fit_plot", 
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "id": "#main/nucl_occ_fit_plot"
                }, 
                {
                    "outputSource": "#main/nucl_position_calling/nucl_occ_fit_txt", 
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "id": "#main/nucl_occ_fit_txt"
                }, 
                {
                    "outputSource": "#main/nucl_position_calling/nucl_occ_lower_bound_tracks", 
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "id": "#main/nucl_occ_lower_bound_tracks"
                }, 
                {
                    "outputSource": "#main/nucl_position_calling/nucl_occ_peaks_bed", 
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "id": "#main/nucl_occ_peaks_bed"
                }, 
                {
                    "outputSource": "#main/nucl_position_calling/nucl_occ_tracks", 
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "id": "#main/nucl_occ_tracks"
                }, 
                {
                    "outputSource": "#main/nucl_position_calling/nucl_occ_upper_bound_tracks", 
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "id": "#main/nucl_occ_upper_bound_tracks"
                }, 
                {
                    "outputSource": "#main/nucl_position_calling/nucl_pos_bed", 
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "id": "#main/nucl_pos_bed"
                }, 
                {
                    "outputSource": "#main/nucl_position_calling/nucl_pos_redundant_bed", 
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "id": "#main/nucl_pos_redundant_bed"
                }, 
                {
                    "outputSource": "#main/nucl_position_calling/nucl_vplot_data", 
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "id": "#main/nucl_vplot_data"
                }, 
                {
                    "outputSource": "#main/nucl_position_calling/nucleoatac_stderr", 
                    "type": "File", 
                    "id": "#main/nucleoatac_stderr"
                }, 
                {
                    "outputSource": "#main/nucl_position_calling/nucleoatac_stdout", 
                    "type": "File", 
                    "id": "#main/nucleoatac_stdout"
                }, 
                {
                    "outputSource": "#main/peak_calling_macs2_tn5_bind_region/peaks_xls", 
                    "type": "File", 
                    "id": "#main/peaks_xls_tn5_bind_region_signal"
                }, 
                {
                    "outputSource": "#main/peak_calling_macs2_tn5_center/peaks_xls", 
                    "type": "File", 
                    "id": "#main/peaks_xls_tn5_center_shift_ext_signal"
                }, 
                {
                    "outputSource": "#main/merge_duprem_filter/picard_markdup_stdout", 
                    "type": "File", 
                    "id": "#main/picard_markdup_stdout"
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
                    "outputSource": "#main/qc_phantompeakqualtools/qc_crosscorr_plot", 
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "id": "#main/qc_crosscorr_plot"
                }, 
                {
                    "outputSource": "#main/qc_phantompeakqualtools/qc_crosscorr_summary", 
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "id": "#main/qc_crosscorr_summary"
                }, 
                {
                    "outputSource": "#main/qc_phantompeakqualtools/qc_phantompeakqualtools_stderr", 
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "id": "#main/qc_phantompeakqualtools_stderr"
                }, 
                {
                    "outputSource": "#main/qc_plot_fingerprint/qc_plot_fingerprint_plot", 
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "id": "#main/qc_plot_fingerprint_plot"
                }, 
                {
                    "outputSource": "#main/qc_plot_fingerprint/qc_plot_fingerprint_stderr", 
                    "type": "File", 
                    "id": "#main/qc_plot_fingerprint_stderr"
                }, 
                {
                    "outputSource": "#main/qc_plot_fingerprint/qc_plot_fingerprint_tsv", 
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "id": "#main/qc_plot_fingerprint_tsv"
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
                    "doc": "clips features exceeding the chromosome boundaries", 
                    "out": [
                        "#main/clip_shift_ext_bedgraph/bed_clipped"
                    ], 
                    "run": "#bedtools_slop_clip.cwl", 
                    "id": "#main/clip_shift_ext_bedgraph", 
                    "in": [
                        {
                            "source": "#main/peak_calling_macs2_tn5_center/treat_pileup_bdg", 
                            "id": "#main/clip_shift_ext_bedgraph/bed"
                        }, 
                        {
                            "source": "#main/reference_info", 
                            "id": "#main/clip_shift_ext_bedgraph/reference_info"
                        }
                    ]
                }, 
                {
                    "doc": "bedtools bamtobed", 
                    "out": [
                        "#main/converting_bam_to_bedpe/bedpe"
                    ], 
                    "run": "#bedtools_bamtobed_pe.cwl", 
                    "id": "#main/converting_bam_to_bedpe", 
                    "in": [
                        {
                            "source": "#main/name_sorting_filtered_bam/bam_sorted", 
                            "id": "#main/converting_bam_to_bedpe/bam"
                        }
                    ]
                }, 
                {
                    "doc": "bedGraphToBigWig (kentUtils)", 
                    "out": [
                        "#main/converting_shift_ext_bedgraph_to_bigwig/bigwig"
                    ], 
                    "run": "#kentutils_bedGraphToBigWig.cwl", 
                    "id": "#main/converting_shift_ext_bedgraph_to_bigwig", 
                    "in": [
                        {
                            "source": "#main/sorting_shift_ext_bedgraph/bedgraph_sorted", 
                            "id": "#main/converting_shift_ext_bedgraph_to_bigwig/bedgraph_sorted"
                        }, 
                        {
                            "source": "#main/reference_info", 
                            "id": "#main/converting_shift_ext_bedgraph_to_bigwig/reference_info"
                        }
                    ]
                }, 
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
                                "#main/merge_duprem_filter/post_filter_fastqc_html", 
                                "#main/generating_atac_signal_tags/frag_size_stats_tsv", 
                                "#main/generating_atac_signal_tags/fragment_sizes_tsv", 
                                "#main/generating_atac_signal_tags/filtering_stats_tsv", 
                                "#main/plot_fragment_size_distribution/frag_size_distr_tsv", 
                                "#main/peak_calling_macs2_tn5_bind_region/peaks_xls", 
                                "#main/peak_calling_macs2_tn5_center/peaks_xls", 
                                "#main/qc_plot_coverage_fragments_tn5_excl/qc_plot_coverage_tsv", 
                                "#main/qc_plot_fingerprint/qc_plot_fingerprint_tsv", 
                                "#main/qc_phantompeakqualtools/qc_phantompeakqualtools_stdout", 
                                "#main/qc_phantompeakqualtools/qc_crosscorr_summary", 
                                "#main/merge_duprem_filter/picard_markdup_stdout"
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
                    "doc": null, 
                    "out": [
                        "#main/generating_atac_signal_tags/bed_tn5_bind_region_signal", 
                        "#main/generating_atac_signal_tags/bed_tn5_center_1bp_signal", 
                        "#main/generating_atac_signal_tags/bed_nucl_free_signal", 
                        "#main/generating_atac_signal_tags/bed_nucl_bound_signal", 
                        "#main/generating_atac_signal_tags/bed_fragments_tn5_incl_signal", 
                        "#main/generating_atac_signal_tags/fragment_sizes_tsv", 
                        "#main/generating_atac_signal_tags/filtering_stats_tsv", 
                        "#main/generating_atac_signal_tags/frag_size_stats_tsv", 
                        "#main/generating_atac_signal_tags/irreg_mappings_bedpe"
                    ], 
                    "run": "#generate_atac_signal_tags.cwl", 
                    "id": "#main/generating_atac_signal_tags", 
                    "in": [
                        {
                            "source": "#main/converting_bam_to_bedpe/bedpe", 
                            "id": "#main/generating_atac_signal_tags/bedpe_alignm"
                        }, 
                        {
                            "source": "#main/sample_id", 
                            "id": "#main/generating_atac_signal_tags/output_basename"
                        }
                    ]
                }, 
                {
                    "doc": null, 
                    "out": [
                        "#main/generating_fragments_tn5_excl_signal_tracks/bigwig", 
                        "#main/generating_fragments_tn5_excl_signal_tracks/bam"
                    ], 
                    "run": "#bed_to_coverage_track.cwl", 
                    "id": "#main/generating_fragments_tn5_excl_signal_tracks", 
                    "in": [
                        {
                            "source": "#main/generating_atac_signal_tags/bed_fragments_tn5_incl_signal", 
                            "id": "#main/generating_fragments_tn5_excl_signal_tracks/bed"
                        }, 
                        {
                            "source": "#main/reference_info", 
                            "id": "#main/generating_fragments_tn5_excl_signal_tracks/reference_info"
                        }
                    ]
                }, 
                {
                    "doc": null, 
                    "out": [
                        "#main/generating_nucl_bound_signal_tracks/bigwig", 
                        "#main/generating_nucl_bound_signal_tracks/bam"
                    ], 
                    "run": "#bed_to_coverage_track.cwl", 
                    "id": "#main/generating_nucl_bound_signal_tracks", 
                    "in": [
                        {
                            "source": "#main/generating_atac_signal_tags/bed_nucl_bound_signal", 
                            "id": "#main/generating_nucl_bound_signal_tracks/bed"
                        }, 
                        {
                            "source": "#main/reference_info", 
                            "id": "#main/generating_nucl_bound_signal_tracks/reference_info"
                        }
                    ]
                }, 
                {
                    "doc": null, 
                    "out": [
                        "#main/generating_nucl_free_signal_tracks/bigwig", 
                        "#main/generating_nucl_free_signal_tracks/bam"
                    ], 
                    "run": "#bed_to_coverage_track.cwl", 
                    "id": "#main/generating_nucl_free_signal_tracks", 
                    "in": [
                        {
                            "source": "#main/generating_atac_signal_tags/bed_nucl_free_signal", 
                            "id": "#main/generating_nucl_free_signal_tracks/bed"
                        }, 
                        {
                            "source": "#main/reference_info", 
                            "id": "#main/generating_nucl_free_signal_tracks/reference_info"
                        }
                    ]
                }, 
                {
                    "doc": null, 
                    "out": [
                        "#main/generating_tn5_bind_region_signal_tracks/bigwig", 
                        "#main/generating_tn5_bind_region_signal_tracks/bam"
                    ], 
                    "run": "#bed_to_coverage_track.cwl", 
                    "id": "#main/generating_tn5_bind_region_signal_tracks", 
                    "in": [
                        {
                            "source": "#main/generating_atac_signal_tags/bed_tn5_bind_region_signal", 
                            "id": "#main/generating_tn5_bind_region_signal_tracks/bed"
                        }, 
                        {
                            "source": "#main/reference_info", 
                            "id": "#main/generating_tn5_bind_region_signal_tracks/reference_info"
                        }
                    ]
                }, 
                {
                    "doc": null, 
                    "out": [
                        "#main/generating_tn5_center_1bp_signal_tracks/bigwig", 
                        "#main/generating_tn5_center_1bp_signal_tracks/bam"
                    ], 
                    "run": "#bed_to_coverage_track.cwl", 
                    "id": "#main/generating_tn5_center_1bp_signal_tracks", 
                    "in": [
                        {
                            "source": "#main/generating_atac_signal_tags/bed_tn5_center_1bp_signal", 
                            "id": "#main/generating_tn5_center_1bp_signal_tracks/bed"
                        }, 
                        {
                            "source": "#main/reference_info", 
                            "id": "#main/generating_tn5_center_1bp_signal_tracks/reference_info"
                        }
                    ]
                }, 
                {
                    "out": [
                        "#main/merge_duprem_filter/post_filter_fastqc_zip", 
                        "#main/merge_duprem_filter/post_filter_fastqc_html", 
                        "#main/merge_duprem_filter/picard_markdup_stdout", 
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
                            "default": true, 
                            "id": "#main/merge_duprem_filter/is_paired_end"
                        }, 
                        {
                            "source": "#main/sample_id", 
                            "id": "#main/merge_duprem_filter/sample_id"
                        }
                    ]
                }, 
                {
                    "doc": "samtools sort - sorting of filtered bam file by read name", 
                    "out": [
                        "#main/name_sorting_filtered_bam/bam_sorted"
                    ], 
                    "run": "#samtools_sort_name.cwl", 
                    "id": "#main/name_sorting_filtered_bam", 
                    "in": [
                        {
                            "source": "#main/merge_duprem_filter/bam", 
                            "id": "#main/name_sorting_filtered_bam/bam_unsorted"
                        }
                    ]
                }, 
                {
                    "doc": "NucleoATAC", 
                    "out": [
                        "#main/nucl_position_calling/nucl_occ_tracks", 
                        "#main/nucl_position_calling/nucl_occ_lower_bound_tracks", 
                        "#main/nucl_position_calling/nucl_occ_upper_bound_tracks", 
                        "#main/nucl_position_calling/nucl_dist_txt", 
                        "#main/nucl_position_calling/nucl_dist_plot", 
                        "#main/nucl_position_calling/fragsize_in_peaks_txt", 
                        "#main/nucl_position_calling/nucl_occ_fit_txt", 
                        "#main/nucl_position_calling/nucl_occ_fit_plot", 
                        "#main/nucl_position_calling/nucl_occ_peaks_bed", 
                        "#main/nucl_position_calling/nucl_vplot_data", 
                        "#main/nucl_position_calling/nucl_pos_bed", 
                        "#main/nucl_position_calling/nucl_pos_redundant_bed", 
                        "#main/nucl_position_calling/nucl_norm_crosscor_tracks", 
                        "#main/nucl_position_calling/nucl_norm_smooth_crosscor_tracks", 
                        "#main/nucl_position_calling/combined_nucl_pos_bed", 
                        "#main/nucl_position_calling/nfr_pos_bed", 
                        "#main/nucl_position_calling/nucleoatac_stderr", 
                        "#main/nucl_position_calling/nucleoatac_stdout"
                    ], 
                    "run": "#nucleoatac.cwl", 
                    "id": "#main/nucl_position_calling", 
                    "in": [
                        {
                            "source": "#main/merge_duprem_filter/bam", 
                            "id": "#main/nucl_position_calling/bam"
                        }, 
                        {
                            "source": "#main/peak_calling_macs2_tn5_center/broad_peaks_bed", 
                            "id": "#main/nucl_position_calling/bed"
                        }, 
                        {
                            "source": "#main/reference", 
                            "id": "#main/nucl_position_calling/fasta"
                        }, 
                        {
                            "source": "#main/sample_id", 
                            "valueFrom": "$(self + \".nucleoatac\")", 
                            "id": "#main/nucl_position_calling/output_basename"
                        }
                    ]
                }, 
                {
                    "doc": "peak calling using macs2", 
                    "out": [
                        "#main/peak_calling_macs2_tn5_bind_region/broad_peaks_bed", 
                        "#main/peak_calling_macs2_tn5_bind_region/gapped_peaks_bed", 
                        "#main/peak_calling_macs2_tn5_bind_region/peaks_xls"
                    ], 
                    "run": "#macs2_callpeak_atac_tn5_binding_region.cwl", 
                    "id": "#main/peak_calling_macs2_tn5_bind_region", 
                    "in": [
                        {
                            "source": "#main/macs2_genome_size", 
                            "id": "#main/peak_calling_macs2_tn5_bind_region/genome_size"
                        }, 
                        {
                            "source": "#main/sample_id", 
                            "valueFrom": "$(self + \"_tn5_bind_region.macs2\")", 
                            "id": "#main/peak_calling_macs2_tn5_bind_region/output_basename"
                        }, 
                        {
                            "source": "#main/generating_atac_signal_tags/bed_tn5_bind_region_signal", 
                            "id": "#main/peak_calling_macs2_tn5_bind_region/treatment_bed"
                        }
                    ]
                }, 
                {
                    "doc": "peak calling using macs2", 
                    "out": [
                        "#main/peak_calling_macs2_tn5_center/broad_peaks_bed", 
                        "#main/peak_calling_macs2_tn5_center/gapped_peaks_bed", 
                        "#main/peak_calling_macs2_tn5_center/peaks_xls", 
                        "#main/peak_calling_macs2_tn5_center/treat_pileup_bdg"
                    ], 
                    "run": "#macs2_callpeak_atac_tn5_center.cwl", 
                    "id": "#main/peak_calling_macs2_tn5_center", 
                    "in": [
                        {
                            "source": "#main/macs2_genome_size", 
                            "id": "#main/peak_calling_macs2_tn5_center/genome_size"
                        }, 
                        {
                            "source": "#main/sample_id", 
                            "valueFrom": "$(self + \"_tn5_center_shift_ext.macs2\")", 
                            "id": "#main/peak_calling_macs2_tn5_center/output_basename"
                        }, 
                        {
                            "source": "#main/generating_atac_signal_tags/bed_tn5_center_1bp_signal", 
                            "id": "#main/peak_calling_macs2_tn5_center/treatment_bed"
                        }
                    ]
                }, 
                {
                    "out": [
                        "#main/plot_fragment_size_distribution/frag_size_distr_plot", 
                        "#main/plot_fragment_size_distribution/frag_size_distr_tsv"
                    ], 
                    "run": "#plot_frag_size_distr.cwl", 
                    "id": "#main/plot_fragment_size_distribution", 
                    "in": [
                        {
                            "source": "#main/generating_atac_signal_tags/fragment_sizes_tsv", 
                            "id": "#main/plot_fragment_size_distribution/fragment_sizes_tsv"
                        }, 
                        {
                            "source": "#main/sample_id", 
                            "id": "#main/plot_fragment_size_distribution/output_basename"
                        }
                    ]
                }, 
                {
                    "out": [
                        "#main/qc_phantompeakqualtools/qc_crosscorr_summary", 
                        "#main/qc_phantompeakqualtools/qc_crosscorr_plot", 
                        "#main/qc_phantompeakqualtools/qc_phantompeakqualtools_stderr", 
                        "#main/qc_phantompeakqualtools/qc_phantompeakqualtools_stdout"
                    ], 
                    "run": "#phantompeakqualtools.cwl", 
                    "id": "#main/qc_phantompeakqualtools", 
                    "in": [
                        {
                            "source": "#main/merge_duprem_filter/bam", 
                            "id": "#main/qc_phantompeakqualtools/bam"
                        }
                    ]
                }, 
                {
                    "doc": "deeptools plotCoverage - plots how many times a certain fraction of the \ngenome was covered (consideres the complete fragment between a reads pair).\n", 
                    "out": [
                        "#main/qc_plot_coverage_fragments_tn5_excl/qc_plot_coverage_plot", 
                        "#main/qc_plot_coverage_fragments_tn5_excl/qc_plot_coverage_tsv"
                    ], 
                    "run": "#deeptools_plotCoverage.cwl", 
                    "id": "#main/qc_plot_coverage_fragments_tn5_excl", 
                    "in": [
                        {
                            "source": "#main/generating_fragments_tn5_excl_signal_tracks/bam", 
                            "id": "#main/qc_plot_coverage_fragments_tn5_excl/bam"
                        }, 
                        {
                            "source": "#main/sample_id", 
                            "id": "#main/qc_plot_coverage_fragments_tn5_excl/sample_id"
                        }
                    ]
                }, 
                {
                    "out": [
                        "#main/qc_plot_fingerprint/qc_plot_fingerprint_plot", 
                        "#main/qc_plot_fingerprint/qc_plot_fingerprint_tsv", 
                        "#main/qc_plot_fingerprint/qc_plot_fingerprint_stderr"
                    ], 
                    "run": "#deeptools_plotFingerprint.cwl", 
                    "id": "#main/qc_plot_fingerprint", 
                    "in": [
                        {
                            "source": "#main/generating_tn5_bind_region_signal_tracks/bam", 
                            "id": "#main/qc_plot_fingerprint/bam"
                        }, 
                        {
                            "default": true, 
                            "id": "#main/qc_plot_fingerprint/is_paired_end"
                        }, 
                        {
                            "source": "#main/sample_id", 
                            "id": "#main/qc_plot_fingerprint/sample_id"
                        }
                    ]
                }, 
                {
                    "doc": "LC_COLLATE=C sort -k1,1 -k2,2n", 
                    "out": [
                        "#main/sorting_shift_ext_bedgraph/bedgraph_sorted"
                    ], 
                    "run": "#bedgraph_sort.cwl", 
                    "id": "#main/sorting_shift_ext_bedgraph", 
                    "in": [
                        {
                            "source": "#main/clip_shift_ext_bedgraph/bed_clipped", 
                            "id": "#main/sorting_shift_ext_bedgraph/bedgraph"
                        }, 
                        {
                            "source": "#main/sample_id", 
                            "valueFrom": "$(self + \"_tn5_center_shift_ext.bedgraph\")", 
                            "id": "#main/sorting_shift_ext_bedgraph/output_name"
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
                            "default": true, 
                            "id": "#main/trim_and_map/is_paired_end"
                        }, 
                        {
                            "source": "#main/max_mapping_insert_length", 
                            "id": "#main/trim_and_map/max_mapping_insert_length"
                        }, 
                        {
                            "source": "#main/reference", 
                            "id": "#main/trim_and_map/reference"
                        }, 
                        {
                            "source": "#main/sample_id", 
                            "id": "#main/trim_and_map/sample_id"
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