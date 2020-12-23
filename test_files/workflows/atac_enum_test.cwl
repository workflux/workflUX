{
    "$graph": [
        {
            "class": "CommandLineTool",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
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
            ],
            "baseCommand": [
                "bedtools",
                "bamtobed"
            ],
            "arguments": [
                {
                    "valueFrom": "-bedpe",
                    "position": 1
                }
            ],
            "stdout": "$(inputs.bam.nameroot + \".bedpe\")",
            "inputs": [
                {
                    "type": "File",
                    "inputBinding": {
                        "prefix": "-i",
                        "position": 10
                    },
                    "id": "#bedtools_bamtobed_pe.cwl/bam"
                }
            ],
            "outputs": [
                {
                    "type": "stdout",
                    "id": "#bedtools_bamtobed_pe.cwl/bedpe"
                }
            ],
            "id": "#bedtools_bamtobed_pe.cwl"
        },
        {
            "class": "CommandLineTool",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
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
            ],
            "baseCommand": [
                "bedtools",
                "bedtobam"
            ],
            "stdout": "$(inputs.bed.nameroot + \".bam\")",
            "inputs": [
                {
                    "type": "File",
                    "inputBinding": {
                        "position": 1,
                        "prefix": "-i"
                    },
                    "id": "#bedtools_bedtobam.cwl/bed"
                },
                {
                    "type": "File",
                    "inputBinding": {
                        "position": 2,
                        "prefix": "-g"
                    },
                    "id": "#bedtools_bedtobam.cwl/genome_info"
                }
            ],
            "outputs": [
                {
                    "type": "stdout",
                    "id": "#bedtools_bedtobam.cwl/bam"
                }
            ],
            "id": "#bedtools_bedtobam.cwl"
        },
        {
            "class": "CommandLineTool",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "class": "StepInputExpressionRequirement"
                }
            ],
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
            ],
            "baseCommand": [
                "bowtie2"
            ],
            "arguments": [
                {
                    "valueFrom": "--very-sensitive",
                    "position": 1
                },
                {
                    "valueFrom": "$(runtime.cores)",
                    "prefix": "-p",
                    "position": 1
                },
                {
                    "position": 10,
                    "valueFrom": "${\n  if ( inputs.is_paired_end ){\n     return \"-1\";\n  }\n  else {\n    return \"-U\";\n  }\n}\n"
                },
                {
                    "valueFrom": "$(inputs.fastq1.nameroot + \".sam\")",
                    "prefix": "-S",
                    "position": 6
                }
            ],
            "stderr": "$( inputs.fastq1.nameroot + \".bowtie2_stderr\")",
            "inputs": [
                {
                    "type": "File",
                    "inputBinding": {
                        "position": 11
                    },
                    "id": "#bowtie2.cwl/fastq1"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "valueFrom": "${\n    if ( inputs.is_paired_end ){\n        return self;\n    }\n    else {\n      return null;\n    }\n}  \n",
                        "position": 12,
                        "prefix": "-2"
                    },
                    "id": "#bowtie2.cwl/fastq2"
                },
                {
                    "doc": "path to the FM-index files for the chosen genome genome",
                    "type": "File",
                    "secondaryFiles": [
                        ".fai",
                        "^.1.bt2",
                        "^.2.bt2",
                        "^.3.bt2",
                        "^.4.bt2",
                        "^.rev.1.bt2",
                        "^.rev.2.bt2"
                    ],
                    "inputBinding": {
                        "position": 2,
                        "prefix": "-x",
                        "valueFrom": "$(self.path.replace(/\\.fa/i,\"\"))"
                    },
                    "id": "#bowtie2.cwl/genome_index"
                },
                {
                    "type": "boolean",
                    "id": "#bowtie2.cwl/is_paired_end"
                },
                {
                    "doc": "usefull for very long fragments, as expected for ATAC",
                    "type": [
                        "null",
                        "long"
                    ],
                    "inputBinding": {
                        "prefix": "--maxins",
                        "position": 1
                    },
                    "id": "#bowtie2.cwl/max_mapping_insert_length"
                }
            ],
            "outputs": [
                {
                    "type": "stderr",
                    "id": "#bowtie2.cwl/bowtie2_log"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*.sam"
                    },
                    "id": "#bowtie2.cwl/sam"
                }
            ],
            "id": "#bowtie2.cwl"
        },
        {
            "class": "CommandLineTool",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "class": "StepInputExpressionRequirement"
                }
            ],
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
            ],
            "baseCommand": [
                "bamCoverage"
            ],
            "arguments": [
                {
                    "valueFrom": "$(inputs.bam.nameroot + \".bigwig\")",
                    "prefix": "--outFileName",
                    "position": 10
                },
                {
                    "valueFrom": "bigwig",
                    "prefix": "--outFileFormat",
                    "position": 10
                },
                {
                    "valueFrom": "${ \n  if( inputs.spike_in_count == null ){\n    return \"RPGC\"\n  }\n  else{\n    return null \n  }\n}\n",
                    "prefix": "--normalizeUsing",
                    "position": 10
                }
            ],
            "inputs": [
                {
                    "doc": "bam file as input; needs bai index file in the same directory",
                    "type": "File",
                    "secondaryFiles": ".bai",
                    "inputBinding": {
                        "position": 100,
                        "prefix": "--bam"
                    },
                    "id": "#deeptools_bamCoverage.cwl/bam"
                },
                {
                    "type": "int",
                    "default": 10,
                    "inputBinding": {
                        "prefix": "--binSize",
                        "position": 10
                    },
                    "id": "#deeptools_bamCoverage.cwl/bin_size"
                },
                {
                    "doc": "the effectively mappable genome size, \nsee: https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html\n",
                    "type": "long",
                    "inputBinding": {
                        "position": 10,
                        "prefix": "--effectiveGenomeSize"
                    },
                    "id": "#deeptools_bamCoverage.cwl/effective_genome_size"
                },
                {
                    "doc": "List of space-delimited chromosome names that shall be ignored \nwhen calculating the scaling factor. \nSpecify as space-delimited string. \nDefault: \"chrX chrY chrM\"\n",
                    "type": [
                        "null",
                        "string"
                    ],
                    "default": "chrX chrY chrM",
                    "inputBinding": {
                        "prefix": "--ignoreForNormalization",
                        "position": 10
                    },
                    "id": "#deeptools_bamCoverage.cwl/ignoreForNormalization"
                },
                {
                    "doc": "number of reads aligned to the spike in genome, optional",
                    "type": [
                        "null",
                        "long"
                    ],
                    "inputBinding": {
                        "position": 10,
                        "prefix": "--scaleFactor",
                        "valueFrom": "${ \n  if( self == null ){\n    return null\n  }\n  else{\n    return (1.0 / parseFloat(self)) \n  }\n}\n"
                    },
                    "id": "#deeptools_bamCoverage.cwl/spike_in_count"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.bam.nameroot + \".bigwig\")"
                    },
                    "id": "#deeptools_bamCoverage.cwl/bigwig"
                }
            ],
            "id": "#deeptools_bamCoverage.cwl"
        },
        {
            "class": "CommandLineTool",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "class": "StepInputExpressionRequirement"
                }
            ],
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
            ],
            "baseCommand": [
                "plotFingerprint"
            ],
            "arguments": [
                {
                    "valueFrom": "${\n  if ( inputs.is_paired_end ){\n     return null;\n  }\n  else {\n    return \"--extendReads\";\n  }\n}\n",
                    "position": 1
                },
                {
                    "valueFrom": "${\n  if ( inputs.is_paired_end ){\n     return null;\n  }\n  else {\n    return inputs.fragment_size;\n  }\n}\n",
                    "position": 2
                },
                {
                    "valueFrom": "$(inputs.sample_id)",
                    "prefix": "--labels",
                    "position": 10
                },
                {
                    "valueFrom": "$(inputs.sample_id + \".plot_fingerp.png\")",
                    "prefix": "--plotFile",
                    "position": 10
                },
                {
                    "valueFrom": "$(inputs.sample_id + \".plot_fingerp.tsv\")",
                    "prefix": "--outRawCounts",
                    "position": 10
                }
            ],
            "stderr": "$( inputs.sample_id + \".plot_fingerp.stderr\")",
            "inputs": [
                {
                    "doc": "must be indexed",
                    "type": "File",
                    "secondaryFiles": ".bai",
                    "inputBinding": {
                        "position": 100,
                        "prefix": "--bamfiles"
                    },
                    "id": "#deeptools_plotFingerprint.cwl/bam"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#deeptools_plotFingerprint.cwl/fragment_size"
                },
                {
                    "doc": "if paired end, reads are extended",
                    "type": "boolean",
                    "default": true,
                    "id": "#deeptools_plotFingerprint.cwl/is_paired_end"
                },
                {
                    "type": "string",
                    "id": "#deeptools_plotFingerprint.cwl/sample_id"
                }
            ],
            "outputs": [
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "$(inputs.sample_id + \".plot_fingerp.png\")"
                    },
                    "id": "#deeptools_plotFingerprint.cwl/qc_plot_fingerprint_plot"
                },
                {
                    "type": "stderr",
                    "id": "#deeptools_plotFingerprint.cwl/qc_plot_fingerprint_stderr"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "$(inputs.sample_id + \".plot_fingerp.tsv\")"
                    },
                    "id": "#deeptools_plotFingerprint.cwl/qc_plot_fingerprint_tsv"
                }
            ],
            "successCodes": [
                0,
                1,
                2
            ],
            "temporaryFailCodes": [],
            "permanentFailCodes": [],
            "id": "#deeptools_plotFingerprint.cwl"
        },
        {
            "class": "CommandLineTool",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "class": "StepInputExpressionRequirement"
                }
            ],
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
            ],
            "baseCommand": "fastqc",
            "arguments": [
                {
                    "valueFrom": "$(runtime.outdir)",
                    "prefix": "-o"
                },
                {
                    "valueFrom": "--noextract"
                }
            ],
            "inputs": [
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "position": 1
                    },
                    "id": "#fastqc.cwl/bam"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "position": 1
                    },
                    "id": "#fastqc.cwl/fastq1"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "position": 2
                    },
                    "id": "#fastqc.cwl/fastq2"
                }
            ],
            "outputs": [
                {
                    "doc": "html report showing results from zip",
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "outputBinding": {
                        "glob": "*_fastqc.html"
                    },
                    "id": "#fastqc.cwl/fastqc_html"
                },
                {
                    "doc": "all data e.g. figures",
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "outputBinding": {
                        "glob": "*_fastqc.zip"
                    },
                    "id": "#fastqc.cwl/fastqc_zip"
                }
            ],
            "id": "#fastqc.cwl"
        },
        {
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
            ],
            "requirements": [
                {
                    "listing": [
                        {
                            "entryname": "generate_atac_signal_tag.sh",
                            "entry": "BEDPE=\"$1\"\nOUTPUT_BASENAME=\"$2\"\ntouch irreg_mappings.bedpe \ntouch fragment_sizes.txt\ntouch tn5_center_29bp_unsorted.bed \ntouch tn5_center_73bp_unsorted.bed\ntouch tn5_center_200bp_unsorted.bed \ntouch tn5_center_fragment_unsorted.bed\ntouch tn5_center_1bp_unsorted.bed\ncat \"\\${BEDPE}\" | awk -f generate_atac_signal_tag.awk\nsort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n tn5_center_29bp_unsorted.bed > \"\\${OUTPUT_BASENAME}_tn5_center_29bp.bed\"\nsort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n tn5_center_73bp_unsorted.bed > \"\\${OUTPUT_BASENAME}_tn5_center_73bp.bed\"\nsort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n tn5_center_200bp_unsorted.bed > \"\\${OUTPUT_BASENAME}_tn5_center_200bp.bed\"\nsort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n tn5_center_fragment_unsorted.bed > \"\\${OUTPUT_BASENAME}_tn5_center_fragment.bed\"\nsort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n tn5_center_1bp_unsorted.bed > \"\\${OUTPUT_BASENAME}_tn5_center_1bp.bed\"\nrm tn5_center_29bp_unsorted.bed \nrm tn5_center_73bp_unsorted.bed   \nrm tn5_center_200bp_unsorted.bed  \nrm tn5_center_fragment_unsorted.bed \nrm tn5_center_1bp_unsorted.bed \nmv irreg_mappings.bedpe \"\\${OUTPUT_BASENAME}_irreg_mappings.bedpe\"\nmv fragment_sizes.txt \"\\${OUTPUT_BASENAME}_fragment_sizes.txt\"\nmv filtering_stats_mqc.tsv \"\\${OUTPUT_BASENAME}_filtering_stats_mqc.tsv\"\nmv frag_size_classification_mqc.tsv \"\\${OUTPUT_BASENAME}_frag_size_classification_mqc.tsv\"\n"
                        },
                        {
                            "entryname": "generate_atac_signal_tag.awk",
                            "entry": "function tn5_center_ext(seqname, name, score, bp, center1, center2, output_file) { \n  if( (bp%2) == 0 ) {\n    start_shift = bp/2\n    end_shift = bp/2\n  }\n  else {\n    start_shift = (bp-1)/2\n    end_shift = (bp-1)/2+1\n  }\n  tag_start = (center1-start_shift);\n  tag_end = (center1+end_shift+1);\n  if(tag_start < 0){tag_start=0};\n  if(tag_end < 0){tag_end=0};\n  print seqname, tag_start, tag_end, name, score > output_file;\n  tag_start = (center2-start_shift);\n  tag_end = (center2+end_shift+1);\n  if(tag_start < 0){tag_start=0};\n  if(tag_end < 0){tag_end=0};\n  print seqname, tag_start, tag_end, name\"_mate\", score > output_file;\n}\nBEGIN {\n    OFS=\"\\t\";\n    chrM_read_count=0;\n    interchrom_map_read_count=0;\n    regular_read_count=0;\n    irregular_read_count=0;\n    too_small_fragment_count=0;\n    nucl_free_fragment_count=0;\n    nucl_bound_fragment_count=0;\n    wrong_strand_orient_count=0;\n}\n{\n  if ( $1==\"chrM\" || $4==\"chrM\") {\n        irregular_read_count += 2;\n        if ( $1==$4 ) {\n            chrM_read_count += 2;\n            print $0, \"chrM\" > \"irreg_mappings.bedpe\" ;\n        }\n        else {\n            chrM_read_count += 1;\n            interchrom_map_read_count += 2;\n            print $0, \"interchrom_map\" > \"irreg_mappings.bedpe\";\n        }\n    }\n    else if ( $1!=$4 ) {\n        irregular_read_count += 2;\n        interchrom_map_read_count += 2;\n        print $0, \"interchrom_map\" > \"irreg_mappings.bedpe\";\n    }\n    else if ( $9==$10 ) {\n        irregular_read_count += 2;\n        wrong_strand_orient_count += 2;\n        print $0, \"wrong_read_pair_orientation\" > \"irreg_mappings.bedpe\";\n    }\n    else {\n        right_orientation=0;\n        if ( $9==\"+\" && $2<=$5 && $3<=$6 ){\n            right_orientation=1;\n            fragment_size=$6-$2;\n            center1=$2+4;\n            center2=$6-5;\n        }\n        else if ( $9==\"-\" && $2>=$5 && $3>=$6 ){\n            right_orientation=1;\n            fragment_size=$3-$5;\n            center1=$5+4;\n            center2=$3-5;\n        }\n        else {\n            wrong_strand_orient_count += 2;\n            print $0, \"wrong_read_pair_orientation\" > \"irreg_mappings.bedpe\";\n        }\n        print fragment_size >  \"fragment_sizes.txt\";\n        if ( fragment_size>=38 && right_orientation ) {\n            regular_read_count += 2;\n            tn5_center_ext($1, $7, $8, 29, center1, center2, \"tn5_center_29bp_unsorted.bed\");\n            tn5_center_ext($1, $7, $8, 73, center1, center2, \"tn5_center_73bp_unsorted.bed\");\n            tn5_center_ext($1, $7, $8, 200, center1, center2, \"tn5_center_200bp_unsorted.bed\");\n            tn5_center_ext($1, $7, $8, 1, center1, center2, \"tn5_center_1bp_unsorted.bed\");\n            print $1, center1, center2, $7, $8 >  \"tn5_center_fragment_unsorted.bed\";\n            if (fragment_size>=185) { \n                nucl_bound_fragment_count += 1;\n            }\n            else {\n                nucl_free_fragment_count += 1;\n            }\n        }\n        else {\n            if ( right_orientation ){\n                irregular_read_count += 2;\n                too_small_fragment_count += 1;\n                print $0, \"too_small_frag_size\" > \"irreg_mappings.bedpe\" ;\n            }\n        }\n  }\n}\nEND {\n    print \"\\# id: \\\"Filtering Statistics\\\"\" >  \"filtering_stats_mqc.tsv\";\n    print \"\\# description: \\\"- This section shows statistics on read filtering: (1) reads pairs that are mapping to different chromosomes as well as (2) read pairs located on ChrM are filtered out; (3) read pairs that have a wrong orientation towards each other (for instance both reads on same strand, or reads pointing to different direction) are removed, too; (4) Only the remaining regular reads are used for the fragment size analysis and the generation of atac signal tracks; in addition to the filtering shown here, reads were also selected for high mapping quality and to be mapped in a proper pair\\\"\" >  \"filtering_stats_mqc.tsv\";\n    print \"\\# plot_type: \\\"bargraph\\\"\" >  \"filtering_stats_mqc.tsv\";\n    print \"chrM_reads\", chrM_read_count >  \"filtering_stats_mqc.tsv\";\n    print \"interchrom_map_reads\", interchrom_map_read_count >  \"filtering_stats_mqc.tsv\";\n    print \"wrong_read_pair_orientation\", wrong_strand_orient_count >  \"filtering_stats_mqc.tsv\";\n    print \"regular_reads\", regular_read_count >  \"filtering_stats_mqc.tsv\";\n    print \"\\# id: \\\"Fragment Length Classification\\\"\" >  \"frag_size_classification_mqc.tsv\";\n    print \"\\# description: \\\"Fragments are classified by their size: (1) nucleosome free fragements are smaller than a typical a nucleosome binding region while (2)) potentially nucleosome bound fragments are larger;(3) fragments that are classified as too small are shorter than expected for a Tn5 digestion; for these calculations, the DNA span that is covered by the Tn5 enzyme during the trasposition is taken into account\\\"\" >  \"frag_size_classification_mqc.tsv\";\n    print \"\\# plot_type: \\\"bargraph\\\"\" >  \"frag_size_classification_mqc.tsv\";\n    print \"too_small_fragments\", too_small_fragment_count >  \"frag_size_classification_mqc.tsv\";\n    print \"nucl_free_fragments\", nucl_free_fragment_count >  \"frag_size_classification_mqc.tsv\";\n    print \"pot_nucl_bound_fragments\", nucl_bound_fragment_count >  \"frag_size_classification_mqc.tsv\";\n}\n"
                        }
                    ],
                    "class": "InitialWorkDirRequirement"
                },
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "class": "StepInputExpressionRequirement"
                }
            ],
            "baseCommand": [
                "bash",
                "generate_atac_signal_tag.sh"
            ],
            "inputs": [
                {
                    "type": "File",
                    "inputBinding": {
                        "position": 1
                    },
                    "id": "#generate_atac_signal_tags.cwl/bedpe_alignm"
                },
                {
                    "type": "string",
                    "inputBinding": {
                        "position": 2
                    },
                    "id": "#generate_atac_signal_tags.cwl/output_basename"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.output_basename + \"_tn5_center_1bp.bed\")"
                    },
                    "id": "#generate_atac_signal_tags.cwl/bed_tn5_center_1bp"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.output_basename + \"_tn5_center_200bp.bed\")"
                    },
                    "id": "#generate_atac_signal_tags.cwl/bed_tn5_center_200bp"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.output_basename + \"_tn5_center_29bp.bed\")"
                    },
                    "id": "#generate_atac_signal_tags.cwl/bed_tn5_center_29bp"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.output_basename + \"_tn5_center_73bp.bed\")"
                    },
                    "id": "#generate_atac_signal_tags.cwl/bed_tn5_center_73bp"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.output_basename + \"_tn5_center_fragment.bed\")"
                    },
                    "id": "#generate_atac_signal_tags.cwl/bed_tn5_center_fragment"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.output_basename + \"_filtering_stats_mqc.tsv\")"
                    },
                    "id": "#generate_atac_signal_tags.cwl/filtering_stats_tsv"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.output_basename + \"_frag_size_classification_mqc.tsv\")"
                    },
                    "id": "#generate_atac_signal_tags.cwl/frag_size_stats_tsv"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.output_basename + \"_fragment_sizes.txt\")"
                    },
                    "id": "#generate_atac_signal_tags.cwl/fragment_sizes_tsv"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.output_basename + \"_irreg_mappings.bedpe\")"
                    },
                    "id": "#generate_atac_signal_tags.cwl/irreg_mappings_bedpe"
                }
            ],
            "id": "#generate_atac_signal_tags.cwl"
        },
        {
            "class": "CommandLineTool",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "class": "StepInputExpressionRequirement"
                }
            ],
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
            ],
            "baseCommand": [
                "macs2",
                "callpeak"
            ],
            "arguments": [
                {
                    "valueFrom": "BED",
                    "prefix": "--format",
                    "position": 1
                },
                {
                    "valueFrom": "--nomodel",
                    "position": 2
                },
                {
                    "valueFrom": "all",
                    "prefix": "--keep-dup",
                    "position": 2
                },
                {
                    "valueFrom": "$(inputs.treatment_bed.nameroot + \".macs2\")",
                    "prefix": "--name",
                    "position": 100
                }
            ],
            "inputs": [
                {
                    "type": "boolean",
                    "inputBinding": {
                        "position": 3,
                        "prefix": "--broad"
                    },
                    "id": "#macs2_callpeak_atac.cwl/broad"
                },
                {
                    "type": "long",
                    "inputBinding": {
                        "position": 3,
                        "prefix": "--gsize"
                    },
                    "id": "#macs2_callpeak_atac.cwl/genome_size"
                },
                {
                    "type": "float",
                    "inputBinding": {
                        "position": 3,
                        "prefix": "--qvalue"
                    },
                    "id": "#macs2_callpeak_atac.cwl/qvalue"
                },
                {
                    "type": "File",
                    "inputBinding": {
                        "position": 101,
                        "prefix": "--treatment"
                    },
                    "id": "#macs2_callpeak_atac.cwl/treatment_bed"
                }
            ],
            "outputs": [
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "outputBinding": {
                        "glob": "*Peak"
                    },
                    "id": "#macs2_callpeak_atac.cwl/peaks_bed"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*_peaks.xls"
                    },
                    "id": "#macs2_callpeak_atac.cwl/peaks_xls"
                }
            ],
            "id": "#macs2_callpeak_atac.cwl"
        },
        {
            "class": "CommandLineTool",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "class": "StepInputExpressionRequirement"
                }
            ],
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
            ],
            "baseCommand": [
                "bash",
                "-c"
            ],
            "arguments": [
                {
                    "valueFrom": "${\n    var qc_files_array = inputs.qc_files_array;\n    var qc_files_array_of_array = inputs.qc_files_array_of_array;\n    var cmdline = \"echo 'copying input file ...'\";\n\n    if ( qc_files_array != null ){\n      for (var i=0; i<qc_files_array.length; i++){\n        if( qc_files_array[i] != null ){\n          cmdline += \"; cp \" + qc_files_array[i].path + \" .\";\n        }\n      }\n    }\n\n    if ( qc_files_array_of_array != null ){\n      for (var i=0; i<qc_files_array_of_array.length; i++){ \n        for (var ii=0; ii<qc_files_array_of_array[i].length; ii++){\n          if( qc_files_array_of_array[i][ii] != null ){\n            cmdline += \"; cp \" + qc_files_array_of_array[i][ii].path + \" .\";\n          }\n        }\n      }\n    }\n    \n    cmdline += \"; echo \\'copying done\\'\" +\n        \"; multiqc --zip-data-dir --cl_config \\'log_filesize_limit: 100000000\\' \" +\n        \"--outdir \" + runtime.outdir +\n        \" --filename \" + inputs.report_name + \"_report .\";\n\n    return cmdline\n  }\n"
                }
            ],
            "inputs": [
                {
                    "doc": "qc files which shall be part of the multiqc summary;\noptional, only one of qc_files_array or qc_files_array_of_array \nmust be provided\n",
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": [
                                "File",
                                "null"
                            ]
                        }
                    ],
                    "id": "#multiqc_hack.cwl/qc_files_array"
                },
                {
                    "doc": "qc files which shall be part of the multiqc summary;\noptional, only one of qc_files_array or qc_files_array_of_array \nmust be provided\n",
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": {
                                "type": "array",
                                "items": [
                                    "File",
                                    "null"
                                ]
                            }
                        }
                    ],
                    "id": "#multiqc_hack.cwl/qc_files_array_of_array"
                },
                {
                    "doc": "name used for the html report and the corresponding zip file",
                    "type": "string",
                    "default": "multiqc",
                    "id": "#multiqc_hack.cwl/report_name"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.report_name + \"_report.html\")"
                    },
                    "id": "#multiqc_hack.cwl/multiqc_html"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.report_name + \"_report_data.zip\")"
                    },
                    "id": "#multiqc_hack.cwl/multiqc_zip"
                }
            ],
            "id": "#multiqc_hack.cwl"
        },
        {
            "class": "CommandLineTool",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "class": "StepInputExpressionRequirement"
                }
            ],
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
            ],
            "baseCommand": [
                "Rscript",
                "--verbose",
                "--max-ppsize=500000",
                "/usr/bin/phantompeakqualtools-1.2/run_spp.R"
            ],
            "arguments": [
                {
                    "valueFrom": "$(runtime.tmpdir)",
                    "prefix": "-tmpdir=",
                    "separate": false,
                    "position": 10
                },
                {
                    "valueFrom": "$(runtime.outdir)",
                    "prefix": "-odir=",
                    "separate": false,
                    "position": 10
                },
                {
                    "valueFrom": "$(inputs.bam.nameroot + \".crosscor.pdf\")",
                    "prefix": "-savp=",
                    "separate": false,
                    "position": 100
                },
                {
                    "valueFrom": "$(inputs.bam.nameroot + \".spp.out\")",
                    "prefix": "-out=",
                    "separate": false,
                    "position": 100
                }
            ],
            "stderr": "$(inputs.bam.nameroot + \".phantompeakqualtools_stderr\")",
            "stdout": "$(inputs.bam.nameroot + \".phantompeakqualtools_stdout\")",
            "inputs": [
                {
                    "type": "File",
                    "inputBinding": {
                        "prefix": "-c=",
                        "separate": false,
                        "position": 10
                    },
                    "id": "#phantompeakqualtools.cwl/bam"
                }
            ],
            "outputs": [
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "*.pdf"
                    },
                    "id": "#phantompeakqualtools.cwl/qc_crosscorr_plot"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "*.spp.out"
                    },
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
            "successCodes": [
                0,
                1,
                2
            ],
            "temporaryFailCodes": [],
            "permanentFailCodes": [],
            "id": "#phantompeakqualtools.cwl"
        },
        {
            "class": "CommandLineTool",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "class": "StepInputExpressionRequirement"
                }
            ],
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
            ],
            "baseCommand": [
                "java",
                "-jar"
            ],
            "arguments": [
                {
                    "valueFrom": "MarkDuplicates",
                    "position": 2
                },
                {
                    "valueFrom": "$(inputs.bam_sorted.nameroot + \"_duprem.bam\")",
                    "prefix": "OUTPUT=",
                    "separate": false,
                    "position": 13
                },
                {
                    "valueFrom": "$(inputs.bam_sorted.nameroot + \"_duprem.log\")",
                    "prefix": "METRICS_FILE=",
                    "separate": false,
                    "position": 13
                },
                {
                    "valueFrom": "REMOVE_DUPLICATES=TRUE",
                    "position": 14
                },
                {
                    "valueFrom": "ASSUME_SORTED=TRUE",
                    "position": 15
                },
                {
                    "valueFrom": "VALIDATION_STRINGENCY=SILENT",
                    "position": 16
                },
                {
                    "valueFrom": "VERBOSITY=INFO",
                    "position": 17
                },
                {
                    "valueFrom": "QUIET=false",
                    "position": 17
                }
            ],
            "stderr": "$(inputs.bam_sorted.nameroot + \".picard_markdup.log\")",
            "inputs": [
                {
                    "doc": "sorted bam input file",
                    "type": "File",
                    "inputBinding": {
                        "prefix": "INPUT=",
                        "separate": false,
                        "position": 11
                    },
                    "id": "#picard_markdup.cwl/bam_sorted"
                },
                {
                    "type": "string",
                    "default": "/bin/picard.jar",
                    "inputBinding": {
                        "position": 1
                    },
                    "id": "#picard_markdup.cwl/path_to_picards"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.bam_sorted.nameroot + \"_duprem.bam\")"
                    },
                    "id": "#picard_markdup.cwl/bam_duprem"
                },
                {
                    "type": "stderr",
                    "id": "#picard_markdup.cwl/picard_markdup_log"
                }
            ],
            "id": "#picard_markdup.cwl"
        },
        {
            "class": "CommandLineTool",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
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
            ],
            "baseCommand": [
                "Rscript",
                "-e"
            ],
            "arguments": [
                {
                    "valueFrom": "${\n  var r_cmd = \"library(tidyverse); \" +\n            \"library(gridExtra); \" +\n            \"size_freqs <- read.table('\" + \n                inputs.fragment_sizes_tsv.path +\n                  \"', header=F)[,1] %>% as.integer() %>% table(); \" +\n            \"size_freq_tb <- tibble(frag_size=1:max(as.integer(names(size_freqs)))) %>% \" +\n              \"full_join( tibble( frag_size=as.integer(names(size_freqs)), freq=size_freqs), by='frag_size') %>% \" +\n              \"mutate(freq=sapply(freq,function(f) ifelse(is.na(f), 0, f))); \" +\n            \"plot_ <- ggplot(size_freq_tb) + \" +\n              \"geom_line(mapping=aes(x=frag_size,y=freq)) + \" +\n              \"ggtitle('Fragment Size Distribution') + \" +\n              \"xlab('') + ylab('frequency'); \" +\n            \"plot_log <- ggplot(size_freq_tb) + \" +\n              \"geom_line(mapping=aes(x=frag_size,y=freq)) + \" +\n              \"scale_y_log10() + \" +\n              \"xlab('fragment size [bp]') + ylab('frequency'); \" +\n            \"png( '\" + inputs.output_basename + \"_frag_size_distr.png', \" +\n                \"width = 850, height = 850); \" +\n            \"grid.arrange(plot_, plot_log, nrow=2);\" +\n            \"dev.off();\"+\n            \"mqc_file_headers <- c(\"+\n              \"'# id: \\\\\\\"Fragment Size Distribution\\\\\\\"', \" +\n              \"'# description: \\\\\\\"Distribution of fragment sizes (<1000) which typically shows the DNA pitch and multimers of nucleosomes. \" + \n                              \"For a better visualization, which also includes higher fragment sizes, please see the files ending with _frag_size_distr.png\\\\\\\"', \" +\n              \"'# pconfig:', \" +\n              \"'#    xmax: 1000'\" +\n              \");\" +\n            \"mqc_file <- file('\" + inputs.output_basename + \"_fragm_sizes_mqc.tsv');\" +\n            \"writeLines(mqc_file_headers, mqc_file);\" +\n            \"close(mqc_file);\" +\n            \"write.table(size_freq_tb, file=\\'\" + inputs.output_basename + \"_fragm_sizes_mqc.tsv\\', \" +\n              \"sep='\\\\t', row.names=F, col.names=F, append=T);\" +\n            \"print('done')\"\n  return r_cmd\n}\n"
                }
            ],
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
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.output_basename + \"_frag_size_distr.png\")"
                    },
                    "id": "#plot_frag_size_distr.cwl/frag_size_distr_plot"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.output_basename + \"_fragm_sizes_mqc.tsv\")"
                    },
                    "id": "#plot_frag_size_distr.cwl/frag_size_distr_tsv"
                }
            ],
            "id": "#plot_frag_size_distr.cwl"
        },
        {
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
            ],
            "baseCommand": [
                "samtools",
                "flagstat"
            ],
            "stdout": "$(inputs.bam.nameroot + inputs.output_suffix)",
            "inputs": [
                {
                    "type": "File",
                    "inputBinding": {
                        "position": 2
                    },
                    "id": "#samtools_flagstat.cwl/bam"
                },
                {
                    "type": "string",
                    "default": ".flagStat",
                    "id": "#samtools_flagstat.cwl/output_suffix"
                }
            ],
            "outputs": [
                {
                    "type": "stdout",
                    "id": "#samtools_flagstat.cwl/flagstat_output"
                }
            ],
            "id": "#samtools_flagstat.cwl"
        },
        {
            "class": "CommandLineTool",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
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
            ],
            "baseCommand": [
                "bash",
                "-c"
            ],
            "arguments": [
                {
                    "valueFrom": "$(\"cp \" + inputs.bam_sorted.path + \" . && samtools index -b \" + inputs.bam_sorted.basename )"
                }
            ],
            "inputs": [
                {
                    "doc": "sorted bam input file",
                    "type": "File",
                    "id": "#samtools_index_hack.cwl/bam_sorted"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "secondaryFiles": ".bai",
                    "outputBinding": {
                        "glob": "$(inputs.bam_sorted.basename)"
                    },
                    "id": "#samtools_index_hack.cwl/bam_sorted_indexed"
                }
            ],
            "id": "#samtools_index_hack.cwl"
        },
        {
            "class": "CommandLineTool",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
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
            ],
            "baseCommand": [
                "samtools",
                "merge"
            ],
            "inputs": [
                {
                    "id": "#samtools_merge.cwl/output_name",
                    "doc": "name of merged bam file",
                    "type": "string",
                    "inputBinding": {
                        "position": 1
                    }
                },
                {
                    "id": "#samtools_merge.cwl/bams",
                    "doc": "bam files to be merged",
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "inputBinding": {
                        "position": 2
                    }
                }
            ],
            "outputs": [
                {
                    "id": "#samtools_merge.cwl/bam_merged",
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.output_name)"
                    }
                }
            ],
            "id": "#samtools_merge.cwl"
        },
        {
            "doc": "Sort a bam file by read names.",
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
            ],
            "baseCommand": [
                "samtools",
                "sort"
            ],
            "arguments": [
                {
                    "valueFrom": "$(runtime.cores)",
                    "prefix": "-@"
                }
            ],
            "inputs": [
                {
                    "doc": "aligned reads to be checked in sam or bam format",
                    "type": "File",
                    "inputBinding": {
                        "position": 2
                    },
                    "id": "#samtools_sort.cwl/bam_unsorted"
                }
            ],
            "stdout": "$(inputs.bam_unsorted.basename)",
            "outputs": [
                {
                    "type": "stdout",
                    "id": "#samtools_sort.cwl/bam_sorted"
                }
            ],
            "id": "#samtools_sort.cwl"
        },
        {
            "doc": "Sort a bam file by read names.",
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
            ],
            "baseCommand": [
                "samtools",
                "sort"
            ],
            "arguments": [
                {
                    "valueFrom": "$(runtime.cores)",
                    "prefix": "-@",
                    "position": 1
                },
                {
                    "valueFrom": "-n",
                    "position": 1
                }
            ],
            "inputs": [
                {
                    "doc": "aligned reads to be checked in sam or bam format",
                    "type": "File",
                    "inputBinding": {
                        "position": 2
                    },
                    "id": "#samtools_sort_name.cwl/bam_unsorted"
                }
            ],
            "stdout": "$(inputs.bam_unsorted.basename)",
            "outputs": [
                {
                    "type": "stdout",
                    "id": "#samtools_sort_name.cwl/bam_sorted"
                }
            ],
            "id": "#samtools_sort_name.cwl"
        },
        {
            "class": "CommandLineTool",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
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
            ],
            "baseCommand": [
                "samtools",
                "view"
            ],
            "inputs": [
                {
                    "doc": "aligned reads to be checked in bam format",
                    "type": "File",
                    "inputBinding": {
                        "position": 10
                    },
                    "id": "#samtools_view_filter.cwl/bam"
                },
                {
                    "doc": "if paired end, only properly paired reads pass",
                    "type": "boolean",
                    "default": true,
                    "id": "#samtools_view_filter.cwl/is_paired_end"
                }
            ],
            "arguments": [
                {
                    "valueFrom": "-h",
                    "position": 1
                },
                {
                    "valueFrom": "-b",
                    "position": 1
                },
                {
                    "valueFrom": "4",
                    "prefix": "-F",
                    "position": 1
                },
                {
                    "valueFrom": "20",
                    "prefix": "-q",
                    "position": 1
                },
                {
                    "valueFrom": "${\n  if ( inputs.is_paired_end ){\n     return \"-f\";\n  }\n  else {\n    return null;\n  }\n}\n",
                    "position": 2
                },
                {
                    "valueFrom": "${\n  if ( inputs.is_paired_end ){\n     return \"3\";\n  }\n  else {\n    return null;\n  }\n}\n",
                    "position": 3
                }
            ],
            "stdout": "$(inputs.bam.nameroot)_filt.bam",
            "outputs": [
                {
                    "type": "stdout",
                    "id": "#samtools_view_filter.cwl/bam_filtered"
                }
            ],
            "id": "#samtools_view_filter.cwl"
        },
        {
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
            ],
            "baseCommand": [
                "samtools",
                "view"
            ],
            "inputs": [
                {
                    "doc": "aligned reads to be checked in sam or bam format",
                    "type": "File",
                    "inputBinding": {
                        "position": 2
                    },
                    "id": "#samtools_view_sam2bam.cwl/sam"
                }
            ],
            "arguments": [
                {
                    "valueFrom": "-h",
                    "position": 1
                },
                {
                    "valueFrom": "-b",
                    "position": 1
                }
            ],
            "stdout": "$(inputs.sam.nameroot).bam",
            "outputs": [
                {
                    "type": "stdout",
                    "id": "#samtools_view_sam2bam.cwl/bam_unsorted"
                }
            ],
            "id": "#samtools_view_sam2bam.cwl"
        },
        {
            "class": "CommandLineTool",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
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
            ],
            "baseCommand": "trim_galore",
            "inputs": [
                {
                    "doc": "Adapter to be trimmed from first reads. Cab be one of the following: \\n\n- \"nextera\" for the Nextera adapter (CTGTCTCTTATA)\\n\n- \"illumina\" for the Illumina universal adapter (AGATCGGAAGAGC)\\n\n- \"small_rna\" for the Illumina Small RNA 3-prime Adapter (TGGAATTCTCGG)\\n\n- \"auto\" to automatically detect the write setting\n",
                    "type": {
                        "type": "enum",
                        "symbols": [
                            "#trim_galore.cwl/adapter1/nextera",
                            "#trim_galore.cwl/adapter1/illumina",
                            "#trim_galore.cwl/adapter1/small_rna",
                            "#trim_galore.cwl/adapter1/auto"
                        ]
                    },
                    "id": "#trim_galore.cwl/adapter1"
                },
                {
                    "doc": "Adapters to be trimmed from second read. Cab be one of the following: \\n\n- \"nextera\" for the Nextera adapter (CTGTCTCTTATA)\\n\n- \"illumina\" for the Illumina universal adapter (AGATCGGAAGAGC)\\n\n- \"small_rna\" for the Illumina Small RNA 3-prime Adapter (TGGAATTCTCGG)\\n\n- \"auto\" to automatically detect the write setting\n",
                    "type": {
                        "type": "enum",
                        "symbols": [
                            "#trim_galore.cwl/adapter2/nextera",
                            "#trim_galore.cwl/adapter2/illumina",
                            "#trim_galore.cwl/adapter2/small_rna",
                            "#trim_galore.cwl/adapter2/auto"
                        ]
                    },
                    "id": "#trim_galore.cwl/adapter2"
                },
                {
                    "doc": "raw reads in fastq format; can be gzipped;\nif paired end, the file contains the first reads;\nif single end, the file contains all reads\n",
                    "type": "File",
                    "inputBinding": {
                        "position": 10
                    },
                    "id": "#trim_galore.cwl/fastq1"
                },
                {
                    "doc": "(optional) raw reads in fastq format; can be gzipped;\nif paired end, the file contains the second reads;\nif single end, the file does not exist\n",
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "position": 11
                    },
                    "id": "#trim_galore.cwl/fastq2"
                },
                {
                    "doc": "minimum overlap with adapter seq in bp needed to trim",
                    "type": "int",
                    "default": 1,
                    "inputBinding": {
                        "prefix": "--stringency",
                        "position": 1
                    },
                    "id": "#trim_galore.cwl/min_adapter_overlap"
                },
                {
                    "doc": "discard reads that get shorter than this value",
                    "type": "int",
                    "default": 20,
                    "inputBinding": {
                        "prefix": "--length",
                        "position": 1
                    },
                    "id": "#trim_galore.cwl/min_read_length"
                },
                {
                    "doc": "if only one read of a pair passes the qc and adapter trimming,\nit needs at least this length to be rescued\n",
                    "type": "int",
                    "default": 35,
                    "id": "#trim_galore.cwl/min_unpaired_read_rescue_length"
                },
                {
                    "doc": "trim all base with a phred score lower than this valueFrom",
                    "type": "int",
                    "default": 20,
                    "inputBinding": {
                        "prefix": "--quality",
                        "position": 1
                    },
                    "id": "#trim_galore.cwl/qual_trim_cutoff"
                }
            ],
            "arguments": [
                {
                    "prefix": "--fastqc_args",
                    "valueFrom": "\"--noextract\"",
                    "position": 1
                },
                {
                    "prefix": "--gzip",
                    "position": 1
                },
                {
                    "valueFrom": "${\n  if ( inputs.adapter1 == \"illumina\" ){ return \"--illumina\" }\n  else if ( inputs.adapter1 == \"nextera\" ){ return \"--nextera\" }\n  else if ( inputs.adapter1 == \"small_rna\" ){ return \"--small_rna\" }\n  else { return null }\n}\n",
                    "position": 1
                },
                {
                    "prefix": "--adapter",
                    "valueFrom": "${\n  if ( inputs.apdater1 != null && inputs.adapter1 != \"illumina\" && inputs.adapter1 != \"nextera\" && inputs.adapter1 != \"small_rna\" ){\n    return inputs.adapter1\n  } else {\n    return null\n  }\n}\n",
                    "position": 1
                },
                {
                    "prefix": "--adapter2",
                    "valueFrom": "${\n  if ( inputs.fastq2 != null && inputs.apdater2 != null && inputs.adapter1 != \"illumina\" && inputs.adapter1 != \"nextera\" && inputs.adapter1 != \"small_rna\" ){\n    return inputs.adapter2\n  } else {\n    return null\n  }\n}\n",
                    "position": 1
                },
                {
                    "valueFrom": "${\n  if ( inputs.fastq2 == null ){ return null }\n  else { return \"--paired\" }\n}\n",
                    "position": 1
                },
                {
                    "valueFrom": "${\n  if ( inputs.fastq2 == null ){ return null }\n  else { return \"--retain_unpaired\" }\n}\n",
                    "position": 1
                },
                {
                    "prefix": "--length_1",
                    "valueFrom": "${\n  if ( inputs.fastq2 == null ){ return null }\n  else { return inputs.min_unpaired_read_rescue_length }\n}\n",
                    "position": 1
                },
                {
                    "prefix": "--length_2",
                    "valueFrom": "${\n  if ( inputs.fastq2 == null ){ return null }\n  else { return inputs.min_unpaired_read_rescue_length }\n}\n",
                    "position": 1
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "${\n    if ( inputs.fastq2 == null  ){ return \"*trimmed.fq*\" }\n    else { return \"*val_1.fq*\" }\n}\n"
                    },
                    "id": "#trim_galore.cwl/fastq1_trimmed"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "*unpaired_1.fq*"
                    },
                    "id": "#trim_galore.cwl/fastq1_trimmed_unpaired"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "*val_2.fq*"
                    },
                    "id": "#trim_galore.cwl/fastq2_trimmed"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "*unpaired_2.fq*"
                    },
                    "id": "#trim_galore.cwl/fastq2_trimmed_unpaired"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "outputBinding": {
                        "glob": "*trimming_report.txt"
                    },
                    "id": "#trim_galore.cwl/trim_galore_log"
                },
                {
                    "doc": "html report of post-trimming fastqc",
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "outputBinding": {
                        "glob": "*fastqc.html"
                    },
                    "id": "#trim_galore.cwl/trimmed_fastqc_html"
                },
                {
                    "doc": "all data of post-trimming fastqc e.g. figures",
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "outputBinding": {
                        "glob": "*fastqc.zip"
                    },
                    "id": "#trim_galore.cwl/trimmed_fastqc_zip"
                }
            ],
            "id": "#trim_galore.cwl"
        },
        {
            "class": "Workflow",
            "inputs": [
                {
                    "type": "File",
                    "id": "#bed_to_coverage_track.cwl/bed"
                },
                {
                    "type": "int",
                    "id": "#bed_to_coverage_track.cwl/bin_size"
                },
                {
                    "type": "long",
                    "id": "#bed_to_coverage_track.cwl/effective_genome_size"
                },
                {
                    "type": "File",
                    "id": "#bed_to_coverage_track.cwl/genome_info"
                },
                {
                    "doc": "List of space-delimited chromosome names that shall be ignored \nwhen calculating the scaling factor. \nSpecify as space-delimited string. \nDefault: \"chrX chrY chrM\"\n",
                    "type": [
                        "null",
                        "string"
                    ],
                    "default": "chrX chrY chrM",
                    "id": "#bed_to_coverage_track.cwl/ignoreForNormalization"
                }
            ],
            "steps": [
                {
                    "doc": "deeptools bamCoverage\n",
                    "run": "#deeptools_bamCoverage.cwl",
                    "in": [
                        {
                            "source": "#bed_to_coverage_track.cwl/indexing_bam/bam_sorted_indexed",
                            "id": "#bed_to_coverage_track.cwl/converting_bam_to_bigwig/bam"
                        },
                        {
                            "source": "#bed_to_coverage_track.cwl/bin_size",
                            "id": "#bed_to_coverage_track.cwl/converting_bam_to_bigwig/bin_size"
                        },
                        {
                            "source": "#bed_to_coverage_track.cwl/effective_genome_size",
                            "id": "#bed_to_coverage_track.cwl/converting_bam_to_bigwig/effective_genome_size"
                        },
                        {
                            "source": "#bed_to_coverage_track.cwl/ignoreForNormalization",
                            "id": "#bed_to_coverage_track.cwl/converting_bam_to_bigwig/ignoreForNormalization"
                        }
                    ],
                    "out": [
                        "#bed_to_coverage_track.cwl/converting_bam_to_bigwig/bigwig"
                    ],
                    "id": "#bed_to_coverage_track.cwl/converting_bam_to_bigwig"
                },
                {
                    "doc": "bedtools bedtobam - converts bed to bam;\nas most tools can handle bam as input but not always\nthe bed format(e.g. deeptools), moreover, bam is compressed;\ntherefore, it will be used as final output (instead of the bed file)\n",
                    "run": "#bedtools_bedtobam.cwl",
                    "in": [
                        {
                            "source": "#bed_to_coverage_track.cwl/bed",
                            "id": "#bed_to_coverage_track.cwl/converting_bed_to_bam/bed"
                        },
                        {
                            "source": "#bed_to_coverage_track.cwl/genome_info",
                            "id": "#bed_to_coverage_track.cwl/converting_bed_to_bam/genome_info"
                        }
                    ],
                    "out": [
                        "#bed_to_coverage_track.cwl/converting_bed_to_bam/bam"
                    ],
                    "id": "#bed_to_coverage_track.cwl/converting_bed_to_bam"
                },
                {
                    "doc": "samtools index - indexes sorted bam\n",
                    "run": "#samtools_index_hack.cwl",
                    "in": [
                        {
                            "source": "#bed_to_coverage_track.cwl/sorting_bam/bam_sorted",
                            "id": "#bed_to_coverage_track.cwl/indexing_bam/bam_sorted"
                        }
                    ],
                    "out": [
                        "#bed_to_coverage_track.cwl/indexing_bam/bam_sorted_indexed"
                    ],
                    "id": "#bed_to_coverage_track.cwl/indexing_bam"
                },
                {
                    "doc": "samtools sort - sorting of merged bam",
                    "run": "#samtools_sort.cwl",
                    "in": [
                        {
                            "source": "#bed_to_coverage_track.cwl/converting_bed_to_bam/bam",
                            "id": "#bed_to_coverage_track.cwl/sorting_bam/bam_unsorted"
                        }
                    ],
                    "out": [
                        "#bed_to_coverage_track.cwl/sorting_bam/bam_sorted"
                    ],
                    "id": "#bed_to_coverage_track.cwl/sorting_bam"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputSource": "#bed_to_coverage_track.cwl/indexing_bam/bam_sorted_indexed",
                    "id": "#bed_to_coverage_track.cwl/bam"
                },
                {
                    "type": "File",
                    "outputSource": "#bed_to_coverage_track.cwl/converting_bam_to_bigwig/bigwig",
                    "id": "#bed_to_coverage_track.cwl/bigwig"
                }
            ],
            "id": "#bed_to_coverage_track.cwl"
        },
        {
            "class": "Workflow",
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
            "inputs": [
                {
                    "type": {
                        "type": "array",
                        "items": "File"
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
            "steps": [
                {
                    "doc": "samtools view",
                    "run": "#samtools_view_filter.cwl",
                    "in": [
                        {
                            "source": "#merge_duprem_filter.cwl/sorting_merged_bam/bam_sorted",
                            "id": "#merge_duprem_filter.cwl/filter_by_mapq/bam"
                        },
                        {
                            "source": "#merge_duprem_filter.cwl/is_paired_end",
                            "id": "#merge_duprem_filter.cwl/filter_by_mapq/is_paired_end"
                        }
                    ],
                    "out": [
                        "#merge_duprem_filter.cwl/filter_by_mapq/bam_filtered"
                    ],
                    "id": "#merge_duprem_filter.cwl/filter_by_mapq"
                },
                {
                    "doc": "samtools flagstat on merged bams",
                    "run": "#samtools_flagstat.cwl",
                    "in": [
                        {
                            "source": "#merge_duprem_filter.cwl/remove_duplicates/bam_duprem",
                            "id": "#merge_duprem_filter.cwl/flagstat_duprem/bam"
                        }
                    ],
                    "out": [
                        "#merge_duprem_filter.cwl/flagstat_duprem/flagstat_output"
                    ],
                    "id": "#merge_duprem_filter.cwl/flagstat_duprem"
                },
                {
                    "doc": "samtools flagstat on quality filtered bams",
                    "run": "#samtools_flagstat.cwl",
                    "in": [
                        {
                            "source": "#merge_duprem_filter.cwl/filter_by_mapq/bam_filtered",
                            "id": "#merge_duprem_filter.cwl/flagstat_filtered/bam"
                        }
                    ],
                    "out": [
                        "#merge_duprem_filter.cwl/flagstat_filtered/flagstat_output"
                    ],
                    "id": "#merge_duprem_filter.cwl/flagstat_filtered"
                },
                {
                    "doc": "samtools flagstat on merged bams",
                    "run": "#samtools_flagstat.cwl",
                    "in": [
                        {
                            "source": "#merge_duprem_filter.cwl/sorting_merged_bam/bam_sorted",
                            "id": "#merge_duprem_filter.cwl/flagstat_merged/bam"
                        }
                    ],
                    "out": [
                        "#merge_duprem_filter.cwl/flagstat_merged/flagstat_output"
                    ],
                    "id": "#merge_duprem_filter.cwl/flagstat_merged"
                },
                {
                    "doc": "samtools index - indexes sorted bam\n",
                    "run": "#samtools_index_hack.cwl",
                    "in": [
                        {
                            "source": "#merge_duprem_filter.cwl/remove_duplicates/bam_duprem",
                            "id": "#merge_duprem_filter.cwl/indexing_duprem_bam/bam_sorted"
                        }
                    ],
                    "out": [
                        "#merge_duprem_filter.cwl/indexing_duprem_bam/bam_sorted_indexed"
                    ],
                    "id": "#merge_duprem_filter.cwl/indexing_duprem_bam"
                },
                {
                    "doc": "samtools merge - merging bam files of lane replicates",
                    "run": "#samtools_merge.cwl",
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
                    ],
                    "out": [
                        "#merge_duprem_filter.cwl/lane_replicate_merging/bam_merged"
                    ],
                    "id": "#merge_duprem_filter.cwl/lane_replicate_merging"
                },
                {
                    "doc": "fastqc - quality control for reads directly after mapping",
                    "run": "#fastqc.cwl",
                    "in": [
                        {
                            "source": "#merge_duprem_filter.cwl/indexing_duprem_bam/bam_sorted_indexed",
                            "id": "#merge_duprem_filter.cwl/qc_duprem/bam"
                        }
                    ],
                    "out": [
                        "#merge_duprem_filter.cwl/qc_duprem/fastqc_zip",
                        "#merge_duprem_filter.cwl/qc_duprem/fastqc_html"
                    ],
                    "id": "#merge_duprem_filter.cwl/qc_duprem"
                },
                {
                    "doc": "picard markdup - emoves duplicates from a single sorted bam file.",
                    "run": "#picard_markdup.cwl",
                    "in": [
                        {
                            "source": "#merge_duprem_filter.cwl/filter_by_mapq/bam_filtered",
                            "id": "#merge_duprem_filter.cwl/remove_duplicates/bam_sorted"
                        }
                    ],
                    "out": [
                        "#merge_duprem_filter.cwl/remove_duplicates/bam_duprem",
                        "#merge_duprem_filter.cwl/remove_duplicates/picard_markdup_log"
                    ],
                    "id": "#merge_duprem_filter.cwl/remove_duplicates"
                },
                {
                    "doc": "samtools sort - sorting of merged bam",
                    "run": "#samtools_sort.cwl",
                    "in": [
                        {
                            "source": "#merge_duprem_filter.cwl/lane_replicate_merging/bam_merged",
                            "id": "#merge_duprem_filter.cwl/sorting_merged_bam/bam_unsorted"
                        }
                    ],
                    "out": [
                        "#merge_duprem_filter.cwl/sorting_merged_bam/bam_sorted"
                    ],
                    "id": "#merge_duprem_filter.cwl/sorting_merged_bam"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "secondaryFiles": ".bai",
                    "outputSource": "#merge_duprem_filter.cwl/indexing_duprem_bam/bam_sorted_indexed",
                    "id": "#merge_duprem_filter.cwl/bam"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "outputSource": "#merge_duprem_filter.cwl/qc_duprem/fastqc_html",
                    "id": "#merge_duprem_filter.cwl/duprem_fastqc_html"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "outputSource": "#merge_duprem_filter.cwl/qc_duprem/fastqc_zip",
                    "id": "#merge_duprem_filter.cwl/duprem_fastqc_zip"
                },
                {
                    "type": "File",
                    "outputSource": "#merge_duprem_filter.cwl/flagstat_duprem/flagstat_output",
                    "id": "#merge_duprem_filter.cwl/duprem_flagstat_output"
                },
                {
                    "type": "File",
                    "outputSource": "#merge_duprem_filter.cwl/flagstat_filtered/flagstat_output",
                    "id": "#merge_duprem_filter.cwl/filtered_flagstat_output"
                },
                {
                    "type": "File",
                    "outputSource": "#merge_duprem_filter.cwl/flagstat_merged/flagstat_output",
                    "id": "#merge_duprem_filter.cwl/merged_flagstat_output"
                },
                {
                    "type": "File",
                    "outputSource": "#merge_duprem_filter.cwl/remove_duplicates/picard_markdup_log",
                    "id": "#merge_duprem_filter.cwl/picard_markdup_log"
                }
            ],
            "id": "#merge_duprem_filter.cwl"
        },
        {
            "class": "Workflow",
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
            "inputs": [
                {
                    "type": {
                        "type": "enum",
                        "symbols": [
                            "#trim_and_map.cwl/adapter1/nextera",
                            "#trim_and_map.cwl/adapter1/illumina",
                            "#trim_and_map.cwl/adapter1/small_rna",
                            "#trim_and_map.cwl/adapter1/auto"
                        ]
                    },
                    "id": "#trim_and_map.cwl/adapter1"
                },
                {
                    "type": {
                        "type": "enum",
                        "symbols": [
                            "#trim_and_map.cwl/adapter2/nextera",
                            "#trim_and_map.cwl/adapter2/illumina",
                            "#trim_and_map.cwl/adapter2/small_rna",
                            "#trim_and_map.cwl/adapter2/auto"
                        ]
                    },
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
                    "type": "File",
                    "secondaryFiles": [
                        ".fai",
                        "^.1.bt2",
                        "^.2.bt2",
                        "^.3.bt2",
                        "^.4.bt2",
                        "^.rev.1.bt2",
                        "^.rev.2.bt2"
                    ],
                    "id": "#trim_and_map.cwl/genome"
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
                }
            ],
            "steps": [
                {
                    "doc": "trim galore - adapter trimming using trim_galore",
                    "run": "#trim_galore.cwl",
                    "in": [
                        {
                            "source": "#trim_and_map.cwl/adapter1",
                            "id": "#trim_and_map.cwl/adaptor_trimming_and_qc_trimmed/adapter1"
                        },
                        {
                            "source": "#trim_and_map.cwl/adapter2",
                            "id": "#trim_and_map.cwl/adaptor_trimming_and_qc_trimmed/adapter2"
                        },
                        {
                            "source": "#trim_and_map.cwl/fastq1",
                            "id": "#trim_and_map.cwl/adaptor_trimming_and_qc_trimmed/fastq1"
                        },
                        {
                            "source": "#trim_and_map.cwl/fastq2",
                            "id": "#trim_and_map.cwl/adaptor_trimming_and_qc_trimmed/fastq2"
                        }
                    ],
                    "out": [
                        "#trim_and_map.cwl/adaptor_trimming_and_qc_trimmed/fastq1_trimmed",
                        "#trim_and_map.cwl/adaptor_trimming_and_qc_trimmed/fastq2_trimmed",
                        "#trim_and_map.cwl/adaptor_trimming_and_qc_trimmed/fastq1_trimmed_unpaired",
                        "#trim_and_map.cwl/adaptor_trimming_and_qc_trimmed/fastq2_trimmed_unpaired",
                        "#trim_and_map.cwl/adaptor_trimming_and_qc_trimmed/trim_galore_log",
                        "#trim_and_map.cwl/adaptor_trimming_and_qc_trimmed/trimmed_fastqc_html",
                        "#trim_and_map.cwl/adaptor_trimming_and_qc_trimmed/trimmed_fastqc_zip"
                    ],
                    "id": "#trim_and_map.cwl/adaptor_trimming_and_qc_trimmed"
                },
                {
                    "doc": "bowite2 - mapper, produces sam file",
                    "run": "#bowtie2.cwl",
                    "in": [
                        {
                            "source": "#trim_and_map.cwl/adaptor_trimming_and_qc_trimmed/fastq1_trimmed",
                            "id": "#trim_and_map.cwl/mapping/fastq1"
                        },
                        {
                            "source": "#trim_and_map.cwl/adaptor_trimming_and_qc_trimmed/fastq2_trimmed",
                            "id": "#trim_and_map.cwl/mapping/fastq2"
                        },
                        {
                            "source": "#trim_and_map.cwl/genome",
                            "id": "#trim_and_map.cwl/mapping/genome_index"
                        },
                        {
                            "source": "#trim_and_map.cwl/is_paired_end",
                            "id": "#trim_and_map.cwl/mapping/is_paired_end"
                        },
                        {
                            "source": "#trim_and_map.cwl/max_mapping_insert_length",
                            "id": "#trim_and_map.cwl/mapping/max_mapping_insert_length"
                        }
                    ],
                    "out": [
                        "#trim_and_map.cwl/mapping/sam",
                        "#trim_and_map.cwl/mapping/bowtie2_log"
                    ],
                    "id": "#trim_and_map.cwl/mapping"
                },
                {
                    "doc": "fastqc - quality control for trimmed fastq",
                    "run": "#fastqc.cwl",
                    "in": [
                        {
                            "source": "#trim_and_map.cwl/fastq1",
                            "id": "#trim_and_map.cwl/qc_raw/fastq1"
                        },
                        {
                            "source": "#trim_and_map.cwl/fastq2",
                            "id": "#trim_and_map.cwl/qc_raw/fastq2"
                        }
                    ],
                    "out": [
                        "#trim_and_map.cwl/qc_raw/fastqc_zip",
                        "#trim_and_map.cwl/qc_raw/fastqc_html"
                    ],
                    "id": "#trim_and_map.cwl/qc_raw"
                },
                {
                    "doc": "samtools view - convert sam to bam",
                    "run": "#samtools_view_sam2bam.cwl",
                    "in": [
                        {
                            "source": "#trim_and_map.cwl/mapping/sam",
                            "id": "#trim_and_map.cwl/sam2bam/sam"
                        }
                    ],
                    "out": [
                        "#trim_and_map.cwl/sam2bam/bam_unsorted"
                    ],
                    "id": "#trim_and_map.cwl/sam2bam"
                },
                {
                    "doc": "samtools sort - sorts unsorted bam file by coordinates.",
                    "run": "#samtools_sort.cwl",
                    "in": [
                        {
                            "source": "#trim_and_map.cwl/sam2bam/bam_unsorted",
                            "id": "#trim_and_map.cwl/sort_bam/bam_unsorted"
                        }
                    ],
                    "out": [
                        "#trim_and_map.cwl/sort_bam/bam_sorted"
                    ],
                    "id": "#trim_and_map.cwl/sort_bam"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputSource": "#trim_and_map.cwl/sort_bam/bam_sorted",
                    "id": "#trim_and_map.cwl/bam"
                },
                {
                    "type": "File",
                    "outputSource": "#trim_and_map.cwl/mapping/bowtie2_log",
                    "id": "#trim_and_map.cwl/bowtie2_log"
                },
                {
                    "type": "File",
                    "outputSource": "#trim_and_map.cwl/adaptor_trimming_and_qc_trimmed/fastq1_trimmed",
                    "id": "#trim_and_map.cwl/fastq1_trimmed"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputSource": "#trim_and_map.cwl/adaptor_trimming_and_qc_trimmed/fastq2_trimmed",
                    "id": "#trim_and_map.cwl/fastq2_trimmed"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "outputSource": "#trim_and_map.cwl/qc_raw/fastqc_html",
                    "id": "#trim_and_map.cwl/raw_fastqc_html"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "outputSource": "#trim_and_map.cwl/qc_raw/fastqc_zip",
                    "id": "#trim_and_map.cwl/raw_fastqc_zip"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "outputSource": "#trim_and_map.cwl/adaptor_trimming_and_qc_trimmed/trim_galore_log",
                    "id": "#trim_and_map.cwl/trim_galore_log"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "outputSource": "#trim_and_map.cwl/adaptor_trimming_and_qc_trimmed/trimmed_fastqc_html",
                    "id": "#trim_and_map.cwl/trimmed_fastqc_html"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "outputSource": "#trim_and_map.cwl/adaptor_trimming_and_qc_trimmed/trimmed_fastqc_zip",
                    "id": "#trim_and_map.cwl/trimmed_fastqc_zip"
                }
            ],
            "id": "#trim_and_map.cwl"
        },
        {
            "class": "Workflow",
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
            "inputs": [
                {
                    "doc": "Adapter to be trimmed from first reads. Cab be one of the following: \\n\n- \"nextera\" for the Nextera adapter (CTGTCTCTTATA)\\n\n- \"illumina\" for the Illumina universal adapter (AGATCGGAAGAGC)\\n\n- \"small_rna\" for the Illumina Small RNA 3-prime Adapter (TGGAATTCTCGG)\\n\n- \"auto\" to automatically detect the write setting\n",
                    "type": {
                        "type": "enum",
                        "symbols": [
                            "#main/adapter1/nextera",
                            "#main/adapter1/illumina",
                            "#main/adapter1/small_rna",
                            "#main/adapter1/auto"
                        ]
                    },
                    "id": "#main/adapter1"
                },
                {
                    "doc": "Adapters to be trimmed from second read. Cab be one of the following: \\n\n- \"nextera\" for the Nextera adapter (CTGTCTCTTATA)\\n\n- \"illumina\" for the Illumina universal adapter (AGATCGGAAGAGC)\\n\n- \"small_rna\" for the Illumina Small RNA 3-prime Adapter (TGGAATTCTCGG)\\n\n- \"auto\" to automatically detect the write setting\n",
                    "type": {
                        "type": "enum",
                        "symbols": [
                            "#main/adapter2/nextera",
                            "#main/adapter2/illumina",
                            "#main/adapter2/small_rna",
                            "#main/adapter2/auto"
                        ]
                    },
                    "id": "#main/adapter2"
                },
                {
                    "doc": "Bin size used for generation of coverage tracks. \nThe larger the bin size the smaller are the coverage tracks, however, \nthe less precise is the signal. For single bp resolution set to 1.\n",
                    "type": "int",
                    "default": 10,
                    "id": "#main/bin_size"
                },
                {
                    "doc": "The effectively mappable genome size, please see: \nhttps://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html\n",
                    "type": "long",
                    "id": "#main/effective_genome_size"
                },
                {
                    "doc": "List of fastq files containing the first mate of raw reads. \nMuliple files are provided if  multiplexing of the same library has been done \non multiple lanes. The reads comming from different fastq files are pooled \nafter alignment. Also see parameter \"fastq2\". \n",
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "id": "#main/fastq1"
                },
                {
                    "doc": "List of fastq files containing the second mate of raw reads. \nImportant: this list has to be of same length as parameter \"fastq1\".\n",
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "id": "#main/fastq2"
                },
                {
                    "doc": "Path to reference genome in fasta format. \nBowtie2 index files (\".1.bt2\", \".2.bt2\", ...) as well as a samtools index (\".fai\") \nhas to be located in the same directory.\\n\nAll of these files can be downloaded for the most common genome builds at  \nhttps://support.illumina.com/sequencing/sequencing_software/igenome.html. \nAlternatively, you can use \"bowtie2-build\" or \"samtools index\" to create them yourself.\n",
                    "type": "File",
                    "secondaryFiles": [
                        ".fai",
                        "^.1.bt2",
                        "^.2.bt2",
                        "^.3.bt2",
                        "^.4.bt2",
                        "^.rev.1.bt2",
                        "^.rev.2.bt2"
                    ],
                    "id": "#main/genome"
                },
                {
                    "doc": "Path to a tab-delimited file listing chromosome sizes in following fashion:\\n\n\"chromosome_name<tab>total_number_of_bp\".\\n\nFor the most common UCSC genome build, you can find corresponding files at: \nhttps://github.com/CompEpigen/ATACseq_workflows/tree/master/chrom_sizes. \nOr you can generate them yourself using UCSC script fetchChromSizes \n(http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/fetchChromSizes) in following fashion:\\n\n\"fetchChromSizes hg38 > hg38.chrom.sizes\".\\n\nIf you are dealing with a non-UCSC build, you can generate such a file from a samtools index using:\\n\n\"awk -v OFS='\\t' {'print $1,$2'} hg38.fa.fai > hg38.chrom.sizes\".\n",
                    "type": "File",
                    "id": "#main/genome_info"
                },
                {
                    "doc": "List of space-delimited chromosome names that shall be ignored \nwhen calculating the scaling factor. \nSpecify as space-delimited string. \nDefault: \"chrX chrY chrM\"\n",
                    "type": [
                        "null",
                        "string"
                    ],
                    "default": "chrX chrY chrM",
                    "id": "#main/ignoreForNormalization"
                },
                {
                    "doc": "Q-value cutoff used for peak calling by MACS2. \nThe default is 0.05.\n",
                    "type": "float",
                    "default": 0.05,
                    "id": "#main/macs2_qvalue"
                },
                {
                    "doc": "Maximum insert length between two reads of a pair. In case of ATACseq, \nvery long insert sizes are possible. So it is recommended to use at least \na value of 1500. However, please note that alignment will take significantly \nlonger for higher insert sizes. The default is 2500.\n",
                    "type": "long",
                    "default": 2500,
                    "id": "#main/max_mapping_insert_length"
                },
                {
                    "doc": "Sample ID used for naming the output files.\n",
                    "type": "string",
                    "id": "#main/sample_id"
                }
            ],
            "steps": [
                {
                    "doc": "bedtools bamtobed",
                    "run": "#bedtools_bamtobed_pe.cwl",
                    "in": [
                        {
                            "source": "#main/name_sorting_filtered_bam/bam_sorted",
                            "id": "#main/converting_bam_to_bedpe/bam"
                        }
                    ],
                    "out": [
                        "#main/converting_bam_to_bedpe/bedpe"
                    ],
                    "id": "#main/converting_bam_to_bedpe"
                },
                {
                    "doc": "multiqc summarizes the qc results from fastqc \nand other tools\n",
                    "run": "#multiqc_hack.cwl",
                    "in": [
                        {
                            "source": [
                                "#main/trim_and_map/bowtie2_log",
                                "#main/merge_duprem_filter/duprem_fastqc_zip",
                                "#main/merge_duprem_filter/duprem_fastqc_html",
                                "#main/generating_atac_signal_tags/frag_size_stats_tsv",
                                "#main/generating_atac_signal_tags/fragment_sizes_tsv",
                                "#main/generating_atac_signal_tags/filtering_stats_tsv",
                                "#main/peak_calling_macs2_broad/peaks_xls",
                                "#main/peak_calling_macs2_narrow/peaks_bed",
                                "#main/peak_calling_macs2_narrow/peaks_xls",
                                "#main/plot_fragment_size_distribution/frag_size_distr_tsv",
                                "#main/qc_plot_fingerprint/qc_plot_fingerprint_tsv",
                                "#main/qc_phantompeakqualtools/qc_phantompeakqualtools_stdout",
                                "#main/qc_phantompeakqualtools/qc_crosscorr_summary",
                                "#main/merge_duprem_filter/picard_markdup_log"
                            ],
                            "linkMerge": "merge_flattened",
                            "id": "#main/create_summary_qc_report/qc_files_array"
                        },
                        {
                            "source": [
                                "#main/trim_and_map/raw_fastqc_zip",
                                "#main/trim_and_map/raw_fastqc_html",
                                "#main/trim_and_map/trimmed_fastqc_html",
                                "#main/trim_and_map/trimmed_fastqc_zip",
                                "#main/trim_and_map/trim_galore_log",
                                "#main/peak_calling_macs2_broad/peaks_bed"
                            ],
                            "linkMerge": "merge_flattened",
                            "id": "#main/create_summary_qc_report/qc_files_array_of_array"
                        },
                        {
                            "source": "#main/sample_id",
                            "id": "#main/create_summary_qc_report/report_name"
                        }
                    ],
                    "out": [
                        "#main/create_summary_qc_report/multiqc_zip",
                        "#main/create_summary_qc_report/multiqc_html"
                    ],
                    "id": "#main/create_summary_qc_report"
                },
                {
                    "doc": null,
                    "run": "#generate_atac_signal_tags.cwl",
                    "in": [
                        {
                            "source": "#main/converting_bam_to_bedpe/bedpe",
                            "id": "#main/generating_atac_signal_tags/bedpe_alignm"
                        },
                        {
                            "source": "#main/sample_id",
                            "id": "#main/generating_atac_signal_tags/output_basename"
                        }
                    ],
                    "out": [
                        "#main/generating_atac_signal_tags/bed_tn5_center_29bp",
                        "#main/generating_atac_signal_tags/bed_tn5_center_73bp",
                        "#main/generating_atac_signal_tags/bed_tn5_center_200bp",
                        "#main/generating_atac_signal_tags/bed_tn5_center_1bp",
                        "#main/generating_atac_signal_tags/bed_tn5_center_fragment",
                        "#main/generating_atac_signal_tags/fragment_sizes_tsv",
                        "#main/generating_atac_signal_tags/filtering_stats_tsv",
                        "#main/generating_atac_signal_tags/frag_size_stats_tsv",
                        "#main/generating_atac_signal_tags/irreg_mappings_bedpe"
                    ],
                    "id": "#main/generating_atac_signal_tags"
                },
                {
                    "doc": null,
                    "run": "#bed_to_coverage_track.cwl",
                    "scatter": [
                        "#main/generating_coverage_tracks/bed"
                    ],
                    "scatterMethod": "dotproduct",
                    "in": [
                        {
                            "source": [
                                "#main/generating_atac_signal_tags/bed_tn5_center_29bp",
                                "#main/generating_atac_signal_tags/bed_tn5_center_73bp",
                                "#main/generating_atac_signal_tags/bed_tn5_center_200bp",
                                "#main/generating_atac_signal_tags/bed_tn5_center_1bp",
                                "#main/generating_atac_signal_tags/bed_tn5_center_fragment"
                            ],
                            "id": "#main/generating_coverage_tracks/bed"
                        },
                        {
                            "source": "#main/bin_size",
                            "id": "#main/generating_coverage_tracks/bin_size"
                        },
                        {
                            "source": "#main/effective_genome_size",
                            "id": "#main/generating_coverage_tracks/effective_genome_size"
                        },
                        {
                            "source": "#main/genome_info",
                            "id": "#main/generating_coverage_tracks/genome_info"
                        },
                        {
                            "source": "#main/ignoreForNormalization",
                            "id": "#main/generating_coverage_tracks/ignoreForNormalization"
                        }
                    ],
                    "out": [
                        "#main/generating_coverage_tracks/bigwig",
                        "#main/generating_coverage_tracks/bam"
                    ],
                    "id": "#main/generating_coverage_tracks"
                },
                {
                    "run": "#merge_duprem_filter.cwl",
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
                    ],
                    "out": [
                        "#main/merge_duprem_filter/duprem_fastqc_zip",
                        "#main/merge_duprem_filter/duprem_fastqc_html",
                        "#main/merge_duprem_filter/merged_flagstat_output",
                        "#main/merge_duprem_filter/filtered_flagstat_output",
                        "#main/merge_duprem_filter/duprem_flagstat_output",
                        "#main/merge_duprem_filter/picard_markdup_log",
                        "#main/merge_duprem_filter/bam"
                    ],
                    "id": "#main/merge_duprem_filter"
                },
                {
                    "doc": "samtools sort - sorting of filtered bam file by read name",
                    "run": "#samtools_sort_name.cwl",
                    "in": [
                        {
                            "source": "#main/merge_duprem_filter/bam",
                            "id": "#main/name_sorting_filtered_bam/bam_unsorted"
                        }
                    ],
                    "out": [
                        "#main/name_sorting_filtered_bam/bam_sorted"
                    ],
                    "id": "#main/name_sorting_filtered_bam"
                },
                {
                    "doc": "peak calling using macs2",
                    "run": "#macs2_callpeak_atac.cwl",
                    "scatter": [
                        "#main/peak_calling_macs2_broad/treatment_bed"
                    ],
                    "scatterMethod": "dotproduct",
                    "in": [
                        {
                            "valueFrom": "${return(true)}",
                            "id": "#main/peak_calling_macs2_broad/broad"
                        },
                        {
                            "source": "#main/effective_genome_size",
                            "id": "#main/peak_calling_macs2_broad/genome_size"
                        },
                        {
                            "source": "#main/macs2_qvalue",
                            "id": "#main/peak_calling_macs2_broad/qvalue"
                        },
                        {
                            "source": [
                                "#main/generating_atac_signal_tags/bed_tn5_center_73bp",
                                "#main/generating_atac_signal_tags/bed_tn5_center_200bp",
                                "#main/generating_atac_signal_tags/bed_tn5_center_fragment"
                            ],
                            "linkMerge": "merge_flattened",
                            "id": "#main/peak_calling_macs2_broad/treatment_bed"
                        }
                    ],
                    "out": [
                        "#main/peak_calling_macs2_broad/peaks_bed",
                        "#main/peak_calling_macs2_broad/peaks_xls"
                    ],
                    "id": "#main/peak_calling_macs2_broad"
                },
                {
                    "doc": "peak calling using macs2",
                    "run": "#macs2_callpeak_atac.cwl",
                    "in": [
                        {
                            "valueFrom": "${return(false)}",
                            "id": "#main/peak_calling_macs2_narrow/broad"
                        },
                        {
                            "source": "#main/effective_genome_size",
                            "id": "#main/peak_calling_macs2_narrow/genome_size"
                        },
                        {
                            "source": "#main/macs2_qvalue",
                            "id": "#main/peak_calling_macs2_narrow/qvalue"
                        },
                        {
                            "source": "#main/generating_atac_signal_tags/bed_tn5_center_29bp",
                            "id": "#main/peak_calling_macs2_narrow/treatment_bed"
                        }
                    ],
                    "out": [
                        "#main/peak_calling_macs2_narrow/peaks_bed",
                        "#main/peak_calling_macs2_narrow/peaks_xls"
                    ],
                    "id": "#main/peak_calling_macs2_narrow"
                },
                {
                    "run": "#plot_frag_size_distr.cwl",
                    "in": [
                        {
                            "source": "#main/generating_atac_signal_tags/fragment_sizes_tsv",
                            "id": "#main/plot_fragment_size_distribution/fragment_sizes_tsv"
                        },
                        {
                            "source": "#main/sample_id",
                            "id": "#main/plot_fragment_size_distribution/output_basename"
                        }
                    ],
                    "out": [
                        "#main/plot_fragment_size_distribution/frag_size_distr_plot",
                        "#main/plot_fragment_size_distribution/frag_size_distr_tsv"
                    ],
                    "id": "#main/plot_fragment_size_distribution"
                },
                {
                    "run": "#phantompeakqualtools.cwl",
                    "in": [
                        {
                            "source": "#main/merge_duprem_filter/bam",
                            "id": "#main/qc_phantompeakqualtools/bam"
                        }
                    ],
                    "out": [
                        "#main/qc_phantompeakqualtools/qc_crosscorr_summary",
                        "#main/qc_phantompeakqualtools/qc_crosscorr_plot",
                        "#main/qc_phantompeakqualtools/qc_phantompeakqualtools_stderr",
                        "#main/qc_phantompeakqualtools/qc_phantompeakqualtools_stdout"
                    ],
                    "id": "#main/qc_phantompeakqualtools"
                },
                {
                    "run": "#deeptools_plotFingerprint.cwl",
                    "in": [
                        {
                            "source": "#main/merge_duprem_filter/bam",
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
                    ],
                    "out": [
                        "#main/qc_plot_fingerprint/qc_plot_fingerprint_plot",
                        "#main/qc_plot_fingerprint/qc_plot_fingerprint_tsv",
                        "#main/qc_plot_fingerprint/qc_plot_fingerprint_stderr"
                    ],
                    "id": "#main/qc_plot_fingerprint"
                },
                {
                    "run": "#trim_and_map.cwl",
                    "scatter": [
                        "#main/trim_and_map/fastq1",
                        "#main/trim_and_map/fastq2"
                    ],
                    "scatterMethod": "dotproduct",
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
                            "source": "#main/genome",
                            "id": "#main/trim_and_map/genome"
                        },
                        {
                            "default": true,
                            "id": "#main/trim_and_map/is_paired_end"
                        },
                        {
                            "source": "#main/max_mapping_insert_length",
                            "id": "#main/trim_and_map/max_mapping_insert_length"
                        }
                    ],
                    "out": [
                        "#main/trim_and_map/raw_fastqc_zip",
                        "#main/trim_and_map/raw_fastqc_html",
                        "#main/trim_and_map/fastq1_trimmed",
                        "#main/trim_and_map/fastq2_trimmed",
                        "#main/trim_and_map/trim_galore_log",
                        "#main/trim_and_map/trimmed_fastqc_html",
                        "#main/trim_and_map/trimmed_fastqc_zip",
                        "#main/trim_and_map/bam",
                        "#main/trim_and_map/bowtie2_log"
                    ],
                    "id": "#main/trim_and_map"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "secondaryFiles": ".bai",
                    "outputSource": "#main/merge_duprem_filter/bam",
                    "id": "#main/bam"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "outputSource": "#main/generating_coverage_tracks/bam",
                    "id": "#main/bam_signal_tags"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "outputSource": "#main/generating_coverage_tracks/bigwig",
                    "id": "#main/bigwig_signal_tags"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "outputSource": "#main/trim_and_map/bowtie2_log",
                    "id": "#main/bowtie2_log"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "outputSource": "#main/merge_duprem_filter/duprem_fastqc_html",
                    "id": "#main/duprem_fastqc_html"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "outputSource": "#main/merge_duprem_filter/duprem_fastqc_zip",
                    "id": "#main/duprem_fastqc_zip"
                },
                {
                    "type": "File",
                    "outputSource": "#main/merge_duprem_filter/duprem_flagstat_output",
                    "id": "#main/duprem_flagstat_output"
                },
                {
                    "type": "File",
                    "outputSource": "#main/merge_duprem_filter/filtered_flagstat_output",
                    "id": "#main/filtered_flagstat_output"
                },
                {
                    "type": "File",
                    "outputSource": "#main/generating_atac_signal_tags/filtering_stats_tsv",
                    "id": "#main/filtering_stats_tsv"
                },
                {
                    "type": "File",
                    "outputSource": "#main/plot_fragment_size_distribution/frag_size_distr_plot",
                    "id": "#main/frag_size_distr_plot"
                },
                {
                    "type": "File",
                    "outputSource": "#main/plot_fragment_size_distribution/frag_size_distr_tsv",
                    "id": "#main/frag_size_distr_tsv"
                },
                {
                    "type": "File",
                    "outputSource": "#main/generating_atac_signal_tags/frag_size_stats_tsv",
                    "id": "#main/frag_size_stats_tsv"
                },
                {
                    "type": "File",
                    "outputSource": "#main/generating_atac_signal_tags/fragment_sizes_tsv",
                    "id": "#main/fragment_sizes_tsv"
                },
                {
                    "type": "File",
                    "outputSource": "#main/generating_atac_signal_tags/irreg_mappings_bedpe",
                    "id": "#main/irreg_mappings_bedpe"
                },
                {
                    "type": "File",
                    "outputSource": "#main/merge_duprem_filter/merged_flagstat_output",
                    "id": "#main/merged_flagstat_output"
                },
                {
                    "type": "File",
                    "outputSource": "#main/create_summary_qc_report/multiqc_html",
                    "id": "#main/multiqc_html"
                },
                {
                    "type": "File",
                    "outputSource": "#main/create_summary_qc_report/multiqc_zip",
                    "id": "#main/multiqc_zip"
                },
                {
                    "type": {
                        "type": "array",
                        "items": {
                            "type": "array",
                            "items": "File"
                        }
                    },
                    "outputSource": "#main/peak_calling_macs2_broad/peaks_bed",
                    "id": "#main/peaks_bed_macs2_broad"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "outputSource": "#main/peak_calling_macs2_narrow/peaks_bed",
                    "id": "#main/peaks_bed_macs2_narrow"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "outputSource": "#main/peak_calling_macs2_broad/peaks_xls",
                    "id": "#main/peaks_xls_macs2_broad"
                },
                {
                    "type": "File",
                    "outputSource": "#main/peak_calling_macs2_narrow/peaks_xls",
                    "id": "#main/peaks_xls_macs2_narrow"
                },
                {
                    "type": "File",
                    "outputSource": "#main/merge_duprem_filter/picard_markdup_log",
                    "id": "#main/picard_markdup_log"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputSource": "#main/qc_phantompeakqualtools/qc_crosscorr_plot",
                    "id": "#main/qc_crosscorr_plot"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputSource": "#main/qc_phantompeakqualtools/qc_crosscorr_summary",
                    "id": "#main/qc_crosscorr_summary"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputSource": "#main/qc_phantompeakqualtools/qc_phantompeakqualtools_stderr",
                    "id": "#main/qc_phantompeakqualtools_stderr"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputSource": "#main/qc_plot_fingerprint/qc_plot_fingerprint_plot",
                    "id": "#main/qc_plot_fingerprint_plot"
                },
                {
                    "type": "File",
                    "outputSource": "#main/qc_plot_fingerprint/qc_plot_fingerprint_stderr",
                    "id": "#main/qc_plot_fingerprint_stderr"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputSource": "#main/qc_plot_fingerprint/qc_plot_fingerprint_tsv",
                    "id": "#main/qc_plot_fingerprint_tsv"
                },
                {
                    "type": {
                        "type": "array",
                        "items": {
                            "type": "array",
                            "items": "File"
                        }
                    },
                    "outputSource": "#main/trim_and_map/raw_fastqc_html",
                    "id": "#main/raw_fastqc_html"
                },
                {
                    "type": {
                        "type": "array",
                        "items": {
                            "type": "array",
                            "items": "File"
                        }
                    },
                    "outputSource": "#main/trim_and_map/raw_fastqc_zip",
                    "id": "#main/raw_fastqc_zip"
                },
                {
                    "type": {
                        "type": "array",
                        "items": {
                            "type": "array",
                            "items": "File"
                        }
                    },
                    "outputSource": "#main/trim_and_map/raw_fastqc_zip",
                    "id": "#main/trim_galore_log"
                },
                {
                    "type": {
                        "type": "array",
                        "items": {
                            "type": "array",
                            "items": "File"
                        }
                    },
                    "outputSource": "#main/trim_and_map/trimmed_fastqc_html",
                    "id": "#main/trimmed_fastqc_html"
                },
                {
                    "type": {
                        "type": "array",
                        "items": {
                            "type": "array",
                            "items": "File"
                        }
                    },
                    "outputSource": "#main/trim_and_map/trimmed_fastqc_zip",
                    "id": "#main/trimmed_fastqc_zip"
                }
            ],
            "id": "#main"
        }
    ],
    "cwlVersion": "v1.0"
}