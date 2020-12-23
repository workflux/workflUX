{
    "$graph": [
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
            "class": "Workflow",
            "doc": "bla\n",
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
                    "type": "File",
                    "doc": "test",
                    "id": "#main/fastq1"
                },
                {
                    "type": "File",
                    "id": "#main/fastq2"
                }
            ],
            "steps": [
                {
                    "doc": "fastqc - quality control for trimmed fastq",
                    "run": "#fastqc.cwl",
                    "in": [
                        {
                            "source": "#main/fastq1",
                            "id": "#main/fastqc/fastq1"
                        },
                        {
                            "source": "#main/fastq2",
                            "id": "#main/fastqc/fastq2"
                        }
                    ],
                    "out": [
                        "#main/fastqc/fastqc_zip",
                        "#main/fastqc/fastqc_html"
                    ],
                    "id": "#main/fastqc"
                }
            ],
            "outputs": [
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "outputSource": "#main/fastqc/fastqc_html",
                    "id": "#main/fastqc_html"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "outputSource": "#main/fastqc/fastqc_zip",
                    "id": "#main/fastqc_zip"
                }
            ],
            "id": "#main"
        }
    ],
    "cwlVersion": "v1.0"
}