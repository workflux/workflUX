{
    "$graph": [
        {
            "class": "CommandLineTool",
            "hints": {
                "DockerRequirement": {
                    "dockerPull": "ubuntu:latest"
                },
                "ResourceRequirement": {
                    "coresMin": 1,
                    "ramMin": 5000,
                    "tmpdirMin": 1000
                }
            },
            "baseCommand": [
                "touch"
            ],
            "inputs": {
                "filename": {
                    "type": "string",
                    "inputBinding": {
                        "position": 1
                    }
                }
            },
            "outputs": {
                "test_file": {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.filename)"
                    }
                }
            },
            "id": "#touch.cwl"
        },
        {
            "class": "Workflow",
            "requirements": [
                {
                    "class": "ScatterFeatureRequirement"
                }
            ],
            "inputs": [
                {
                    "type": {
                        "type": "array",
                        "items": "string"
                    },
                    "id": "#main/filename"
                }
            ],
            "steps": [
                {
                    "scatter": [
                        "#main/touch/filename"
                    ],
                    "scatterMethod": "dotproduct",
                    "run": "#touch.cwl",
                    "in": [
                        {
                            "source": "#main/filename",
                            "id": "#main/touch/filename"
                        }
                    ],
                    "out": [
                        "#main/touch/test_file"
                    ],
                    "id": "#main/touch"
                }
            ],
            "outputs": [
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "outputSource": "#main/touch/test_file",
                    "id": "#main/test_file"
                }
            ],
            "id": "#main"
        }
    ],
    "cwlVersion": "v1.0"
}