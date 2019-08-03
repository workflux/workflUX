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
                "sleep"
            ],
            "inputs": {
                "sleep_time": {
                    "type": "int",
                    "inputBinding": {
                        "position": 1
                    }
                }
            },
            "outputs": [],
            "id": "#sleep.cwl"
        },
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
                },
                {
                    "type": {
                        "type": "array",
                        "items": "int"
                    },
                    "id": "#main/sleep_time"
                }
            ],
            "steps": [
                {
                    "scatter": [
                        "#main/sleep/sleep_time"
                    ],
                    "scatterMethod": "dotproduct",
                    "run": "#sleep.cwl",
                    "in": [
                        {
                            "source": "#main/sleep_time",
                            "id": "#main/sleep/sleep_time"
                        }
                    ],
                    "out": [],
                    "id": "#main/sleep"
                },
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