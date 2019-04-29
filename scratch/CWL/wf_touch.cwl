{
    "$graph": [
        {
            "class": "CommandLineTool",
            "hints": [
                {
                    "coresMin": 1,
                    "ramMin": 5000,
                    "tmpdirMin": 1000,
                    "class": "ResourceRequirement"
                }
            ],
            "baseCommand": [
                "touch"
            ],
            "inputs": [
                {
                    "type": "string",
                    "inputBinding": {
                        "position": 1
                    },
                    "id": "#touch_no_container.cwl/filename"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.filename)"
                    },
                    "id": "#touch_no_container.cwl/test_file"
                }
            ],
            "id": "#touch_no_container.cwl"
        },
        {
            "class": "Workflow",
            "inputs": [
                {
                    "type": "string",
                    "id": "#main/filename"
                }
            ],
            "steps": [
                {
                    "run": "#touch_no_container.cwl",
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
                    "type": "File",
                    "outputSource": "#main/touch/test_file",
                    "id": "#main/test_file"
                }
            ],
            "id": "#main"
        }
    ],
    "cwlVersion": "v1.0"
}