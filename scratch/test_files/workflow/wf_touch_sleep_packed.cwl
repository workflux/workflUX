{
    "$graph": [
        {
            "class": "CommandLineTool",
            "requirements": {
                "InlineJavascriptRequirement": {}
            },
            "hints": {
                "ResourceRequirement": {
                    "coresMin": 1,
                    "ramMin": 100
                }
            },
            "baseCommand": [
                "sleep"
            ],
            "arguments": [
                {
                    "valueFrom": "${return(Math.random()*10)}",
                    "position": 1
                }
            ],
            "inputs": {
                "file": {
                    "type": "File"
                }
            },
            "outputs": [],
            "id": "#sleep.cwl"
        },
        {
            "class": "CommandLineTool",
            "hints": {
                "ResourceRequirement": {
                    "coresMin": 1,
                    "ramMin": 100
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
                "file": {
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
            "inputs": [
                {
                    "type": "string",
                    "id": "#main/filename"
                }
            ],
            "steps": [
                {
                    "run": "#sleep.cwl",
                    "in": [
                        {
                            "source": "#main/touch/file",
                            "id": "#main/sleep/file"
                        }
                    ],
                    "out": [],
                    "id": "#main/sleep"
                },
                {
                    "run": "#touch.cwl",
                    "in": [
                        {
                            "source": "#main/filename",
                            "id": "#main/touch/filename"
                        }
                    ],
                    "out": [
                        "#main/touch/file"
                    ],
                    "id": "#main/touch"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputSource": "#main/touch/file",
                    "id": "#main/file"
                }
            ],
            "id": "#main"
        }
    ],
    "cwlVersion": "v1.0"
}