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
            "class": "Workflow",
            "inputs": [
                {
                    "type": "int",
                    "id": "#main/sleep_time"
                }
            ],
            "steps": [
                {
                    "run": "#sleep.cwl",
                    "in": [
                        {
                            "source": "#main/sleep_time",
                            "id": "#main/sleep/sleep_time"
                        }
                    ],
                    "out": [],
                    "id": "#main/sleep"
                }
            ],
            "outputs": [],
            "id": "#main"
        }
    ],
    "cwlVersion": "v1.0"
}