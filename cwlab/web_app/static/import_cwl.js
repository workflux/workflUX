// class ImportButton extends React.Component {
//     constructor(props){
//         super(props);
//         // props.file
//         // props.pathOrUrl
//         // props.isUrl
//         // props.onStart
//         // props.onCompletion
//         // props.route
//     }



// }


class ImportCWLZip extends React.Component{
    constructor(props){
        super(props);
        this.state = {
            importName: "",
            actionStatus: "none",
            importMessages: [],
            extractedTmpDir: "",
            cwlPath: ""
        }

        this.changeInputField = this.changeInputField.bind(this);
        this.onUploadCompletion = this.onUploadCompletion.bind(this);
        this.onCWLSelection = this.onCWLSelection.bind(this);
        this.importCWLPath = this.importCWLPath.bind(this);
        this.browseContent = this.browseContent.bind(this);
        this.ajaxRequest = ajaxRequest.bind(this);
    }

    changeInputField(event){
        this.setState({
            [event.currentTarget.name]: event.currentTarget.value
        })
    }

    onUploadCompletion(isSuccess, data){
        if (isSuccess){
            this.setState({
                extractedTmpDir: data["temp_dir"],
                cwlPath: ""
            })
        }
        else {
            this.setState({
                extractedTmpDir: "",
                cwlPath: ""
            })
        }
    }

    onCWLSelection(changes, selectedItem){
        if (changes){
            this.setState({
                cwlPath: selectedItem,
                importName: selectedItem.split('\\').pop().split('/').pop().replace(".cwl", "").replace(".CWL",""),
                actionStatus: "none"
            })
        }
        else {
            this.setState({
                actionStatus: "none"
            })
        }
    }

    browseContent(){
        this.setState({actionStatus: "browse"})
    }

    importCWLPath(){
        this.ajaxRequest({
            statusVar: "actionStatus",
            statusValueDuringRequest: "import",
            messageVar: "importMessages",
            sendData: {
                cwl_path: this.state.cwlPath,
                import_name: this.state.importName
            },
            route: routeImportCwlByPath
        })
    }

    render(){
        const cwlBasename = this.state.cwlPath.split('\\').pop().split('/').pop()
        return(
            <div className="w3-panel">
                <p>
                    Import CWL document(s) via a zip file
                </p>
                <span className="w3-text-green">1. Choose a zip file containing CWL documents:</span>&nbsp;
                <FileUploadComponent
                    requestRoute={routeUploadCwlZip}
                    buttonLabel="extract content"
                    onUploadCompletion={this.onUploadCompletion}
                />
                {this.state.extractedTmpDir && (
                    <span>
                        
                        <span className="w3-text-green">
                            2. Browse the zip content and select a CWL document to import:
                        </span>
                        <br/>
                        <ActionButton
                            name="browse"
                            value="browse"
                            onAction={this.browseContent}
                            label="browse and select"
                            loading={this.state.actionStatus == "browse"}
                            disabled={this.state.actionStatus != "none"}
                        />&nbsp;
                        <IneditableValueField>
                            {this.state.cwlPath ? (
                                    cwlBasename
                                ) : (
                                    "no CWL document selected"
                                )
                            }
                        </IneditableValueField>
                        {this.state.actionStatus == "browse" && (
                            <BrowseDir
                                path={this.state.cwlPath}
                                fileExts={["cwl", "CWL"]}
                                allowInput={true}
                                fixedBaseDir={this.state.extractedTmpDir}
                                fixedBaseDirName={"ZIP_CONTENT"}
                                defaultBaseDir={this.state.extractedTmpDir}
                                selectDir={false}
                                includeTmpDir={true}
                                showCancelButton={true}
                                terminateBrowseDialog={this.onCWLSelection}
                            />
                        )}
                        {this.state.cwlPath && (
                            <span>
                                <br/>
                                <div className="w3-text-green" style={ {paddingTop: "10px"} }>3. Choose a name and import:</div>
                                <input type="text"
                                    className="w3-input w3-border"
                                    name="importName"
                                    style={ {width: "50%", display: "inline-block"} }
                                    value={this.state.importName}
                                    onChange={this.changeInputField}
                                />
                                <ActionButton
                                    name="import"
                                    value="import"
                                    onAction={this.importCWLPath}
                                    label="import using selected name"
                                    loading={this.state.actionStatus == "import"}
                                    disabled={this.state.actionStatus != "none"}
                                />
                                <DisplayServerMessages messages={this.state.importMessages} />
                            </span>
                        )}
                    </span>

                )}
            </div>
        )
    }

}

class ImportSingleCWL extends React.Component{
    constructor(props){
        super(props);
        this.state = {
            importName: "",
            actionStatus: "none",
            serverMessages: []
        }

        this.changeInputField = this.changeInputField.bind(this);
        // this.importCWL = this.importCWL.bind(this);
        // this.ajaxRequest = ajaxRequest.bind(this)
    }

    changeInputField(event){
        this.setState({
            [event.currentTarget.name]: event.currentTarget.value
        })
    }

    // importCWL(){
    //     this.ajaxRequest({
    //         statusVar: "actionStatus",
    //         statusValueDuringRequest: "importing",
    //         messageVar: "serverMessages",
    //         sendViaFormData=true,
    //         sendData: {
    //             job_id: this.props.jobId,
    //             run_id: this.props.runId
    //         },
    //         route: routeGetRunDetails,
    //         onSuccess: (data, messages) => {
    //             return({
    //                 logContent: data.log,
    //                 yamlContent: data.yaml,
    //                 doNotUpdate: !this.mounted
    //             })
    //         },
    //         onError: (messages) => {
    //             return({
    //                 doNotUpdate: !this.mounted
    //             })
    //         }
    //     })
    // }

    render(){
        return(
            <div className="w3-panel">
                <p>
                    Import a CWL-wrapped tool or a packed CWL workflow.
                </p>
                <span className="w3-text-green">1. Choose a name:</span>&nbsp;
                <input type="text"
                    className="w3-input w3-border"
                    name="importName"
                    style={ {width: "50%"} }
                    value={this.state.importName}
                    onChange={this.changeInputField}
                />
                <br/>
                <span className="w3-text-green">2. Choose and import a CWL file:</span>&nbsp;
                <Message type="hint">
                    Please Note: Currently, CWL workflows are only supported in packed format. 
                    Please see <a href="https://github.com/common-workflow-language/cwltool#combining-parts-of-a-workflow-into-a-single-document">the documentation of cwltool</a> for details.
                </Message>
                <FileUploadComponent
                    requestRoute={routeImportPackedCwl}
                    metaData={ {"import_name": this.state.importName} }
                />
            </div>
        )
    }

}

class ImportCWLRoot extends React.Component {
    constructor(props){
        super(props);
        this.state = {
            importMethod: "singleCWL"
        }

        this.importMethods = {
            file: {
                descr: "upload a single CWL document (CWL-wrapped tool or a packed CWL Workflow)",
                component: <ImportSingleCWL />
            },
            CWLZip: {
                descr: "upload a zip file containing CWL documents (e.g. a CWL workflow with its dependencies)",
                component: <ImportCWLZip />
            },
            public: {
                descr: "URL to public CWL document or CWL-containing zip ",
                component: <ImportCWLZip />
            }
        }

        this.changeInputField = this.changeInputField.bind(this);
    }

    changeInputField(event){
        this.setState({
            [event.currentTarget.name]: event.currentTarget.value
        })
    }

    render() {
        return(
            <div className="w3-panel">
                <h3>Import a Tool or Workflow</h3>
                <select className="w3-button w3-white w3-border w3-padding-small" 
                    name="importMethod"
                    onChange={this.changeInputField}
                    value={this.state.importMethod}
                >
                    {
                        Object.keys(this.importMethods).map((importMethod) =>
                            <option key={importMethod} value={importMethod}>
                                {this.importMethods[importMethod].descr}
                            </option>
                        )
                    }
                </select> 
                <hr/>
                {this.importMethods[this.state.importMethod].component}
            </div>
        );
    }
}
