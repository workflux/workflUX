
class ImportJanisFile extends React.Component{
    constructor(props){
        super(props);
        
        this.default_states = {
            wfFile: null,
            availableWfNames: [],
            wfNameSelectedForImport: null,
            translateToCWL: true,
            translateToWDL: false,
            importName: "",
            actionStatus: "none",
            serverMessages: []
        }

        this.state = this.default_states

        this.changeInputField = this.changeInputField.bind(this);
        this.upload = this.upload.bind(this);
        this.handleFileChange = this.handleFileChange.bind(this);
        this.handleBooleanChange = this.handleBooleanChange.bind(this);
    }

    changeInputField(event){
        let newState = {
            [event.currentTarget.name]: event.currentTarget.value,
            serverMessages: []
        }
        if (event.currentTarget.name == "wfNameSelectedForImport"){
            newState["importName"] = event.currentTarget.value
        }
        this.setState(newState)
    }

    handleFileChange(event){
        let states = this.default_states
        states[event.currentTarget.name] = event.currentTarget.files[0]
        this.setState(states)
    }

    handleBooleanChange(value, isSet){
        this.setState({
            [value]: isSet
        })
    }

    async upload(event){ // ajax request to server
        if (this.state.wfFile == null){
            let newState = this.default_states
            newState["serverMessages"] = [{type:"error", text:"No file selected."}]
            this.setState(newState)
        }
        else{
            let newState = event.currentTarget.value == "list_contained_workflows" ?
                ( 
                    this.default_states
                ) : (
                    {}
                )
            newState["actionStatus"] = event.currentTarget.value
            this.setState(newState)
            let formData = new FormData()
            formData.append("wf_file", this.state.wfFile)
            formData.append("meta", JSON.stringify({
                "import_name": this.state.importName,
                "wf_type": "janis",
                "translate_to_cwl": this.state.translateToCWL,
                "translate_to_wdl": this.state.translateToWDL,
                "wf_name_in_script": this.state.wfNameSelectedForImport,
                "access_token": await getUserInfo("accessToken")
            }))
            fetch(
                event.currentTarget.value == "list_contained_workflows" ? (
                    routeListAvailWfsInJanis
                ) :(
                    routeUploadWf
                ), 
                {
                    method: "POST",
                    body: formData,
                    cache: "no-cache"
                }
            ).then(res => res.json())
            .then(
                (result) => {
                    var messages = result.messages;
                    var data = result.data;
                    let errorOccured = false;
                    for( let i=0;  i<messages.length; i++){
                        if(messages[i].type == "error"){
                            errorOccured = true;
                            break;
                        }
                    }
                    if (errorOccured){
                        this.setState({actionStatus: "none", error: true, serverMessages: messages})
                    } 
                    else {
                        let newState = {
                            actionStatus: "none", 
                            error: false, 
                            serverMessages: messages
                        }
                        if (data.hasOwnProperty('avail_wfs')){
                            newState["availableWfNames"] = data.avail_wfs
                            newState["wfNameSelectedForImport"] = data.avail_wfs[0]
                            newState["importName"] = data.avail_wfs[0]
                        }
                        this.setState(newState);
                    }
                },
                (error) => {
                    // server could not be reached
                    var messages = [{
                        time: get_time_str(),
                        type: "error", 
                        text: serverNotReachableError
                    }];
                    this.setState({actionStatus: "none", error: true, serverMessages: messages});
                }
            )
        }
    }

    render(){
        return(
            <div className="w3-panel">
                <p>
                    Import a Janis workflow script.
                </p>
                <ExperimentalTag />
                <span className="w3-text-green">1. Choose a Janis script:</span>&nbsp;
                <input 
                    className="w3-button w3-border w3-border-grey"
                    style={ {display: "block"} }
                    type="file" 
                    name="wfFile" 
                    onChange={this.handleFileChange}
                    disabled={this.state.actionStatus != "none"}
                />
                <ActionButton 
                    name="list_contained_workflows"
                    value="list_contained_workflows"
                    label="list contained workflows"
                    loading={this.state.actionStatus == "list_contained_workflows"} 
                    onAction={this.upload}
                    forwardEvent={true}
                    disabled={this.state.actionStatus != "none" || demo}
                />
                {this.state.availableWfNames.length > 0 && (
                    <div>
                        <br/>
                        <span className="w3-text-green">2. Choose the workflow to import:</span><br/>
                        <select className="w3-button w3-white w3-border w3-padding-small" 
                            name="wfNameSelectedForImport"
                            onChange={this.changeInputField}
                            value={this.state.wfNameSelectedForImport}
                        >
                            {
                                this.state.availableWfNames.map((availWfName) => (
                                    <option 
                                        key={availWfName}
                                        value={availWfName}
                                    >
                                        {availWfName}
                                    </option>
                                ))
                            }
                        </select>
                        <br/><br/>
                        {/* <span className="w3-text-green">3. Choose at least one target format:</span><br/>
                        <Checkbox
                            name="translateToCWL"
                            value="translateToCWL"
                            checked={this.state.translateToCWL}
                            onChange={this.handleBooleanChange}
                        />
                        &nbsp; translate to CWL
                        <br/>
                        <Checkbox
                            name="translateToWDL"
                            value="translateToWDL"
                            checked={this.state.translateToWDL}
                            onChange={this.handleBooleanChange}
                        />
                        &nbsp; translate to WDL
                        { (! this.state.translateToCWL) && (! this.state.translateToWDL) ? (
                                <Message type="hint">
                                    <b>Attention:</b> Please select at least one of the two formats.
                                </Message>
                            ) : (
                                <span><br/><br/></span>
                            )
                        } */}
                        <span className="w3-text-green">3. Choose a name and import:</span><br/>
                        <input type="text"
                            className="w3-input w3-border"
                            name="importName"
                            style={ {width: "50%"} }
                            value={this.state.importName}
                            onChange={this.changeInputField}
                        />
                        <ActionButton 
                            name="import_wf"
                            value="import_wf"
                            label="import"
                            loading={this.state.actionStatus == "import_wf"} 
                            onAction={this.upload}
                            forwardEvent={true}
                            disabled={
                                this.state.actionStatus != "none" ||
                                ((! this.state.translateToCWL) && (! this.state.translateToWDL)) ||
                                demo
                            }
                        />
                    </div>
                )}
                <DisplayServerMessages messages={this.state.serverMessages} />
            </div>
        )
    }
}

class ImportCwlZip extends React.Component{
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

    async importCWLPath(){
        this.ajaxRequest({
            statusVar: "actionStatus",
            statusValueDuringRequest: "import",
            messageVar: "importMessages",
            sendData: {
                wf_path: this.state.cwlPath,
                import_name: this.state.importName
            },
            route: routeImportCwlByPathOrUrl
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
                    disabled={demo}
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
                            disabled={this.state.actionStatus != "none" || demo}
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
                                    disabled={this.state.actionStatus != "none" || demo}
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

class ImportWfFile extends React.Component{
    constructor(props){
        super(props);
        // props.wfType cwl or wdl

        this.state = {
            importName: "",
            actionStatus: "none",
            serverMessages: [],
            wfFile: null,
            wfImportsZip: null,
            provideWfImportsZip: false
        }

        this.changeInputField = this.changeInputField.bind(this);
        this.upload = this.upload.bind(this);
        this.handleFileChange = this.handleFileChange.bind(this);
        this.handleBooleanChange = this.handleBooleanChange.bind(this);
    }

    changeInputField(event){
        this.setState({
            [event.currentTarget.name]: event.currentTarget.value
        })
    }

    handleFileChange(event){
        this.setState({
            [event.currentTarget.name]: event.currentTarget.files[0],
            serverMessages: []
        })
    }

    handleBooleanChange(value, isSet){
        this.setState({
            [value]: isSet
        })
    }

    async upload(){ // ajax request to server
        if (this.state.wfFile == ""){
            this.state.serverMessages = [{type:"error", text:"No file selected."}]
            this.setState({actionStatus: "none", error: true})
        }
        else{
            this.setState({actionStatus:"uploading"})
            let formData = new FormData()
            formData.append("wf_file", this.state.wfFile)
            if (this.state.provideWfImportsZip){
                formData.append("wf_imports_zip", this.state.wfImportsZip)
            }
            formData.append("meta", JSON.stringify({
                "import_name": this.state.importName,
                "wf_type": this.props.wfType,
                "access_token": await getUserInfo("accessToken")
            }))
            fetch(routeUploadWf, {
                method: "POST",
                body: formData,
                cache: "no-cache"
            }).then(res => res.json())
            .then(
                (result) => {
                    var messages = result.messages;
                    var data = result.data;
                    let errorOccured = false;
                    for( let i=0;  i<messages.length; i++){
                        if(messages[i].type == "error"){
                            errorOccured = true;
                            break;
                        }
                    }
                    if (errorOccured){
                        this.setState({actionStatus: "none", error: true, serverMessages: messages})
                    } 
                    else {
                        this.setState({actionStatus: "none", error: false, serverMessages: messages});
                    }
                },
                (error) => {
                    // server could not be reached
                    var messages = [{
                        time: get_time_str(),
                        type: "error", 
                        text: serverNotReachableError
                    }];
                    this.setState({actionStatus: "none", error: true, serverMessages: messages});
                }
            )
        }
    }

    render(){
        return(
            <div className="w3-panel">
                <p>
                    {this.props.wfType == "CWL" ? (
                        "Import a CWL-wrapped tool or a packed CWL workflow."
                        ):(
                        <span>
                            Import a WDL workflow.   
                        </span>
                    )}
                </p>
                {this.props.wfType == "WDL" && (
                    <ExperimentalTag />
                )}
                <span className="w3-text-green">
                    1. Choose a {this.props.wfType} file:
                </span>&nbsp;   
                {this.props.wfType == "CWL" && (             
                    <Message type="hint">
                        <b>Please Note: CWL workflows are only supported in packed format (workflow with all contained tools).</b>&nbsp;
                        You may provide a ZIP file containing non-packed CWL workflows with all it's dependencies (see above "from ZIP file") or see&nbsp;
                        <a href="https://github.com/common-workflow-language/cwltool#combining-parts-of-a-workflow-into-a-single-document">
                            the documentation of cwltool
                        </a> for details on how to pack a workflow.
                    </Message>  
                )}  
                <input 
                    className="w3-button w3-border w3-border-grey"
                    style={ {display: "block"} }
                    type="file" 
                    name="wfFile" 
                    onChange={this.handleFileChange}
                    disabled={this.state.actionStatus != "none"}
                />
                
                {this.props.wfType == "WDL" && (
                    <div className="w3-container">
                        <br/>
                        <Checkbox
                            name="provideWfImportsZip"
                            value="provideWfImportsZip"
                            checked={this.state.provideWfImportsZip}
                            onChange={this.handleBooleanChange}
                        />
                        &nbsp; include WDL imports via a zip file
                        <br/>
                        {this.state.provideWfImportsZip && (
                            <span>
                                <span className="w3-text-green">
                                    Choose an imports.zip:
                                </span>&nbsp;   
                                <input 
                                    className="w3-button w3-border w3-border-grey"
                                    style={ {display: "block"} }
                                    type="file" 
                                    name="wfImportsZip" 
                                    onChange={this.handleFileChange}
                                    disabled={this.state.actionStatus != "none"}
                                />
                            </span>
                        )}
                    </div>
                )}
                <br/>
                
                <span className="w3-text-green">2. Choose a name:</span>&nbsp;
                <input type="text"
                    className="w3-input w3-border"
                    name="importName"
                    style={ {width: "50%"} }
                    value={this.state.importName}
                    onChange={this.changeInputField}
                />
                <ActionButton 
                    name="import"
                    value="import"
                    label="import"
                    loading={this.state.actionStatus == "uploading"} 
                    onAction={this.upload}
                    disabled={this.state.actionStatus != "none" || demo}
                />
                <DisplayServerMessages messages={this.state.serverMessages} />
            </div>
        )
    }
}

class ImportCwlUrl extends React.Component{
    constructor(props){
        super(props);
        this.state = {
            importName: "",
            cwlUrl: "",
            actionStatus: "none",
            importMessages: []
        }

        this.changeInputField = this.changeInputField.bind(this);
        this.importCWLUrl = this.importCWLUrl.bind(this);
        this.ajaxRequest = ajaxRequest.bind(this);
    }

    changeInputField(event){
        this.setState({
            [event.currentTarget.name]: event.currentTarget.value
        })
    }

    async importCWLUrl(){
        this.ajaxRequest({
            statusVar: "actionStatus",
            statusValueDuringRequest: "import",
            messageVar: "importMessages",
            sendData: {
                wf_path: this.state.cwlUrl,
                is_url: true,
                import_name: this.state.importName
            },
            route: routeImportCwlByPathOrUrl
        })
    }

    render(){
        return(
            <div className="w3-panel">
                <p>
                    Import a CWL-wrapped tool or CWL workflow (packed or not) from a public URL.
                </p>
                <span className="w3-text-green">1. Provide a URL to a CWL document:</span>&nbsp;
                <Message type="hint">
                    <table>
                        <tbody>
                            <tr>
                                <td>
                                    <i className="fab fa-github" style={ {fontSize:"48px"} }/>
                                </td>
                                <td>
                                    You can directly use a URL from github. However, please make sure to provide the URL to the raw file, e.g.:<br/>
                                    <a href="https://raw.githubusercontent.com/CompEpigen/ChIPseq_workflows/master/CWL/workflows/ACTseq.cwl">
                                    https://raw.githubusercontent.com/CompEpigen/ChIPseq_workflows/master/CWL/workflows/ACTseq.cwl
                                    </a>
                                </td>
                            </tr>
                        </tbody>
                    </table>
                </Message>
                <input type="text"
                    className="w3-input w3-border"
                    name="cwlUrl"
                    style={ {width: "50%"} }
                    value={this.state.cwlUrl}
                    onChange={this.changeInputField}
                />
                <br/>
                <span className="w3-text-green">2. Choose a name:</span>&nbsp;
                <input type="text"
                    className="w3-input w3-border"
                    name="importName"
                    style={ {width: "50%"} }
                    value={this.state.importName}
                    onChange={this.changeInputField}
                />
                <br/>
                <ActionButton
                    name="import"
                    value="import"
                    onAction={this.importCWLUrl}
                    label="import using selected name"
                    loading={this.state.actionStatus == "import"}
                    disabled={this.state.actionStatus != "none" || demo}
                />
                <DisplayServerMessages messages={this.state.importMessages} />
            </div>
        )
    }
}


class ImportTrsUri extends React.Component{
    constructor(props){
        super(props);
        this.state = {
            importName: "",
            trsUri: "",
            actionStatus: "none",
            importMessages: []
        }

        this.changeInputField = this.changeInputField.bind(this);
        this.importTrsUri = this.importTrsUri.bind(this);
        this.ajaxRequest = ajaxRequest.bind(this);
    }

    changeInputField(event){
        this.setState({
            [event.currentTarget.name]: event.currentTarget.value
        })
    }

    async importTrsUri(){
        this.ajaxRequest({
            statusVar: "actionStatus",
            statusValueDuringRequest: "import",
            messageVar: "importMessages",
            sendData: {
                trs_uri: this.state.trsUri,
                import_name: this.state.importName
            },
            route: routeImportCwlByTrsUri
        })
    }

    render(){
        return(
            <div className="w3-panel">
                <p>
                    Import a workflow via a Tool Repository Service (TRS) URI.
                </p>
                <span className="w3-text-green">1. Provide a TRS URI:</span>&nbsp;
                <Message type="hint">
                    <span>
                        TRS URIs look for instance like this:
                        "trs://trs-filer-test.c03.k8s-popup.csc.fi/gwas-fasp/versions/e214485"
                    </span>
                </Message>
                <input type="text"
                    className="w3-input w3-border"
                    name="trsUri"
                    style={ {width: "50%"} }
                    value={this.state.trsUri}
                    onChange={this.changeInputField}
                />
                <br/>
                <span className="w3-text-green">2. Choose a name:</span>&nbsp;
                <input type="text"
                    className="w3-input w3-border"
                    name="importName"
                    style={ {width: "50%"} }
                    value={this.state.importName}
                    onChange={this.changeInputField}
                />
                <br/>
                <ActionButton
                    name="import"
                    value="import"
                    onAction={this.importTrsUri}
                    label="import using selected name"
                    loading={this.state.actionStatus == "import"}
                    disabled={this.state.actionStatus != "none" || demo}
                />
                <DisplayServerMessages messages={this.state.importMessages} />
            </div>
        )
    }
}


class ImportCWLRoot extends React.Component {
    constructor(props){
        super(props);
        this.state = {
            importMethod: "trsUrl"
        }

        this.importMethods = {
            trsUrl: {
                descr: "Tool Repository Service (TRS) URI",
                component: <ImportTrsUri />
            },
            cwlUrl: {
                descr: "URL to public CWL document (e.g. from github)",
                component: <ImportCwlUrl />
            },
            cwlFile: {
                descr: "from (packed) CWL file",
                component: <ImportWfFile wfType="CWL"/>
            },
            // wdlFile: {
            //     descr: "from WDL file",
            //     component: <ImportWfFile wfType="WDL"/>
            // },
            cwlZip: {
                descr: "from ZIP file (e.g. a CWL workflow + its dependencies)",
                component: <ImportCwlZip />
            },
            janisFile: {
                descr: "from Janis file",
                component: <ImportJanisFile />
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
                {demo && (
                    <Message type="info">
                        <span className="w3-text-black">
                            <b>
                                Please note:
                            </b>
                            <p>
                                Since this is a demonstration instance &nbsp;
                                <b>upload of custom workflows is disabled</b>. 
                                Please use the pre-existing workflows.
                                <br/>
                                Therefore, head over to "Create Job" in the top bar.
                            </p>
                        </span>
                    </Message>
                )}
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