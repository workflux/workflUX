
class CreateJobButton extends React.Component {
    constructor(props){
        super(props);
        // props.jobId
        // props.sheet_format
        // props.disabled
        
        this.state = {
            job_creation_status: "none",
            jobCreationMessages: []
        }

        this.createJob = this.createJob.bind(this)
        this.ajaxRequest = ajaxRequest.bind(this)
    }

    createJob() {
        this.ajaxRequest({
            statusVar: "job_creation_status",
            statusValueDuringRequest: "in_progress",
            messageVar: "jobCreationMessages",
            sendData: {
                job_id: this.props.jobId,
                sheet_format: this.props.sheet_format //#! problematic: if format selector is changed after sheet was already submitted
            },
            route: routeCreateJob
        })
    }



    render(){
        return(
            <div>
                <ActionButton
                        name="create job"
                        value="create job"
                        label="create job"
                        disabled={this.props.disabled}
                        loading={this.state.job_creation_status != "none"}
                        onAction={this.createJob}
                />
                <DisplayServerMessages messages={this.state.jobCreationMessages} />
            </div>

        )
    }
}

class ParamField extends React.Component{
    constructor(props){
        super(props);
        // props.type
        // props.nullAllowed
        // props.onChange
        // props.paramValue
        // props.name
        // props.index for array params

        this.key = this.props.index ? (
            this.props.name + "%" + this.props.index.toString()
        ) :(
            this.props.name
        )

        this.disabled = (this.props.paramValue === "null")

        this.inputTypes = {
            "int":"input_number",
            "long":"input_number",
            "string":"input_text",
            "file":"input_text"
        }
        

    }

    handleChange(event){
        const paramValue = event.currentTarget.value
        this.props.onChange(this.props.key, )
    }

    toggleEnabled(event){
        const paramValue = this.disabled ? (
            ""
        ) : (
            "null"
        )
        this.props.onChange(this.props.key, )
    }

    render(){
        const inputType = this.inputTypes.hasOwnPropterty(this.props.type) ? (
            this.inputTypes[this.props.types]
        ) : (
            "input_text"
        )
        
        let input_field
        switch(inputType){
            case "input_number":
                input_field = (
                    <input
                        type="number"
                        name={"input_" + this.key}
                        value={this.paramValue}
                        onChange={this.handleChange}
                        required={true}
                        disabled={this.disabled}
                    />
                )
            case "input_text":
                input_field = (
                    <input
                        type="text"
                        name={"input_" + this.key}
                        value={this.paramValue}
                        onChange={this.handleChange}
                        required={true}
                        disabled={this.disabled}
                    />
                )
        }

        return(
            <span style={ {witeSpace: "nowrap"} }>
                <input
                    type="checkbox"
                    name={"enable_" + this.key}
                    value={"enable_" + this.key}
                    checked={!this.disabled}
                />
                {input_field}
            </span>

        )
    }
}


class JobParamFormHTML extends React.Component {
    constructor(props){
        super(props);
        // props.cwlTarget
        // props.param_modes
        // props.run_mode
        // props.run_names
        // props.jobId

        this.state = {
            action_status: "none",
            form_passed_validation: false,
            serverMessages: [],
            paramConfigs: {},
            paramValues: {},
            paramsByMode: {}
        }

        this.getParamValues = this.getParamValues.bind(this);
        this.ajaxRequest = ajaxRequest.bind(this)
    }

    componentDidMount(){
        // setup timer to automatically update
        this.getParamValues()
    }

    getParamValues(){
            this.ajaxRequest({
                statusVar: "action_status",
                statusValueDuringRequest: "loading",
                messageVar: "serverMessages",
                sendData: {
                    cwl_target: this.props.cwlTarget,
                    param_modes: this.props.param_modes,
                    run_mode: this.props.run_mode, 
                    run_names: this.props.run_names.filter((r) => r != "")
                },
                route: routeGetParamValues,
                onSuccess: (data, messages) => {
                    const paramNames = Object.keys(data.configs).filter((p) => data.configs[p].type != "helper")
                    let paramsByMode = {
                        global_single: [],
                        global_array: [],
                        run_single: [],
                        run_array: []
                    }
                    paramNames.map((p) =>{
                        let run_or_global = this.props.run_mode ? (
                                this.props.param_mode[p] ? ("run") : ("global")
                            ) : (
                                "global"
                            )
                        let single_or_array = data.configs[p].is_array ? ("array") : ("single")
                        let mode = run_or_global + "_" + single_or_array
                        paramsByMode[mode].push(p)
                    })

                    return({
                        paramConfigs: data.configs,
                        paramValues: data.param_values,
                        paramsByMode: paramsByMode
                    })
                }
            })

    }

    render() {
        const globalSingValueForm = {
            
        }


        return(
            <div className="w3-container">
                <DisplayServerMessages messages={this.state.serverMessages} />
                <h4>Configs:</h4>
                <p>
                    {JSON.stringify(this.state.paramConfigs)}
                </p>
                <h4>Param Values:</h4>
                <p>
                    {JSON.stringify(this.state.paramValues)}
                </p>
                <h4>Param by Mode:</h4>
                <p>
                    {JSON.stringify(this.state.paramsByMode)}
                </p>
                <h3>Global Parameters</h3>
                <span className="w3-text-green">Single Values:</span>
                {/* {globalSingValueForm} */}
                <span className="w3-text-green">Lists/arrays:</span>
                <h3>Run-Specific Parameters</h3>
                <span className="w3-text-green">Single Values:</span>
                <span className="w3-text-green">Lists/arrays:</span>

                <CreateJobButton
                    jobId={this.props.jobId}
                    sheet_format="xlsx"
                    disabled={!this.state.form_passed_validation}
                />
            </div>
        )
    }
}


class JobParamFormSpreadsheet extends React.Component {
    constructor(props){
        super(props);
        // props.cwlTarget
        // props.param_modes
        // props.run_mode
        // props.run_names
        // props.jobId

        this.state = {
            sheet_format: "xlsx",
            file_transfer_status: "none",
            form_passed_validation: false,
            sheetFormMessages: [],
        }

        this.changeSheetFormat = this.changeSheetFormat.bind(this);
        this.genFormSheet = this.genFormSheet.bind(this);
        this.handleFormSheetUpload = this.handleFormSheetUpload.bind(this)
        this.ajaxRequest = ajaxRequest.bind(this)
    }

    changeSheetFormat(event){
        this.setState({"sheet_format": event.currentTarget.value})
    }

    genFormSheet(){
        this.ajaxRequest({
            statusVar: "file_transfer_status",
            statusValueDuringRequest: "downloading",
            messageVar: "sheetFormMessages",
            sendData: {
                cwl_target: this.props.cwlTarget,
                param_modes: this.props.param_modes,
                run_mode: this.props.run_mode, 
                run_names: this.props.run_names.filter((r) => r != ""),
                job_id: this.props.jobId,
                sheet_format: this.state.sheet_format
            },
            route: routeGenParamFormSheet,
            onSuccess: (data, messages) => {
                window.location.href = data.get_form_sheet_href
                return({sheetFormMessages: []})
            }
        })
    }

    handleFormSheetUpload(isSuccess){
        this.setState({form_passed_validation: isSuccess})
    }

    render() {
        return(
            <div className="w3-container">
                <span className="w3-text-green">As spreadsheet form:</span>
                <ol>
                    <li>
                        export/download:
                        <select className="w3-button w3-white w3-border" 
                            name="sheet_format"
                            onChange={this.changeSheetFormat}
                            value={this.state.sheet_format}
                            >
                            <option value="xlsx">excel format (xlsx)</option>
                            <option value="xls">excel format (xls)</option>
                            <option value="ods">open office format (ods)</option>
                        </select> 
                        <ActionButton
                            name="export"
                            value="export"
                            label="export"
                            onAction={this.genFormSheet}
                            loading={this.state.file_transfer_status == "downloading"}
                            disabled={this.state.file_transfer_status != "none"}
                        />
                    </li>
                    <li>
                        open in excel or open office and fill in the form
                    </li>
                    <li>
                        <FileUploadComponent
                            requestRoute={routeSendFilledParamFormSheet}
                            instruction="import/upload"
                            oneLine={true}
                            disabled={this.state.file_transfer_status != "none"}
                            meta_data={this.props.jobId}
                            onUploadCompletion={this.handleFormSheetUpload}
                        />
                    </li>
                </ol>
                <DisplayServerMessages messages={this.state.sheetFormMessages} />
                <CreateJobButton
                    jobId={this.props.jobId}
                    sheet_format={this.state.sheet_format}
                    disabled={!this.state.form_passed_validation}
                />
            </div>
        )
    }
}


class JobCreationPrep extends React.Component {
    // Essential information are reqeusted by the user in order to
    // build a job template form. 
    // Following informations are needed:
    //  - which parameter shall be global and which run-specific
    //  - how many jobs (/samples) shall be created
    constructor(props) {
        super(props);
        // Inputs:
        // props.configData
        // props.cwlTarget
        let paramModes={} // container for param mode (is_run_specific true/false)
        this.props.configData.params.map((p) =>
            paramModes[p.param_name] = p.is_run_specific
        )
        this.state = {
            run_mode: false, 
            run_names: ["run1", "run2", "run3"],
            param_modes: paramModes,
            job_name: "new_job",
            display: "prep" // one of prep, form_ssheet, from_html
        }

        // construct job_id:
        const date = new Date()
        let year = date.getFullYear().toString()
        let month = date.getMonth()
        month = (month < 10) ? ("0" + month.toString()) : (month.toString())
        let day = date.getDate()
        day = (day < 10) ? ("0" + day.toString()) : (day.toString())
        const dateString = year + month + day
        let randomNumber = Math.round(Math.random()*1000)
        randomNumber = (randomNumber < 100) ? ("0" + randomNumber.toString()) : (randomNumber.toString())
        this.jobIdNum = dateString + "_" + randomNumber
        

        this.changeJobName = this.changeJobName.bind(this);
        this.changeParamMode = this.changeParamMode.bind(this);
        this.changeRunMode = this.changeRunMode.bind(this);
        this.changeRunNames = this.changeRunNames.bind(this);
        this.toggleParamForm = this.toggleParamForm.bind(this);
    }

    changeJobName(event){
        let jobNameString = event.currentTarget.value.trim().replaceAll(" ", "_")
        this.setState({"job_name": jobNameString})
    }

    changeParamMode(param_name, is_run_specific){
        let update = {}
        update[param_name] = is_run_specific
        this.setState({
            param_modes: Object.assign(this.state.param_modes, update)
        })
    }

    changeRunMode(value, new_bool){
        this.setState({"run_mode": new_bool})
    }

    changeRunNames(event){
        let runNameString = event.currentTarget.value
        runNameString = runNameString.replaceAll("\n", ",").replaceAll("\t", ",").replaceAll(";", ",")
        let runNames = runNameString.split(",")
        runNames = runNames.map((n) => n.trim().replace(" ", "_"))
        this.setState({"run_names": runNames})
    }

    toggleParamForm(value){
        this.setState({display: value})
    }

    render() {
        const paramTable = (
            <div style={ {maxHeight:"50vh", overflowY: "auto"} }>
                <table className="w3-table w3-bordered w3-border">
                    <thead className="w3-text-green">
                        <tr>
                            <th>Parameter</th>
                            <th>Type</th>
                            {this.state.run_mode ? (<th>Mode</th>) : (null)}
                        </tr>
                    </thead>
                    <tbody>
                        {this.props.configData.params.map( (p) => (
                            <tr key={p.param_name}> 
                                <td>{p.param_name}</td>
                                <td>{p.type}</td>
                                {this.state.run_mode ? (
                                        <td>
                                            global &nbsp;
                                            <BooleanSlider
                                                name="param_mode_select"
                                                value={p.param_name}
                                                onChange={this.changeParamMode}
                                                checked={this.state.param_modes[p.param_name]}
                                            />
                                            &nbsp; run-specific
                                        </td>
                                    ) : (null)
                                }
                            </tr>
                        ))}
                    </tbody>
                </table>
            </div>
        )

        const jobIdForm = (
            <div>
                <p>
                    <span className="w3-text-green">Job ID:</span>&nbsp;
                    <IneditableValueField>
                        {this.jobIdNum + "_" + this.state.job_name }
                    </IneditableValueField>
                </p>
                <p>
                    <label className="w3-text-green">Job name:</label><br/>
                    Please enter a job name (no whitespaces allowed).<br/>
                    <input type="text"
                        className="w3-input w3-border"
                        name="job_name"
                        value={this.state.job_name}
                        onChange={this.changeJobName}
                    />
                </p>
            </div>
        )

        const runNumberForm = (
            <div>
                <span className="w3-text-green">Runs per job:</span>&nbsp;
                single run &nbsp;
                <BooleanSlider
                    name="create_multi_run_job"
                    value="create_multi_run_job"
                    onChange={this.changeRunMode}
                    checked={this.run_mode}
                />
                &nbsp; multiple runs
                {this.state.run_mode ? (
                    <span>
                        <Message type="hint">Please choose parameter modes in the above table.</Message>
                        <label className="w3-text-green">Run names/IDs:</label><br/>
                        Please enter unique, comma-seperated IDs.
                        No whitespaces allowed (if inserted will be automatically converted to "_"). 
                        Hint: you may copy&paste cells from an excel file.
                        <textarea className="w3-input w3-border"
                            rows="2"
                            name="create_multi_run_job" 
                            value={this.state.run_names.join(", ").trim()}
                            onChange={this.changeRunNames}
                        />
                    </span>
                    ) : (null)
                }
            </div>
        )
        
        if (this.state.display == "prep"){
            return(
                <div>
                    <h3>Input parameters:</h3>
                    {paramTable}
                    <hr/>
                    <h3>Job ID:</h3>
                    {jobIdForm}
                    <hr/>
                    <h3>Number of runs:</h3>
                    {runNumberForm}
                    <hr/>
                    <h3>Provide Parameter Values:</h3>
                    <p>
                        To provide parameter values, the convenient HTML form can be used (recommended for most use cases). 
                        Alternatively, a spreadsheet can be used (Excel or OpenOffice) which is mainly for very large datasets (over 50 runs) or
                        when you would like to make use of the advanced parameter validation and manipulation options.
                    </p>
                    <p>
                        <span className="w3-text-green">Provide parameter using:</span>&nbsp;
                        <ActionButton
                            name="form_html"
                            value="form_html"
                            label={<span><i className="fas fa-list-alt"></i>&nbsp;HTML form</span>}
                            onAction={this.toggleParamForm}
                        />&nbsp; or &nbsp;
                        <ActionButton
                            name="form_ssheet"
                            value="form_ssheet"
                            label={<span><i className="fas fa-file-excel"></i>&nbsp;Spreadsheet</span>}
                            onAction={this.toggleParamForm}
                        />
                    </p>
                </div>
            )
        }
        else {
            return(
                <div>
                    <ActionButton
                        name="prep"
                        value="prep"
                        label={<span><i className="fas fa-caret-left"></i>&nbsp;back to job overview</span>}
                        onAction={this.toggleParamForm}
                    />
                    {this.state.display == "form_ssheet" ? (
                            <div>
                                <h2>Generate Parameter Form:</h2>
                                <JobParamFormSpreadsheet
                                    cwlTarget={this.props.cwlTarget}
                                    param_modes={this.state.param_modes}
                                    run_mode={this.state.run_mode}
                                    run_names={this.state.run_names}
                                    jobId={this.jobIdNum + "_" + this.state.job_name}
                                />
                            </div>
                        ) : (
                            <div>
                                <h2>Parameter Form:</h2>
                                <JobParamFormHTML
                                    cwlTarget={this.props.cwlTarget}
                                    param_modes={this.state.param_modes}
                                    run_mode={this.state.run_mode}
                                    run_names={this.state.run_names}
                                    jobId={this.jobIdNum + "_" + this.state.job_name}
                                />
                            </div>
                        )
                    }
                </div>
            )
        }
    }
}


class JobTemplConfigInfoAjax extends React.Component {
    // Request information on the job template configuration
    constructor(props) {
        super(props);
        // Inputs:
        // props.cwlTarget
        this.buildContentOnSuccess = this.buildContentOnSuccess.bind(this);
    }

    buildContentOnSuccess(data){ // when AJAX request succeeds
        return (<JobCreationPrep configData={data} cwlTarget={this.props.cwlTarget} />);
    }

    render() {
        return (
            <AjaxComponent
                key={this.props.cwlTarget}
                requestRoute={routeGetJobTemplConfigInfo}
                sendData={ {cwl_target: this.props.cwlTarget} }
                buildContentOnSuccess={this.buildContentOnSuccess}
                loaderSize="large"
                loaderMessage="Loading template infos."
            />
        );
    }
}

class JobTemplList extends React.Component {
    constructor(props) {
        super(props);
        this.state = {whichFocus: ""}; // no list item is focues by default
        this.changeFocus = this.changeFocus.bind(this);
    }

    changeFocus(newFocusValue){
        this.setState({whichFocus: newFocusValue});
    }

    render() {
        const itemValues = this.props.templFilesInfo.map( (tf) => tf.cwl_target);
        const itemNames = itemValues; // here no distinction between values and names neccessary
        let itemContent = (
            <div>
                <DisplayServerMessages messages={this.props.initMessages}/> 
                <p>
                    <i className="fas fa-arrow-left"></i>
                    Select a CWL document to make a new job.
                </p>
            </div>
        )
        if(this.state.whichFocus != "") {
            itemContent = <JobTemplConfigInfoAjax cwlTarget={this.state.whichFocus} />
        }

        return (
            <SideBarPanel
                label="CWL documents:"
                itemValues={itemValues}
                itemNames={itemNames}
                whichFocus={this.state.whichFocus}
                itemContent={itemContent}
                onChange={this.changeFocus}
            />
        );
    }
}

class CeateJobRoot extends React.Component {
    constructor(props) {
        super(props);
        this.buildContentOnSuccess = this.buildContentOnSuccess.bind(this);
    }

    buildContentOnSuccess(data, messages){ // when AJAX request succeeds
        return (
            <div style={ {height:"100%"} }>
                { data.length > 0 ? (
                        <JobTemplList templFilesInfo={data} initMessages={messages}/>
                    ) : (
                        <Message type="info">
                            No job templates found. Please import a CWL document.
                        </Message>
                    )
                }
            </div>
        );
    }

    render() {
        return (
            <AjaxComponent
                requestRoute={routeGetJobTemplList}
                sendData={ {} }
                buildContentOnSuccess={this.buildContentOnSuccess}
                loaderSize="large"
                loaderMessage="Loading available job templates"
            />
        );
    }
}
