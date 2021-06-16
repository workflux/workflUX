class RunDetailsLog extends React.Component {
    constructor(props) {
        super(props);
        // props.logContent
        // props.toggleAutoRefresh
        // props.autoRefresh
        this.scrollToBottom = this.scrollToBottom.bind(this);
        this.scrollToTop = this.scrollToTop.bind(this);
    }

    componentDidMount() {
        this.scrollToTop("auto");
        this.scrollToBottom();
    }
  
    componentDidUpdate() {
        this.scrollToBottom();
    }
  
    scrollToTop(behavior="smooth") {
        this.logContentContainer.scrollTop = 0
    }
    
    scrollToBottom(behavior="smooth") {
        // const {logContentContainer} = this.refs
        this.logContentContainer.scrollTop = this.logContentContainer.scrollHeight - this.logContentContainer.clientHeight
        // this.logEnd.scrollIntoView({ behavior: behavior });
    }

    render(){
        return(
            <div>
                <div className="w3-cell-row">
                    <div className="w3-cell w3-right">
                        auto refresh: off &nbsp;
                        <BooleanSlider
                            name="toggle_auto_refresh"
                            value={"toggle_auto_refresh"}
                            onChange={this.props.toggleAutoRefresh}
                            checked={this.props.autoRefresh}
                        />
                        &nbsp; on

                    </div>
                </div>
                <div>
                    (Maximally the last {readMaxCharsFromFile} characters are shown. To scroll up disable auto refresh.)
                </div>
                <div 
                    className="w3-metro-darken w3-panel"
                    style={ {whiteSpace: "pre-wrap", maxHeight:"60vh", overflowY: "auto"} }  
                    ref={(el) => { this.logContentContainer = el; }}
                >
                    {this.props.logContent}
                </div>
            </div>
        )
    }

}

class RunDetails extends React.Component {
    constructor(props) {
        super(props);
        // props.jobName
        // props.runName
        // props.handleBack
        this.state = {
            actionStatus: "none",
            logContent: "Loading",
            yamlContent: "Loading",
            serverMessages: [],
            autoRefresh: true,
            whichFocus: "Input Parameters"
        }
        this.mounted = false
        this.getRunDetails = this.getRunDetails.bind(this);
        this.toggleAutoRefresh = this.toggleAutoRefresh.bind(this);
        this.handleBackButton = this.handleBackButton.bind(this);
        this.downloadFileOrFolder = this.downloadFileOrFolder.bind(this);
        this.changeFocus = this.changeFocus.bind(this);
        this.ajaxRequest = ajaxRequest.bind(this)
    }

    componentDidMount(){
        this.mounted = true
        // setup timer to automatically update
        this.getRunDetails()
        this.timerId = setInterval(
            () => {
                if(this.state.autoRefresh){
                    this.getRunDetails()
                }
            },
            autoRefreshInterval
          );
    }

    componentWillUnmount() {
        clearInterval(this.timerId);
        this.mounted = false
    }
  
    
    async getRunDetails(){
        this.ajaxRequest({
            statusVar: "actionStatus",
            statusValueDuringRequest: "updating",
            messageVar: "serverMessages",
            sendData: {
                job_name: this.props.jobName,
                run_name: this.props.runName
            },
            route: routeGetRunDetails,
            onSuccess: (data, messages) => {
                return({
                    logContent: data.log,
                    yamlContent: data.yaml,
                    doNotUpdate: !this.mounted
                })
            },
            onError: (messages) => {
                return({
                    doNotUpdate: !this.mounted
                })
            }
        })
    }

    downloadFileOrFolder(changes, selectedItem){

    }

    toggleAutoRefresh(dummy, autoRefresh){
        this.setState({autoRefresh: autoRefresh})
    }

    handleBackButton(event){
        this.props.handleBack()
    }

    changeFocus(tab){
        this.setState({whichFocus: tab})
    }

    render() {
        let content
        if(this.state.whichFocus == "Input Parameters"){
            content=(
                <div 
                    className="w3-metro-darken w3-panel"
                    style={ {whiteSpace: "pre-wrap", maxHeight:"60vh", overflowY: "auto"} } 
                >
                    {this.state.yamlContent}
                </div>
            )
        }
        else if(this.state.whichFocus == "Execution Log"){
            content=(
                <RunDetailsLog 
                    logContent={this.state.logContent}
                    toggleAutoRefresh={this.toggleAutoRefresh}
                    autoRefresh={this.state.autoRefresh}
                />
            )
        }
        else {
            content=( 
                <BrowseDir
                    allowDownload={true}
                    jobName={this.props.jobName}
                    runName={this.props.runName}
                    defaultBaseDir="OUTPUT_DIR_CURRENT_RUN"
                    terminateBrowseDialog={this.downloadFileOrFolder}
                    disableOnTop={true}
                />
            )
        }
        return (
            <div>
                <ActionButton
                    name="back"
                    value="back"
                    label={
                        <span><i className="fas fa-caret-left"></i>&nbsp;back to job overview</span>
                    }
                    onAction={this.handleBackButton}
                />
                <div
                    style={ {paddingTop: "10px", paddingBottom: "10px"} }
                >
                    <TabPanel
                        title=""
                        tabs={["Input Parameters", "Execution Log", "Output Files"]}
                        changeFocus={this.changeFocus}
                        whichFocus={this.state.whichFocus}
                    >
                    {content}
                </TabPanel>
                </div>
                <DisplayServerMessages messages={this.state.serverMessages} />
            </div>
        );
    }
}

class RunList extends React.Component {
    constructor(props){
        super(props)
        // props.runNames
        // props.jobName
        // props.changeRunSelection
        // props.toggelRunSelectionAll
        // props.runSelection
        // props.showRunDetails
        this.initRunInfo = {
            status: "Loading", 
            duration: "-", 
                exec_profile: "-", 
                retry_count: 0
        }
        this.errorRunInfo = {
            status: "could not connect to database", 
            duration: "-", 
                exec_profile: "-", 
                retry_count: 0
        }
        this.mounted = false;
        let runInfo = {}
        this.props.runNames.map( (r) => runInfo[r] = this.initRunInfo)
        this.state = {
            actionStatus: "none", 
            serverMessages: [],
            runInfo: runInfo
        }
        this.columnNames = {
            runName: "Run ID",
            status: "Status",
            duration: "duration",
            execProfile: "Exec. Profile",
            details: "Details"
        }

        this.statusClassName = {
            "Loading": "w3-grey",
            "not started yet": "w3-grey",
            "waiting to queue": "w3-grey",
            "queued": "w3-grey",
            "submitting": "w3-grey",
            "preparing": "w3-amber",
            "executing": "w3-amber",
            "evaluating": "w3-amber",
            "finalizing": "w3-amber",
            "finishing": "w3-amber",
            "finished": "w3-green",
        }

        this.nonProblematicStatusList = [
            "Loading",
            "not started yet",
            "waiting to queue",
            "queued",
            "preparing",
            "executing",
            "evaluating",
            "finalizing",
            "finished"
        ]

        this.getRunInfo = this.getRunInfo.bind(this)
        this.getDurationString = this.getDurationString.bind(this)
        this.getStatus = this.getStatus.bind(this)
        this.ajaxRequest = ajaxRequest.bind(this)
    }

    async getRunInfo(){
        this.ajaxRequest({
            statusVar: "actionStatus",
            statusValueDuringRequest: "updating",
            messageVar: "serverMessages",
            sendData: {
                job_name: this.props.jobName,
                run_names: this.props.runNames
            },
            route: routeGetRunStatus,
            onSuccess: (data, messages) => {
                return({runInfo: data, doNotUpdate: !this.mounted})
            },
            onError: (messages) => {
                let runInfo = {}
                this.props.runNames.map( (r) => runInfo[r] = this.errorRunInfo)
                return({runInfo: runInfo, doNotUpdate: !this.mounted})
            }
        })
    }

    async componentDidMount(){
        this.mounted = true
        await this.getRunInfo()
        // setup timer to automatically update
        this.timerId = setInterval(
            () => {
                this.getRunInfo()
            },
            autoRefreshInterval
          );
    }

    componentWillUnmount() {
        clearInterval(this.timerId);
        this.mounted = false
    }

    getStatus(runInfo, returnColorClass=false){
        if (runInfo.custom_status && this.nonProblematicStatusList.includes(runInfo.status)){
            if (returnColorClass){
                return("w3-" + runInfo.custom_status_color)
            }
            else {
                return(runInfo.custom_status)
            }
        }
        else {
            if (returnColorClass){
                if (Object.keys(this.statusClassName).includes(runInfo.status)){
                    return(this.statusClassName[runInfo.status])
                } else{
                    return("w3-red")
                }
            }
            else {
                return(runInfo.status)
            }
        }
    }

    getDurationString(duration){
        let durationString = ""
        if ( !duration || duration == "-" ){
            durationString = "-"
        }
        else{
            if (duration[0] > 0){
                durationString += duration[0].toString() + "d "
            }
            if (duration[1] > 0){
                durationString += duration[1].toString() + "h "
            }
            durationString += duration[2].toString() + "m "
        }
        return(durationString)
    }

    render(){
        let runInfo = {}
        this.props.runNames.map( (r) => (
            Object.keys(this.state.runInfo).includes(r) ? (
                runInfo[r] = this.state.runInfo[r]
            ) : (
                runInfo[r] = this.initRunInfo
            )
        ))

        const rowData = this.props.runNames.map( (r) => (
            {
                runName: {value: r},
                status: {
                    value: (
                        <span>
                            {this.getStatus(runInfo[r])}
                            {
                                runInfo[r].retry_count > 0 && 
                                    ' (retry: ' + runInfo[r].retry_count.toString() + ')'
                            }
                        </span>
                    ),
                    className: this.getStatus(runInfo[r], true)
                },
                duration: {value: this.getDurationString(runInfo[r].duration)},
                execProfile: {value: runInfo[r].exec_profile},
                details: {
                    value: (
                        <a
                            onClick={(event) => this.props.showRunDetails(r)}
                        >
                            <i className="fas fa-eye w3-button w3-text-green" />
                        </a>
                    )
                }    
            }
        ))

        return(
            <div>
                <DisplayServerMessages messages={this.state.serverMessages} />
                <Table
                    columnKeys={Object.keys(this.columnNames)}
                    columnNames={this.columnNames}
                    selectionEnabled={true}
                    handleSelectionChange={this.props.changeRunSelection}
                    selection={this.props.runSelection}
                    rowData={rowData}
                    selectionKey="runName"
                />
            </div>
        )
    }
}


class JobContent extends React.Component {
    // Inputs:
    // props.runs list of run ids
    // props.jobName
    // props.cwlTarget
    // props.execProfiles
    // props.execProfileParams
    // props.whichRunDetails
    // props.showRunDetails
    // props.showRunList
    // props.triggerJobListReload
    constructor(props){
        super(props)
        this.initState = {
            runNames: [],
            runSelection: [],
            actionStatus: "get_run_list",
            execProfile: this.props.execProfiles[0],
            parallelExec: this.props.execProfileParams[this.props.execProfiles[0]]["max_parallel_exec"],
            maxRetries: this.props.execProfileParams[this.props.execProfiles[0]]["max_retries"],
            maxParallelExec: this.props.execProfileParams[this.props.execProfiles[0]]["max_parallel_exec"],
            enableQueueing: this.props.execProfileParams[this.props.execProfiles[0]]["enable_queueing"],
            allowUserDecreaseMaxParallelExec: this.props.execProfileParams[this.props.execProfiles[0]]["allow_user_decrease_max_parallel_exec"],
            runDangerZoneUnlocked: false,
            globalDangerZoneUnlocked: false,
            serverMessages: [],
            actionRunExecMessages: [],
            actionRunDangerMessages: [],
            actionGlobalDangerMessages: []
        }
        this.state = this.initState
        this.actionMessages = []

        this.changeRunSelection = this.changeRunSelection.bind(this)
        this.execRuns = this.execRuns.bind(this)
        this.changeExecProfile = this.changeExecProfile.bind(this)
        this.toggleDangerZoneLock = this.toggleDangerZoneLock.bind(this)
        this.terminateRuns = this.terminateRuns.bind(this)
        this.ajaxRequest = ajaxRequest.bind(this)
        this.getRunList = this.getRunList.bind(this)
        this.deleteJob = this.deleteJob.bind(this)
    }

    componentDidMount(){
        this.getRunList()
    }

    componentDidUpdate(prevProps){
        if (prevProps.jobName != this.props.jobName){
            this.setState(this.initState)
            this.getRunList()
        }
    }

    toggleDangerZoneLock(value, unlocked){
        if (value == "unlock run danger zone"){
            this.setState({
                runDangerZoneUnlocked: unlocked
            })
        }
        else if (value == "unlock global danger zone"){
            this.setState({
                globalDangerZoneUnlocked: unlocked
            })
        }
    }

    changeRunSelection(newSelection){
        this.setState({
            runSelection: newSelection
        })
    }

    changeExecProfile(event){
        if(event.currentTarget.name == "exec_profile"){
            this.setState({
                execProfile: event.currentTarget.value,
                parallelExec: this.props.execProfileParams[event.currentTarget.value]["max_parallel_exec"],
                maxRetries: this.props.execProfileParams[event.currentTarget.value]["max_retries"],
                maxParallelExec: this.props.execProfileParams[event.currentTarget.value]["max_parallel_exec"],
                enableQueueing: this.props.execProfileParams[event.currentTarget.value]["enable_queueing"],
                allowUserDecreaseMaxParallelExec: this.props.execProfileParams[event.currentTarget.value]["allow_user_decrease_max_parallel_exec"]
            })
        }
        else if(event.currentTarget.name == "parallel_exec"){
            this.setState({
                parallelExec: event.currentTarget.value
            })
        }
    }

    async getRunList(){
        this.ajaxRequest({
            statusVar: "actionStatus",
            statusValueDuringRequest: "get_run_list",
            messageVar: "serverMessages",
            sendData: {
                job_name: this.props.jobName
            },
            route: routeGetRunList,
            onSuccess: (data, messages) => {
                return({
                    runNames: data.run_names, 
                    runSelection: []
                }) 
            }
        })
    }

    async execRuns(){
        const userInfo = await getUserInfo("all")
        if (loginEnabled && userInfo.expiresIn <= min_remaining_access_token_time_for_exec) {
            const response = confirm(`
                Your access token is about to expire.
                
                Please login in again to continue.
            `)

            if (response) {
                refreshAccessToken()
            }
        }
        else {
            this.ajaxRequest({
                statusVar: "actionStatus",
                statusValueDuringRequest: "starting",
                messageVar: "actionRunExecMessages",
                sendData: {
                    job_name: this.props.jobName,
                    run_names: this.state.runSelection,
                    exec_profile: this.state.execProfile,
                    parallel_exec: this.state.parallelExec
                },
                route: routeStartExec
            })
        }
    }

    async terminateRuns(mode="terminate"){
        this.ajaxRequest({
            statusVar: "actionStatus",
            statusValueDuringRequest: mode,
            messageVar: "actionRunDangerMessages",
            sendData: {
                job_name: this.props.jobName,
                run_names: this.state.runSelection,
                mode: mode
            },
            route: routeTerminateRuns,
            onSuccess: (data, messages) => {
                this.getRunList()       
                return({runDangerZoneUnlocked: false})
            },
            onError: (data, messages) => {
                this.getRunList()       
                return({runDangerZoneUnlocked: false})
            }
        })
    }

    async deleteJob(){
        this.ajaxRequest({
            statusVar: "actionStatus",
            statusValueDuringRequest: "delete_job",
            messageVar: "actionGlobalDangerMessages",
            sendData: {
                job_name: this.props.jobName
            },
            route: routeDeleteJob,
            onSuccess: (data, messages) => {
                this.props.triggerJobListReload()
                return({doNotUpdate: true}) 
            }
        })
    }

    render(){
        if (this.state.actionStatus == "get_run_list"){
            return(
                <LoadingIndicator 
                    message="Loading list of runs." 
                    size="large" 
                />
            )
        }
        else if (this.props.whichRunDetails == null){
            const is_run_selected= this.state.runSelection.length > 0
            const disable_run_actions = (! is_run_selected) || (this.state.actionStatus != "none")
            const disable_danger_run_actions = disable_run_actions || (! this.state.runDangerZoneUnlocked)
            return(
                <div>
                    <DisplayServerMessages messages={this.state.serverMessages} />
                    <h3>List of Runs:</h3>
                    <RunList 
                        jobName={this.props.jobName}
                        runNames={this.state.runNames}
                        changeRunSelection={this.changeRunSelection}
                        toggelRunSelectionAll={this.toggelRunSelectionAll}
                        runSelection={this.state.runSelection}
                        showRunDetails={this.props.showRunDetails}
                    />
                    <h3>Actions on Selected Runs:</h3>
                    {disable_run_actions &&
                        <Message type="hint">
                            Please select one or multiple runs in the above list to enable following actions.
                        </Message>
                    }
                    <div
                        style={
                            disable_run_actions ? ({opacity: 0.4}) : ({})
                        }
                    >  
                        <span className="w3-text-green">Start Execution:</span>
                        <div className="w3-panel">
                            <p>
                                <label>
                                    Select execution profile:
                                    <select className="w3-button w3-white w3-border w3-padding-small" 
                                        name="exec_profile"
                                        onChange={this.changeExecProfile}
                                        value={this.state.execProfile}
                                        disabled={disable_run_actions}
                                    >
                                        {
                                            this.props.execProfiles.map((execProfile) =>
                                                <option key={execProfile} value={execProfile}>{execProfile}</option>
                                            )
                                        }
                                    </select> 
                                </label>
                            </p>
                            {this.state.enableQueueing && this.state.allowUserDecreaseMaxParallelExec && (
                                <p>
                                    <label>
                                        Select how many runs may execute in parallel: &nbsp;
                                        <select className="w3-button w3-white w3-border w3-padding-small" 
                                            name="parallel_exec"
                                            onChange={this.changeExecProfile}
                                            value={this.state.parallelExec}
                                        >
                                            {
                                                [...Array(this.state.maxParallelExec).keys()].map((key) =>
                                                    <option key={key} value={key+1}>{key+1}</option>
                                                )
                                            }
                                        </select>
                                    </label>
                                </p>
                            )}
                            <p>
                                <ActionButton
                                    name="start"
                                    value="start"
                                    onAction={this.execRuns}
                                    label={<span><i className="fas fa-rocket w3-text-green"/>&nbsp;start</span>}
                                    disabled={disable_run_actions}
                                    loading={this.state.actionStatus == "starting"} 
                                />
                            </p>
                        </div >
                        <DisplayServerMessages messages={this.state.actionRunExecMessages} />
                        <span className="w3-text-red">Danger Zone:</span>
                        <div 
                            className="w3-panel"
                            style={ 
                                disable_danger_run_actions ? (
                                        {backgroundColor: "hsl(0, 20%, 50%)"}
                                    ) : (
                                        {backgroundColor: "hsl(0, 40%, 50%)"}
                                    )
                            }
                        >
                            <div className="w3-padding-16">
                                <BooleanSlider
                                    name="unlock run danger zone"
                                    value="unlock run danger zone"
                                    onChange={this.toggleDangerZoneLock}
                                    disabled={disable_run_actions}
                                    checked={this.state.runDangerZoneUnlocked}
                                /> &nbsp; unlock
                            </div>
                            <div
                                className="w3-padding-16"
                                style={
                                    disable_danger_run_actions && !disable_run_actions ? ({opacity: 0.4}) : ({})
                                }
                            >
                                <table className="w3-table">
                                    <tbody>
                                        <tr>
                                            <td>
                                                <ActionButton
                                                    name="terminate"
                                                    value="terminate"
                                                    onAction={this.terminateRuns}
                                                    label={<span><i className="fas fa-stop-circle w3-text-red"/>&nbsp;terminate</span>}
                                                    disabled={disable_danger_run_actions}
                                                    loading={this.state.actionStatus == "terminate"} 
                                                />
                                            </td>
                                            <td>
                                                Stop execution of selected runs. Intermediate results will be maintained.
                                                Runs will be marked as "terminated by user".
                                            </td>
                                        </tr>
                                        <tr>
                                            <td>
                                                <ActionButton
                                                    name="reset"
                                                    value="reset"
                                                    onAction={this.terminateRuns}
                                                    label={<span><i className="fas fa-undo w3-text-red"/>&nbsp;reset</span>}
                                                    disabled={disable_danger_run_actions}
                                                    loading={this.state.actionStatus == "reset"} 
                                                />
                                            </td>
                                            <td>
                                                Stop execution of selected runs and clear their (intermediate) results.
                                                Runs will be appear as "not started yet".
                                            </td>
                                        </tr>
                                        <tr>
                                            <td>
                                                <ActionButton
                                                    name="delete"
                                                    value="delete"
                                                    onAction={this.terminateRuns}
                                                    label={<span><i className="fas fa-trash-alt w3-text-red"/>&nbsp;delete</span>}
                                                    disabled={disable_danger_run_actions}
                                                    loading={this.state.actionStatus == "delete"} 
                                                />
                                            </td>
                                            <td>
                                                Stop execution of selected runs and deleted them entirely.
                                                They will no longer show up in the list of runs.
                                            </td>
                                        </tr>
                                    </tbody>
                                </table>
                            </div>
                        </div>
                    </div>
                    <DisplayServerMessages messages={this.state.actionRunDangerMessages} />
                    <h3>Global Actions:</h3><div>
                        <span className="w3-text-red">Danger Zone:</span>
                        <div 
                            className="w3-panel"
                            style={ 
                                this.state.globalDangerZoneUnlocked ? (
                                        {backgroundColor: "hsl(0, 40%, 50%)"}
                                    ) : (
                                        {backgroundColor: "hsl(0, 20%, 50%)"}
                                    )
                            }
                        >
                            <div className="w3-padding-16">
                                <BooleanSlider
                                    name="unlock global danger zone"
                                    value="unlock global danger zone"
                                    onChange={this.toggleDangerZoneLock}
                                    checked={this.state.globalDangerZoneUnlocked}
                                /> &nbsp; unlock
                            </div>
                            <div
                                className="w3-padding-16"
                            >
                                <table className="w3-table">
                                    <tbody>
                                        <tr>
                                            <td>
                                                <ActionButton
                                                    name="delete_job"
                                                    value="delete_job"
                                                    onAction={this.deleteJob}
                                                    label={<span><i className="fas fa-trash-alt w3-text-red"/>&nbsp;delete job</span>}
                                                    disabled={!this.state.globalDangerZoneUnlocked}
                                                    loading={this.state.actionStatus == "delete_job"} 
                                                />
                                            </td>
                                            <td>
                                                Delete the entire job. All runs will be stopped and deleted.
                                            </td>
                                        </tr>
                                    </tbody>
                                </table>
                            </div>
                        </div>
                    </div>
                    <DisplayServerMessages messages={this.state.actionGlobalDangerMessages} />
                </div>
            )
        } else {
            return(
                <RunDetails 
                    jobName={this.props.jobName}
                    runName={this.props.whichRunDetails}
                    handleBack={this.props.showRunList}
                />
            )
        }
    }
}


class JobList extends React.Component {
    constructor(props) {
        super(props);
        // props.jobInfo
        // props.execProfiles
        // props.execProfileParams
        // props.triggerReload
        this.state = {
            whichFocus: "",
            whichRunDetails: null
        }; // no list item is focues by default
        this.changeFocus = this.changeFocus.bind(this);
        this.showRunDetails = this.showRunDetails.bind(this);
        this.showRunList = this.showRunList.bind(this);
    }

    changeFocus(newFocusValue){
        this.setState({
            whichFocus: newFocusValue,
            whichRunDetails: null
        });
    }

    showRunDetails(runName){
        this.setState({whichRunDetails: runName})
    }

    showRunList(){
        this.setState({whichRunDetails: null})
    }


    render() {
        const itemValues = this.props.jobInfo.map( (job) => job.job_name);
        const itemNames = this.props.jobInfo.map( (job) => (
            <p key={job.job_name}>
                <i style={ {fontFamily: "courier"} }>
                    {(
                        job.job_name.substring(0,4) + "." +
                        job.job_name.substring(4,6) + "." + 
                        job.job_name.substring(6,8) + "/" + 
                        job.job_name.substring(9,12)
                    )}
                </i><br/>
                {
                    job.job_name.substring(13,job.job_name.length)
                }<br/>
                <IneditableValueField
                    backclassName="w3-theme"
                >
                    {job.wf_name}
                </IneditableValueField>
            </p>
        ))

        let itemContent
        if(this.state.whichFocus && this.state.whichFocus != "") {
            let jobInfo={}
            this.props.jobInfo.map((job) =>
                job.job_name == this.state.whichFocus && (jobInfo = job)
            )
            itemContent = <JobContent 
                jobName={this.state.whichFocus} 
                cwlTarget={jobInfo.wf_target} 
                execProfiles={this.props.execProfiles}
                execProfileParams={this.props.execProfileParams}
                whichRunDetails={this.state.whichRunDetails}
                showRunDetails={this.showRunDetails}
                showRunList={this.showRunList}
                triggerJobListReload={this.props.triggerReload}
            />
        }
        else {
            itemContent = (
                <div>
                    <DisplayServerMessages messages={this.props.initMessages}/> 
                    <p>
                        <i className="fas fa-arrow-left"></i>
                        Select a job.
                    </p>
                </div>
            );
        }

        return (
            <SideBarPanel
                label="Jobs:"
                itemValues={itemValues}
                itemNames={itemNames}
                whichFocus={this.state.whichFocus}
                itemContent={itemContent}
                onChange={this.changeFocus}
            />
        );
    }
}

class JobExecRoot extends React.Component {
    constructor(props) {
        super(props);
        this.buildContentOnSuccess = this.buildContentOnSuccess.bind(this);
    }

    buildContentOnSuccess(data, messages, triggerReload){ // when AJAX request succeeds
        return (
            <div>
                { data.jobs.length > 0 ? (
                        <JobList 
                            jobInfo={data.jobs} 
                            execProfiles={data.exec_profiles} 
                            execProfileParams={data.exec_profile_params} 
                            initMessages={messages}
                            triggerReload={triggerReload}
                        />
                    ) : (
                        <Message type="info">
                            No jobs found.
                        </Message>
                    )
                }
            </div>
        );
    }

    render() {
        return (
            <AjaxComponent
                requestRoute={routeGetJobList}
                sendData={ {} }
                buildContentOnSuccess={this.buildContentOnSuccess}
                loaderSize="large"
                loaderMessage="Loading available job templates"
                suppressMessages={false}
            />
        );
    }
}
