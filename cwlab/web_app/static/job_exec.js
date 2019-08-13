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
        this.logStart.scrollIntoView({ behavior: behavior });
    }
    
    scrollToBottom(behavior="smooth") {
        this.logEnd.scrollIntoView({ behavior: behavior });
    }

    render(){
        return(
            <div>
                <div className="w3-cell-row">
                    <div className="w3-cell">
                        <h3>Execution Log:</h3>
                    </div>
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
                    style={ {whiteSpace: "pre-wrap", maxHeight:"50vh", overflowY: "auto"} }  
                >   
                    <div style={{ float:"left", clear: "both" }}
                        ref={(el) => { this.logStart = el; }}>
                    </div>
                    {this.props.logContent}
                    <div style={{ float:"left", clear: "both" }}
                        ref={(el) => { this.logEnd = el; }}>
                    </div>
                </div>
            </div>
        )
    }

}

class RunDetails extends React.Component {
    constructor(props) {
        super(props);
        // props.jobId
        // props.runId
        this.state = {
            actionStatus: "none",
            logEndPos: 0,
            logContent: "Loading",
            yamlContent: "Loading",
            actionMessages: [],
            autoRefresh: true
        }
        this.getRunDetails = this.getRunDetails.bind(this);
        this.toggleAutoRefresh = this.toggleAutoRefresh.bind(this);
    }

    componentDidMount(){
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
    }
  
    
    getRunDetails(){
        this.setState({actionStatus: "updating"})
        const sendData = {
            job_id: this.props.jobId,
            run_id: this.props.runId,
            log_end_pos: this.state.logEndPos
        }
        fetch(routeGetRunDetails, {
            method: "POST",
            body: JSON.stringify(sendData),
            headers: new Headers({
                'Content-Type': 'application/json'
            }),
            cache: "no-cache"
        }).then(res => res.json())
        .then(
            (result) => {
                this.messages = result.messages;
                let errorOccured = false;
                for( let i=0;  i<this.messages.length; i++){
                    if(this.messages[i].type == "error"){
                        errorOccured = true;
                        break;
                    }
                }
                if (! errorOccured){
                    // nothing just display messages
                    this.setState({
                        actionStatus: "none", 
                        logContent: result.data.log,
                        yamlContent: result.data.yaml,
                        logEndPos: result.data.log_end_pos, 
                        actionMessages: this.messages
                    }) 
                }
                else{
                    this.setState({actionStatus: "none", actionMessages: this.messages}) 
                }       
            },
            (error) => {
                // server could not be reached
                this.setState({actionStatus: "none", actionMessages: [{type: "error", text: serverNotReachableError}]}) 
            }
        )

    }

    toggleAutoRefresh(dummy, autoRefresh){
        console.log("peep")
        this.setState({autoRefresh: autoRefresh})
    }

    render() {
        return (
            <div>
                <h3>Input Parameters:</h3>
                <div 
                    className="w3-metro-darken w3-panel"
                    style={ {whiteSpace: "pre-wrap", maxHeight:"50vh", overflowY: "auto"} } 
                >
                    {this.state.yamlContent}
                </div>
                <RunDetailsLog 
                    logContent={this.state.logContent}
                    toggleAutoRefresh={this.toggleAutoRefresh}
                    autoRefresh={this.state.autoRefresh}
                />
            </div>
        );
    }
}

class RunListElement extends React.Component {
    // Inputs:
    // props.runId
    // props.checked
    // props.status
    // props.duration
    // props.execProfile
    // props.retryCount
    // props.onSelectionChange function to handle change
    //  takes 2 arguments: runId, is_checked
    // props.showRunDetails
    constructor(props){
        super(props)
        this.handleSelectionChange = this.handleSelectionChange.bind(this)
        this.get_status_color = this.get_status_color.bind(this)
        this.handleShowDetails = this.handleShowDetails.bind(this)
    }

    get_status_color(status){
        const statusColorClass = {
            "Loading": "w3-grey",
            "not started yet": "w3-grey",
            "queued": "w3-grey",
            "preparing for execution": "w3-amber",
            "executing": "w3-amber",
            "evaluating results": "w3-amber",
            "finishing": "w3-amber",
            "finished": "w3-green",
        }
        if (Object.keys(statusColorClass).includes(status)){
            return(statusColorClass[status])
        } else{
            return("w3-red")
        }
    }
    
    handleSelectionChange(event){
        this.props.onSelectionChange(this.props.runId, event.target.checked)
    }

    handleShowDetails(){
        this.props.showRunDetails(this.props.runId)
    }

    render(){
        let duration = ""
        if ( this.props.duration == "-" ){
            duration = "-"
        }
        else{
            if (this.props.duration[0] > 0){
                duration += this.props.duration[0].toString() + "d "
            }
            if (this.props.duration[1] > 0){
                duration += this.props.duration[1].toString() + "h "
            }
            duration += this.props.duration[2].toString() + "m "
        }

        return (
            <tr>
                <td>
                    <input 
                        type="checkbox" 
                        name="runs"
                        value={this.props.runId}
                        checked={this.props.checked}
                        onChange={this.handleSelectionChange}
                    />
                </td>
                <td>{this.props.runId}</td>
                <td className={this.get_status_color(this.props.status)}>
                    {this.props.status}
                    {
                        this.props.retryCount > 0 && 
                            "(retry: " + this.props.retryCount.toString() + "\")"
                    }
                </td>
                <td>
                    {duration}
                </td>
                <td>{this.props.execProfile}</td>
                <td>
                    <a onClick={this.handleShowDetails}><i className="fas fa-eye w3-button w3-text-green"></i></a>
                </td>
            </tr>
        )
    }
}

class RunList extends React.Component {
    constructor(props){
        super(props)
        // props.runsIds
        // props.jobId
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
        let runInfo = {}
        this.props.runIds.map( (r) => runInfo[r] = this.initRunInfo)
        this.messages = []
        this.state = {
            actionStatus: "none", 
            runInfo: runInfo,
            mirroredJobId: this.props.jobId
        }
        this.getRunInfo = this.getRunInfo.bind(this)
    }

    getRunInfo(){
        this.setState({actionStatus: "updating"})
        const sendData = {
            job_id: this.props.jobId,
            run_ids: this.props.runIds
        }
        fetch(routeGetRunStatus, {
            method: "POST",
            body: JSON.stringify(sendData),
            headers: new Headers({
                'Content-Type': 'application/json'
            }),
            cache: "no-cache"
        }).then(res => res.json())
        .then(
            (result) => {
                this.messages = result.messages;
                let errorOccured = false;
                for( let i=0;  i<this.messages.length; i++){
                    if(this.messages[i].type == "error"){
                        errorOccured = true;
                        break;
                    }
                }
                if (! errorOccured){
                    // nothing just display messages
                    this.setState({actionStatus: "none", runInfo: result.data}) 
                }
                else{
                    let runInfo = {}
                    this.props.runIds.map( (r) => runInfo[r] = this.errorRunInfo)
                    this.setState({actionStatus: "none", runInfo: runInfo}) 
                }       
            },
            (error) => {
                // server could not be reached
                this.actionMessages = [{type: "error", text: serverNotReachableError}];
                let runInfo = {}
                this.props.runIds.map( (r) => runInfo[r] = this.errorRunInfo)
                this.setState({actionStatus: "none", runInfo: runInfo}) 
            }
        )

    }

    componentDidMount(){
        this.getRunInfo()
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
    }

    static getDerivedStateFromProps(nextProps, prevState) {
        if(nextProps.jobId !== prevState.mirroredJobId){
            let runInfo = {}
            nextProps.runIds.map( (r) => runInfo[r] = {
                    status: "Loading", 
                    duration: "-", 
                        exec_profile: "-", 
                        retry_count: 0
                }
            )
            return({runInfo: runInfo, mirroredJobId: nextProps.jobId, actionStatus: "none"})
        }
        return(null)
    }

    // componentDidUpdate(){
    //     if(this.state.runInfo === {}){
    //         this.getRunInfo()
    //     }
    // }
    
    render(){
        return(
            <div style={ {maxHeight:"50vh", overflowY: "auto"} }>
                <table className="w3-table w3-bordered w3-border">
                    <thead className="w3-text-green">
                        <tr>
                            <th>
                                <i 
                                    className="fas fa-arrow-down" 
                                    style={ {paddingRight:"10px"} }
                                />
                                <ActionButton 
                                    name="select all"
                                    value="select all"
                                    label="all"
                                    onAction={this.props.toggelRunSelectionAll}
                                    colorClass="w3-black w3-text-green"
                                    smallPadding={true}
                                />
                            </th>
                            <th>Run Id</th>
                            <th>Status</th>
                            <th>Duration</th>
                            <th>Exec Profile</th>
                            <th>Details</th>
                        </tr>
                    </thead>
                    <tbody>
                        {this.props.runIds.map( (r) => (
                            <RunListElement 
                                key={r}
                                runId={r}
                                status={this.state.runInfo[r].status}
                                duration={this.state.runInfo[r].duration}
                                execProfile={this.state.runInfo[r].exec_profile}
                                retryCount={this.state.runInfo[r].retry_count}
                                checked={this.props.runSelection[r]}
                                onSelectionChange={this.props.changeRunSelection}
                                showRunDetails={this.props.showRunDetails}
                            />
                        ))}
                    </tbody>
                </table>
                <DisplayServerMessages messages={this.messages} />
            </div>
        )
    }
}


class JobContent extends React.Component {
    // Inputs:
    // props.runs list of run ids
    // props.jobId
    // props.cwlTarget
    // props.execProfiles
    constructor(props){
        super(props)
        let runSelection = {}
        this.props.runs.map((r) =>
            (runSelection[r] = false)
        )
        this.state = {
            runSelection: runSelection,
            actionStatus: "none",
            execProfile: this.props.execProfiles[0],
            showRunDetails: null
        }
        this.actionMessages = []

        this.changeRunSelection = this.changeRunSelection.bind(this)
        this.toggelRunSelectionAll = this.toggelRunSelectionAll.bind(this)
        this.execRuns = this.execRuns.bind(this)
        this.changeExecProfile = this.changeExecProfile.bind(this)
        this.showRunDetails = this.showRunDetails.bind(this)
        this.showRunList = this.showRunList.bind(this)
    }

    
    toggelRunSelectionAll(){
        // if all runs are selected, deselect all:
        const select = !Object.values(this.state.runSelection).every(Boolean) 
        let update = {}
        this.props.runs.map((r) =>
            (update[r] = select)
        )
        this.setState({
            runSelection: update
        })
    }

    changeRunSelection(runId, is_checked){
        let update = {}
        update[runId] = is_checked
        this.setState({
            runSelection: Object.assign(this.state.runSelection, update)
        })
    }

    changeExecProfile(event){
        this.setState({execProfile: event.currentTarget.value})
    }

    execRuns(){
        this.setState({actionStatus: "starting"})
        const runSelection = this.state.runSelection
        let selectedRuns = []
        this.props.runs.map((run) =>
            runSelection[run] && (selectedRuns.push(run))
        )
        const sendData = {
            cwl_target: this.props.cwlTarget,
            job_id: this.props.jobId,
            run_ids: selectedRuns,
            exec_profile: this.state.execProfile
        }
        fetch(routeStartExec, {
            method: "POST",
            body: JSON.stringify(sendData),
            headers: new Headers({
                'Content-Type': 'application/json'
            }),
            cache: "no-cache"
        }).then(res => res.json())
        .then(
            (result) => {
                this.actionMessages = result.messages;
                let errorOccured = false;
                for( let i=0;  i<this.actionMessages.length; i++){
                    if(this.actionMessages[i].type == "error"){
                        errorOccured = true;
                        break;
                    }
                }
                if (! errorOccured){
                    // nothing just display messages
                }    
                this.setState({actionStatus: "none"})        
            },
            (error) => {
                // server could not be reached
                this.actionMessages = [{type: "error", text: serverNotReachableError}];
                this.setState({actionStatus: "none"}) 
            }
        )
    }

    showRunDetails(runId){
        this.setState({showRunDetails:runId})
    }

    showRunList(){
        this.setState({showRunDetails:null})
    }

    render(){
        const is_run_selected=Object.values(this.state.runSelection).includes(true)
        if (this.state.showRunDetails == null){
            return(
                <div>
                    <h3>List of Runs:</h3>
                    <RunList 
                        jobId={this.props.jobId}
                        runIds={this.props.runs}
                        changeRunSelection={this.changeRunSelection}
                        toggelRunSelectionAll={this.toggelRunSelectionAll}
                        runSelection={this.state.runSelection}
                        showRunDetails={this.showRunDetails}
                    />
                    <i 
                        className="fas fa-arrow-up" 
                        style={ {paddingLeft:"20px", paddingRight:"10px"} }
                    />
                    select one or multiple jobs for following actions:
                    <h3>Actions:</h3>
                    <span className="w3-text-green">Start Execution:</span>
                    <div className="w3-container">
                        <label>
                            Select execution profile:
                            <select className="w3-button w3-white w3-border" 
                                name="exec_profile"
                                onChange={this.changeExecProfile}
                                value={this.state.execProfile}
                                >
                                {
                                    this.props.execProfiles.map((execProfile) =>
                                        <option key={execProfile} value={execProfile}>{execProfile}</option>
                                    )
                                }
                            </select> 
                        </label>
                        <br/>
                        <ActionButton
                            name="start"
                            value="start"
                            onAction={this.execRuns}
                            label={<span><i className="fas fa-rocket w3-text-green"/>&nbsp;start</span>}
                            disabled={!is_run_selected}
                        />
                    </div>
                    <DisplayServerMessages messages={this.actionMessages} />
                </div>
            )
        } else {
            return(
                <RunDetails 
                    jobId={this.props.jobId}
                    runId={this.state.showRunDetails}
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
        this.state = {whichFocus: ""}; // no list item is focues by default
        this.changeFocus = this.changeFocus.bind(this);
    }

    changeFocus(newFocusValue){
        this.setState({whichFocus: newFocusValue});
    }

    render() {
        const itemValues = this.props.jobInfo.map( (job) => job.job_id);
        const itemNames = this.props.jobInfo.map( (job) => (
            <p key={job.job_id}>
                <i style={ {fontFamily: "courier"} }>
                    {(
                        job.job_id.substring(0,4) + "." +
                        job.job_id.substring(4,6) + "." + 
                        job.job_id.substring(6,8) + "/" + 
                        job.job_id.substring(9,12)
                    )}
                </i><br/>
                {
                    job.job_id.substring(13,job.job_id.length)
                }<br/>
                <IneditableValueField
                    backColorClass="w3-theme"
                >
                    {job.cwl_target}
                </IneditableValueField>
            </p>
        ))
        let itemContent = (
            <div>
                <DisplayServerMessages messages={this.props.initMessages}/> 
                <p>
                    <i className="fas fa-arrow-left"></i>
                    Select a job.
                </p>
            </div>
        );
        if(this.state.whichFocus && this.state.whichFocus != "") {
            let jobInfo={}
            this.props.jobInfo.map((job) =>
                job.job_id == this.state.whichFocus && (jobInfo = job)
            )
            itemContent = <JobContent 
                jobId={this.state.whichFocus} 
                runs={jobInfo.runs} 
                cwlTarget={jobInfo.cwl_target} 
                execProfiles={this.props.execProfiles} 
            />
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

    buildContentOnSuccess(data, messages){ // when AJAX request succeeds
        return (
            <div>
                { data.jobs.length > 0 ? (
                        <JobList jobInfo={data.jobs} execProfiles={data.exec_profiles} initMessages={messages}/>
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
            />
        );
    }
}
