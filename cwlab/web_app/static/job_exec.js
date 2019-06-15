class RunListElement extends React.Component {
    // Inputs:
    // props.runID
    // props.checked
    // props.onSelectionChange function to handle change
    //  takes 2 arguments: runID, is_checked
    constructor(props){
        super(props)
        this.state = {
            checked: false
        }
        this.handleSelectionChange = this.handleSelectionChange.bind(this)
        this.statusColorClass = {
            "not started yet": "w3-grey",
            "running": "w3-amber",
            "finished": "w3-green",
            "failed": "w3-failed"
        }
    }

    handleSelectionChange(event){
        this.props.onSelectionChange(this.props.runID, event.target.checked)
    }

    render(){
        return (
            <tr>
                <td>
                    <input 
                            type="checkbox" 
                            name="runs"
                            value={this.props.runID}
                            checked={this.props.checked}
                            onChange={this.handleSelectionChange}
                    />
                </td>
                <td>{this.props.runID}</td>
                <td className={this.statusColorClass[this.props.status]}>
                    {this.props.status}
                </td>
                <td>{this.props.duration}</td>
                <td>{this.props.execType}</td>
                <td>
                    <a><i className="fas fa-eye w3-button w3-text-green"></i></a>
                </td>
            </tr>
        )
    }
}

class RunList extends React.Component {
    constructor(props){
        super(props)
        // props.runsIDs
        // props.changeRunSelection
    }

    render(){
        return(
            <div style={ {maxHeight:"50vh", overflowY: "auto"} }>
                <table className="w3-table w3-bordered w3-border">
                    <thead className="w3-text-green">
                        <tr>
                            <th></th>
                            <th>Run ID</th>
                            <th>Status</th>
                            <th>Duration</th>
                            <th>Exec Type</th>
                            <th>Details</th>
                        </tr>
                    </thead>
                    <tbody>
                        {this.props.runIDs.map( (r) => (
                            <RunListElement 
                                key={r}
                                runID={r}
                                status="not started yet"
                                duration="-"
                                execType="-"
                                onSelectionChange={this.props.changeRunSelection}
                            />
                        ))}
                    </tbody>
                </table>
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
            execProfile: this.props.execProfiles[0]
        }
        this.actionMessages = []

        this.changeRunSelection = this.changeRunSelection.bind(this)
        this.execRuns = this.execRuns.bind(this)
        this.changeExecProfile = this.changeExecProfile.bind(this)
    }

    changeRunSelection(runID, is_checked){
        let update = {}
        update[runID] = is_checked
        this.setState({
            runSelection: Object.assign(this.state.runSelection, update)
        })
    }

    changeExecProfile(execProfile){
        this.setState({execProfile: execProfile})
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

    render(){
        const is_job_selected=Object.values(this.state.runSelection).includes(true)
        
        return(
            <div>
                <h3>List of Runs:</h3>
                <RunList 
                    runIDs={this.props.runs}
                    changeRunSelection={this.changeRunSelection}
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
                        disabled={!is_job_selected}
                    />
                </div>
                <DisplayServerMessages messages={this.actionMessages} />
            </div>
        )
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
