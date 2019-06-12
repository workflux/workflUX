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
    constructor(props){
        super(props)
        let runSelection = {}
        this.props.runs.map((r) =>
            runSelection[r] = false
        )
        this.state = {
            runSelection: runSelection
        }

        this.changeRunSelection = this.changeRunSelection.bind(this)
    }

    changeRunSelection(runID, is_checked){
        let update = {}
        update[runID] = is_checked
        this.setState({
            runSelection: Object.assign(this.state.runSelection, update)
        })
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
                <ActionButton
                    name="start"
                    value="start"
                    onAction={console.log}
                    label={<span><i className="fas fa-rocket w3-text-green"/>&nbsp;start</span>}
                    disabled={!is_job_selected}
                />
            </div>
        )
    }
}


class JobList extends React.Component {
    constructor(props) {
        super(props);
        this.state = {whichFocus: ""}; // no list item is focues by default
        this.changeFocus = this.changeFocus.bind(this);
    }

    changeFocus(newFocusValue){
        this.setState({whichFocus: newFocusValue});
    }

    render() {
        const itemValues = this.props.jobInfo.map( (job) => job.job_id);
        const itemNames = this.props.jobInfo.map( (job) => (
            <span key={job.job_id}>
                <i><span style={ {fontFamily: "courier"} }>
                    {(
                        job.job_id.substring(0,4) + "." +
                        job.job_id.substring(4,6) + "." + 
                        job.job_id.substring(6,8) + "/" + 
                        job.job_id.substring(9,12)
                    )}
                </span></i><br/>
                {
                    job.job_id.substring(13,job.job_id.length)
                }<br/>
                <IneditableValueField
                    backColorClass="w3-theme"
                >
                    {job.cwl_target}
                </IneditableValueField>
            </span>
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
        if(this.state.whichFocus != "") {
            let runs=[]
            this.props.jobInfo.map((job) =>
                job.job_id == this.state.whichFocus && (runs = job.runs)
            )
            itemContent = <JobContent jobId={this.state.whichFocus} runs={runs} />
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
                { data.length > 0 ? (
                        <JobList jobInfo={data} initMessages={messages}/>
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
