class RunListElement extends React.Component {
    // Inputs:
    // props.runId
    // props.checked
    // props.onChange function to handle change
    //  takes 2 arguments: runId, is_checked
    constructor(props){
        super(props)
        this.state = {
            checked: false
        }
        this.handleChange = this.handleChange.bind(this)
    }

    handleChange(event){
        this.props.onChange(this.props.runId, event.target.checked)
    }

    render(){
        return (
            <div>
                <label>
                    <input 
                        type="checkbox" 
                        name="runs"
                        value={this.props.runId}
                        checked={this.props.checked}
                        onChange={this.handleChange}
                    />
                    {this.props.runId}

                </label>
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

    changeRunSelection(runId, is_checked){
        let update = {}
        update[runId] = is_checked
        this.setState({
            runSelection: Object.assign(this.state.runSelection, update)
        })
    }

    render(){
        return(
            <div className="w3-panel">
                {this.props.runs.map((r) =>(
                    <RunListElement
                        key={r}
                        runId={r}
                        onChange={this.changeRunSelection}
                        checked={this.state.runSelection[r]}
                    />
                ))}
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
