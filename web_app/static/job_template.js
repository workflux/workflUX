
class JobTemplateFormPreparationPanel extends React.Component {
    // Essential information are reqeusted by the user in order to
    // build a job template form. 
    // Following informations are needed:
    //  - which parameter shall be global and which run-specific
    //  - how many jobs (/samples) shall be created
    constructor(props) {
        super(props);
        // Inputs:
        // props.templFileRelpath
        this.buildContentOnSuccess = this.buildContentOnSuccess.bind(this);
    }

    buildContentOnSuccess(data){ // when AJAX request succeeds
        console.log(data)
        return (
            <div>
                Targeted CWL document: <IneditableValueField>{data.templ_attr.CWL}</IneditableValueField>
                <div style={ {maxHeight:"50vh", overflowY: "auto"} }>
                    <table className="w3-table w3-bordered">
                        <thead>
                            <tr>
                                <th>Parameter</th>
                                <th>Mode</th>
                            </tr>
                        </thead>
                        <tbody>
                            {data.params.map( (param) => (
                                <tr key={param.param_name}> 
                                    <td>{param.param_name}</td>
                                    <td>
                                        global &nbsp;
                                        <label className="switch">
                                            <input 
                                                type="checkbox" 
                                                name="param_mode_select"
                                                value={param.param_name}>
                                            </input>
                                            <span className="slider round"></span>
                                        </label>
                                        &nbsp; run-specific
                                    </td>
                                </tr>
                            ))}
                        </tbody>
                    </table>
                </div>
            </div>
        );
    }

    render() {
        return (
            <AjaxComponent
                requestRoute={routeGetTemplConfigInfo}
                sendData={ {file_relpath: this.props.templFileRelpath} }
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
        const itemValues = this.props.templFilesInfo.map( (tf) => tf.file_relpath);
        const itemNames = this.props.templFilesInfo.map( (tf) => tf.file_nameroot);
        let itemContent = "";
        if(this.state.whichFocus != "") {
            itemContent = <JobTemplateFormPreparationPanel templFileRelpath={this.state.whichFocus} />
        }

        return (
            <div>
                <p>Select one of the following job templates:</p>
                <CollapsibleList
                    itemValues={itemValues}
                    itemNames={itemNames}
                    whichFocus={this.state.whichFocus}
                    itemContent={itemContent}
                    onChange={this.changeFocus}
                />
            </div>
        );
    }
}

class JobTemplRoot extends React.Component {
    constructor(props) {
        super(props);
        this.buildContentOnSuccess = this.buildContentOnSuccess.bind(this);
    }

    buildContentOnSuccess(data){ // when AJAX request succeeds
        return (
            <div>
                <Title>Job templates</Title>
                { data.length > 0 ? (
                        <JobTemplList templFilesInfo={data}>  </JobTemplList>
                    ) : (
                        <Message type="info">There are no Job templates.</Message>
                    )
                }
            </div>
        );
    }

    render() {
        return (
            <AjaxComponent
                requestRoute={routeGetTmplInfo}
                sendData={ {} }
                buildContentOnSuccess={this.buildContentOnSuccess}
                loaderSize="large"
                loaderMessage="Loading available job templates"
            />
        );
    }
}
