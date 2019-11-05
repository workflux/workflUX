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

class ImportSingleCWL extends React.Component{
    constructor(props){
        super(props);
        this.state = {
            importName: ""
        }

        this.changeInputField = this.changeInputField.bind(this);
    }

    changeInputField(event){
        this.setState({
            [event.currentTarget.name]: event.currentTarget.value
        })
    }

    render(){
        return(
            <div className="w3-panel">
                <p>
                    Import a CWL-wrapped tool or a packed CWL workflow.
                </p>
                <span className="w3-text-green">Choose a name:</span>&nbsp;
                <input type="text"
                    className="w3-input w3-border"
                    name="importName"
                    style={ {width: "50%"} }
                    value={this.state.importName}
                    onChange={this.changeInputField}
                />
                <br/>
                <span className="w3-text-green">Choose a CWL file:</span>&nbsp;
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
            importName: "",
            importMethod: "singleCWL"
        }

        this.importMethods = {
            singleCWL: {
                descr: "a single CWL file (CWL-wrapped tool or a packed CWL Workflow)",
                component: <ImportSingleCWL />
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
                <span className="w3-text-green">Choose how and what you like to import:</span>
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
