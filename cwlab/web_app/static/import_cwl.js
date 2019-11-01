
class ImportCWLRoot extends React.Component {
    render() {
        return(
            <div className="w3-panel">
                <span className="w3-text-green">Choose a name:</span>&nbsp;
                <input type="text"
                    className="w3-input w3-border"
                    name="job_name"
                    value={this.state.job_name}
                    onChange={this.changeJobName}
                />
                <FileUploadComponent
                    requestRoute={routeImportPackedCwl}
                    instruction="Please select a CWL file:"
                />
            </div>
        );
    }
}
