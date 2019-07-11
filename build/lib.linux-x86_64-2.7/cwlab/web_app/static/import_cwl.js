
class ImportCWLRoot extends React.Component {
    render() {
        return(
            <div className="w3-panel">
                <FileUploadComponent
                    requestRoute={routeCwlImport}
                    instruction="Please select a CWL file:"
                />
            </div>
        );
    }
}
