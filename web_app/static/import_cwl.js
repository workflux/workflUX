
class ImportCWLRoot extends React.Component {
    render() {
        return(
            <FileUploadComponent
                requestRoute={routeCwlImport}
                instruction="Please select a CWL file:"
            />
        );
    }
}
