

const modules = {   // each entry corresponds to one wep app module,
                    // each module has a dedicated button in the top bar
                    // and specific content which is displayed on click of that button
    home: {
        text: (<span><font color="white">CW</font><font color="lightgreen">Lab</font></span>),
        icon: "",
        content: "Welcome to CWLab"
    },
    import_cwl:  {
        text: "Import CWL Workflow/Tool",
        icon: "fas fa-file-import",
        content: "Under construction"
    },
    create_job:  {
        text: "Create New Job",
        icon: "fas  fa-plus",
        content: "Under construction"
    },
    jobs: {
        text:  "Job Execution & Results",
        icon: "fas fa-rocket",
        content: "Under construction"
    },
    help: {
        text: "Help",
        icon: "fas fa-info-circle",
        content: "Under construction"
    }
};


class TopBar extends React.Component { // controlled by Root Component
    constructor(props){
        super(props);
        this.handleClick = this.handleClick.bind(this);
    }

    handleClick(key) {
        this.props.handleModuleChange(key)
    }
    
    render() {
        return (
            <div className="w3-bar w3-metro-darken">
                {
                    Object.keys(modules).map( (key) => (
                            <a className="w3-bar-item w3-button" 
                                key={key} 
                                onClick={this.handleClick.bind(this, key)}>
                                <i className={modules[key].icon} style={ {paddingLeft: "10px", paddingRight: "10px"} }></i>
                                {modules[key].text}
                            </a>
                        )
                    )
                }
            </div>
        );
    }
}

class MainContent extends React.Component {
    render() {
        return (
            <div className="w3-panel">
                {this.props.children}
            </div>
        );
    }
}

class Root extends React.Component {
    constructor(props) {
        super(props);
        this.state = {module: "home"};  // determines which web app module is displayed 
        this.changeModule = this.changeModule.bind(this)
    }

    changeModule(target_module){
        this.setState({module:target_module})
    }

    render() {
        return (
            <div className="w3-theme-d3 w3-medium">
                <TopBar handleModuleChange={this.changeModule} />
                <MainContent> {modules[this.state.module].content} </MainContent>
            </div>
        );
    }
}