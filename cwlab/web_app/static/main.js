
class Welcome extends React.Component {
    render() {
        return(
            <div className="w3-panel">
                <p style={ {fontSize: 40} } className="w3-center">
                    Welcome to CW<span className="w3-text-green">Lab</span>
                </p>
                <p style={ {fontSize: 20} } className="w3-center">
                    An open-source framework for simplified deployment of the 
                    Common Workflow Language using a graphical web interface
                </p>
            </div>
        )
    }
}

const modules = {   // each entry corresponds to one wep app module,
                    // each module has a dedicated button in the top bar
                    // and specific content which is displayed on click of that button
    home: {
        text: (<span>CW<span className="w3-text-green">Lab</span></span>),
        icon: "",
        content: (<Welcome />)
    },
    import_cwl:  {
        text: "Import CWL Workflow/Tool",
        icon: "fas fa-file-import",
        content: (<ImportCWLRoot />)
    },
    create_job:  {
        text: "Create New Job",
        icon: "fas  fa-plus",
        content: (<CeateJobRoot />)
    },
    jobs: {
        text:  "Job Execution & Results",
        icon: "fas fa-rocket",
        content: (<JobExecRoot />)
    }
    // ,
    // help: {
    //     text: "Help",
    //     icon: "fas fa-info-circle",
    //     content: "Under construction"
    // }
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
            <div 
                className="w3-card w3-bar w3-metro-darken"
                style={ {position: "fixed", top: "0", height: "39px", zIndex: "999"} }

            >
                {
                    Object.keys(modules).map( (key) => (
                            <a
                                className={this.props.whichFocus == key ?
                                    "w3-bar-item w3-button w3-theme-d3" : "w3-bar-item w3-button"}
                                key={key} 
                                onClick={this.handleClick.bind(this, key)}>
                                {modules[key].icon != "" &&
                                    <i className={modules[key].icon} style={ {paddingLeft: "10px", paddingRight: "10px"} }></i>
                                }
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
            <div style={ {paddingTop: "39px"} }>
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
                <TopBar handleModuleChange={this.changeModule} whichFocus={this.state.module} />
                <MainContent> {modules[this.state.module].content} </MainContent>
            </div>
        );
    }
}