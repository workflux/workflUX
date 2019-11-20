
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
    import_wf:  {
        text: "Import CWL Workflow/Tool",
        icon: "fas fa-file-import",
        content: (<ImportCWLRoot />),
        disabled: loginEnabled && ! loggedIn
    },
    create_job:  {
        text: "Create Job",
        icon: "fas  fa-plus",
        content: (<CreateJobRoot />),
        disabled: loginEnabled && ! loggedIn
    },
    jobs: {
        text:  "Job Execution & Results",
        icon: "fas fa-rocket",
        content: (<JobExecRoot />),
        disabled: loginEnabled && ! loggedIn
    },
    // ,
    // help: {
    //     text: "Help",
    //     icon: "fas fa-info-circle",
    //     content: "Under construction"
    // }
    users: {
        text:  loggedIn ? username : "login / register",
        icon: "fas fa-user",
        content: (<UserRoot />),
        align: "right",
        disabled: ! loginEnabled
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
            <div 
                className="w3-card w3-metro-darken"
                style={ this.props.fixed ? (
                        {position: "fixed", top: "0", width:"100%", zIndex:"1000"}
                    ):(
                        {width:"100%"}
                )}
            >
                {
                    Object.keys(modules).map( (key) => (
                            modules[key].hasOwnProperty("disabled") && modules[key].disabled ?
                                (
                                    <span key={key}></span>
                                ) : (
                                    <a
                                        key={key}
                                        className={
                                            "w3-button" + 
                                            (
                                                this.props.whichFocus == key ?
                                                    " w3-theme-d3" : ""
                                            ) +
                                            (
                                                modules[key].hasOwnProperty("align") && modules[key].align == "right" ?
                                                    " w3-right" : ""
                                            )
                                            
                                        }
                                        style={ {display: "inline-block"} }
                                        key={key} 
                                        onClick={this.handleClick.bind(this, key)}>
                                        {modules[key].icon != "" &&
                                            <i className={modules[key].icon} style={ {paddingLeft: "10px", paddingRight: "10px"} }></i>
                                        }
                                        {modules[key].text}
                                    </a>
                                )
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
            <div>
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
                <TopBar handleModuleChange={this.changeModule} whichFocus={this.state.module} fixed={true}/>
                <MainContent> {modules[this.state.module].content} </MainContent>
            </div>
        );
    }
}