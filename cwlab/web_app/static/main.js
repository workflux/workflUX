
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
                {demo && (
                    <Message type="info">
                        <span className="w3-text-black w3-center">
                            <p>
                                <b>Please note: </b>
                            </p> 
                            <p>
                                This instance is for demonstration purposes only and some functionalites are locked.
                                Please do not submit large-scale workflow executions.
                            </p>
                        </span>
                    </Message>
                )}
            </div>
        )
    }
}



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
                    Object.keys(this.props.modules).map( (key) => (
                            this.props.modules[key].hasOwnProperty("disabled") && this.props.modules[key].disabled ?
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
                                                this.props.modules[key].hasOwnProperty("align") && this.props.modules[key].align == "right" ?
                                                    " w3-right" : ""
                                            )
                                            
                                        }
                                        style={ {display: "inline-block"} }
                                        key={key} 
                                        onClick={this.handleClick.bind(this, key)}>
                                        {this.props.modules[key].hasOwnProperty("customIconHTML") && this.props.modules[key].customIconHTML ? (
                                                <span dangerouslySetInnerHTML={ {__html: this.props.modules[key].customIconHTML} } style={ {height: "20px"} } />
                                            ) : (
                                                <span>    
                                                    {this.props.modules[key].icon != "" ? (
                                                            <i className={this.props.modules[key].icon} style={ {paddingLeft: "10px", paddingRight: "10px"} }></i>
                                                        ) : (
                                                            <span/>
                                                        )
                                                    }
                                                </span>
                                            )
                                        }
                                        {this.props.modules[key].text}
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
        this.state = {
            userInfoLoaded: false,
            module: "home"
        };  // determines which web app module is displayed 
        this.changeModule = this.changeModule.bind(this);
    }

    async componentDidMount() {
        const userInfo = await getUserInfo()

        this.modules = {   // each entry corresponds to one wep app module,
                            // each module has a dedicated button in the top bar
                            // and specific content which is displayed on click of that button
            home: {
                text: (
                    <span>
                        CW<span className="w3-text-green">Lab</span>
                        {buildNumber != "none" && (<span className="w3-text-orange">&nbsp;build {buildNumber}</span>)}
                    </span>
                ),
                icon: "",
                content: (<Welcome />)
            },
            import_wf:  {
                text: "Import CWL Workflow/Tool",
                icon: "fas fa-file-import",
                content: (<ImportCWLRoot />),
                disabled: loginEnabled && !userInfo.isLoggedIn
            },
            create_job:  {
                text: "Create Job",
                icon: "fas  fa-plus",
                content: (<CreateJobRoot />),
                disabled: loginEnabled && !userInfo.isLoggedIn
            },
            jobs: {
                text:  "Job Execution & Results",
                icon: "fas fa-rocket",
                content: (<JobExecRoot />),
                disabled: loginEnabled && !userInfo.isLoggedIn
            },
            // ,
            // help: {
            //     text: "Help",
            //     icon: "fas fa-info-circle",
            //     content: "Under construction"
            // }
            user: {
                text:  userInfo.isLoggedIn ? (
                    userInfo.username
                ) : (
                    "login / register"
                ),
                icon: "fas fa-user",
                customIconHTML: customLoginIconHTML,
                content: (<UserRoot />),
                align: "right",
                disabled: !loginEnabled
            }
        }

        this.setState({
            userInfoLoaded: true
        })
    }

    changeModule(target_module){
        this.setState({module: target_module})
    }

    render() {
        return (
            <div className="w3-theme-d3 w3-medium">
                {this.state.userInfoLoaded ? (
                        <div>
                            <TopBar 
                                modules={this.modules} 
                                handleModuleChange={this.changeModule} 
                                whichFocus={this.state.module} 
                            />
                            <TopBar 
                                modules={this.modules} 
                                handleModuleChange={this.changeModule} 
                                whichFocus={this.state.module} 
                                fixed={true}
                            />
                            <MainContent> {this.modules[this.state.module].content} </MainContent>
                        </div>
                    ) : (
                        <LoadingIndicator
                            size="large"
                            message="Please wait."
                        />
                    )
                }
            </div>
        );
    }
}

