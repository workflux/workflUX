// contains general utilities important for multiple modules
// import styled from "styled-components"

function get_time_str(){
    const date = new Date()
    return(date.toLocaleTimeString())
}

function storeSessionInfo(param_obj){
    Object.keys(param_obj).map( (param) => {
        sessionStorage.setItem(param, param_obj[param])
        localStorage.setItem(param, param_obj[param])
    })
}

function cleanSessionInfo(){
    sessionStorage.clear()
    localStorage.clear()
}

function getSessionInfo(param_list){
    let param_obj = {}
    param_list.map( (param) => {
        let value = sessionStorage.getItem(param)
        if (!value){
            value = localStorage.getItem(param)
        }
        if (value == "true"){
            value = true
        }
        else if (value == "false"){
            value = false
        }
        param_obj[param] = value
    })
    return(param_obj)
}

function seconds_to_duration_str(s){
    const h = Math.floor( s / 3600 )
    let s_remain = s % 3600
    const m = Math.floor( s_remain / 60 )
    s_remain = s_remain % 60
    let duration_str = h == 0 ? "" : h.toString() + " h "
    duration_str += m == 0 ? "" : m.toString() + " m "
    duration_str += s_remain == 0 ? "" : s_remain.toString() + " s"
    return(duration_str.trim())
}

async function ajaxRequest({
    // in a component bind this function: this.ajaxRequest = ajaxRequest.bind(this)
    statusVar="actionStatus",
    statusValueDuringRequest="action",
    statusValueAfterRequest="none",
    messageVar="serverMessages",
    sendData={},
    sendViaFormData=false,
    route,
    onSuccess= (data, messages) => {
                return({})
            }, //function taking arguments data and messages, return state update
    onError= (messages) => { //function taking argument messages, return state update
        return({})
    } 
}){
    let setDataAT = sendData
    const userInfo = await getUserInfo("all")
    setDataAT.access_token = userInfo.accessToken
    if (! setDataAT.hasOwnProperty("username")){
        setDataAT.username = userInfo.username
    }
    let formData
    if (sendViaFormData){
        formData = new FormData()
        formData.append("meta", JSON.stringify(setDataAT))
    }
    this.setState({
        [statusVar]: statusValueDuringRequest
    })
    fetch(route, {
        method: "POST",
        body: sendViaFormData ? formData : JSON.stringify(setDataAT),
        headers: new Headers(sendViaFormData ? (
                {}
            ) : (
                {
                    'Content-Type': 'application/json'
                }
            )
        ),
        cache: "no-cache"
    }).then(res => res.json())
    .then(
        (result) => {
            const messages = result.messages;
            const data = result.data
            let errorOccured = false;
            for( let i=0;  i<messages.length; i++){
                if(messages[i].type == "error"){
                    errorOccured = true;
                    break;
                }
            }
            let stateUpdate = {
                [statusVar]: statusValueAfterRequest,
                [messageVar]: messages
            }
            if (! errorOccured){
                // success
                Object.assign(stateUpdate, onSuccess(data, messages))
            }
            else{
                // server returned error
                Object.assign(stateUpdate, onError(messages))
            } 
            if (!(stateUpdate.hasOwnProperty("doNotUpdate") && stateUpdate.doNotUpdate)){
                this.setState(stateUpdate)
            }
        },
        (error) => {
            // server could not be reached
            const messages = [{
                time: get_time_str(),
                type: "error", 
                text: serverNotReachableError
            }];
            let stateUpdate = {
                [statusVar]: statusValueAfterRequest,
                [messageVar]: messages
            }
            Object.assign(stateUpdate, onError(messages))
            if (!(stateUpdate.hasOwnProperty("doNotUpdate") && stateUpdate.doNotUpdate)){
                this.setState(stateUpdate)
            }
        }
    )

}

function getAllowedDirs({
    statusVar="actionStatus",
    statusValueDuringRequest="action",
    statusValueAfterRequest="none",
    messageVar="serverMessages",
    jobName=null,
    runName=null,
    onSuccess= (data, messages) => {
                return({})
            },
    onError= (messages) => {
        return({})
    } 
}){
    let sendData = {}
    if (jobName){
        sendData["job_name"] = jobName
        if (runName){
            sendData["run_name"] = runName
        }
    }

}

String.prototype.replaceAll = function(search, replacement) {
    // in contrast to replace() replaces not only the first occurence
    var target = this;
    return target.split(search).join(replacement);
};

class IneditableValueField extends React.Component {
    // Input:
    // props.backclassName e.g. w3-metro-darken w3-theme-dark
    // props.textclassName e.g. w3-text-green
    render() {
        const style = {
            fontFamily: "courier",
            padding: "2px"
        }
        let backclassName = "w3-metro-darken"
        let textclassName = ""
        if (this.props.backclassName){
            backclassName = this.props.backclassName
        }
        if (this.props.textclassName){
            textclassName = this.props.textclassName
        }
        return( 
            <span 
                className={backclassName + " " + textclassName} 
                style={style}> 
                {this.props.children} 
            </span> 
        )
    }
}

class TabPanel extends React.Component { // controlled by Root Component
    constructor(props){
        super(props);
        //props.title
        // props.tabs
        // props.changeFocus
        // props.whichFocus
        // props.children
        this.handleClick = this.handleClick.bind(this);
    }

    handleClick(tab) {
        this.props.changeFocus(tab)
    }
    
    render() {
        return (
            <div>
                <div 
                    className="w3-card w3-metro-darken"
                    style={ {width:"100%"} }
                >
                    {this.props.title && (
                        <div 
                            className="w3-padding-small"
                            style={ {display: "inline-block"} }
                        >
                            {this.props.title}
                        </div>
                    )}
                    {
                        this.props.tabs.map( (tab) => (
                                <a
                                    className={this.props.whichFocus == tab ?
                                        "w3-button w3-theme-d4" : "w3-button"}
                                    style={ {display: "inline-block"} }
                                    key={tab} 
                                    onClick={this.handleClick.bind(this, tab)}>
                                    {tab}
                                </a>
                            )
                        )
                    }
                </div>
                <div 
                    className="w3-container w3-padding w3-theme-d4"
                >
                    {this.props.children}
                </div>
            </div>
        );
    }
}

class Checkbox extends React.Component{
    constructor(props){
        super(props);
        // props.name
        // props.value
        // props.onChange function executed on change
        //  takes two arguments: (1) value, (2) whether set or not
        // props.checked
        // props.disabled
        // props.size
        this.handleChange = this.handleChange.bind(this)
    }

    handleChange(event){
        const value=event.currentTarget.value
        const isSet=event.currentTarget.checked
        this.props.onChange(value, isSet)
    }

    render(){
        return(
            <input
                className="w3-check"
                style={ {height: this.props.size ? (this.props.size) : ("20px")} }
                type="checkbox"
                name={this.props.name}
                value={this.props.name}
                checked={this.props.checked}
                onChange={this.handleChange}
            />
        )
    }
}

class BooleanSlider extends React.Component {
    constructor(props) {
        super(props);
        // props.name
        // props.value
        // props.onChange function executed on change
        //  takes two arguments: (1) value, (2) whether set or not
        // props.checked
        // props.disabled
        // props.doNotSendValue
        // props.forwardEvent
        this.handleChange = this.handleChange.bind(this)
    }
    
    handleChange(event){
        if (this.props.forwardEvent){
            const event_ = {
                currentTarget: {
                    name: this.props.name,
                    value: event.currentTarget.checked
                }
            }
            this.props.onChange(event_)
        }
        else{
            const value=event.currentTarget.value
            const isSet=event.currentTarget.checked
            if (this.props.doNotSendValue){
                this.props.onChange(isSet)
            }
            else {
                this.props.onChange(value, isSet)
            }
        }
    }

    render() {
        return( 
            <label className="switch">
                <input 
                    type="checkbox" 
                    name={this.props.name}
                    value={this.props.value}
                    onChange={this.handleChange}
                    disabled={this.props.disabled ? true : false}
                    checked={this.props.checked}
                />
                <span className="slider round"></span>
            </label>
        );
    }
}

class ActionButton extends React.Component {
    constructor(props) {
        super(props);
        // props.name
        // props.value
        // props.className set color class default is w3-black
        // props.onAction function to execute on click
        // props.label labeling
        // props.loading if true, button will be disabled and loading indicator is shown
        // props.disabled if true, disable without loading indicator
        // props.smallPadding true/false
        // props.style
        // props.forwardEvent
        this.handleAction = this.handleAction.bind(this)
    }
    
    handleAction(event){
        if (this.props.forwardEvent){
            this.props.onAction({
                currentTarget:{
                    value: this.props.value,
                    name: this.props.name
                }
            })
        }
        else{
            this.props.onAction(this.props.value)
        }
    }

    render() {
        const style = this.props.style ? (
                this.props.style
            ) : (
                {marginTop: "2px", marginBottom: "2px", width: this.props.width}
            )
        return( 
            <button style={style}
                name={this.props.name}
                value={this.props.value}
                className={
                    "w3-button" +
                    (this.props.className ? (" " + this.props.className) : " w3-black") +
                    (this.props.smallPadding ? (" w3-padding-small") : "")
                }
                onClick={this.handleAction}
                disabled={this.props.loading || this.props.disabled ? true : false}
                > 
                {this.props.loading ? (<LoadingIndicator message="" size="tiny" />) : null}
                {this.props.label}
            </button> )
    }
}

class Tooltip extends React.Component {
    constructor(props) {
        super(props);
        // props.hoverButton
        // props.children
        // props.preview
        // props.title
        // props.disable
        this.state = {
            showTopPanel: false
        }

        this.showTopPanel = this.showTopPanel.bind(this);
    }

    showTopPanel(event){
        this.setState({
            showTopPanel: event.currentTarget.name == "show_top_panel" ? true : false
        })
    }

    render(){
        let preview = this.props.preview ? this.props.preview : this.props.children
        preview = preview.length > 40 ? (
            preview.substring(0, 28).concat(" ... [click]")
        ) : (
            preview
        )
        let width = preview.length*8+20
        if (width < 140) {
            width = 140
        }
        else if (width > 400){
            width = 400
        }
        return(
            <div 
                className={this.props.disabled ? "" : "tooltip"}
                style = { {display: "inline-block"} }
            >
                {this.state.showTopPanel && !this.props.disabled ? (
                        <TopPanel>
                            <div 
                                className="w3-bar"
                                style={ {height: "40px"} }
                            >
                                <div 
                                    className="w3-bar-item w3-display-topleft"
                                    style={ {padding: "8px 16px"} }
                                >
                                    {this.props.title}
                                </div>
                                <ActionButton
                                    name="exit_top_panel"
                                    value="exit_top_panel"
                                    className="w3-black w3-bar-item w3-display-topright"
                                    label={<i className="fas fa-times"/>}
                                    onAction={this.showTopPanel}
                                    forwardEvent={true}
                                />
                            </div>
                            <div className="w3-container" style={ {overflowY: "auto", whiteSpace: "pre-wrap"} }>
                                {this.props.children}
                            </div>
                        </TopPanel>
                    ) : (
                        <span 
                            className="tooltiptext w3-metro-darken"
                            style={ {width: width.toString().concat("px"), whiteSpace: "nowrap", zIndex: "100"} }
                        >
                            {preview}
                        </span>   
                    )
                }
                {this.props.hoverButton && !this.props.disabled ? (
                        this.props.hoverButton
                    ) : (
                        <a
                            name="show_top_panel"
                            onClick={this.showTopPanel}
                        >
                            <i className="fas fa-info-circle w3-text-blue"/>
                        </a>
                    )
                }
            </div>
        )
    }

}

class CollapsibleListItem extends React.Component {
    constructor(props) {
        super(props);
        this.handleChange = this.handleChange.bind(this);
    }

    handleChange(event) {
        this.props.onChange(event.currentTarget.value)
    }

    render() {
        // Inputs:
        // props.value
        // props.name
        // props.itemContent contains content; "" for empty
        
        const headerStyle = {
            marginTop: "10px",
            marginBottom: "0px",
            paddingTop: "10px",
            paddingBottom: "10px",
            paddingLeft: "16px",
            paddingRight: "16px"
        }
        const contentStyle = {
            marginTop: "0px",
            marginBottom: "0px",
            paddingTop: "10px",
            paddingBottom: "10px",
            paddingLeft: "16px",
            paddingRight: "16px"
        }

        return (
            <div>
                <label>
                    <div key={this.props.value} style={headerStyle}
                        className="w3-metro-darken">
                        <input 
                            type="radio" name="templ_select"
                            value={this.props.value}
                            onChange={this.handleChange}>
                        </input>
                        &nbsp; {this.props.name}
                    </div>
                </label>
                {this.props.itemContent != "" &&
                    <div style={contentStyle} className="w3-theme-d5">
                        {this.props.itemContent}
                    </div>
                }
            </div>
        );
    }
}

class Table extends React.Component {
    constructor(props){
        super(props);
        // props.columnKeys
        // props.columnNames
        // props.selectionEnabled
        // props.handleSelectionChange
        // props.selection
        // props.rowData
        // props.selectionKey

        this.handleSelectionChange = this.handleSelectionChange.bind(this)
        this.toggelRunSelectionAll = this.toggelRunSelectionAll.bind(this)
    }

    handleSelectionChange(event){
        let newSelection
        if (this.props.selection.includes(event.currentTarget.value)){
            newSelection = this.props.selection.filter(item => item != event.currentTarget.value)
        }
        else {
            newSelection = this.props.selection
            newSelection.push(event.currentTarget.value)
        }
        this.props.handleSelectionChange(newSelection)
    }

    toggelRunSelectionAll(){
        // if all runs are selected, deselect all:
        let newSelection
        if (this.props.selection.length == 0 ){
            newSelection = this.props.rowData.map((r) => (r[this.props.selectionKey].hasOwnProperty("value") ? 
                r[this.props.selectionKey].value : r[this.props.selectionKey]
            ))
        }
        else {
            newSelection = []
        }
        this.props.handleSelectionChange(newSelection)
    }


    render(){
        return(
            <div style={ {maxHeight:"50vh", overflowY: "auto"} }>
                <table className="w3-table w3-bordered w3-border">
                    <thead className="w3-text-green">
                        <tr>
                            {this.props.selectionEnabled && (
                                <th>
                                    <i 
                                        className="fas fa-arrow-down" 
                                        style={ {paddingRight:"10px"} }
                                    />
                                    <ActionButton 
                                        name="select all"
                                        value="select all"
                                        label="all"
                                        onAction={this.toggelRunSelectionAll}
                                        className="w3-black w3-text-green"
                                        smallPadding={true}
                                    />
                                </th>
                            )}
                            {this.props.columnKeys.map( (c) => (
                                <th key={c}>{this.props.columnNames[c]}</th>
                            ))}
                        </tr>
                    </thead>
                    <tbody>
                        {this.props.rowData.map( (r) => {
                            let key = r[this.props.selectionKey].hasOwnProperty("value") ? r[this.props.selectionKey].value : r[this.props.selectionKey];
                            return(
                                <tr key={key}>
                                    {this.props.selectionEnabled && (
                                        <td>
                                            <input 
                                                type="checkbox" 
                                                name="select"
                                                value={key}
                                                checked={this.props.selection.includes(key)}
                                                onChange={this.handleSelectionChange}
                                            />
                                        </td>
                                    )}
                                    {this.props.columnKeys.map( (c) => {
                                        let value
                                        let className
                                        if (r[c]){
                                            value = r[c].hasOwnProperty("value") ? (
                                                        r[c].value ? r[c].value : "-"
                                                    ) :(
                                                        r[c]
                                                    )
                                            className =  r[c].hasOwnProperty("className") ? r[c].className : "";
                                        }
                                        else {
                                            value = "-"
                                            className = ""
                                        }
                                        return(
                                            <td 
                                                key={key + c}
                                                className={className}
                                            >
                                                {value} 
                                            </td>
                                        )
                                    })}
                                </tr>
                            )
                        })}
                    </tbody>
                </table>
            </div>
        )
    }
}

class CollapsibleList extends React.Component {
    render () {
        // Inputs:
        // props.itemValues contains values of list items
        // props.itemNames contains names of list items
        // props.itemContent contains content for the
        //  list item which has focus; maximally one item can have focus
        // props.whichFocus contains the value of the
        //  the list items that has focus; "" if none has focus
        // props.onChange function executed once the focus changes
        
        return (
            <div>
                {
                    this.props.itemValues.map( (itemValue, index) => (
                        <CollapsibleListItem 
                            key={itemValue}
                            value={itemValue}
                            name={this.props.itemNames[index]}
                            itemContent={this.props.whichFocus == itemValue ? this.props.itemContent : ''}
                            onChange={this.props.onChange}
                        />
                    ))
                }
            </div>
        );

    }
}


class SideBarPanel extends React.Component {
    // Inputs:
    // props.itemValues contains values of list items
    // props.itemNames contains names of list items
    // props.itemContent contains content for the
    //  list item which has focus; maximally one item can have focus
    // props.whichFocus contains the value of the
    //  the list items that has focus; "" if none has focus
    // props.onChange function executed once the focus changes
    // props.label

    constructor(props){
        super(props)
        this.handleClick = this.handleClick.bind(this)
    }

    handleClick(event) {
        this.props.onChange(event.currentTarget.value)
    }

    render () {
        return (
            <div className="w3-row">
                <div 
                    className="w3-card w3-sidebar w3-bar-block w3-col l2 m4 s4 w3-metro-darken"
                    style={ {width:"280px", overflowY: "auto"} }
                >
                    <div className="w3-bar-item w3-text-green">
                        {this.props.label}
                    </div>
                    {
                        this.props.itemValues.map( (itemValue, index) => (
                            <button
                                className={this.props.whichFocus == itemValue ?
                                    "w3-bar-item w3-button w3-theme-d3" : "w3-bar-item w3-button"}
                                key={itemValue}
                                value={itemValue}
                                name={this.props.itemNames[index]}
                                onClick={this.handleClick}
                            >
                                {this.props.itemNames[index]}
                            </button>
                        ))
                    }
                </div>
                <div 
                    className="w3-col s10 m8 s8 w3-panel"  
                    style={ {marginLeft:"280px", width: "calc(100% - 280px)"} }
                >
                    {this.props.itemContent}
                </div>
            </div>
        );

    }
}


class LoadingIndicator extends React.Component {
    render() {
        let loaderClass
        if( this.props.size == 'large' ){
            loaderClass='large_loader';
        } else if( this.props.size == 'small' ){
            loaderClass='small_loader';
        } else{
            loaderClass='tiny_loader';
        }
        return (
            <div className="w3-container w3-center">
                    <div className={loaderClass}></div>
                    {this.props.message ? (
                        <div className="w3-center-align">{this.props.message}</div>
                        ):( null )
                    }
            </div>
        );
    }
}


class Title extends React.Component {
    render() {
        return ( <h1> {this.props.children} </h1> )
    }
}


class Message extends React.Component {
    render() {
        const color = {
            error: "w3-red",
            warning: "w3-orange",
            info: "w3-blue",
            success: "w3-green",
            hint: "w3-khaki"
        };
        return ( <div className={"w3-panel " + color[this.props.type]}> {this.props.children} </div> );
    }
}


class DisplayServerMessages extends React.Component {
    // Inputs:
    // props.messages
    render() {
        if( this.props.messages.length > 0) {
            return (
            <div>
                {
                    this.props.messages.map( (message, index) => (
                        <div key={index}>
                            <Message type={message.type}>
                                {message.hasOwnProperty("time") &&(
                                    <span>{message.time} - </span>
                                )}
                                {message.text}
                            </Message>
                        </div>
                    ))
                }
            </div>
            );
        }
        else {
            return null;
        }
    }
}

class ExperimentalTag extends React.Component {
    render() {
        return(
            <Message type="warning">
                <b>Please Note: This feature is experimental.</b><br/>
                If you experience any issues, 
                please report them to&nbsp;
                <a href="https://github.com/CompEpigen/CWLab/issues">
                    https://github.com/CompEpigen/CWLab/issues
                </a> or contact k.breuer@dkfz.de.
            </Message>
        )
    }
}

function sleep(secs) {
    secs = (+new Date) + secs * 1000;
    while ((+new Date) < secs);
}


class AjaxComponent extends React.Component {
    // diplays content based on server information
    constructor(props) {
        super(props);
        // Inputs:
        // props.requestRoute route to perform the post request on
        // props.buildContentOnSuccess function which builds content on successful request
        //  function receives a single argument (data)
        // props.loaderSize size of loading indicator (large, small, tiny)
        // props.loaderMessage message displayed on loading indicator
        // props.suppressMessages

        this.state = {
            loading: true,
            serverMessages: [],
            data: [],
            success: false
        };  
                                    
        this.request = this.request.bind(this);
        this.ajaxRequest = ajaxRequest.bind(this)
    }

    async request(){ // ajax request to server
        let sendData = this.props.sendData
        sendData["access_token"] = await getUserInfo("accessToken")
        this.ajaxRequest({
            statusVar: "loading",
            statusValueDuringRequest: true,
            statusValueAfterRequest: false,
            messageVar: "serverMessages",
            sendData: sendData,
            route: this.props.requestRoute,
            onSuccess: (data, messages) => {
                return({
                    data: data,
                    success: true
                })
            },
            onError: (messages) => {
                return({
                    data: [],
                    success: false
                })
            }
        })
    }

    componentDidMount() {
        if(this.state.loading) {
            this.request();
        }
    }

    render() {
        if (this.state.loading) {
            return(
                <LoadingIndicator message={this.props.loaderMessage} size={this.props.loaderSize} />
            )
        } else{
            return (
                <div>
                    {!this.props.suppressMessages &&
                        <DisplayServerMessages messages={this.state.serverMessages} />
                    }
                    {this.state.success &&
                        this.props.buildContentOnSuccess(this.state.data, this.state.serverMessages, this.request)
                    }
                </div>
            )
        }
    }
}

class TopPanel extends React.Component {
    constructor(props){
        super(props);
        // props.children
    }
    render(){
        return(
            <div>
                <div
                    style={ {
                        position: "fixed",  
                        zIndex: "100002",
                        top: "0",
                        left: "0",
                        width: "100vw",
                        height: "100vh"
                    } }
                >
                    <div
                        className="w3-card w3-metro-darken w3-display-middle"
                        style={ {
                            width: "80%"
                        } }
                    >
                        {this.props.children}
                    </div>
                </div>
                <div
                    className="w3-white"
                    style={ {
                        position: "fixed",  
                        zIndex: "100000",
                        top: "0",
                        left: "0",
                        width: "100vw",
                        height: "100vh",
                        opacity: "0.5"
                    } }
                />
            </div>
            
        )
    }
}

class BrowseDir extends React.Component {
    constructor(props){
        super(props);
        // props.path
        // props.prevPath
        // props.ignoreFiles
        // props.fileExts
        // props.showOnlyHits
        // props.selectDir
        // props.allowInput
        // props.allowUpload
        // props.allowDownload
        // props.jobName
        // props.runName
        // props.defaultBaseDir
        // props.showCancelButton
        // props.terminateBrowseDialog
        // props.fixedBaseDir
        // props.fixedBaseDirName
        // props.includeTmpDir
        // props.disableOnTop

        this.baseDirInfo = {
            input: "select only",
            upload: (
                (this.props.allowInput ? "select, " : "") + 
                "upload"
            ),
            download: (
                "download"
            ),
        }

        this.allowInput = this.props.allowInput
        this.allowUpload = this.props.allowUpload
        this.allowDownload = this.props.allowInput ? false : this.props.allowDownload
        this.selectDir = this.props.allowInput ? this.props.selectDir : false

        this.state = {
            dirPath: this.props.path ? this.props.path : "",
            address: this.props.path ? this.props.path : "",
            actionStatus: "init",
            serverMessages: [],
            downloadMessages: [],
            items: [],
            selectedItem: this.props.path && this.props.path != "Please fill" && this.props.path != "" ? (
                    this.props.path
                ) : (
                    null
                ),
            baseDir: "error",
            allowedDirs: {error: {path: "error", mode: "error"}}
        };  

        this.ajaxRequest = ajaxRequest.bind(this);
        this.getItemsInDir = this.getItemsInDir.bind(this);
        this.handleAction = this.handleAction.bind(this);
        this.changeStateVar = this.changeStateVar.bind(this);
        this.changeBaseDir = this.changeBaseDir.bind(this);
        this.terminateBrowseDialog = this.terminateBrowseDialog.bind(this);
        this.handleFileUpload = this.handleFileUpload.bind(this);
        this.triggerFileSend = this.triggerFileSend.bind(this);
        this.downloadFileOrFolder = this.downloadFileOrFolder.bind(this);
    }

    componentDidMount(){
        let path = ""
        if (this.state.selectedItem){
            path = this.state.selectedItem
        }
        else if (this.props.prevPath){
            path = this.props.prevPath
        }
        this.getItemsInDir(false, path, true)
    }

    async getItemsInDir(getParentDir, targetDir, init){
        this.ajaxRequest({
            statusVar: "actionStatus",
            statusValueDuringRequest: init ? "init" : "loading",
            messageVar: "serverMessages",
            sendData: {
                path: targetDir ? targetDir : this.state.dirPath,
                ignore_files: this.props.ignoreFiles ? this.props.ignoreFiles : false,
                file_exts: this.props.fileExts ? this.props.fileExts : [], 
                show_only_hits: this.props.showOnlyHits ? this.props.showOnlyHits : false,
                get_parent_dir: getParentDir ? true : false,
                allow_input: this.allowInput ? true : false,
                allow_upload: this.allowUpload ? true : false,
                allow_download: this.allowDownload ? true : false,
                job_name: this.props.jobName ? this.props.jobName : null,
                run_name: this.props.runName ? this.props.runName : null,
                on_error_return_base_dir_items: init ? true : false,
                default_base_dir: this.props.defaultBaseDir ? this.props.defaultBaseDir : null,
                fixed_base_dir: this.props.fixedBaseDir ? this.props.fixedBaseDir : null,
                fixed_base_dir_name: this.props.fixedBaseDirName ? this.props.fixedBaseDirName: "FIXED_BASE_DIR",
                include_tmp_dir: this.props.includeTmpDir ? true : false
            },
            route: routeBrowseDir,
            onSuccess: (data, messages) => {
                return({
                    items: data.items,
                    dirPath: data.dir,
                    address: data.dir,
                    baseDir: data.base_dir,
                    allowedDirs: data.allowed_dirs
                })
            },
            onError: (message) => {
                if(this.state.address != this.state.dirPath){
                    return({items:[]})
                }
                else{
                    return({})
                }
            }
        })
    }

    handleAction(event){
        if (event.currentTarget.name == "up"){
            this.getItemsInDir(true, this.state.dirPath)
        }
        else if (event.currentTarget.name == "refresh"){
            this.getItemsInDir(false, this.state.dirPath)
        }
        else if (event.currentTarget.name == "open"){
            this.getItemsInDir(false, event.currentTarget.value)
        }
        else if (event.currentTarget.name == "address"){
            if (event.key == "Enter"){
                this.getItemsInDir(false, event.currentTarget.value)
            }
        }
    }

    changeStateVar(event){
        this.setState({
            [event.currentTarget.name]: event.currentTarget.value
        })
    }

    changeBaseDir(event){
        this.setState({
            baseDir: event.currentTarget.value,
            address: this.state.allowedDirs[event.currentTarget.value]
        })
        this.getItemsInDir(false, this.state.allowedDirs[event.currentTarget.value].path)
    }

    handleFileUpload(isSuccess, data){
        if(isSuccess && this.allowInput){
            this.setState({
                selectedItem: data.file_path
            })
            this.getItemsInDir(false, this.state.dirPath)
        }
    }

    terminateBrowseDialog(event){
        let changes
        let selectedItem
        if (event.currentTarget.nam == "cancel"){
            changes = false
            selectedItem = null
        }
        else if (event.currentTarget.name == "select_file"){
            changes = true
            selectedItem = this.state.selectedItem
        }
        else if (event.currentTarget.name == "select_dir"){
            changes = true
            selectedItem = this.state.dirPath
        }
        this.props.terminateBrowseDialog(changes, selectedItem)
    }
    
    async triggerFileSend(event, data, path){
        let form = document.createElement('form');
        form.method = 'post';
        form.action = routeDownload;
        let input = document.createElement('input');
        input.type = 'hidden';
        input.name = "meta";
        input.value = JSON.stringify({
            path: event.currentTarget.name == "download_dir" ? data["zip_path"] : path,
            job_name: this.props.jobName ? this.props.jobName : null,
            run_name: this.props.runName ? this.props.runName : null,
            send_file: true,
            access_token: await getUserInfo("accessToken")
        });
        form.appendChild(input);
        document.body.appendChild(form);
        form.submit();
        return({downloadMessages: []})
    }

    async downloadFileOrFolder(event){
        if (event.currentTarget.name == "download_dir"){
            this.setState({
                downloadMessages: {
                    time: get_time_str(),
                    type: "warning",
                    text: "Creating zip. This might take several minutes. Please wait."
                }
            })
        }
        const path = event.currentTarget.name == "download_dir" ? this.state.dirPath : this.state.selectedItem
        this.ajaxRequest({
            statusVar: "actionStatus",
            statusValueDuringRequest: event.currentTarget.name,
            messageVar: "downloadMessages",
            sendData: {
                path: path,
                job_name: this.props.jobName ? this.props.jobName : null,
                run_name: this.props.runName ? this.props.runName : null,
                send_file: false
            },
            sendViaFormData: true,
            route: routeDownload,
            onSuccess: (data, messages) => this.triggerFileSend(event, data, path)
        })
    }

    render(){
        const allowInput = this.allowInput && ["upload", "input"].includes(this.state.allowedDirs[this.state.baseDir].mode)
        const allowUpload = this.allowUpload && this.state.allowedDirs[this.state.baseDir].mode == "upload"
        const allowDownload = this.allowDownload && this.state.allowedDirs[this.state.baseDir].mode == "download"
        const content=(
            <div>
                {this.state.actionStatus == "init" ? (
                        <LoadingIndicator size="large" message="Loading directory. Please wait." />
                    ) : (
                        <span>
                            <div>
                                <div className="w3-bar">
                                    <ActionButton
                                        name="up"
                                        value="up"
                                        className="w3-bar-item"
                                        label={<span><i className="fas fa-arrow-up" />&nbsp;up</span>}
                                        onAction={this.handleAction}
                                        disabled={this.state.actionStatus != "none"}
                                        forwardEvent={true}
                                    />
                                    <ActionButton
                                        name="refresh"
                                        value="refresh"
                                        className="w3-bar-item"
                                        label={<span><i className="fas fa-sync" />&nbsp;refresh</span>}
                                        onAction={this.handleAction}
                                        disabled={this.state.actionStatus != "none"}
                                        forwardEvent={true}
                                    />
                                    <div className="w3-bar-item">
                                        Change base directory:&nbsp;
                                        <select
                                            className="w3-button w3-white w3-border w3-padding-small" 
                                            name="inputBaseDir"
                                            onChange={this.changeBaseDir}
                                            value={this.state.baseDir}
                                        >
                                            {Object.keys(this.state.allowedDirs).map((d) =>(
                                                <option
                                                    key={d}
                                                    value={d}
                                                >
                                                    {d + " (" + this.baseDirInfo[this.state.allowedDirs[d].mode] + ")"}
                                                </option>
                                            ))}
                                        </select>
                                    </div>
                                </div>
                                <span className="w3-text-green">
                                    Current Dir:
                                </span>
                                <input
                                    className="w3-input w3-bar-item"
                                    type="text"
                                    name="address"
                                    value={this.state.address}
                                    disabled={this.state.actionStatus != "none"}
                                    onChange={this.changeStateVar}
                                    onKeyPress={this.handleAction}
                                />

                                <DisplayServerMessages messages={this.state.serverMessages} />
                            </div>
                            <div 
                                className="w3-panel w3-theme-d3"
                                style={ {
                                    overflowY: "auto",
                                    height: "50vh"
                                } }
                            >
                                {["none", "download_file", "download_dir"].includes(this.state.actionStatus) ? (
                                        this.state.items.length == 0 ? (
                                                <span>No items found.</span>
                                            ) : (
                                                this.state.items.map((i) =>(
                                                    <button 
                                                        key={i.abs_path}
                                                        className={this.state.selectedItem == i.abs_path ? (
                                                                "w3-button w3-block w3-left-align w3-green"
                                                            ) : (
                                                                "w3-button w3-block w3-left-align"
                                                            )
                                                        }
                                                        name={i.is_dir ? "open" : "selectedItem"}
                                                        value={i.abs_path}
                                                        disabled={!i.is_dir && (!i.hit || this.selectDir)}
                                                        onClick={i.is_dir ? (
                                                                this.handleAction
                                                            ) : (
                                                                this.changeStateVar
                                                            )
                                                        }
                                                    >
                                                            {i.is_dir ? (
                                                                    <i className="fas fa-folder" />
                                                                ) : (
                                                                    <i className="fas fa-file" />
                                                                )
                                                            }&nbsp;
                                                            {i.name}
                                                    </button>
                                                ))
                                            )
                                    ) : (
                                        <LoadingIndicator
                                            size="large"
                                            message="Loading directory. Please wait."
                                        />
                                    )
                                }
                            </div>
                            {allowUpload && (
                                <div className="w3-container">
                                    <FileUploadComponent
                                        requestRoute={routeUploadFile}
                                        instruction="Upload file:"
                                        buttonLabel="upload"
                                        oneLine={true}
                                        disabled={this.state.actionStatus != "none"}
                                        metaData={ 
                                            {
                                                job_name: this.props.jobName,
                                                dir_path: this.state.dirPath
                                            }
                                        }
                                        showProgress={true}
                                        onUploadCompletion={this.handleFileUpload}
                                    />
                                </div>
                            )}
                            {(allowDownload || (!this.selectDir && allowInput)) && (
                                <div className="w3-container">
                                    <span className="w3-text-green">Selected file:</span>&nbsp;
                                    <IneditableValueField>
                                        {this.state.selectedItem ? this.state.selectedItem :"No file selected."}
                                    </IneditableValueField>
                                </div>
                            )}
                            <div className="w3-bar">
                                {this.props.showCancelButton && (
                                    <ActionButton
                                        name="cancel"
                                        value="cancel"
                                        className="w3-bar-item w3-right"
                                        label={<span><i className="fas fa-times" />&nbsp;cancel</span>}
                                        onAction={this.terminateBrowseDialog}
                                        forwardEvent={true}
                                    />
                                )}
                                {allowInput && (
                                    this.selectDir ? (
                                            <ActionButton
                                                name="select_dir"
                                                value="select_dir"
                                                className="w3-bar-item w3-right"
                                                label={<span><i className="fas fa-check" />&nbsp;choose current dir</span>}
                                                onAction={this.terminateBrowseDialog}
                                                disabled={this.state.actionStatus != "none"}
                                                forwardEvent={true}
                                            />
                                        ) : (
                                            <ActionButton
                                                name="select_file"
                                                value="select_file"
                                                className="w3-bar-item w3-right"
                                                label={<span><i className="fas fa-check" />&nbsp;choose selected file</span>}
                                                disabled={!this.state.selectedItem || this.state.actionStatus != "none"}
                                                onAction={this.terminateBrowseDialog}
                                                forwardEvent={true}
                                            />
                                        )
                                )}
                                {allowDownload && (
                                    <span>
                                        <ActionButton
                                            name="download_dir"
                                            value="download_dir"
                                            className="w3-bar-item w3-right"
                                            label={<span><i className="fas fa-check" />&nbsp;download current dir as zip</span>}
                                            onAction={this.downloadFileOrFolder}
                                            disabled={this.state.actionStatus != "none"}
                                            loading={this.state.actionStatus == "download_dir"}
                                            forwardEvent={true}
                                        />
                                        <ActionButton
                                            name="download_file"
                                            value="download_file"
                                            className="w3-bar-item w3-right"
                                            label={<span><i className="fas fa-check" />&nbsp;download selected file</span>}
                                            disabled={!this.state.selectedItem || this.state.actionStatus != "none"}
                                            onAction={this.downloadFileOrFolder}
                                            loading={this.state.actionStatus == "download_file"}
                                            forwardEvent={true}
                                        />
                                    </span>
                                )}
                            </div>
                            {allowDownload && (
                                <DisplayServerMessages messages={this.state.downloadMessages} />
                            )}
                        </span>
                    )
                }
            </div>
        )
        return( this.props.disableOnTop ? (
                content
            ) : (
                <TopPanel>
                    {content}
                </TopPanel>
            )
        )
    }
}

class BrowseDirTextField extends React.Component {
    constructor(props){
        super(props);
        // props.ignoreFiles
        // props.fileExts
        // props.showOnlyHits
        // props.onChange
        // props.name
        // props.value
        // props.selectDir
        // props.disabled
        // props.allowInput
        // props.allowUpload
        // props.allowDownload
        // props.jobName
        // props.defaultBaseDir
        // props.prevPath
        // props.changePrevPath
        // props.smallSize
        // props.placeholder


        this.state = {
            actionStatus: "none",
            serverMessages: []
        };  

        this.browse = this.browse.bind(this);
        this.terminateBrowseDialog = this.terminateBrowseDialog.bind(this);
    }

    browse(){
        this.setState({
            actionStatus: "browse"
        })
    }

    terminateBrowseDialog(changes, selectedItem){
        this.setState({
            actionStatus: "none"
        })
        if(changes){
            const event = {
                currentTarget: {
                    name: this.props.name,
                    value: selectedItem
                }
            }
            this.props.onChange(event)
            this.props.changePrevPath(selectedItem)
        }
    }

    render(){
        return(
            <span>
                <input
                    className={this.props.smallSize ? "param-input-browse" : "w3-input"}
                    style={ {width: "calc(100% - 28px)", display: "inline-block"} }
                    type="text"
                    name={this.props.name}
                    value={this.props.value}
                    onChange={this.props.onChange}
                    required={true}
                    disabled={this.props.disabled}
                    placeholder={this.props.placeholder}
                />
                <ActionButton
                    name="browse"
                    value="browse"
                    style={ this.props.smallSize ? ({
                            width: "14px", 
                            display: "inline-block", 
                            padding: "1px",
                            fontSize: "80%"
                        }) : ({
                            width: "24px", 
                            display: "inline-block", 
                            paddingLeft: "2px",
                            paddingRight: "2px"
                        })
                    }
                    label={<i className="fas fa-folder-open" />}
                    onAction={this.browse}
                    disabled={this.props.disabled}
                />
                {this.state.actionStatus == "browse" && (
                    <BrowseDir
                        path={this.props.value}
                        ignoreFiles={this.props.ignoreFiles}
                        fileExts={this.props.fileExts}
                        showOnlyHits={this.props.showOnlyHits}
                        selectDir={this.props.selectDir}
                        allowInput={this.props.allowInput}
                        allowUpload={this.props.allowUpload}
                        allowDownload={this.props.allowDownload}
                        jobName={this.props.jobName}
                        defaultBaseDir={this.props.defaultBaseDir}
                        prevPath={this.props.prevPath}
                        changePrevPath={this.props.changePrevPath}
                        showCancelButton={true}
                        terminateBrowseDialog={this.terminateBrowseDialog}
                    />
                )}
            </span>
        )
    }

} 

class FileUploadComponent extends React.Component {
    // diplays content based on server information
    constructor(props) {
        super(props);
        // Inputs:
        // props.requestRoute the route to send the post request to
        // props.instruction short instructive statement
        // props.oneLine if true, everything will be in one line
        // props.diabled if true, upload button disabled
        // props.metaData meta data send together with the file
        // props.onUploadCompletion function to exectute on completion
        // props.buttonLabel
        // props.showProgress

        this.state = {
            status: "wait_for_upload", // can be "wait_for_upload"/"uploading"/"done"
            error: false,
            file: null, 
            serverMessages: []
        };  
        
        this.upload = this.upload.bind(this);
        this.handleFileChange = this.handleFileChange.bind(this);
        this.handleCompletion = this.handleCompletion.bind(this)
    }

    handleFileChange(event){
        this.setState({
            file: event.currentTarget.files[0],
            error: false,
            serverMessages: []
        })
    }

    handleCompletion(isSuccess, data){
        if(this.props.onUploadCompletion){
            this.props.onUploadCompletion(isSuccess, data)
        }
    }

    async upload(value){ // ajax request to server
        const fileToUpload = this.state.file
        if (this.state.fileToUpload == ""){
            this.state.serverMessages = [{type:"error", text:"No file selected."}]
            this.setState({status: "wait_for_upload", error: true})
        }
        else{
            let metaData = this.props.metaData ? this.props.metaData : {}
            const userInfo = await getUserInfo("all")
            metaData["access_token"] = userInfo.accessToken
            metaData["username"] = userInfo.username
            this.setState({status:"uploading"})
            let formData = new FormData()
            formData.append("file", fileToUpload)
            formData.append("meta", JSON.stringify(metaData))

            if (this.props.showProgress){let request = new XMLHttpRequest()
                request.upload.addEventListener("progress", event => {
                    if (event.lengthComputable) {
                        this.setState({
                            serverMessages: [{
                                time: get_time_str(),
                                type: "warning",
                                text: (
                                    "Uploading file: " + 
                                    parseInt((event.loaded / event.total) * 100).toString() + "%"
                                )
                            }]
                        })
                    }
                })
                
                request.open("POST", this.props.requestRoute)
    
                request.onreadystatechange = () => {
                    var status;
                    if (request.readyState == 4) {
                        status = request.status;
                        if (status == 200) {
                            var result = JSON.parse(request.responseText);
                            var messages = result.messages;
                            var data = result.data;
                            let errorOccured = false;
                            for( let i=0;  i<messages.length; i++){
                                if(messages[i].type == "error"){
                                    errorOccured = true;
                                    break;
                                }
                            }
                            if (errorOccured){
                                this.setState({status: "wait_for_upload", error: true, serverMessages: messages})
                                this.handleCompletion(false, data)
                            } 
                            else {
                                this.setState({status: "done", error: false, serverMessages: messages});
                                this.handleCompletion(true, data)
                            }
                        } 
                        else {
                            // server could not be reached
                            var messages = [{
                                time: get_time_str(),
                                type: "error", 
                                text: serverNotReachableError
                            }];
                            this.handleCompletion(false, {})
                            this.setState({status: "wait_for_upload", error: true, serverMessages: messages});
                        }
                    }
                };
    
                request.send(formData)    
            }
            else {
                fetch(this.props.requestRoute, {
                    method: "POST",
                    body: formData,
                    cache: "no-cache"
                }).then(res => res.json())
                .then(
                    (result) => {
                        var messages = result.messages;
                        var data = result.data;
                        let errorOccured = false;
                        for( let i=0;  i<messages.length; i++){
                            if(messages[i].type == "error"){
                                errorOccured = true;
                                break;
                            }
                        }
                        if (errorOccured){
                            this.setState({status: "wait_for_upload", error: true, serverMessages: messages})
                            this.handleCompletion(false, data)
                        } 
                        else {
                            this.setState({status: "done", error: false, serverMessages: messages});
                            this.handleCompletion(true, data)
                        }
                    },
                    (error) => {
                        // server could not be reached
                        var messages = [{
                            time: get_time_str(),
                            type: "error", 
                            text: serverNotReachableError
                        }];
                        this.handleCompletion(false, {})
                        this.setState({status: "wait_for_upload", error: true, serverMessages: messages});
                    }
                )
            }
            
        }
    }

    render() {
        const messages = <DisplayServerMessages messages={this.state.serverMessages} />
        const instruction = this.props.instruction
        const upload_selector = (
            <input 
                className="w3-button w3-border w3-border-grey"
                type="file" 
                name="file" 
                onChange={this.handleFileChange}
                disabled={(this.props.disabled || this.state.status == "uploading") ? true : false}
            />
        )
        const action_button = (
            <ActionButton 
                name="import"
                value="import"
                label={this.props.buttonLabel ? (this.props.buttonLabel) : ("import")}
                loading={this.state.status == "uploading"} 
                onAction={this.upload}
                disabled={this.props.disabled ? true : false}
            />
        )
        if (this.props.oneLine){
            return (
                <table>
                    <tbody>
                        <tr>
                            {this.props.instruction && (
                                <td>
                                    {this.props.instruction}
                                </td>
                            )}
                            <td>
                                {upload_selector}
                            </td>
                            <td>
                                {action_button}
                            </td>
                            <td>
                                {messages}
                            </td>
                        </tr>
                    </tbody>
                </table>
            )
        }
        else {
            return (
                <div>
                    <div>
                        {this.props.instruction && (
                            <span>
                                {this.props.instruction}<br/>
                            </span>
                        )}
                        {upload_selector}{action_button}
                    </div>
                    <div>
                        {messages}
                    </div>
                </div>
            )

        }

    }
}

class AjaxButton extends React.Component {
    // diplays content based on server information
    constructor(props) {
        super(props);
        // Inputs:
        // props.name
        // props.route the route to send the post request to
        // props.instruction short instructive statement
        // props.oneLine if true, everything will be in one line
        // props.diabled if true, upload button disabled
        // props.sendData
        // props.onSuccess 
        // props.onError 
        // props.label
        // props.disabled
        // props.className

        this.state = {
            status: "none",
            serverMessages: []
        };  
                                    
        this.handleOnClick = this.handleOnClick.bind(this);
        this.ajaxRequest = ajaxRequest.bind(this);
    }

    async handleOnClick(value){
        this.ajaxRequest({
            statusVar: "status",
            statusValueDuringRequest: "loading",
            messageVar: "serverMessages",
            sendData: this.props.sendData,
            route: this.props.route,
            onSuccess: this.props.onSuccess,
            onError: this.props.onError
        })        
    }

    render() {
        return(
            <div>
                <ActionButton
                    name={this.props.name}
                    value={this.props.name}
                    className={this.props.className}
                    label={this.props.label}
                    loading={this.state.status == "loading"}
                    disabled={this.props.disabled || this.state.status == "loading"}
                    onAction={this.handleOnClick}
                />
                <DisplayServerMessages messages={this.state.serverMessages} />
            </div>
        )
    }
}