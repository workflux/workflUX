// contains general utilities important for multiple modules
// import styled from "styled-components"


function ajaxRequest({
    // in a component bind this function: this.ajaxRequest = ajaxRequest.bind(this)
    statusVar="actionStatus",
    statusValueDuringRequest="action",
    statusValueAfterRequest="none",
    messageVar="serverMessages",
    sendData={},
    route,
    onSuccess= (data, messages) => {
                return({})
            }, //function taking arguments data and messages, return state update
    onError= (messages) => { //function taking argument messages, return state update
        return({})
    } 
}){
    this.setState({
        [statusVar]: statusValueDuringRequest
    })
    fetch(route, {
        method: "POST",
        body: JSON.stringify(sendData),
        headers: new Headers({
            'Content-Type': 'application/json'
        }),
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
            const messages = [{type: "error", text: serverNotReachableError}];
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
    jobId=null,
    runId=null,
    onSuccess= (data, messages) => {
                return({})
            },
    onError= (messages) => {
        return({})
    } 
}){
    let sendData = {}
    if (jobId){
        sendData["job_id"] = jobId
        if (runId){
            sendData["run_id"] = runId
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
                    <div 
                        className="w3-padding-small"
                        style={ {display: "inline-block"} }
                    >
                        {this.props.title}
                    </div>
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
        const isSet=event.target.checked
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
        this.handleChange = this.handleChange.bind(this)
    }
    
    handleChange(event){
        const value=event.currentTarget.value
        const isSet=event.target.checked
        if (this.props.doNotSendValue){
            this.props.onChange(isSet)
        }
        else {
            this.props.onChange(value, isSet)
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
                target:{
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
        let newSelection = this.props.selection
        if (this.props.selection.includes(event.target.value)){
            delete newSelection[newSelection.indexOf(event.target.value)]
        }
        else {
            newSelection.push(event.target.value)
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
                            <Message type={message.type}>{message.text}</Message>
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

    request(){ // ajax request to server
        this.ajaxRequest({
            statusVar: "loading",
            statusValueDuringRequest: true,
            statusValueAfterRequest: false,
            messageVar: "serverMessages",
            sendData: this.props.sendData,
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
        // props.jobId
        // props.defaultBaseDir

        this.baseDirInfo = {
            input: "select only",
            upload: (
                (this.props.allowInput ? "select, " : "") + 
                "upload"
            ),
            download: (
                (this.props.allowInput ? "select, " : "") + 
                (this.props.allowUpload ? "upload, " : "") + 
                "download"
            ),
        }

        this.state = {
            dirPath: this.props.path,
            address: this.props.path,
            actionStatus: "init",
            serverMessages: [],
            items: [],
            selectedItem: this.props.path,
            baseDir: "error",
            allowedDirs: {error: {path: "error", mode: "error"}},
        };  

        this.ajaxRequest = ajaxRequest.bind(this);
        this.getItemsInDir = this.getItemsInDir.bind(this);
        this.handleAction = this.handleAction.bind(this);
        this.changeStateVar = this.changeStateVar.bind(this);
        this.changeBaseDir = this.changeBaseDir.bind(this);
    }

    componentDidMount(){
        let path = ""

        if (this.props.path && this.props.path != "Please fill" && this.props.path != ""){
            path = this.props.path
        }
        else if (this.props.prevPath){
            path = this.props.prevPath
        }
        console.log(path)
        this.getItemsInDir(false, path, true)
    }

    getItemsInDir(getParentDir, targetDir, init){
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
                allow_input: this.props.allowInput ? true : false,
                allow_upload: this.props.allowUpload ? true : false,
                allow_download: this.props.allowDownload ? true : false,
                job_id: this.props.jobId ? this.props.jobId : null,
                on_error_return_base_dir_items: init ? true : false,
                default_base_dir: this.props.defaultBaseDir ? this.props.defaultBaseDir : null
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
        if (event.target.name == "up"){
            this.getItemsInDir(true, this.state.dirPath)
        }
        else if (event.target.name == "refresh"){
            this.getItemsInDir(false, this.state.dirPath)
        }
        else if (event.target.name == "open"){
            this.getItemsInDir(false, event.target.value)
        }
        else if (event.target.name == "address"){
            if (event.key == "Enter"){
                this.getItemsInDir(false, event.target.value)
            }
        }
    }

    changeStateVar(event){
        this.setState({
            [event.target.name]: event.target.value
        })
    }

    changeBaseDir(event){
        this.setState({
            baseDir: event.target.value,
            address: this.state.allowedDirs[event.target.value]
        })
        this.getItemsInDir(false, this.state.allowedDirs[event.target.value].path)
    }

    render(){
        return(
            <div>
                <div
                    style={ {
                        position: "fixed",  
                        zIndex: "10002",
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
                                        {this.state.actionStatus == "none" ? (
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
                                                                disabled={!i.is_dir && (!i.hit || this.props.selectDir)}
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
                                </span>
                            )
                        }
                        
                    </div>
                </div>
                <div
                    className="w3-white"
                    style={ {
                        position: "fixed",  
                        zIndex: "10000",
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
        // props.jobId
        // props.defaultBaseDir

        this.state = {
            actionStatus: "none",
            serverMessages: []
        };  

        this.handleOnChange = this.handleOnChange.bind(this);
        this.browse = this.browse.bind(this);
    }

    handleOnChange(event){
        this.props.onChange(event.target.name, event.target.value)
    }

    browse(){
        this.setState({
            actionStatus: "browse"
        })
    }

    render(){
        return(
            <span>
                <input
                    className="w3-input"
                    style={ {width: "calc(100% - 28px)", display: "inline-block"} }
                    type="text"
                    name={this.props.name}
                    value={this.props.value}
                    onChange={this.handleOnChange}
                    required={true}
                    disabled={this.props.disabled}
                />
                <ActionButton
                    name="browse"
                    value="browse"
                    style={ {
                        width: "24px", 
                        display: "inline-block", 
                        paddingLeft: "2px",
                        paddingRight: "2px"
                    } }
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
                        jobId={this.props.jobId}
                        defaultBaseDir={this.props.defaultBaseDir}
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
        // props.meta_data meta data send together with the file
        // props.onUploadCompletion function to exectute on completion, 
        //  takes one argument: true (on success)/ false (on error)
        //props.buttonLabel

        this.state = {
            status: "wait_for_upload", // can be "wait_for_upload"/"uploading"/"done"
            error: false
        };  

        this.serverMessages = [];    // container for errors/warnings/infos 
                                    //from server upon ajax request
                                    
        this.file = "" // container to store file object
        
        this.upload = this.upload.bind(this);
        this.handleFileChange = this.handleFileChange.bind(this);
        this.handleCompletion = this.handleCompletion.bind(this)
    }

    handleFileChange(event){
        this.file = event.target.files[0]
    }

    handleCompletion(isSuccess){
        if(this.props.onUploadCompletion){
            this.props.onUploadCompletion(isSuccess)
        }
    }

    upload(value){ // ajax request to server
        const fileToUpload = this.file
        if (this.fileToUpload == ""){
            this.serverMessages = [{type:"error", text:"No file selected."}]
            this.setState({status: "wait_for_upload", error: true})
        }
        else{
            this.setState({status:"uploading"})
            let formData = new FormData()
            formData.append("file", fileToUpload)
            if (this.props.meta_data){
                formData.append("meta", JSON.stringify(this.props.meta_data))
            }
            fetch(this.props.requestRoute, {
                method: "POST",
                body: formData,
                cache: "no-cache"
            }).then(res => res.json())
            .then(
                (result) => {
                    this.serverMessages = result.messages;
                    let errorOccured = false;
                    for( let i=0;  i<this.serverMessages.length; i++){
                        if(this.serverMessages[i].type == "error"){
                            errorOccured = true;
                            break;
                        }
                    }
                    if (errorOccured){
                        this.setState({status: "wait_for_upload", error: true})
                        this.handleCompletion(false)
                    } 
                    else {
                        this.setState({status: "done", error: false});
                        this.handleCompletion(true)
                    }
                },
                (error) => {
                    // server could not be reached
                    this.serverMessages = [{type: "error", text: serverNotReachableError}];
                    this.handleCompletion(false)
                    this.setState({status: "wait_for_upload", error: true});
                }
            )
        }

        
    }

    render() {
        status = this.state.status
        const messages = <DisplayServerMessages messages={this.serverMessages} />
        const instruction = this.props.instruction
        const upload_selector = (
            <input 
            className="w3-button w3-border w3-border-grey"
            type="file" 
            name="file" 
            onChange={this.handleFileChange}
            disabled={this.props.disabled ? true : false}
            />
        )
        const action_button = (
            <ActionButton 
                name="import"
                value="import"
                label={this.props.buttonLabel ? (this.props.buttonLabel) : ("import")}
                loading={status == "uploading"} 
                onAction={this.upload}
                disabled={this.props.disabled ? true : false}
            />
        )
        if (this.props.oneLine){
            return (
                <span>
                    {instruction}&nbsp;
                    {upload_selector}&nbsp;
                    {action_button}
                    {messages}
                </span>
            )
        }
        else {
            return (
                <div>
                    <div>
                        <p>
                            {instruction}<br/>
                            {upload_selector}
                        </p>
                    </div>
                    <div>
                        {action_button}
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

    handleOnClick(value){
        this.ajaxRequest({
            statusVar: "status",
            statusValueDuringRequest: "loading",
            messageVar: "serverMessages",
            sendData: this.props.sendData,
            route: this.props.route,
            onSuccess: this.props.onSuccess,
            onError: this.props.onError,
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