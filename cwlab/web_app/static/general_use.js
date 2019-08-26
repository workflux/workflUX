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

String.prototype.replaceAll = function(search, replacement) {
    // in contrast to replace() replaces not only the first occurence
    var target = this;
    return target.split(search).join(replacement);
};

class IneditableValueField extends React.Component {
    // Input:
    // props.backColorClass e.g. w3-metro-darken w3-theme-dark
    // props.textColorClass e.g. w3-text-green
    render() {
        const style = {
            fontFamily: "courier",
            padding: "2px"
        }
        let backColorClass = "w3-metro-darken"
        let textColorClass = ""
        if (this.props.backColorClass){
            backColorClass = this.props.backColorClass
        }
        if (this.props.textColorClass){
            textColorClass = this.props.textColorClass
        }
        return( 
            <span 
                className={backColorClass + " " + textColorClass} 
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
        const is_set=event.target.checked
        this.props.onChange(value, is_set)
    }

    render(){
        return(
            <input
                className="w3-check"
                style={ {height: this.props.size ? (this.props.size) : ("20px")} }
                type="checkbox"
                name={this.props.name}
                value={this.props.name}
                checked={!this.props.isNull}
                onChange={this.handleToggleNull}
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
        this.handleChange = this.handleChange.bind(this)
    }
    
    handleChange(event){
        const value=event.currentTarget.value
        const is_set=event.target.checked
        this.props.onChange(value, is_set)
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
        // props.colorClass set color class default is w3-black
        // props.onAction function to execute on click
        // props.label labeling
        // props.loading if true, button will be disabled and loading indicator is shown
        // props.disabled if true, disable without loading indicator
        // props.smallPadding true/false
        this.handleAction = this.handleAction.bind(this)
    }
    
    handleAction(event){
        this.props.onAction(event.currentTarget.value)
    }

    render() {
        const style = {marginTop: "2px", marginBottom: "2px"}
        return( 
            <button style={style}
                name={this.props.name}
                value={this.props.value}
                className={
                    "w3-button" +
                    (this.props.colorClass ? (" " + this.props.colorClass) : " w3-black") +
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
                    style={ {width:"200px", overflowY: "auto"} }
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
                    className="w3-col s10 m8 s8 w3-container"  
                    style={ {marginLeft:"200px", width: "calc(100% - 200px)"} }
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
            hint: "w3-yellow"
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
                    <DisplayServerMessages messages={this.state.serverMessages} />
                    {this.state.success &&
                        this.props.buildContentOnSuccess(this.state.data, this.state.serverMessages, this.request)
                    }
                </div>
            )
        }
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
                label="import" 
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