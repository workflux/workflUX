
const inputStyle = {width: "100%", maxWidth: "400px", minWidth: "100px"}

function changeInputField(event){
    this.setState({[event.target.name]: event.target.value})
}

class GeneralInfo extends React.Component {
    constructor(props){
        super(props);

        this.labels = {
            username: "Username",
            email: "Email",
            level: "Level"
        }

        this.buildContentOnSuccess = this.buildContentOnSuccess.bind(this)
    }

    buildContentOnSuccess(data, serverMessages, request){
        return(
            <div>
                <h3>General User Account Information:</h3>
                {Object.keys(this.labels).map((key) => (
                    <p key={key}>
                        <span className="w3-text-green">{this.labels[key]}:</span>&nbsp;{data[key]}
                    </p>
                ))}
            </div>
        )
    }

    render(){
        return(
            <AjaxComponent
                requestRoute={routeGetGeneralUserInfo}
                buildContentOnSuccess={this.buildContentOnSuccess}
                loaderMessage="Loading user information."
            />
        )

    }
}

class ChangePassword extends React.Component {
    constructor(props){
        super(props);

        this.state = {
            oldPassword: "",
            newPassword: "",
            repNewPassword: "",
            actionStatus: "none",
            serverMessages: [],
            whichFocus: "login"
        }

        this.inputFieldLabels = {
            oldPassword: "Old password",
            newPassword: "New password",
            repNewPassword: "Repeat new password",
        }

        this.ajaxRequest = ajaxRequest.bind(this)
        this.changeInputField = changeInputField.bind(this)
        this.changePassword = this.changePassword.bind(this)
    }

    changePassword(){
        this.ajaxRequest({
            sendData: {
                old_password: this.state.oldPassword,
                new_password: this.state.newPassword,
                new_rep_password: this.state.repNewPassword
            },
            route: routeChangePassword,
            onSuccess: (data, messages) => {
                if (data.success){
                    window.location.reload(true)
                }
            }
        })
    }

    render(){
        return(
            <div>
                <h3>Change Password</h3>

                {Object.keys(this.inputFieldLabels).map((key) => (
                    <p key={key}>
                        <span className="w3-text-green">{this.inputFieldLabels[key]}:</span>
                        <input type="password"
                            className="w3-input w3-border"
                            name={key}
                            style={inputStyle}
                            value={this.state[key]}
                            onChange={this.changeInputField}
                        />
                    </p>
                ))}
                
                <ActionButton
                    name="change_password"
                    value="change_password"
                    label= "change password and logout"
                    loading={this.state.actionStatus != "none"}
                    onAction={this.changePassword}
                />

                <DisplayServerMessages messages={this.state.serverMessages} />
            </div>
        )
    }
}

class DeleteAccount extends React.Component {
    constructor(props){
        super(props);

        this.state = {
            username: "",
            actionStatus: "none",
            serverMessages: []
        }

        this.deleteAccount = this.deleteAccount.bind(this)
        this.ajaxRequest = ajaxRequest.bind(this)
        this.changeInputField = changeInputField.bind(this)
    }

    deleteAccount(){
        this.ajaxRequest({
            route: routeDeleteAccount,
            sendData: {
                username: this.state.username
            },
            onSuccess: (data, messages) => {
                if (data.success){
                    window.location.reload(true)
                }
            }
        })
    }

    render(){
        return(
            <div>
                <h3>Delete Your User Account</h3>
                <Message type="warning">
                    Attention: This will delete your user account. This action can not be reverted.
                </Message>
                <p>
                    Please enter your username and confirm.
                </p>
                <p>
                    <span className="w3-text-green">Username:</span>
                    <input type="text"
                        className="w3-input w3-border"
                        name="username"
                        style={inputStyle}
                        value={this.state.username}
                        onChange={this.changeInputField}
                    />
                </p>

                <ActionButton
                    name="delete_account"
                    value="delete_account"
                    label= "confirm deletion of account"
                    loading={this.state.actionStatus != "none"}
                    onAction={this.deleteAccount}
                />

                <DisplayServerMessages messages={this.state.serverMessages} />
            </div>
        )
    }
}

class Logout extends React.Component {
    constructor(props){
        super(props);

        this.state = {
            actionStatus: "loading",
            serverMessages: []
        }

        this.logout = this.logout.bind(this)
        this.ajaxRequest = ajaxRequest.bind(this)
    }

    
    componentDidMount() {
        this.logout()
    }

    logout(){
        this.ajaxRequest({
            route: routeLogout,
            onSuccess: (data, messages) => {
                if (data.success){
                    window.location.reload(true)
                }
            }
        })
    }

    render(){
        return(
            <div>
                { this.state.serverMessages.length == 0 && (
                    <LoadingIndicator
                        message="Logging you out."
                        size="large"
                    />
                )}
                <DisplayServerMessages messages={this.state.serverMessages} />
            </div>
        )
    }
}

class UserAccount extends React.Component {
    constructor(props) {
        super(props);
        
        this.itemValues = ["general_info", "change_password", "delete_account", "logout"]
        this.itemNames = [
            <span><i className="fas fa-user"/>&nbsp;General Info</span>,
            <span><i className="fas fa-lock"/>&nbsp;Change Password</span>,
            <span><i className="fas fa-trash-alt"/>&nbsp;Delete Account</span>,
            <span><i className="fas fa-sign-out-alt"/>&nbsp;Logout</span>
        ]
        this.itemContents = {
            general_info: <GeneralInfo />,
            change_password: <ChangePassword />,
            delete_account: <DeleteAccount />,
            logout: <Logout />,
        }

        this.state = {
            whichFocus: "general_info"
        }

        this.changeFocus = this.changeFocus.bind(this)
    }

    changeFocus(newFocus){
        this.setState({
            whichFocus: newFocus
        })
    }

    render() {
        return(
            <SideBarPanel
                itemValues = {this.itemValues}
                itemNames = {this.itemNames}
                itemContent = {this.itemContents[this.state.whichFocus]}
                whichFocus = {this.state.whichFocus}
                onChange = {this.changeFocus}
                label = "Options:"
            />
        )
    }
}

class LoginForm extends React.Component {
    constructor(props) {
        super(props);
        
        this.itemValues = ["login", "register"]
        this.itemNames = [
            <span><i className="fas fa-sign-in-alt"/>&nbsp;Login</span>,
            <span><i className="fas fa-user-plus"/>&nbsp;Register</span>
        ]

        this.state = {
            username: "",
            email: "",
            password: "",
            repPassword: "",
            rememberMe: false,
            actionStatus: "none",
            serverMessages: [],
            whichFocus: "login"
        }

        this.changeInputField = changeInputField.bind(this)
        this.changeRememberMe = this.changeRememberMe.bind(this)
        this.login = this.login.bind(this)
        this.register = this.register.bind(this)
        this.ajaxRequest = ajaxRequest.bind(this)
        this.changeFocus = this.changeFocus.bind(this)
    }

    changeRememberMe(value, isSet){
        this.setState({rememberMe: isSet})
    }

    changeFocus(newFocus){
        this.setState({
            whichFocus: newFocus,
            password: "",
            repPassword: "",
            serverMessages: []
        })
    }

    login(){
        this.ajaxRequest({
            sendData: {
                username: this.state.username,
                password: this.state.password,
                remember_me: this.state.rememberMe
            },
            route: routeLogin,
            onSuccess: (data, messages) => {
                if (data.success){
                    window.location.reload(true)
                }

            }
        })
    }

    register(){
        this.ajaxRequest({
            sendData: {
                username: this.state.username,
                email: this.state.email,
                password: this.state.password,
                rep_password: this.state.repPassword
            },
            route: routeRegister
        })
    }

    render() {
        const instruction = {
            login: loginInstruction,
            register: registrationInstruction
        }
        const itemContent = (
            <div>
                <h3>
                    { {
                        login: "Login",
                        register: "Register"
                    }[this.state.whichFocus] }:
                </h3>
                
                { 
                    instruction[this.state.whichFocus] != "" && (
                        <p>
                            {instruction[this.state.whichFocus]}
                        </p>
                    )
                }

                <p>
                    <span className="w3-text-green">Username:</span>
                    <input type="text"
                        className="w3-input w3-border"
                        name="username"
                        style={inputStyle}
                        value={this.state.username}
                        onChange={this.changeInputField}
                    />
                </p>

                { this.state.whichFocus == "register" && (
                    <p>
                        <span className="w3-text-green">Email:</span>
                        <input type="email"
                            className="w3-input w3-border"
                            name="email"
                            style={inputStyle}
                            value={this.state.email}
                            onChange={this.changeInputField}
                        />
                    </p>
                )}
                
                <p>
                    <span className="w3-text-green">Password:</span>
                    <input type="password"
                        className="w3-input w3-border"
                        name="password"
                        style={inputStyle}
                        value={this.state.password}
                        onChange={this.changeInputField}
                    />
                </p>

                { this.state.whichFocus == "register" && (
                    <p>
                        <span className="w3-text-green">Repeat password:</span>
                        <input type="password"
                            className="w3-input w3-border"
                            name="repPassword"
                            style={inputStyle}
                            value={this.state.repPassword}
                            onChange={this.changeInputField}
                        />
                    </p>
                )}
                
                { this.state.whichFocus == "login" && (
                    <p>
                        <label>
                            <Checkbox
                                name="remember_me"
                                value="remember_me"
                                checked={this.state.rememberMe}
                                onChange={this.changeRememberMe}
                            />
                            Remember me
                        </label>
                    </p>
                )}

                <ActionButton
                    name="login"
                    value="login"
                    label= { 
                        {
                            login: "login",
                            register: "request registration"
                        }[this.state.whichFocus] 
                    }
                    loading={this.state.actionStatus != "none"}
                    onAction={ 
                        {
                            login: this.login,
                            register: this.register
                        }[this.state.whichFocus] 
                    }
                />

                <DisplayServerMessages messages={this.state.serverMessages} />
            </div>
        )

        return(
            <SideBarPanel
                itemValues = {this.itemValues}
                itemNames = {this.itemNames}
                itemContent = {itemContent}
                whichFocus = {this.state.whichFocus}
                onChange = {this.changeFocus}
                label = "Options:"
            />
        )
    }
}

class UserRoot extends React.Component{
    constructor(props){
        super(props);
    }

    render(){
        if (loggedIn){
            return(<UserAccount />)
        }
        else{
            return(<LoginForm />)
        }
    }
}