
const inputStyle = {width: "100%", maxWidth: "400px", minWidth: "100px"}

function changeInputField(event){
    this.setState({[event.target.name]: event.target.value})
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