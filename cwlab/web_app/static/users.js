
const inputStyle = {width: "100%", maxWidth: "400px", minWidth: "100px"}

function changeInputField(event){
    this.setState({[event.currentTarget.name]: event.currentTarget.value})
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

class AdminDashboard extends React.Component {
    constructor(props){
        super(props);

        this.state = {
            actionStatus: "",
            serverMessages: [],
            actionMessages: [],
            userInfo: [],
            userSelection: [],
            userFilter: "all",
            unlockDelete: false,
            selectStatus: "active",
            selectLevel: "user"
        }

        this.columnNames = {
            username: "Username",
            email: "Email",
            status: "Status",
            level: "Level",
            date_register: "Registr. Date",
            date_last_login: "Last Login",
        }

        this.ajaxRequest = ajaxRequest.bind(this)
        this.handleSelectionChange = this.handleSelectionChange.bind(this)
        this.modifyOrDeleteUser = this.modifyOrDeleteUser.bind(this)
        this.changeUserFilter = this.changeUserFilter.bind(this)
        this.toggleUnlockDelete = this.toggleUnlockDelete.bind(this)
        this.changeSelect = this.changeSelect.bind(this)
        this.getUserInfo = this.getUserInfo.bind(this)
    }

    handleSelectionChange(newSelection){
        this.setState({userSelection: newSelection})
    }

    changeSelect(event){
        this.setState({[event.currentTarget.name]: event.currentTarget.value})
    }
    
    changeUserFilter(event){
        this.setState({
            userSelection: [],
            userFilter: event.currentTarget.value
        })
    }

    modifyOrDeleteUser(action){
        this.ajaxRequest({
            route: routeModifyOrDeleteUsers,
            statusVar: "actionStatus",
            statusValueDuringRequest: action,
            messageVar: "actionMessages",
            sendData: {
                action: action,
                user_selection: this.state.userSelection,
                value: action == "set_status" ? this.state.selectStatus : this.state.selectLevel
            },
            onSuccess: (data, messages) => {
                this.getUserInfo()
                return({userSelection: [], unlockDelete: false})
            }
        })
    }

    toggleUnlockDelete(value, isSet){

        this.setState({unlockDelete: isSet})
    }

    componentDidMount(){
        this.getUserInfo()
    }

    getUserInfo(){
        this.ajaxRequest({
            route: routeGetAllUsersInfo,
            statusValueDuringRequest: "loading",
            onSuccess: (data, messages) => {
                return({userInfo: data})
            }
        })
    }


    render() {
        const disableActionButtons = this.state.userSelection.length == 0 || this.state.actionStatus != "none"
        const disabledDelete = disableActionButtons || !this.state.unlockDelete
        let rowData
        if (this.state.userFilter == "all"){
            rowData = this.state.userInfo
        }
        else if (this.state.userFilter == "only_inactive"){
            rowData = this.state.userInfo.filter( (r) => r.status != "active")
        }
        else if (this.state.userFilter == "only_active"){
            rowData = this.state.userInfo.filter( (r) => r.status == "active")
        }
        else if (this.state.userFilter == "only_admins"){
            rowData = this.state.userInfo.filter( (r) => r.level == "admin")
        }
        else if (this.state.userFilter == "only_users"){
            rowData = this.state.userInfo.filter( (r) => r.level == "user")
        }
        return (
            <div>
                <h3>Administrator Dashboard</h3>
                <p>
                    Your are logged in as administrator. 
                    Therefore, you may manage accounts of other users.
                </p>
                Commen tasks are:
                <ul style={ {listStyleType: "circle"} }>
                    <li>Approving the accounts of newly registered users.</li>
                    <li>Deactivating or deleting accounts of users that should no longer have access.</li>
                    <li>Granting privileged admin rights to trusted users.</li>
                </ul>
                
                
                <h3>List of Users:</h3>
                <p>
                    <span className="w3-text-green">Select which users to show:</span>&nbsp;
                    <select className="w3-button w3-white w3-border" 
                        name="userFilter"
                        onChange={this.changeUserFilter}
                        value={this.state.userFilter}
                    >
                        <option value="all">all</option>
                        <option value="only_inactive">only inactive users</option>
                        <option value="only_active">only active users</option>
                        <option value="only_admins">only admins</option>
                        <option value="only_users">only users</option>
                    </select>
                </p>
                {this.state.actionStatus == "updating" ? (
                        <LoadingIndicator
                            size="large"
                            message="Loading user info. Please wait."
                        />
                    ) : (
                        <Table 
                            columnKeys={Object.keys(this.columnNames)}
                            columnNames={this.columnNames}
                            selectionEnabled={true}
                            handleSelectionChange={this.handleSelectionChange}
                            selection={this.state.userSelection}
                            rowData={rowData}
                            selectionKey="username"
                        />
                    )
                }
                
                <DisplayServerMessages messages={this.state.serverMessages} />

                <Message type="hint">
                    <div>
                        <h4>Hints:</h4>
                        <p>
                            <b>Status:</b>&nbsp;
                            Only users with an "active" status can log in. 
                            Users with the status "awaiting approval" have newly registered.
                            Users with the status "inactive" have been deactivated by an administrator.
                        </p>
                        <p>
                            <b>Level:</b>&nbsp;
                            The level can be "user" or "admin" which have privileges like accessing the admin dashboard.
                        </p>
                    </div>
                </Message>

                <h3>Actions on Selected Users:</h3>
                {this.state.userSelection.length == 0 &&
                    <Message type="hint">
                        Please select one or multiple users in the above list to enable following actions.
                    </Message>
                }
                <div style={ {width: "100%"} }>
                    <div
                        className="w3-padding-small"
                        style={ {display: "inline-block"} }
                    >
                        <span className="w3-text-green">Set status to:</span><br/>
                        <select className="w3-button w3-white w3-border" 
                            name="selectStatus"
                            onChange={this.changeSelect}
                            disabled={disableActionButtons}
                            value={this.state.selectStatus}
                        >
                            <option value="active">active</option>
                            <option value="inactive">inactive</option>
                        </select>
                        <ActionButton
                            name="set_status"
                            value="set_status"
                            label="set"
                            disabled={disableActionButtons}
                            loading={this.state.actionMessages == "set_status"}
                            onAction={this.modifyOrDeleteUser}
                        />
                    </div>
                    <div
                        className="w3-padding-small"
                        style={ {display: "inline-block"} }
                    >
                        <span className="w3-text-green">Set level to:</span><br/>
                        <select className="w3-button w3-white w3-border" 
                            name="selectLevel"
                            disabled={disableActionButtons}
                            onChange={this.changeSelect}
                            value={this.state.selectLevel}
                        >
                            <option value="user">user</option>
                            <option value="admin">admin</option>
                        </select>
                        <ActionButton
                            name="set_level"
                            value="set_level"
                            label="set"
                            disabled={disableActionButtons}
                            loading={this.state.actionMessages == "set_level"}
                            onAction={this.modifyOrDeleteUser}
                        />
                    </div>
                    <div
                        className="w3-padding-small"
                        style={ 
                            {
                                display: "inline-block", 
                            }
                        }
                    >   
                        <span className="w3-text-red">Delete users:</span><br/>
                        <div 
                            className="w3-container"
                            style={ {backgroundColor: disabledDelete ? "hsl(0, 20%, 50%)" : "hsl(0, 40%, 50%)"}}
                        >
                            <BooleanSlider
                                name="unlock_delete"
                                value="unlock_delete"
                                onChange={this.toggleUnlockDelete}
                                disabled={disableActionButtons}
                                checked={this.state.unlockDelete}
                            /> &nbsp; unlock&nbsp;
                            <ActionButton
                                name="delete"
                                value="delete"
                                label="delete"
                                disabled={disabledDelete}
                                loading={this.state.actionMessages == "delete"}
                                onAction={this.modifyOrDeleteUser}
                            />
                        </div>
                    </div>
                </div>                
                
                <DisplayServerMessages messages={this.state.actionMessages} />
                
            </div>
        );
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
        
        this.itemValues = userLevel == "admin" ? (
                ["general_info", "admin_dashboard", "change_password", "delete_account", "logout"]
            ):(
                ["general_info", "change_password", "delete_account", "logout"]
            )
        this.itemNames = userLevel == "admin" ? (
                [
                    <span><i className="fas fa-user"/>&nbsp;General Info</span>,
                    <span><i className="fas fa-users"/>&nbsp;Admin Dashboard</span>,
                    <span><i className="fas fa-lock"/>&nbsp;Change Password</span>,
                    <span><i className="fas fa-trash-alt"/>&nbsp;Delete Account</span>,
                    <span><i className="fas fa-sign-out-alt"/>&nbsp;Logout</span>
                ]
            ) : (
                [
                    <span><i className="fas fa-user"/>&nbsp;General Info</span>,
                    <span><i className="fas fa-lock"/>&nbsp;Change Password</span>,
                    <span><i className="fas fa-trash-alt"/>&nbsp;Delete Account</span>,
                    <span><i className="fas fa-sign-out-alt"/>&nbsp;Logout</span>
                ]
            )
            
        this.itemContents = {
            general_info: <GeneralInfo />,
            admin_dashboard: <AdminDashboard />,
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