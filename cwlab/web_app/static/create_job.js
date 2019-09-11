
class CreateJobButton extends React.Component {
    constructor(props){
        super(props);
        // props.jobId
        // props.sheetFormat
        // props.disabled
        // props.validatePaths
        // props.searchPaths
        // props.searchDir
        // props.includeSubbDirsForSearching
    }

    render(){
        return(
            <AjaxButton
                name="create job"
                label="create job"
                disabled={this.props.disabled}
                sendData={
                    {
                        job_id: this.props.jobId,
                        sheet_format: this.props.sheetFormat,
                        validate_paths: this.props.validatePaths,
                        search_paths: this.props.searchPaths,
                        search_dir: this.props.searchDir,
                        include_subdirs_for_searching: this.props.includeSubbDirsForSearching
                    }
                }
                route={routeCreateJob}
            />
        )
    }
}

class PathValAndSearch extends React.Component{
    constructor(props){
        super(props);
        // props.validatePaths
        // props.searchPaths
        // props.searchDir
        // props.changeSearchDir
        // props.includeSubbDirsForSearching
        // props.changeIncludeSubDirsForSearching
        // props.changePathValAndSearch
        this.handleChangeSearchDir = this.handleChangeSearchDir.bind(this);
        this.handleChangeValAndSearch = this.handleChangeValAndSearch.bind(this);
    }

    handleChangeValAndSearch(event){
        const searchPaths = event.currentTarget.value == "val_search"
        const validatePaths = event.currentTarget.value == "val_search" || event.currentTarget.value == "val"
        this.props.changePathValAndSearch(searchPaths, validatePaths)
    }

    handleChangeSearchDir(event){
        this.props.changeSearchDir(event.currentTarget.value)
    }

    render(){
        const validate_or_search = this.props.validatePaths ? (
            this.props.searchPaths ? (
                "val_search"
            ) : (
                "val"
            )
        ) : (
            "no"
        )

        const info_and_advanced_options = {
            "val_search": (
                <span>
                    <br/>
                    With this option, you may:
                    <ul>
                        <li>
                            specify absolute paths for a file or directory parameter (paths are validated)
                        </li>
                        <li>
                            or provide a character string that uniquely identifies 
                            the file/directory name in the following search directory 
                            (globbing patterns are allowed)
                        </li>
                    </ul>
                    <span className="w3-text-green">Search directory: </span>
                    <input
                        className="w3-input"
                        type="text"
                        name={"select_input_dir"}
                        value={this.props.searchDir}
                        onChange={this.handleChangeSearchDir}
                        required={true}
                    />
                    <br/>
                    <span className="w3-text-green">Include sub-directories for searching: </span>
                    no &nbsp;
                    <BooleanSlider
                        name="include_subdirs_for_searching"
                        value="include_subdirs_for_searching"
                        onChange={this.props.changeIncludeSubDirsForSearching}
                        checked={this.props.includeSubbDirsForSearching}
                        doNotSendValue={true}
                    />
                    &nbsp; yes
                </span>

            ),
            "val": (
                <p>
                    With this option, you have to specify absolute paths for file or directory parameters.
                    The paths will be validated.
                </p>
            ),
            "no": (
                <p>
                    With this option, you have to specify absolute paths for file or directory parameters.
                    The paths will not be validated.
                </p>
            )

        }
        return(
            <div>
                <h3>Path Validation and Searching:</h3>
                <p>
                    Please select how paths of file or directory parameters should be treated.
                </p>
                <span className="w3-text-green">Options: </span>
                <select className="w3-button w3-white w3-border" 
                    name="path_validation_and_search"
                    onChange={this.handleChangeValAndSearch}
                    value={validate_or_search}
                    >
                    <option value="val_search">validate and search</option>
                    <option value="val">only validate</option>
                    <option value="no">do not validate or search</option>
                </select>
                <div className="w3-container">
                    {info_and_advanced_options[validate_or_search]}
                </div>
            </div>
        )
    }
}

class ParamName extends React.Component{
    constructor(props){
        super(props);
        // props.name
    }

    render(){
        return(
            <span className="w3-text-green" style={ {display: "inline-block"} }>
                {this.props.name}:
            </span>
        )
    }
}

class ParamField extends React.Component{
    constructor(props){
        super(props);
        // props.type
        // props.itemNullAllowed
        // props.onChange
        // props.paramValue
        // props.name
        // props.index for array params
        // props.isNull

        this.key = this.props.index ? (
            this.props.name + "%" + this.props.index.toString()
        ) :(
            this.props.name
        )

        this.inputTypes = {
            // "int":"input_number",
            // "long":"input_number",
            "string":"input_text",
            "File":"input_text"
        }
        
        this.handleChange = this.handleChange.bind(this);
    }

    handleChange(event){
        const paramValue = event.currentTarget.value
        this.props.onChange(this.props.name, this.props.index ? (this.props.index) : (0), paramValue)
    }

    render(){
        const isItemNull = (this.props.paramValue == "itemNull" && this.props.itemNullAllowed)
        const disableInput = isItemNull || this.props.isNull 

        const inputType = this.inputTypes.hasOwnProperty(this.props.type) ? (
            this.inputTypes[this.props.type]
        ) : (
            "input_text"
        )
        
        let input_field
        switch(inputType){
            case "input_number":
                return(
                    <input
                        className="param-input"
                        type="number"
                        name={"input_" + this.key}
                        value={this.props.paramValue}
                        onChange={this.handleChange}
                        required={true}
                        disabled={disableInput}
                    />
                )
                break;
            case "input_text":
                return(
                    <input
                        className="param-input"
                        type="text"
                        name={"input_" + this.key}
                        value={this.props.paramValue}
                        onChange={this.handleChange}
                        required={true}
                        disabled={disableInput}
                    />
                )
                break;
        }
    }
}


class ParamNullCheckbox extends React.Component{
    constructor(props){
        super(props);
        // props.name
        // props.isNull
        // props.nullValue
        // props.refersTo
        // props.indexOrRunId
        // props.mode
        // props.toggleNull
        // props.size
        // props.disabled
        this.handleToggleNull = this.handleToggleNull.bind(this)
    }

    handleToggleNull(event){
        this.props.toggleNull(
            this.props.mode,
            this.props.name, 
            !event.target.checked,
            this.props.nullValue,
            this.props.refersTo,
            typeof this.props.indexOrRunId === "undefined" ? (null) : (this.props.indexOrRunId)
        )
    }

    render(){
        return(
            <input
                className="w3-check"
                style={ {
                    height: this.props.size ? (this.props.size) : ("20px"), 
                    width: this.props.size ? (this.props.size) : ("20px"), 
                    verticalAlign: "top",
                    display: "inline-block"
                } }
                type="checkbox"
                name={"null_checkbox_" + this.props.name}
                value={"null_checkbox_" + this.props.name}
                checked={!this.props.isNull}
                disabled={this.props.disabled ? (true) : (false)}
                onChange={this.handleToggleNull}
            />
        )
    }
}

class ParamAddOrRemove extends React.Component{
    constructor(props){
        super(props);
        // props.mode
        // props.name
        // props.runId only for run array form
        // props.handleAddOrRemove
        // props.disabled
        // props.disabledRemove
        this.handleAddOrRemove = this.handleAddOrRemove.bind(this);
    }

    handleAddOrRemove(event){
        this.props.handleAddOrRemove(
            event.currentTarget.value == "add", 
            this.props.mode, 
            this.props.name,
            this.props.runId ? (this.props.runId) : (null)
        )
    }

    render(){
        return(
            <div>
                <button 
                    className="w3-button"
                    style={ {width: "50%"} }
                    name={"add_" + this.props.name}
                    value="add"
                    onClick={this.handleAddOrRemove}
                    disabled={this.props.disabled}
                >
                    <i className="fas fa-plus"/> &nbsp;
                    add
                </button>
                <button 
                    className="w3-button"
                    style={ {width: "50%"} }
                    name={"remove_" + this.props.name}
                    value="remove"
                    onClick={this.handleAddOrRemove}
                    disabled={this.props.disabled || this.props.disabledRemove}
                >
                    <i className="fas fa-minus"/>
                    &nbsp;remove
                </button>
            </div>
        )
    }
}

class ParamForm extends React.Component{
    constructor(props){
        super(props);
        // props.paramValues
        // props.paramsConfigs
        // props.paramsHelperValues only for run array
        // props.runIds
        // props.changeParamValue
        // props.toggleNull
        
        this.state = {
            whichRunIdFocus: null
        }

        this.columnWidth = "250px"
        this.rowHeight = "25px"
        this.headerHeight = "35px"

        this.checkIfNull = this.checkIfNull.bind(this);
        this.fieldBackgroundColorClass = this.fieldBackgroundColorClass.bind(this);
        this.changeRunIdFocus = this.changeRunIdFocus.bind(this);
    }

    checkIfNull(indexByRunId){
        let isNull = {}
        Object.keys(this.props.paramValues).map( (p) => {
            let paramValues = indexByRunId ? (
                this.props.paramValues[p].filter( (v, i) => indexByRunId[p].indexOf(i) != -1)
                ) : (
                this.props.paramValues[p]
            )
            isNull[p] = paramValues == "null" && this.props.paramConfigs[p].null_allowed
        })
        return isNull
    }

    fieldBackgroundColorClass(isNull){
        return(
            isNull ? ("param-field-isnull") : ("param-field-notnull")
        )
    }

    changeRunIdFocus(newFocus){
        this.setState({whichRunIdFocus: newFocus})
    }
}



class ParamFormGlobalSingle extends ParamForm{
    render(){
        const isNull = this.checkIfNull()

        return(
            <div style={ {overflow:"auto"} }>
                <h4>Gobally-defined (Non-list) Parameters:</h4>
                <table style={ {borderSpacing: "0px 8px"} }><tbody>
                    {Object.keys(this.props.paramValues).map( (p) => (
                            <tr 
                                key={p} 
                                className={this.fieldBackgroundColorClass(isNull[p])} 
                                style={ {height: this.headerHeight} }
                            >
                                <td style={ {padding: "8px", width: "auto"} }>
                                    <ParamName name={p} />
                                </td>
                                <td style={ {padding: "8px", width: "auto"} }>
                                    {this.props.paramConfigs[p].null_allowed &&
                                        <ParamNullCheckbox
                                            name={p}
                                            isNull={isNull[p]}
                                            nullValue="null"
                                            refersTo="all"
                                            mode="global_single"
                                            toggleNull={this.props.toggleNull}
                                            size="20px"
                                        />
                                    }
                                </td>
                                <td style={ {padding: "8px", width: "100%", minWidth: this.columnWidth } }>
                                    
                                    <ParamField
                                        name={p}
                                        type={this.props.paramConfigs[p].type}
                                        paramValue={this.props.paramValues[p][0]}
                                        onChange={this.props.changeParamValue}
                                        isNull={isNull[p]}
                                    />
                                </td>
                            </tr>
                        ))
                    }
                </tbody></table>
            </div>

        )
    }
}

class ParamFormGlobalArray extends ParamForm{
    render(){
        const isNull = this.checkIfNull()

        return(
                <div style={ {overflow:"auto"} }>
                    <h4>Gobally-defined List Parameters:</h4>
                    <table style={ {borderSpacing: "8px 0px"} }><tbody>
                        <tr>
                            {Object.keys(this.props.paramValues).map( (p) => (
                                    <td 
                                        key={p} 
                                        className={this.fieldBackgroundColorClass(isNull[p]) + " w3-cell-top"}
                                        style={ {padding: "8px", minWidth: this.columnWidth} }
                                    >
                                        <table><tbody>
                                            <tr style={ {height: this.headerHeight} }>
                                                <td>#</td>
                                                <td style={ {minWidth: this.columnWidth} }>
                                                    <ParamName
                                                        name={p}
                                                    />
                                                    {this.props.paramConfigs[p].null_allowed &&
                                                        <span className="w3-right">
                                                            <ParamNullCheckbox
                                                                name={p}
                                                                isNull={isNull[p]}
                                                                nullValue="null"
                                                                refersTo="all"
                                                                mode="global_array"
                                                                toggleNull={this.props.toggleNull}
                                                                size="20px"
                                                            />
                                                        </span>
                                                    }
                                                </td>
                                            </tr>
                                            {[...Array(this.props.paramValues[p].length).keys()].map( (index) => (
                                                <tr style={ {height: this.rowHeight} } key={index}>
                                                    <td>
                                                        {index+1}
                                                    </td>
                                                    <td style={ {minWidth: this.columnWidth} }>
                                                        {this.props.paramConfigs[p].null_items_allowed &&
                                                            <ParamNullCheckbox
                                                                name={p}
                                                                isNull={this.props.paramValues[p][index]=="itemNull"}
                                                                refersTo="item"
                                                                nullValue="itemNull"
                                                                indexOrRunId={index}
                                                                mode="global_array"
                                                                toggleNull={this.props.toggleNull}
                                                                size="10px"
                                                            />
                                                        }
                                                        <ParamField
                                                            name={p}
                                                            type={this.props.paramConfigs[p].type}
                                                            paramValue={this.props.paramValues[p][index]}
                                                            onChange={this.props.changeParamValue}
                                                            isNull={isNull[p]}
                                                            itemNullAllowed={this.props.paramConfigs[p].null_items_allowed}
                                                            index={index}
                                                        />
                                                    </td>
                                                </tr>
                                            ))}
                                        </tbody></table>
                                        <ParamAddOrRemove
                                            name={p}
                                            mode="global_array"
                                            handleAddOrRemove={this.props.addOrRemoveItem}
                                            disabled={isNull[p]}
                                            disabledRemove={this.props.paramValues[p].length <= 1}
                                        />
                                    </td>
                                ))
                            }
                        </tr>
                    </tbody></table>
                </div>
        )
    }
}

class ParamFormRunSingle extends ParamForm{
    render(){
        return(
                <div style={ {overflow:"auto"} }>
                    <h4>Run-specific (Non-list) Parameters:</h4>
                    <table style={ {borderSpacing: "8px 0px"} }><tbody>
                        <tr>
                            {Object.keys(this.props.paramValues).map( (p) => (
                                <td 
                                    key={p} 
                                    className={this.fieldBackgroundColorClass(false) + " w3-cell-top"}
                                    style={ {padding: "8px", minWidth: this.columnWidth} }
                                >
                                    <table><tbody>
                                        <tr>
                                            <td>Run ID</td>
                                            <td style={ {minWidth: this.columnWidth} }>
                                                <ParamName
                                                    name={p}
                                                />
                                            </td>
                                        </tr>
                                        {[...Array(this.props.runIds.length).keys()].map( (index) => (
                                            <tr key={index}>
                                                <td>
                                                    {this.props.runIds[index]}
                                                </td>
                                                <td style={ {minWidth: this.columnWidth} }>
                                                    {this.props.paramConfigs[p].null_allowed &&
                                                        <ParamNullCheckbox
                                                            name={p}
                                                            isNull={this.props.paramValues[p][index]=="null"}
                                                            refersTo="item"
                                                            nullValue="null"
                                                            indexOrRunId={index}
                                                            mode="run_single"
                                                            toggleNull={this.props.toggleNull}
                                                            size="10px"
                                                        />
                                                    }
                                                    <ParamField
                                                        name={p}
                                                        type={this.props.paramConfigs[p].type}
                                                        paramValue={this.props.paramValues[p][index]}
                                                        onChange={this.props.changeParamValue}
                                                        isNull={this.props.paramValues[p][index]=="null"}
                                                        itemNullAllowed={this.props.paramConfigs[p].null_allowed}
                                                        index={index}
                                                    />
                                                </td>
                                            </tr>
                                        ))}
                                    </tbody></table>
                                </td>
                            ))}
                        </tr>
                    </tbody></table>
                </div>
        )
    }
}

class ParamFormRunArray extends ParamForm{
    render(){
        const whichRunIdFocus = this.state.whichRunIdFocus ? (
                this.state.whichRunIdFocus
            ) : (
                this.props.runIds[0]
            )

        let indexByRunId = {}
        Object.keys(this.props.paramValues).forEach( (p) => {
            indexByRunId[p] = []
            let runIdHelper = this.props.paramHelperValues[this.props.paramConfigs[p].split_into_runs_by[0]]
            this.props.paramValues[p].forEach( (value, index) => {
                if(runIdHelper[index] == whichRunIdFocus){
                    indexByRunId[p].push(index)
                }
            })
        })

        const isNull = this.checkIfNull(indexByRunId)

        return(
                <div style={ {overflow:"auto"} }>
                    <h4>Run-specific (Non-list) Parameters:</h4>
                    <TabPanel
                        title="Run IDs:"
                        tabs={this.props.runIds}
                        whichFocus={whichRunIdFocus}
                        changeFocus={this.changeRunIdFocus}
                    >
                        <table style={ {borderSpacing: "8px 0px"} }><tbody>
                            <tr>
                                {Object.keys(this.props.paramValues).map( (p) => (
                                    <td 
                                        key={p} 
                                        className={this.fieldBackgroundColorClass(isNull[p]) + " w3-cell-top"}
                                        style={ {padding: "8px", minWidth: this.columnWidth} }
                                    >
                                        <table><tbody>
                                            <tr>
                                                <td>#</td>
                                                <td style={ {minWidth: this.columnWidth} }>
                                                    <ParamName
                                                        name={p}
                                                    />
                                                    {this.props.paramConfigs[p].null_allowed &&
                                                        <span className="w3-right">
                                                            <ParamNullCheckbox
                                                                name={p}
                                                                isNull={isNull[p]}
                                                                nullValue="null"
                                                                refersTo="runId"
                                                                mode="run_array"
                                                                indexOrRunId={whichRunIdFocus}
                                                                toggleNull={this.props.toggleNull}
                                                                size="20px"
                                                            />
                                                        </span>
                                                    }
                                                </td>
                                            </tr>
                                            {[...Array(indexByRunId[p].length).keys()].map( (index) => (
                                                <tr key={index}>
                                                    <td>
                                                        {index}
                                                    </td>
                                                    <td style={ {minWidth: this.columnWidth} }>
                                                        {this.props.paramConfigs[p].null_items_allowed &&
                                                            <ParamNullCheckbox
                                                                name={p}
                                                                isNull={this.props.paramValues[p][
                                                                    indexByRunId[p][index]
                                                                ]=="itemNull"}
                                                                refersTo="item"
                                                                nullValue="itemNull"
                                                                indexOrRunId={indexByRunId[p][index]}
                                                                mode="run_array"
                                                                toggleNull={this.props.toggleNull}
                                                                size="10px"
                                                                disabled={isNull[p]}
                                                            />
                                                        }
                                                        <ParamField
                                                            name={p}
                                                            type={this.props.paramConfigs[p].type}
                                                            paramValue={this.props.paramValues[p][
                                                                indexByRunId[p][index]
                                                            ]}
                                                            onChange={this.props.changeParamValue}
                                                            isNull={isNull[p]}
                                                            itemNullAllowed={this.props.paramConfigs[p].null_items_allowed}
                                                        />
                                                    </td>
                                                </tr>
                                            ))}
                                        </tbody></table>
                                        <ParamAddOrRemove
                                            name={p}
                                            mode="run_array"
                                            runId={whichRunIdFocus}
                                            handleAddOrRemove={this.props.addOrRemoveItem}
                                            disabled={isNull[p]}
                                            disabledRemove={indexByRunId[p].length <= 1}
                                        />
                                    </td>
                                ))}
                            </tr>
                        </tbody></table>
                    </TabPanel>
                </div>
        )
    }
}
class JobParamFormHTML extends React.Component {
    constructor(props){
        super(props);
        // props.cwlTarget
        // props.param_modes
        // props.run_mode
        // props.run_names
        // props.jobId
        // props.changeSearchDir
        // props.changeIncludeSubDirsForSearching
        // props.changePathValAndSearch
        // props.validatePaths
        // props.searchPaths
        // props.searchDir
        // props.includeSubbDirsForSearching

        this.state = {
            actionStatus: "loading",
            form_passed_validation: false,
            serverMessages: [],
            validationMessages: [],
            modeExists: {},
            paramConfigs: {},
            paramValuesByMode: {},
            paramHelperValues: {}
        }

        this.getParamValues = this.getParamValues.bind(this);
        this.changeParamValue = this.changeParamValue.bind(this);
        this.toggleNull = this.toggleNull.bind(this);
        this.ajaxRequest = ajaxRequest.bind(this);
        this.addOrRemoveItem = this.addOrRemoveItem.bind(this);
        this.validateParamValues = this.validateParamValues.bind(this);
    }

    componentDidMount(){
        // setup timer to automatically update
        this.getParamValues()
    }

    getParamValues(){
            this.ajaxRequest({
                statusVar: "actionStatus",
                statusValueDuringRequest: "loading",
                messageVar: "serverMessages",
                sendData: {
                    cwl_target: this.props.cwlTarget,
                    param_modes: this.props.param_modes,
                    run_mode: this.props.run_mode, 
                    run_names: this.props.run_names.filter((r) => r != "")
                },
                route: routeGetParamValues,
                onSuccess: (data, messages) => {
                    const paramNames = Object.keys(data.configs).filter((p) => data.configs[p].type != "helper")
                    let paramValuesByMode = {
                        global_single: {},
                        global_array: {},
                        run_single: {},
                        run_array: {}
                    }
                    paramNames.forEach((p) =>{
                        let run_or_global = this.props.run_mode ? (
                                this.props.param_modes[p] ? ("run") : ("global")
                            ) : (
                                "global"
                            )
                        let single_or_array = data.configs[p].is_array ? ("array") : ("single")
                        let mode = run_or_global + "_" + single_or_array
                        paramValuesByMode[mode][p] = data.param_values[p]
                    })
                    let modeExists = {}
                    Object.keys(paramValuesByMode).forEach( (mode) => modeExists[mode] = Object.keys(paramValuesByMode[mode]).length > 0)
                    const paramHelperNames = Object.keys(data.configs).filter((p) => data.configs[p].type == "helper")
                    let paramHelperValues = {}
                    paramHelperNames.forEach( (p) => paramHelperValues[p] = data.param_values[p] )
                    return({
                        paramConfigs: data.configs,
                        paramValuesByMode: paramValuesByMode,
                        paramHelperValues: paramHelperValues,
                        modeExists: modeExists
                    })
                }
            })

    }

    validateParamValues(){
        let paramValues = {}
        Object.keys(this.state.paramValuesByMode).forEach((mode) => {
            Object.assign(paramValues, this.state.paramValuesByMode[mode])
        })
        Object.assign(paramValues, this.state.paramHelperValues)
        
        this.ajaxRequest({
            statusVar: "status",
            statusValueDuringRequest: "validating",
            messageVar: "validationMessages",
            sendData: {
                param_values: paramValues,
                param_configs: this.state.paramConfigs,
                cwl_target: this.props.cwlTarget,
                job_id: this.props.jobId,
                validate_paths: this.props.validatePaths,
                search_paths: this.props.searchPaths,
                search_dir: this.props.searchDir,
                include_subdirs_for_searching: this.props.includeSubbDirsForSearching
            },
            route: routeSendFilledParamValues,
            onSuccess: (data, messages) => {
                return({form_passed_validation: true})
            },
            onError: (messages) => {
                return({form_passed_validation: false})
            },
        })     
    }
    changeParamValue(mode, name, index, newValue){
        let paramValuesByMode = this.state.paramValuesByMode
        paramValuesByMode[mode][name][index] = newValue
        this.setState({paramValuesByMode: paramValuesByMode})
    }

    dissectParamValuesByRunId(mode, name, runID){
        const runIdParamName = this.state.paramConfigs[name].split_into_runs_by[0]
        const runIdIndexes = this.state.paramHelperValues[runIdParamName]
        let dissectedParamValues = {
            before: {
                indexes: [],
                values: []
            },
            match: {
                indexes: [],
                values: []
            },
            after: {
                indexes: [],
                values: []
            },
        }
        const paramValues = this.state.paramValuesByMode[mode][name]
        let before=true
        for (let r=0; r < runIdIndexes.length; r++){
            if(runIdIndexes[r] == runID){
                before=false
                dissectedParamValues.match.indexes.push(runIdIndexes[r])
                dissectedParamValues.match.values.push(paramValues[r])
            }
            else if(before){
                dissectedParamValues.before.indexes.push(runIdIndexes[r])
                dissectedParamValues.before.values.push(paramValues[r])
            }
            else {
                dissectedParamValues.after.indexes.push(runIdIndexes[r])
                dissectedParamValues.after.values.push(paramValues[r])
            }
        }

        return(dissectedParamValues)
    }

    toggleNull(mode, name, setNull, nullValue, refersTo, indexOrRunId){
        let paramValuesByMode = this.state.paramValuesByMode
        let paramHelperValues = this.state.paramHelperValues
        if (refersTo == "all"){
            paramValuesByMode[mode][name] = setNull ? ([nullValue]) : (this.state.paramConfigs[name].default_value)
        }
        else if (refersTo == "item"){
            paramValuesByMode[mode][name][indexOrRunId] = setNull ? (nullValue) : (this.state.paramConfigs[name].default_value[0])
        } 
        else if (refersTo == "runId"){
            const dissectParamValues = this.dissectParamValuesByRunId(mode, name, indexOrRunId)
            const runIdParamName = this.state.paramConfigs[name].split_into_runs_by[0]
            let newParamValues = setNull ? ([nullValue]) : (this.state.paramConfigs[name].default_value)
            paramValuesByMode[mode][name] = dissectParamValues.before.values
                .concat(newParamValues).concat(dissectParamValues.after.values)
            paramHelperValues[runIdParamName] = dissectParamValues.before.indexes.
                concat([indexOrRunId]).concat(dissectParamValues.after.indexes)
        }

        this.setState({
            paramValuesByMode: paramValuesByMode,
            paramHelperValues: paramHelperValues
        })
    }

    addOrRemoveItem(add, mode, name, runId){
        let paramValuesByMode = this.state.paramValuesByMode
        if (runId){
            let paramHelperValues = this.state.paramHelperValues
            const runIdParamName = this.state.paramConfigs[name].split_into_runs_by[0]
            const dissectParamValues = this.dissectParamValuesByRunId(mode, name, runId)
            let newParamValues
            let newParamHelperValues
            if(add){
                newParamValues = dissectParamValues.match.values
                    .concat([this.state.paramConfigs[name].default_value[0]])
                newParamHelperValues = dissectParamValues.match.indexes
                    .concat([runId])
            } 
            else {
                newParamValues = dissectParamValues.match.values.slice(0,-1)
                newParamHelperValues = dissectParamValues.match.indexes.slice(0,-1)
            }
            paramValuesByMode[mode][name] = dissectParamValues.before.values
                .concat(newParamValues).concat(dissectParamValues.after.values)
            paramHelperValues[runIdParamName] = dissectParamValues.before.indexes.
                concat(newParamHelperValues).concat(dissectParamValues.after.indexes)
            this.setState({
                paramValuesByMode: paramValuesByMode,
                paramHelperValues: paramHelperValues
            })
        }
        else {
            paramValuesByMode[mode][name] = add ? (
                paramValuesByMode[mode][name].concat([this.state.paramConfigs[name].default_value[0]])
            ) : (
                paramValuesByMode[mode][name].slice(0,-1)
            )
            this.setState({
                paramValuesByMode: paramValuesByMode
            })
        }
    }

    render() {
        if (this.state.actionStatus == "loading"){
            return(
                <LoadingIndicator
                    size="large"
                    message="Loading parameters."
                />
            )
        } else {          
            return(
                <div>
                    <DisplayServerMessages messages={this.state.serverMessages} />
                    <PathValAndSearch
                        validatePaths={this.props.validatePaths}
                        searchPaths={this.props.searchPaths}
                        changePathValAndSearch={this.props.changePathValAndSearch}
                        searchDir={this.props.searchDir}
                        changeSearchDir={this.props.changeSearchDir}
                        includeSubbDirsForSearching={this.props.includeSubbDirsForSearching}
                        changeIncludeSubDirsForSearching={this.props.changeIncludeSubDirsForSearching}
                    />
                    <hr/>
                    <h3>Provide Parameters:</h3>
                    {this.state.modeExists["global_single"] && 
                        <ParamFormGlobalSingle
                            paramValues={this.state.paramValuesByMode["global_single"]}
                            paramConfigs={this.state.paramConfigs}
                            changeParamValue={(name, index, newValue) => this.changeParamValue("global_single", name, index, newValue)}
                            toggleNull={this.toggleNull}
                        />
                    }
                    {this.state.modeExists["global_array"] && 
                        <ParamFormGlobalArray
                            paramValues={this.state.paramValuesByMode["global_array"]}
                            paramConfigs={this.state.paramConfigs}
                            changeParamValue={(name, index, newValue) => this.changeParamValue("global_array", name, index, newValue)}
                            toggleNull={this.toggleNull}
                            addOrRemoveItem={this.addOrRemoveItem}
                        />
                    }
                    {this.state.modeExists["run_single"] && 
                        <ParamFormRunSingle
                            paramValues={this.state.paramValuesByMode["run_single"]}
                            paramConfigs={this.state.paramConfigs}
                            runIds={this.props.run_names}
                            changeParamValue={(name, index, newValue) => this.changeParamValue("run_single", name, index, newValue)}
                            toggleNull={this.toggleNull}
                        />
                    }
                    {this.state.modeExists["run_array"] && 
                        <ParamFormRunArray
                            paramValues={this.state.paramValuesByMode["run_array"]}
                            paramConfigs={this.state.paramConfigs}
                            paramHelperValues={this.state.paramHelperValues}
                            runIds={this.props.run_names}
                            changeParamValue={(name, index, newValue) => this.changeParamValue("run_array", name, index, newValue)}
                            toggleNull={this.toggleNull}
                            addOrRemoveItem={this.addOrRemoveItem}
                        />
                    }
    
                    <h3>Validate Selection and Create Job</h3>

                    <ActionButton
                        name="validate"
                        value="validate"
                        onAction={this.validateParamValues}
                        label="validate"
                        loading={this.state.actionStatus == "validating"}
                        disabled={this.state.actionStatus != "none"}
                    />
                    <DisplayServerMessages messages={this.state.validationMessages} />

                    <CreateJobButton
                        jobId={this.props.jobId}
                        sheetFormat="xlsx"
                        validatePaths={this.props.validatePaths}
                        searchPaths={this.props.searchPaths}
                        searchDir={this.props.searchDir}
                        includeSubbDirsForSearching={this.props.includeSubbDirsForSearching}
                        disabled={!this.state.form_passed_validation}
                    />

                    <h4>Configs:</h4>
                    <p>
                        {JSON.stringify(this.state.paramConfigs)}
                    </p>
                    <h4>Param Values By Mode:</h4>
                    <p>
                        {JSON.stringify(this.state.paramValuesByMode)}
                    </p>
                    <h4>Param Helper Values:</h4>
                    <p>
                        {JSON.stringify(this.state.paramHelperValues)}
                    </p>
                </div>
            )
        }

    }
}


class JobParamFormSpreadsheet extends React.Component {
    constructor(props){
        super(props);
        // props.cwlTarget
        // props.param_modes
        // props.run_mode
        // props.run_names
        // props.jobId
        // props.changeSearchDir
        // props.changeIncludeSubDirsForSearching
        // props.changePathValAndSearch
        // props.validatePaths
        // props.searchPaths
        // props.searchDir
        // props.includeSubbDirsForSearching

        this.state = {
            sheetFormat: "xlsx",
            file_transfer_status: "none",
            form_passed_validation: false,
            sheetFormMessages: []
        }

        this.changeSheetFormat = this.changeSheetFormat.bind(this);
        this.genFormSheet = this.genFormSheet.bind(this);
        this.handleFormSheetUpload = this.handleFormSheetUpload.bind(this)
        this.ajaxRequest = ajaxRequest.bind(this)
    }

    changeSheetFormat(event){
        this.setState({"sheetFormat": event.currentTarget.value})
    }

    genFormSheet(){
        this.ajaxRequest({
            statusVar: "file_transfer_status",
            statusValueDuringRequest: "downloading",
            messageVar: "sheetFormMessages",
            sendData: {
                cwl_target: this.props.cwlTarget,
                param_modes: this.props.param_modes,
                run_mode: this.props.run_mode, 
                run_names: this.props.run_names.filter((r) => r != ""),
                job_id: this.props.jobId,
                sheetFormat: this.state.sheetFormat
            },
            route: routeGenParamFormSheet,
            onSuccess: (data, messages) => {
                window.location.href = data.get_form_sheet_href
                return({sheetFormMessages: []})
            }
        })
    }

    handleFormSheetUpload(isSuccess){
        this.setState({form_passed_validation: isSuccess})
    }

    render() {
        return(
            <div>
                <PathValAndSearch
                    validatePaths={this.props.validatePaths}
                    searchPaths={this.props.searchPaths}
                    changePathValAndSearch={this.props.changePathValAndSearch}
                    searchDir={this.props.searchDir}
                    changeSearchDir={this.props.changeSearchDir}
                    includeSubbDirsForSearching={this.props.includeSubbDirsForSearching}
                    changeIncludeSubDirsForSearching={this.props.changeIncludeSubDirsForSearching}
                />
                <hr/>
                <h3>Provide Parameters Using a Spreadsheet:</h3>
                <ol>
                    <li>
                        export/download:
                        <select className="w3-button w3-white w3-border" 
                            name="sheet_format"
                            onChange={this.changeSheetFormat}
                            value={this.state.sheetFormat}
                            >
                            <option value="xlsx">excel format (xlsx)</option>
                            <option value="xls">excel format (xls)</option>
                            <option value="ods">open office format (ods)</option>
                        </select> 
                        <ActionButton
                            name="export"
                            value="export"
                            label="export"
                            onAction={this.genFormSheet}
                            loading={this.state.file_transfer_status == "downloading"}
                            disabled={this.state.file_transfer_status != "none"}
                        />
                    </li>
                    <li>
                        open in excel or open office and fill in the form
                    </li>
                    <li>
                        <FileUploadComponent
                            requestRoute={routeSendFilledParamFormSheet}
                            instruction="import/upload"
                            buttonLabel="import & validate"
                            oneLine={true}
                            disabled={this.state.file_transfer_status != "none"}
                            meta_data={ 
                                {
                                    job_id: this.props.jobId,
                                    validate_paths: this.props.validatePaths,
                                    search_paths: this.props.searchPaths,
                                    search_dir: this.props.searchDir,
                                    include_subdirs_for_searching: this.props.includeSubbDirsForSearching
                                }
                            }
                            onUploadCompletion={this.handleFormSheetUpload}
                        />
                    </li>
                </ol>
                <DisplayServerMessages messages={this.state.sheetFormMessages} />
                <CreateJobButton
                    jobId={this.props.jobId}
                    sheetFormat={this.state.sheetFormat}
                    validatePaths={this.props.validatePaths}
                    searchPaths={this.props.searchPaths}
                    searchDir={this.props.searchDir}
                    includeSubbDirsForSearching={this.props.includeSubbDirsForSearching}
                    disabled={!this.state.form_passed_validation}
                />
            </div>
        )
    }
}


class JobCreationPrep extends React.Component {
    // Essential information are reqeusted by the user in order to
    // build a job template form. 
    // Following informations are needed:
    //  - which parameter shall be global and which run-specific
    //  - how many jobs (/samples) shall be created
    constructor(props) {
        super(props);
        // Inputs:
        // props.configData
        // props.cwlTarget
        let paramModes={} // container for param mode (is_run_specific true/false)
        this.props.configData.params.map((p) =>
            paramModes[p.param_name] = p.is_run_specific
        )
        this.state = {
            run_mode: false, 
            run_names: ["run1", "run2", "run3"],
            param_modes: paramModes,
            job_name: "new_job",
            display: "prep", // one of prep, form_ssheet, from_html
            validatePaths: true,
            searchPaths: true,
            searchDir: this.props.configData.default_search_dir,
            includeSubbDirsForSearching: true
        }

        // construct job_id:
        const date = new Date()
        let year = date.getFullYear().toString()
        let month = date.getMonth()
        month = (month < 10) ? ("0" + month.toString()) : (month.toString())
        let day = date.getDate()
        day = (day < 10) ? ("0" + day.toString()) : (day.toString())
        const dateString = year + month + day
        let randomNumber = Math.round(Math.random()*1000)
        randomNumber = (randomNumber < 100) ? ("0" + randomNumber.toString()) : (randomNumber.toString())
        this.jobIdNum = dateString + "_" + randomNumber
        

        this.changeJobName = this.changeJobName.bind(this);
        this.changeParamMode = this.changeParamMode.bind(this);
        this.changeRunMode = this.changeRunMode.bind(this);
        this.changeRunNames = this.changeRunNames.bind(this);
        this.toggleParamForm = this.toggleParamForm.bind(this);
        this.changeSearchDir = this.changeSearchDir.bind(this)
        this.changeIncludeSubDirsForSearching = this.changeIncludeSubDirsForSearching.bind(this)
        this.changePathValAndSearch = this.changePathValAndSearch.bind(this)
    }

    changeJobName(event){
        let jobNameString = event.currentTarget.value.trim().replaceAll(" ", "_")
        this.setState({"job_name": jobNameString})
    }

    changeParamMode(param_name, is_run_specific){
        let update = {}
        update[param_name] = is_run_specific
        this.setState({
            param_modes: Object.assign(this.state.param_modes, update)
        })
    }

    changeRunMode(value, new_bool){
        this.setState({"run_mode": new_bool})
    }

    changeRunNames(event){
        let runNameString = event.currentTarget.value
        runNameString = runNameString.replaceAll("\n", ",").replaceAll("\t", ",").replaceAll(";", ",")
        let runNames = runNameString.split(",")
        runNames = runNames.map((n) => n.trim().replace(" ", "_"))
        this.setState({"run_names": runNames})
    }

    toggleParamForm(value){
        this.setState({display: value})
    }

    changeSearchDir(newSearchDir){
        this.setState({searchDir:newSearchDir})
    }

    changeIncludeSubDirsForSearching(include){
        this.setState({includeSubbDirsForSearching:include})
    }

    changePathValAndSearch(searchPaths, validatePaths){
        this.setState({
            searchPaths: searchPaths,
            validatePaths: validatePaths
        })
    }

    render() {
        const paramTable = (
            <div style={ {maxHeight:"50vh", overflowY: "auto"} }>
                <table className="w3-table w3-bordered w3-border">
                    <thead className="w3-text-green">
                        <tr>
                            <th>Parameter</th>
                            <th>Type</th>
                            {this.state.run_mode ? (<th>Mode</th>) : (null)}
                        </tr>
                    </thead>
                    <tbody>
                        {this.props.configData.params.map( (p) => (
                            <tr key={p.param_name}> 
                                <td>{p.param_name}</td>
                                <td>{p.type}</td>
                                {this.state.run_mode ? (
                                        <td>
                                            global &nbsp;
                                            <BooleanSlider
                                                name="param_mode_select"
                                                value={p.param_name}
                                                onChange={this.changeParamMode}
                                                checked={this.state.param_modes[p.param_name]}
                                            />
                                            &nbsp; run-specific
                                        </td>
                                    ) : (null)
                                }
                            </tr>
                        ))}
                    </tbody>
                </table>
            </div>
        )

        const jobIdForm = (
            <div>
                <p>
                    <span className="w3-text-green">Job ID:</span>&nbsp;
                    <IneditableValueField>
                        {this.jobIdNum + "_" + this.state.job_name }
                    </IneditableValueField>
                </p>
                <p>
                    <label className="w3-text-green">Job name:</label><br/>
                    Please enter a job name (no whitespaces allowed).<br/>
                    <input type="text"
                        className="w3-input w3-border"
                        name="job_name"
                        value={this.state.job_name}
                        onChange={this.changeJobName}
                    />
                </p>
            </div>
        )

        const runNumberForm = (
            <div>
                <span className="w3-text-green">Runs per job:</span>&nbsp;
                single run &nbsp;
                <BooleanSlider
                    name="create_multi_run_job"
                    value="create_multi_run_job"
                    onChange={this.changeRunMode}
                    checked={this.run_mode}
                />
                &nbsp; multiple runs
                {this.state.run_mode ? (
                    <span>
                        <Message type="hint">Please choose parameter modes in the above table.</Message>
                        <label className="w3-text-green">Run names/IDs:</label><br/>
                        Please enter unique, comma-seperated IDs.
                        No whitespaces allowed (if inserted will be automatically converted to "_"). 
                        Hint: you may copy&paste cells from an excel file.
                        <textarea className="w3-input w3-border"
                            rows="2"
                            name="create_multi_run_job" 
                            value={this.state.run_names.join(", ").trim()}
                            onChange={this.changeRunNames}
                        />
                    </span>
                    ) : (null)
                }
            </div>
        )
        
        if (this.state.display == "prep"){
            return(
                <div>
                    <h3>Input parameters:</h3>
                    {paramTable}
                    <hr/>
                    <h3>Job ID:</h3>
                    {jobIdForm}
                    <hr/>
                    <h3>Number of runs:</h3>
                    {runNumberForm}
                    <hr/>
                    <h3>Provide Parameter Values:</h3>
                    <p>
                        To provide parameter values, the convenient HTML form can be used (recommended for most use cases). 
                        Alternatively, a spreadsheet can be used (Excel or OpenOffice) which is mainly for very large datasets (over 50 runs) or
                        when you would like to make use of the advanced parameter validation and manipulation options.
                    </p>
                    <p>
                        <span className="w3-text-green">Provide parameter using:</span>&nbsp;
                        <ActionButton
                            name="form_html"
                            value="form_html"
                            label={<span><i className="fas fa-list-alt"></i>&nbsp;HTML form</span>}
                            onAction={this.toggleParamForm}
                        />&nbsp; or &nbsp;
                        <ActionButton
                            name="form_ssheet"
                            value="form_ssheet"
                            label={<span><i className="fas fa-file-excel"></i>&nbsp;Spreadsheet</span>}
                            onAction={this.toggleParamForm}
                        />
                    </p>
                </div>
            )
        }
        else {
            return(
                <div>
                    <ActionButton
                        name="prep"
                        value="prep"
                        label={<span><i className="fas fa-caret-left"></i>&nbsp;back to job overview</span>}
                        onAction={this.toggleParamForm}
                    />
                    {this.state.display == "form_ssheet" ? (
                            <div>
                                <JobParamFormSpreadsheet
                                    cwlTarget={this.props.cwlTarget}
                                    param_modes={this.state.param_modes}
                                    run_mode={this.state.run_mode}
                                    run_names={this.state.run_names}
                                    jobId={this.jobIdNum + "_" + this.state.job_name}
                                    changeSearchDir={this.changeSearchDir}
                                    changeIncludeSubDirsForSearching={this.changeIncludeSubDirsForSearching}
                                    changePathValAndSearch={this.changePathValAndSearch}
                                    validatePaths={this.state.validatePaths}
                                    searchPaths={this.state.searchPaths}
                                    searchDir={this.state.searchDir}
                                    includeSubbDirsForSearching={this.state.includeSubbDirsForSearching}
                                />
                            </div>
                        ) : (
                            <div>
                                <JobParamFormHTML
                                    cwlTarget={this.props.cwlTarget}
                                    param_modes={this.state.param_modes}
                                    run_mode={this.state.run_mode}
                                    run_names={this.state.run_names}
                                    jobId={this.jobIdNum + "_" + this.state.job_name}
                                    changeSearchDir={this.changeSearchDir}
                                    changeIncludeSubDirsForSearching={this.changeIncludeSubDirsForSearching}
                                    changePathValAndSearch={this.changePathValAndSearch}
                                    validatePaths={this.state.validatePaths}
                                    searchPaths={this.state.searchPaths}
                                    searchDir={this.state.searchDir}
                                    includeSubbDirsForSearching={this.state.includeSubbDirsForSearching}
                                />
                            </div>
                        )
                    }
                </div>
            )
        }
    }
}


class JobTemplConfigInfoAjax extends React.Component {
    // Request information on the job template configuration
    constructor(props) {
        super(props);
        // Inputs:
        // props.cwlTarget
        this.buildContentOnSuccess = this.buildContentOnSuccess.bind(this);
    }

    buildContentOnSuccess(data){ // when AJAX request succeeds
        return (<JobCreationPrep configData={data} cwlTarget={this.props.cwlTarget} />);
    }

    render() {
        return (
            <AjaxComponent
                key={this.props.cwlTarget}
                requestRoute={routeGetJobTemplConfigInfo}
                sendData={ {cwl_target: this.props.cwlTarget} }
                buildContentOnSuccess={this.buildContentOnSuccess}
                loaderSize="large"
                loaderMessage="Loading template infos."
            />
        );
    }
}

class JobTemplList extends React.Component {
    constructor(props) {
        super(props);
        this.state = {whichFocus: ""}; // no list item is focues by default
        this.changeFocus = this.changeFocus.bind(this);
    }

    changeFocus(newFocusValue){
        this.setState({whichFocus: newFocusValue});
    }

    render() {
        const itemValues = this.props.templFilesInfo.map( (tf) => tf.cwl_target);
        const itemNames = itemValues; // here no distinction between values and names neccessary
        let itemContent = (
            <div>
                <DisplayServerMessages messages={this.props.initMessages}/> 
                <p>
                    <i className="fas fa-arrow-left"></i>
                    Select a CWL document to make a new job.
                </p>
            </div>
        )
        if(this.state.whichFocus != "") {
            itemContent = <JobTemplConfigInfoAjax cwlTarget={this.state.whichFocus} />
        }

        return (
            <SideBarPanel
                label="CWL documents:"
                itemValues={itemValues}
                itemNames={itemNames}
                whichFocus={this.state.whichFocus}
                itemContent={itemContent}
                onChange={this.changeFocus}
            />
        );
    }
}

class CreateJobRoot extends React.Component {
    constructor(props) {
        super(props);
        this.buildContentOnSuccess = this.buildContentOnSuccess.bind(this);
    }

    buildContentOnSuccess(data, messages){ // when AJAX request succeeds
        return (
            <div style={ {height:"100%"} }>
                { data.length > 0 ? (
                        <JobTemplList templFilesInfo={data} initMessages={messages}/>
                    ) : (
                        <Message type="info">
                            No job templates found. Please import a CWL document.
                        </Message>
                    )
                }
            </div>
        );
    }

    render() {
        return (
            <AjaxComponent
                requestRoute={routeGetJobTemplList}
                sendData={ {} }
                buildContentOnSuccess={this.buildContentOnSuccess}
                loaderSize="large"
                loaderMessage="Loading available job templates"
            />
        );
    }
}
