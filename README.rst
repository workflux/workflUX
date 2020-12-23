CWLab - An open-source framework for simplified deployment of the Common Workflow Language using a graphical web interface
==========================================================================================================================

Background and Scope:
---------------------

The Common Workflow Language (CWL) allows to wrap and link up
bioinformatic software in a standardized and portable way. However,
setting up and operating a CWL-based workflow management system can be a
labor-intensive challenge for many data-driven laboratories. To this
end, we developed CWLab: a framework for simplified, graphical
deployment of CWL.

CWLab allows life-science researchers with all levels of computational
proficiency to create, execute and monitor jobs for CWL-wrapped tools
and workflows. Input parameters for large sample batches are specified
using a simple HTML form or a spreadsheet and are automatically
validated. The integrated webserver allows to remotely control the
execution on clusters as well as single workstations. Moreover,
automatic infrastructure provisioning and scaling for OpenStack-based
clouds is being implemented. CWLab can also be used as a local desktop
application that supports Linux, MacOS, and Windows by leveraging Docker
containerization. Our Python-based framework is easy to set up and, via
a flexible API, it can be integrated with any CWL runner and adapted to
custom software environments.

With CWLab, we would like to hide the complexity of workflow management
so that scientific users can focus on their data analyses. This might
promote the adoption of CWL in multi-professional life-science
laboratories.

Installation and Quick Start:
-----------------------------

**Attention: CWLab is in alpha state currently and not all features are
available yet. However, the core functionalities are working and we are
happy if you test it. We are working hard to push out a stable version
in the coming weeks. Please press the watch button to not miss it.**

| Installation can be done using pip:
| ``python3 -m pip install cwlab``

Please see the section "Configuration" for a discussion of available
options.

| Start the webserver with your custom configuration (or leave out the
  ``--config`` flag to use the default one):
| ``cwlab up --config config.yaml``

If you like to make use of containers for dependency management, you
need to install `Docker <https://docs.docker.com/install/>`__ or a
Docker-compatible containerization solution like
`singularity <https://singularity.lbl.gov/>`__ or
`udocker <https://github.com/indigo-dc/udocker>`__. To run on Windows or
MacOs, please install the dedicated docker versions: `Docker for
Windows <https://docs.docker.com/docker-for-windows/>`__, `Docker for
Mac <https://docs.docker.com/docker-for-mac/>`__

The usage of the web interface should be self-explanatory with build-in
instruction. The following section gives an overview of the basic usage
scenario.

Supported Systems:
------------------

CWLab is written in platform-agnostic python and can therefore be
executed on:

-  **Linux**
-  **MacOs**
-  **Windows**\ \*

Any CWL runner that has a command-line interface can be integrated into
CWLab in order to execute CWL workflows or tool-wrappers, such as:

-  **cwltool** (the reference implementation) -
   https://github.com/common-workflow-language/cwltool
-  **Toil** (UCSC) - https://github.com/DataBiosphere/toil
-  **Cromwell** (Broad Institute) -
   https://github.com/broadinstitute/cromwell
-  **Reana** (CERN) - https://reana.readthedocs.io/en/latest/index.html
-  **CWLEXEC** (IBM) - https://github.com/IBMSpectrumComputing/cwlexec
   (Please find a constantly updated list at:
   https://www.commonwl.org/#Implementations)

Therefore, CWLab can be used on any infrastructure supported by these
CWL runners, including:

-  **single workstations**
-  **HPC clusters** (PBS, LSF, slurm, ...)
-  **clouds** (AWS, GCP, Azure, OpenStack)

| \*\ **Please Note:**
| Execution on Windows is only supported by cwltool which talks to
  docker for windows. Therefore, CWL-wrapped tools and workflows which
  where originally designed for Linux/MacOs can be executed on Windows
  with a graphical interface provided by CWLab.

Usage:
------

Connect to the web interface:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Open a modern browser of your choice like Chrome, Firefox, Safari, or
Edge (Internet Explorer might be partially incompatible).

Type in the URL of your webserver. The URL depends on your
configuration:

-  If the webserver is running on the same machine and uses port 5000 is
   used (this is the default), type: ``https://localhost:5000/``
-  If CWLab is running on a remote machine in the same network, type in
   the machine's IP address and the used port. For instance, if the IP
   address is 172.22.0.1 and port 5000 is
   used:\ ``https://172.22.0.1:5000/``

| You should see a Welcome page like this:
| |welcome screenshot|

Import a CWL workflow or tool:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

| CWLab can be used to run any workflow or tool that has been wrapped
  using the Common Workflow Language. Of course, you can write workflows
  or tool wrappers yourself (we recommend rabix-composer
  https://rabix.io/), however, for many especially bioinformatic tasks,
  existing CWL solutions are publicly available. Check the CWL website
  as a starting point:
| https://www.commonwl.org/#Repositories\_of\_CWL\_Tools\_and\_Workflows.

| To import a CWL document:
| - Click on the button "Import CWL Workflow/Tool" in the top bar -
  Choose a CWL document (workflow or tool)\* - Press the import button

| The workflow will be automatically validated:
| |import screenshot|

| *\*Please note: Currently, workflows can only be imported in the
  "packed" format. We will add support for the unpacked format soon. To
  pack a CWL workflow, use:*
| ``cwltool --pack my_workflow.cwl > my_workflow_packed.cwl``

Create a new Job:
~~~~~~~~~~~~~~~~~

To run a workflow or tool with your data, you have to create a new job.
One job may contain multiple runs (for instance multiple samples or
conditions). CWLab will automatically present you a list of needed input
parameters. For each parameter, you can choose whether to specify it
globally (all runs of a job will get the same value) or per run.

-  Click on the button "Create New Job" in the top bar and select the
   desired CWL document in the sidebar
-  Specify a descriptive job name (the job ID will be composed of the
   date, time, and the name)
-  If the job shall contain multiple runs toggle the "runs per job"
   switch, then:

   -  Specify run names as a comma-separated list in the dedicated text
      field
   -  In the parameter list, select which parameters should be
      run-specific

-  CWLab will automatically create a parameter form for you to fill in:

   -  Export/download the form in the desired format
   -  Open it in a spreadsheet editor (e.g. Microsoft Excel or Open
      Office)
   -  The file may contain the following sheets (depends on the type of
      input parameters and your selections for "global"/"run-specific"
      specification):

      -  ``global single values``: parameters that take only one value
         and are defined globally (one for all runs)
      -  ``run-specific single values``: parameters that take only one
         value but are specified per run
      -  ``global arrays``: array parameters (takes a list of values)
         that are defined globally
      -  A separate sheet will be created for each run-specific array
         parameter. It will be titled with the parameters name
      -  ``config``: This sheet contains configuration options that only
         need adaption in advanced use cases.

   -  Fill in the sheet and import/upload the edited file to CWLab
      **\***

-  Your parameter settings are automatically validated. (E.g. it is
   checked whether the specified values match the parameter's type and
   whether the paths of specified files or directories exist.)
-  If valid, you can press the "create job" button and head over to "Job
   Execution & Results" in the top bar

**\* Please note:** For specifying file or directory parameters, there
are two options: - Either specify the absolute path - Specify a
character string that can be uniquely matched to a file/directory in the
default input directory (please see the **INPUT\_DIR** parameter in the
config section).

| This is an example screenshot for creating a job for an ATAC-seq
  workflow:
| |create job screenshot|

Job execution:
~~~~~~~~~~~~~~

-  Click on "Job Execution & Results" in the top bar and choose the job
   of interest in the sidebar
-  Select the runs you want to start
-  Select an execution profile (see the "Configuration" for details) and
   press "start"
-  The execution status will be displayed in the run-list
-  Pressing the "DetailsResults" button will show (not implemented yet):

   -  the deployed input parameter
   -  execution logs (from the CWL runner)
   -  a QC report

-  Once finished the output can be found in the "exec" directory (set in
   the configuration) along with the used parameter values, CWL
   document, and log files

| An example screenshot of the execution interface:
| |execution screenshot|

Configuration:
--------------

CWLab is a highly versatile package and makes almost no assumptions on
your hard- and software environment used for the execution of CWL. To
adapt it to your system and use case, a set of configuration options is
available:

-  General configs, including:

   -  webserver (hosting IP address and port, remotely or locally
      available, login protected or not)
   -  paths of working directories

-  Execution profiles:
   This flexible API allows you to adapt CWLab to your local software
   environment and to integrate a CWL runner of your choice (such as
   Cwltool, Toil, or Cromwell).

| All configuration options can be specified in a single YAML file which
  is provided to CWLab upon start:
| ``cwlab up --config my_config.yaml``

| To get an example config file, run the following command:
| ``cwlab print_config > config.yaml`` (or see the example below)

General Configs:
~~~~~~~~~~~~~~~~

-  **WEB\_SERVER\_HOST**:
   Specify the host or IP address on which the webserver shall run. Use
   ``localhost`` for local usage on your machine only. Use ``0.0.0.0``
   to allow remote accessibility by other machines in the same network.
   *Default*: ``localhost``
-  | **WEB\_SERVER\_PORT**:
   | Specify the port used by the webserver.
   | *Default*: 5000

-  **TEMP\_DIR**:
   Directory for temporary files.
   *Default*: a subfolder "cwlab/temp" in the home directory
-  **WORKFLOW\_DIR**:
   Directory for saving CWL documents.
   *Default*: a subfolder "cwlab/temp" in the home directory
-  **EXEC\_DIR**:
   Directory for saving execution data including output files.
   *Default*: a subfolder "cwlab/temp" in the home directory
-  **DEFAULT\_INPUT\_DIR**:
   Default directory where users can search for input files. You may
   specify additional input directories using the "**ADD\_INPUT\_DIRS**"
   parameter. *Default*: a subfolder "cwlab/temp" in the home directory
-  **DB\_DIR**:
   Directory for databases.
   *Default*: a subfolder "cwlab/temp" in the home directory
-  | **ADD\_INPUT\_DIRS**:
   | In addition to "**DEFAULT\_INPUT\_DIR**", these directories can be
     searched by the user for input files.
   | Please specify them in the format "*name: path*" like shown in this
     example:

   ::

       ADD_INPUT_DIRS:
           GENOMES_DIR: '/ngs_share/genomes'
           PUBLIC_GEO_DATA: '/datasets/public/geo'

   *Default*: no additional input dirs.
-  | **ADD\_INPUT\_AND\_UPLOAD\_DIRS**:
   | Users can search these directories for input files (in addition to
     "**DEFAULT\_INPUT\_DIR**") and they may also upload their one
     files.
   | Please specify them in the format "*name: path*" like shown in this
     example:

   ::

       ADD_INPUT_AND_UPLOAD_DIRS:
           UPLOAD_SCRATCH: '/scratch/upload'
           PERMANEN_UPLOAD_STORE: '/datasets/upload'

   *Default*: no additional input dirs.

-  | **DEBUG**:
   | If set to True, the debugging mode is turned on. Do not use on
     production systems.
   | *Default*: False

Exec Profiles:
~~~~~~~~~~~~~~

This is where you configure how to execute cwl jobs on your system. A
profile consists of four steps: pre\_exec, exec, eval, and post\_exec
(only exec required, the rest is optional). For each step, you can
specify commands that are executed in bash or cmd terminal.

You can define multiple execution profile as shown in the config example
below. This allows frontend users to choose between different execution
options (e.g. using different CWL runners, different dependency
management systems, or even choose a between multiple available batch
execution infrastructures like lsf, pbs, ...). For each execution
profile, following configuration parameters are available (but only
**type** and **exec** is required):

-  **shell**:
   Specify which shell to use. For Linux or MacOS use ``bash``. For
   Windows, use ``powershell``.
   *Required*.
-  **max\_retries**: Specify how many times the execution (all steps) is
   retried before marking a run as failed.
-  | **timeout**:
   | For each step in the execution profile, you can set a timeout
     limit.
   | *Default*:

   .. code:: yaml

       prepare: 120
       exec: 86400
       eval: 120
       finalize: 120

-  **pre\_exec**\ \*:
   Shell commands that are executed before the actual CWL execution. For
   instance to load required python/conda environments.
   *Optional*.
-  **exec**\ \*:
   Shell commands to start the CWL execution. Usually, this is only the
   command line to execute the CWL runner. The stdout and stderr of the
   CWL runner should be redirected to the predefined log file.
   *Required*.
-  **eval**\ \*:
   The exit status at the end of the *exec* step is automatically
   checked. Here you can specify shell commands to additionally evaluate
   the content of the execution log to determine if the execution
   succeeded. To communicate failure to CWLab, set the ``SUCCESS``
   variable to ``False``.
   *Optional*.
-  **post\_exec**\ \*: Shell commands that are executed after *exec* and
   *eval*. For instance, this can be used to clean up temporary files.

\* **Additional notes regarding execution profile steps:**

-  In each step following predefined variables are available:

   -  ``JOB_ID``
   -  ``RUN_ID`` (please note: is only unique within a job)
   -  ``WORKFLOW`` (the path to the used CWL document)
   -  ``RUN_INPUT`` (the path to the YAML file containing input
      parameters)
   -  ``OUTPUT_DIR`` (the path of the run-specific output directory)
   -  ``LOG_FILE`` (the path of the log file that should receive the
      stdout and stderr of CWL runner)
   -  ``SUCCESS`` (if set to ``False`` the run will be marked as failed
      and terminated)
   -  ``PYTHON_PATH`` (the path to the python interpreter used to run
      CWLab)

-  The four steps will be executed in the same shell session and
   therefore can be treated as one connected script. (Between the steps,
   CWLab communicates the status to the database allowing the user to
   get status notifications via the front end).
-  Thus you may define your own variables that will also be available in
   all downstream steps.
-  At the end of each step. The exit code is checked. If it is non-zero,
   the run will be marked as failed. Please note, if a step consists of
   multiple commands and an intermediate command fails, this will not be
   recognized by CWLab as long as the final command of the step will
   succeed. To manually communicate failure to CWLab, please set the
   ``SUCCESS`` variable to ``False``.
-  The steps are executed using pexpect
   (https://pexpect.readthedocs.io/en/stable/overview.html), this allows
   you also connect to a remote infrastructure via ssh (recommended to
   use an ssh key). Please be aware that the path of files or
   directories specified in the input parameter YAML will not be adapted
   to the new host. We are working on solutions to achieve an automated
   path correction and/or upload functionality if the execution host is
   not the CWLab server host.
-  On Windows, please be aware that each code block (contained in
   ``{...}``) has to be in one line.

Example configuration files:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Below, you can find example configurations for local execution of CWL
workflows or tools with cwltool.

Linux / MacOs:
^^^^^^^^^^^^^^

.. code:: yaml

    WEB_SERVER_HOST: localhost 
    WEB_SERVER_PORT: 5000

    DEBUG: False  

    TEMP_DIR: '/home/cwlab_user/cwlab/temp'
    WORKFLOW_DIR: '/home/cwlab_user/cwlab/workflows'
    EXEC_DIR: '/datasets/processing_out/'
    DEFAULT_INPUT_DIR: '/home/cwlab_user/cwlab/input'
    DB_DIR: '/home/cwlab_user/cwlab/db'

    ADD_INPUT_DIRS:
        GENOMES_DIR: '/ngs_share/genomes'
        PUBLIC_GEO_DATA: '/datasets/public/geo'

    ADD_INPUT_AND_UPLOAD_DIRS:
        UPLOAD_SCRATCH: '/scratch/upload'
        PERMANEN_UPLOAD_STORE: '/datasets/upload'

    EXEC_PROFILES:
        cwltool_local:
            shell: bash
            max_retries: 2
            timeout:
                prepare: 120
                exec: 86400
                eval: 120
                finalize: 120
            exec: |
                cwltool --outdir "${OUTPUT_DIR}" "${WORKFLOW}" "${RUN_INPUT}" \
                    >> "${LOG_FILE}" 2>&1
            eval: | 
                LAST_LINE=$(tail -n 1 ${LOG_FILE})
                if [[ "${LAST_LINE}" == *"Final process status is success"* ]]
                then
                    SUCCESS=True
                else
                    SUCCESS=False
                    ERR_MESSAGE="cwltool failed - ${LAST_LINE}"
                fi

Windows:
^^^^^^^^

.. code:: yaml

    WEB_SERVER_HOST: localhost
    WEB_SERVER_PORT: 5000

    DEBUG: False  

    TEMP_DIR: 'C:\Users\cwlab_user\cwlab\temp'
    WORKFLOW_DIR: 'C:\Users\cwlab_user\cwlab\workflows'
    EXEC_DIR: 'D:\processing_out\'
    DEFAULT_INPUT_DIR: 'C:\Users\cwlab_user\cwlab\input'
    DB_DIR: 'C:\Users\cwlab_user\cwlab\db'

    ADD_INPUT_DIRS:
        GENOMES_DIR: 'E:\genomes'
        PUBLIC_GEO_DATA: 'D:\public\geo'
        
    ADD_INPUT_AND_UPLOAD_DIRS:
        UPLOAD_SCRATCH: 'E:\upload'
        PERMANEN_UPLOAD_STORE: '\D:\upload'

    EXEC_PROFILES:
        cwltool_windows:
            shell: powershell
            max_retries: 2
            timeout:
                prepare: 120
                exec: 86400
                eval: 120
                finalize: 120
            exec: |
                . "${PYTHON_PATH}" -m cwltool --debug --default-container ubuntu:16.04 --outdir "${OUTPUT_DIR}" "${CWL}" "${RUN_INPUT}" > "${LOG_FILE}" 2>&1

            eval: |
                $LAST_LINES = (Get-Content -Tail 2 "${LOG_FILE}")

                if ($LAST_LINES.Contains("Final process status is success")){$SUCCESS="True"}
                else {$SUCCESS="False"; $ERR_MESSAGE = "cwltool failed - ${LAST_LINE}"}

Documentation:
--------------

**Please note: A much more detailed documentation is on the way. In the
meantime, please notify us if you have any questions (see the "Contact
and Contribution" section). We are happy to help.**

Contact and Contribution:
-------------------------

If you have any question or are experiencing problems with CWLab, please
contact us at ``k.breuer@dkfz.de`` or open an issue in Github.

If you would like to contribute to the development and like to extend
the functionality of CWLab to meet your requirements, you are more than
welcome. We will do our best to support you and your contribution will
be acknowledged.

About Us:
---------

CWLab is developed with love at the Division of Cancer Epigenomics of
the German Cancer Research Center (DKFZ) in the beautiful university
city of Heidelberg. We are an interdisciplinary team with wet-lab
scientists and bioinformaticians working closely together. Our DNA
sequencing-driven methodologies produce challenging amounts of data.
CWLab helps us by giving all members of our team the ability to perform
common bioinformatic analyses autonomously without having to acquire
programming skills. This allows our bioinformatic staff to focus on
method development and interpretation of computationally complex data.

If you like to know more about us, please visit our website
https://www.dkfz.de/en/CanEpi/contact.html.

Licence:
--------

This package is free to use and modify under the Apache 2.0 Licence.

.. |welcome screenshot| image:: https://github.com/CompEpigen/CWLab/blob/master/screenshots/welcome.png?raw=true
.. |import screenshot| image:: https://github.com/CompEpigen/CWLab/blob/master/screenshots/import.png?raw=true
.. |create job screenshot| image:: https://github.com/CompEpigen/CWLab/blob/master/screenshots/create_job.png?raw=true
.. |execution screenshot| image:: https://github.com/CompEpigen/CWLab/blob/master/screenshots/execution.png?raw=true

