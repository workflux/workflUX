# Tutorial: GWAS Workflow with workflUX

TODO general introduction to the workflow
* What does it do?
* How does it work (tools etc.)?
* What are the results?

In order to execute the workflow you need access to a workflux-instance. You can
run an instance on your own computer with docker (see Step 0).

## Step 0: Start the docker container

1. Install docker on your machine. See https://docs.docker.com/get-docker/ for
   detailed instructions
2. Open a terminal on your computer and start a container with 
   `docker run -it --rm --privileged -p 8080:5000 compepigen/workflux:latest`
   This container will print its output to your terminal and will be automatically deleted when it ends.

## Step 1: Import the workflow

1. Open `Import CWL Workflow/Tool`
2. Select `URL to public CWL document (e.g. from github)`
3. Use `https://raw.githubusercontent.com/elixir-cloud-aai/demo-workflows/dev/cwl/gwas_workflow/CWL/gwas_workflow.cwl` as the workflow URL
4. Name the workflow, e.g. "GWAS demo workflow"
5. Click on "Import using selected name"

## Step 2: Prepare the workflow execution

1. Open `Create Job`
2. Choose your imported Workflow, e.g. "GWAS demo workflow"
3. Name your job, e.g. "GWAS demo Job"
4. The workflow requires two parameters, a metadata and variants file. Click on "Provide parameter using: HTML Form"
5. Enter the following public URLs for the fields:
   * metadata: `https://raw.githubusercontent.com/elixir-cloud-aai/demo-workflows/dev/cwl/gwas_workflow/example_job/data/thousand_genomes_meta.csv`
   * variants: `https://raw.githubusercontent.com/elixir-cloud-aai/demo-workflows/blob/dev/cwl/gwas_workflow/example_job/data/ALL.chr21.integrated_phase1_v3.20101123.snps_indels_svs.genotypes_subsample.vcf.gz`
7. Click on "validate and create job"
8. The job is setup to be started

## Step 3: Execute the workflow

1. Open `Job Execution & Results`
2. Choose your created Job on the left side
3. Click on "all" in the List of Runs
4. Click on the "start" button
5. The execution starts
6. Click on the eye in the Details column of the "List of Runs"-table
7. Select `Execution log` to see the workflow logs and the progress
8. Once the workflow finishes (see the status column in the List of Runs table), you can download the output files in the `Output files` section

