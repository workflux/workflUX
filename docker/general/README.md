# Standalone Docker Image

## Concept

This docker image runs the workflux server without authentication enabled. It uses
`cwltool` and `singularity` to execute workflows. The data for workflux (input data,
workflows, database, run data, singularity images) are stored at `/workflux`. This
directory may be mounted from the host to enable persistence of the data. The
configuration is located at `/etc/workflux/config.yaml` and can be replaced on
demand.

The image defines a default non-root user and group (1000:1000). For the
execution of singularity containers it must be run in privileged mode. You may
run the container with different user and group, but therefore you need to
provide an appropriate `/etc/passwd`-file.

### Build the docker image

From the root of the repository call:

```bash
docker build -t <imagename> -f docker/dev/Dockerfile .
```

### Run the container

```bash
# Without persistence
docker run --detach \
           --privileged \
           --name workflux \
           --publish 8080:5000 \
           <imagename>

# With persistence
docker run --detach \
           --privileged \
           --name workflux \
           --publish 8080:5000 \
           --volume <directory>:/workflux \
           <imagename>

# With persistence and custom configuration file
docker run --detach \
           --privileged \
           --name workflux \
           --publish 8080:5000 \
           --volume <directory>:/workflux \
           --volume <config-file>:/etc/workflux/config.yaml:ro \
           <imagename>
```

### Stop and delete the container

```bash
docker stop workflux
docker rm workflux
```
