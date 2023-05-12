# Build context should be set to the root of
# the workflux github repo in the development branch:
# https://github.com/CompEpigen/workflUX/tree/development

FROM quay.io/ucsc_cgl/toil:latest

# Install development utils:
RUN python3 -m pip install --upgrade pip
RUN python3 -m pip install pytest
RUN apt update && apt-get install curl -y
RUN curl -LO https://storage.googleapis.com/kubernetes-release/release/`curl -s https://storage.googleapis.com/kubernetes-release/release/stable.txt`/bin/linux/amd64/kubectl
RUN chmod +x ./kubectl
RUN mv ./kubectl /usr/local/bin/kubectl


# Install workflux:
COPY . /src/workflux
RUN python3 -m pip install -e /src/workflux

# Default command when starting:
CMD ["workflux", "up", "-c", "/src/workflux/docker/toil-wes/config.yaml"]