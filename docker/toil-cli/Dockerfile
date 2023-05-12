# Build context should be set to the root of
# the workflux github repo in the development branch:
# https://github.com/CompEpigen/workflUX/tree/development

FROM quay.io/ucsc_cgl/toil:latest

# Install development utils:
RUN python3 -m pip install pytest

# Install workflux:
COPY . /src/workflux
RUN python3 -m pip install -e /src/workflux

# Default command when starting:
CMD ["workflux", "up", "-c", "/src/workflux/docker/toil-cli/config.yaml"]