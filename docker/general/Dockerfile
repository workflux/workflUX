# Build context should be set to the root of
# the workflux github repo in the development branch:
# https://github.com/CompEpigen/workflUX/tree/development

FROM kerstenbreuer/cwltool:3.0.20201203173111

LABEL author="Kersten Breuer, Cancer Epigenomics, Plass Team, DKFZ 2020" \
    maintainer="kersten-breuer@outlook.com"
    

# Install development utils:
RUN python3 -m pip install pytest

# Install workflux:
COPY . /src/workflux
RUN python3 -m pip install -e /src/workflux


# Default command when starting:
CMD ["workflux", "up", "-c", "/src/workflux/docker/general/config.yaml"]
