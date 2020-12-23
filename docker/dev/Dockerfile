# Build context should be set to the root of
# the cwlab github repo in the development branch:
# https://github.com/CompEpigen/CWLab/tree/development

FROM kerstenbreuer/cwltool:3.0.20201203173111

LABEL author="Kersten Breuer, Cancer Epigenomics, Plass Team, DKFZ 2020" \
    maintainer="kersten-breuer@outlook.com"

# Install cwlab:
COPY . /src/cwlab
RUN python3 -m pip install -e /src/cwlab


# Default command when starting:
CMD ["cwlab", "up", "-c", "/src/cwlab/docker/dev/config.yaml"]
