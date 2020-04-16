FROM ubuntu:18.04

LABEL author="Kersten Breuer, Cancer Epigenomics, Plass Team, DKFZ 2020" \
    maintainer="kersten-breuer@outlook.com"

# Package versions/links:
ENV CWLTOOL_VERSION 1.0.20181012180214

# Install essential dependencies:
RUN apt-get update -qq -y --fix-missing
RUN apt-get install -y build-essential \
    python3-dev \
	python3-pip \
    wget \
    nodejs \
    tar \
	curl \
	libtool \
	libz-dev \
    libssl-dev \
	libcurl4-openssl-dev

# Install singularity:
ENV DEBCONF_FRONTEND noninteractive
RUN echo "console-setup   console-setup/charmap47 select  UTF-8" > encoding.conf && \
	debconf-set-selections encoding.conf && \
	rm encoding.conf
RUN apt-get install -qq -y \
	squashfs-tools \
	uuid-dev \
	libgpgme11-dev \
	libseccomp-dev \
	pkg-config \
	git

# Install cwlab:
COPY . /src/cwltool
RUN python3 -m pip install -e /src/cwltool

# Configure cwlab:
RUN mkdir -p \
        /cwlab/tmp \
        /cwlab/worlflows \
        /cwlab/exec \
        /cwlab/input \
        /cwlab/database && \
    chmod -R 777 /cwlab

COPY docker/dev_elixir/config.yaml /cwlab/config.yaml

# Default command when starting:
CMD ["cwlab", "up", "-c", "/cwlab/config.yaml"]