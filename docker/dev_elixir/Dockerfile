FROM ubuntu:18.04

LABEL author="Kersten Breuer, Cancer Epigenomics, Plass Team, DKFZ 2020" \
    maintainer="kersten-breuer@outlook.com"

# Package versions/links:
ENV CWLTOOL_VERSION 1.0.20181012180214

# Install essential dependencies:
ENV DEBCONF_FRONTEND noninteractive
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
	
RUN python3 -m pip install \
		pymysql \
		trs-cli==0.3.0

# Install cwlab:
COPY . /src/cwlab
RUN python3 -m pip install -e /src/cwlab

# Configure cwlab:
COPY ./docker/dev_elixir/workflows /cwlab/wf_dir
COPY ./docker/dev_elixir/config.yaml /cwlab/config.yaml
COPY ./docker/dev_elixir/wes_exec_profile.py /cwlab/wes_exec_profile.py

# Default command when starting:
CMD ["/src/cwlab/cwlab.py", "up", "-c", "/cwlab/config.yaml"]
