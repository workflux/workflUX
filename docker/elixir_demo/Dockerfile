FROM ubuntu:18.04

LABEL author="Kersten Breuer, Cancer Epigenomics, Plass Team, DKFZ 2020" \
    maintainer="kersten-breuer@outlook.com"

# Package versions/links:
ENV CWLTOOL_VERSION 1.0.20181012180214

# Install essential dependencies:
ENV DEBCONF_FRONTEND noninteractive
RUN apt-get update -qq -y --fix-missing
RUN apt-get install -y build-essential \
    wget \
    nodejs \
    tar \
	curl \
	libtool \
	libz-dev \
    libssl-dev \
	libcurl4-openssl-dev


RUN echo 'debconf debconf/frontend select Noninteractive' | debconf-set-selections

RUN apt-get install -y libreadline-gplv2-dev \
	libncursesw5-dev \
    libsqlite3-dev \
	tk-dev \
	libgdbm-dev \
	libc6-dev \
	libbz2-dev \
	libffi-dev \
	zlib1g-dev

# Install Python 3.7:
RUN cd /usr/src && \
	wget https://www.python.org/ftp/python/3.7.9/Python-3.7.9.tgz && \
	tar xzf Python-3.7.9.tgz && \
	cd Python-3.7.9 && \
	./configure --enable-optimizations && \
	make altinstall	


RUN apt-get install -qq -y \
	squashfs-tools \
	uuid-dev \
	libgpgme11-dev \
	libseccomp-dev \
	pkg-config \
	git \
	python3-pip
	
RUN python3.7 -m pip install \
	pymysql \
	requests \
	pydantic \
	urllib3

# Install workflux:
COPY . /src/workflux
RUN python3.7 -m pip install -e /src/workflux

# Configure workflux:
COPY ./docker/elixir_demo/config.yaml /workflux/config.yaml
COPY ./docker/elixir_demo/wes_exec_profile.py /workflux/wes_exec_profile.py

# Default command when starting:
CMD ["python3.7", "-m", "workflux", "up", "-c", "/workflux/config.yaml"]
