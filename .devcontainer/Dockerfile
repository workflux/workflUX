FROM python:3.8

# install essential dependencies
RUN apt-get update -qq -y
RUN apt-get install -y \
    build-essential \
    nodejs \
    npm \
    less \
    nano

# copy install script to bin
COPY ./dev_install /bin