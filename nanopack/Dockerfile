# Start from official Ubuntu 22.04 image
FROM ubuntu:22.04

LABEL maintainer="BiKC"

# Install the necessary packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    cmake \
    make \
    g++ \
    python3 \
    python3-pip \
    python3-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*


# Install most recent version of NanoComp
RUN pip3 install \
    NanoComp==1.24.0 \
    NanoPlot==1.44.1 

# Stel de werkdirectory in
WORKDIR /data

