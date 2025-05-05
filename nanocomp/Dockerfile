# Use CentOS Stream 9 as base image
# podman pull quay.io/centos/centos:stream9
FROM quay.io/centos/centos:stream9

LABEL maintainer="BiKC"

RUN dnf -y install cmake make gcc gcc-c++ &&\
    dnf -y install python3 python3-devel

RUN pip3 install pyarrow &&\
    pip3 install Nanocomp

WORKDIR /data
