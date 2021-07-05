FROM centos:centos8
LABEL maintainer = "BiKC"

RUN dnf -y install cmake make gcc gcc-c++ &&\
    dnf -y install python3 python3-devel

RUN pip3 install pyarrow==0.9.0 &&\
    pip3 install NanoComp 

WORKDIR /data
