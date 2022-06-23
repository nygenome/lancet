FROM ubuntu:20.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y git g++ cmake liblzma-dev zlib1g-dev libbz2-dev libcurl3-dev libssl-dev

RUN git clone --recursive https://github.com/nygenome/lancet
RUN cd lancet && make all -j$(nproc) && make lancet -j$(nproc)
RUN mkdir /lancet/bin
RUN cp /lancet/lancet /usr/local/bin
RUN PATH=$PATH:/lancet/bin
