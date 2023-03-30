FROM ubuntu:20.04

ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update -y && apt-get install -y make rsync git gcc g++ bzip2 hdf5-tools unzip gfortran curl software-properties-common
WORKDIR /test

RUN mkdir -p /opt/julia-1.6.2 && \
    curl -s -L https://julialang-s3.julialang.org/bin/linux/x64/1.6/julia-1.6.2-linux-x86_64.tar.gz | tar -C /opt/julia-1.6.2 -x -z --strip-components=1 -f -

RUN echo "\nPATH=/opt/julia-1.6.2/bin:\$PATH\n" >> /root/.bashrc
RUN ln -s /opt/julia-1.6.2/bin/julia /usr/local/bin/