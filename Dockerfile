FROM cgr.dev/chainguard/wolfi-base:latest

LABEL version="3.0"
LABEL description="Dockerfile for pipeline environment using Miniconda on Wolfi base image for Flanders"
LABEL maintainer="bruno.ariano@fht.org"

ARG py_v=3.12

RUN apk update && \
    apk add --no-cache \
        bash \
        curl \
        bzip2 \
        ca-certificates \
        openssl \
        posix-libc-utils \
        glibc \
        glibc-dev \
        libgcc \
        libstdc++ \
        git \
        make \
        zlib-dev \
        libxml2-dev \
        openssl-dev \
        coreutils \
        python-${py_v} \
        py${py_v}-pip

# install conda-lock
ENV PIP_ROOT_USER_ACTION=ignore
RUN pip install conda-lock

# Install Miniconda
ENV CONDA_DIR=/opt/conda
RUN curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh" && \
    bash Miniforge3-$(uname)-$(uname -m).sh -b -p $CONDA_DIR && \
    rm Miniforge3-$(uname)-$(uname -m).sh
ENV PATH="$CONDA_DIR/bin:$PATH"

# Create the environment with conda-lock
COPY conda-lock.yml /tmp/conda-lock.yml
RUN conda-lock install -n pipeline_environment /tmp/conda-lock.yml

# Set the environment name
ENV CONDA_ENV=pipeline_environment
ENV PATH="/opt/conda/envs/$CONDA_ENV/bin:$PATH"

# Auto-activation when shell is started
ENTRYPOINT ["conda", "run", "-n", "pipeline_environment"]
