FROM cgr.dev/chainguard/wolfi-base:latest
WORKDIR /app
COPY . /app

RUN apk update && \
    apk add --no-cache \
        bash \
        curl \
        bzip2 \
        ca-certificates \
        openssl \
        libgcc \
        libstdc++ \
        git \
        make \
        zlib-dev \
        libxml2-dev \
        openssl-dev

        # Install Miniconda
ENV CONDA_DIR=/opt/conda
RUN curl -L https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -o miniconda.sh && \
    bash miniconda.sh -b -p $CONDA_DIR && \
    rm miniconda.sh

# Add Conda to PATH
ENV PATH="$CONDA_DIR/bin:$PATH"

# Copy your environment.yml file
COPY coloc_pipe_env_v2.yml /tmp/environment.yml

# Create the Conda environment
RUN conda env create -f /tmp/environment.yml

# Set the default command to activate the environment
SHELL ["conda", "run", "-n", "coloc_pipe_env_v2", "/bin/bash", "-c"]
CMD ["/bin/bash"]
