FROM cgr.dev/chainguard/wolfi-base:latest
WORKDIR /app

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
        openssl-dev \
        coreutils



        # Install Miniconda
ENV CONDA_DIR=/opt/conda
RUN curl -L https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -o miniconda.sh && \
    bash miniconda.sh -b -p $CONDA_DIR && \
    rm miniconda.sh

# Add Conda to PATH
ENV PATH="$CONDA_DIR/bin:$PATH"

# Copy your environment.yml file
COPY pipeline_environment.yml /tmp/environment.yml

# Create the Conda environment
#RUN conda env create -f /tmp/environment.yml

# Set the default command to activate the environment
RUN conda env update -n pipeline_environment -f /tmp/environment.yml
# Set the environment name
ENV CONDA_ENV=pipeline_environment
ENV PATH=/opt/conda/envs/$CONDA_ENV/bin:$PATH

# Optional: add auto-activation when shell is started
ENTRYPOINT ["conda", "run", "-n", "pipeline_environment"]
