# Dockerfile
FROM mambaorg/micromamba:latest

# Let us install into the base env via micromamba
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# Work as the default micromamba user
USER $MAMBA_USER

# Create working directory
WORKDIR /work

# Copy env file
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml

# Install everything into the "base" environment
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean -a -y

# Make sure conda/mamba env binaries are on PATH in shells
SHELL ["bash", "-lc"]
ENV PATH=/opt/conda/bin:$PATH
ENV CONDA_DEFAULT_ENV=base

# Default command: just drop into bash
CMD ["bash"]
