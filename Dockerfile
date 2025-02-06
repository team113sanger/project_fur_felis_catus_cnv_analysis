# syntax=docker/dockerfile:1
# MUTLI-STAGE DOCKERFILE for building a Python image
# STAGE 1: base_stage is the base image for both the development and dev-staging/staging/production images
# STAGE 2: development_only_stage is the image used for development (optimised for VSCode)

########################
# STAGE 1: base_stage #
########################

# IMPORTANT
# If you change the base image, you will need to update the
# PRE_FETCH_BASE_IMAGE variable in the .gitlab-ci.yml file also.
FROM ubuntu:22.04 AS base_stage


# The following ARGs set the versions of Python and R to be installed
#
#  - Changing the Python version will also change installed python
#  - Changeing the R version will only change the checking process;
#    to actually install a different version of R,
#    you will need to pick the correct debian package names, see below.
#
ARG PYTHON_VERSION=3.10 R_VERSION=4.2.2

# Set environment variables.
# 1. Force Python stdout and stderr streams to be unbuffered.
# 2. Do not build if pip can't connect to the internet
# 3. Do not cache pip installs
# 4. Set the poetry version which will handle the rest of the depenedency installation
# 5. The mode in which Django will be executed
USER root
ENV \
    PYTHONUNBUFFERED=1 \
    PIP_DEFAULT_TIMEOUT=100 \
    PIP_NO_CACHE_DIR=1 \
    PIP_VERSION=23.3.1 \
    POETRY_VERSION=1.8.2 \
    DATA_DIRECTORY="/data" \
    OPT_DIRECTORY="/opt" \
    USER_NAME="admin" \
    USER_DIRECTORY="/home/admin" \
    POETRY_VIRTUALENVS_CREATE=false \
    POETRY_VIRTUALENVS_IN_PROJECT=false \
    BIOCONDUCTOR_VERSION="3.16" \
    DNACOPY_VERSION="1.72" \
    COMPLEXHEATMAP_VERSION="2.14"

ENV \
    USER_BASHRC="${USER_DIRECTORY}/.bashrc" \
    USER_BIN_DIRECTORY="${USER_DIRECTORY}/.local/bin" \
    SSH_DIR="${USER_DIRECTORY}/.ssh" \
    VENV_DIRECTORY="${OPT_DIRECTORY}/venv" \
    POETRY_HOME="${OPT_DIRECTORY}/poetry" \
    PIPX_HOME="${OPT_DIRECTORY}/pipx" \
    PIPX_BIN_DIR="${OPT_DIRECTORY}/pipx/bin" \
    POETRY_CACHE_DIR="${OPT_DIRECTORY}/poetry-cache" \
    PROJECT_DIRECTORY="${OPT_DIRECTORY}/repo" \
    LOGGING_DIRECTORY="${DATA_DIRECTORY}/logs"

ENV PATH=${POETRY_HOME}/bin:${PIPX_BIN_DIR}:${PATH}

RUN \
    useradd "${USER_NAME}" --shell /bin/bash --create-home --home-dir "${USER_DIRECTORY}" \
    && mkdir -p "${PROJECT_DIRECTORY}" "${DATA_DIRECTORY}" "${OPT_DIRECTORY}" "${POETRY_CACHE_DIR}" "${PIPX_HOME}" "${POETRY_HOME}" "${VENV_DIRECTORY}" "${PIPX_BIN_DIR}"\
    && chown -R "${USER_NAME}:${USER_NAME}" "${PROJECT_DIRECTORY}" "${DATA_DIRECTORY}" "${USER_DIRECTORY}" "${OPT_DIRECTORY}" \
    && chmod -R 755 "${PROJECT_DIRECTORY}" "${DATA_DIRECTORY}" "${USER_DIRECTORY}" "${OPT_DIRECTORY}"


# Update System
# Install system packages required by Python packages
# We include:
# - curl + wget, for debugging or downloading files
# - nano, for debugging & a minimal editor
# - git, so that pre-commit can be installed and run
# - tree, so that the project structure can be viewed in the container
# - netcat, to do network debugging
# - openssh-client, to be able to do git operations over ssh & also scp for backup-and-restore operations
# - build-essential, meta-packages that are essential to compile software including gcc, g++, make, etc.
# - pkg-config, to manage compile and link flags for libraries
# - software-properties-common/gpg/dirmngr, to add system repositories and install R
# BuildKit logic to cache apt packages
# - https://vsupalov.com/buildkit-cache-mount-dockerfile/
# - https://github.com/moby/buildkit/blob/master/frontend/dockerfile/docs/reference.md#run---mounttypecache
# - https://stackoverflow.com/a/72851168
RUN rm -f /etc/apt/apt.conf.d/docker-clean; echo 'Binary::apt::APT::Keep-Downloaded-Packages "true";' > /etc/apt/apt.conf.d/keep-cache
RUN \
    apt-get update --quiet && \
    apt-get install --yes --quiet software-properties-common && \
    add-apt-repository --yes ppa:deadsnakes/ppa && \
    apt-get update --quiet && \
    apt-get install --yes --quiet --no-install-recommends \
    build-essential \
    pkg-config \
    "python${PYTHON_VERSION:?}" \
    "python${PYTHON_VERSION:?}-venv" \
    vim \
    nano \
    curl \
    wget \
    git \
    tree \
    ncdu \
    openssh-client \
    pipx \
    software-properties-common \
    dirmngr \
    gpg \
    && rm -rf /var/lib/apt/lists/* \
    && pipx ensurepath \
    && pipx install poetry==$POETRY_VERSION \
    && python3 --version | grep ${PYTHON_VERSION} \
    && pipx inject -f poetry poetry-plugin-export \
    && apt-get -s clean

# Install R 4.2.2 from https://cloud.r-project.org/bin/linux/ubuntu/
# The subversions of R-4.* can be found at https://cloud.r-project.org/bin/linux/ubuntu/jammy-cran40/
#
# The system packages from this r-project repository ensure that a specific
# version of R is installed however do not install r-base or r-recommended
# packages, as they will force the installation to the latest R version despite
# specifying different versions of r-base-core and r-base-dev. Yes its mad.
#
RUN \
    wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | gpg --dearmor -o /etc/apt/trusted.gpg.d/cran_ubuntu_key.gpg && \
    add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/' && \
    apt-get update -qq && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    r-base-core=4.2.2-1.2204.0 \
    r-base-dev=4.2.2-1.2204.0 \
    r-base-html=4.2.2-1.2204.0 \
    r-doc-html=4.2.2-1.2204.0 \
    r-doc-info=4.2.2-1.2204.0 \
    r-doc-pdf=4.2.2-1.2204.0 \
    r-mathlib=4.2.2-1.2204.0 \
    && rm -rf /var/lib/apt/lists/* \
    && R --version | grep ${R_VERSION:?}

# Install R packages
#
# To use an older Bioconductor version, we need to specify the version of the
# Bioconductor and the version of the R packages we want to install. Ask is set
# to FALSE to avoid user input.
RUN \
    Rscript -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager', repos='https://cloud.r-project.org')" && \
    Rscript -e "BiocManager::install(version = '${BIOCONDUCTOR_VERSION:?}', ask = FALSE)" && \
    Rscript -e "BiocManager::install('ComplexHeatmap', version = '${BIOCONDUCTOR_VERSION:?}', ask = FALSE)" && \
    Rscript -e "packageVersion('ComplexHeatmap')" | grep -q "${COMPLEXHEATMAP_VERSION}" || (echo \"Got $(Rscript -e \"packageVersion('ComplexHeatmap')\") instead of ${COMPLEXHEATMAP_VERSION}\" && exit 1) && \
    Rscript -e "BiocManager::install('DNAcopy', version = '${BIOCONDUCTOR_VERSION:?}', ask = FALSE)" && \
    Rscript -e "packageVersion('DNAcopy')" | grep -q "${DNACOPY_VERSION}" || (echo \"Got $(Rscript -e \"packageVersion('DNAcopy')\") instead of ${DNACOPY_VERSION}\" && exit 1) && \
    Rscript -e "install.packages('optparse', repos='https://cloud.r-project.org')"


# As the non-root user we install pipx and poetry so that they are available in
# the command line for subsequent package installations.
USER $USER_NAME


# VENV STEP
# There are two motivations:
# 1. To keep the image small and well configured, we want to re-use system
#    installed packages. This is especially important with tensorflow which is a
#    large package.
# 2. We want the extra packages to be installed in a virtual environment so that
#    when working in a Singularity container which is a read-only file system,
#    we can still install packages (by binding and reinstalling into venv directory)

WORKDIR $OPT_DIRECTORY
RUN \
    python3 -m venv ${VENV_DIRECTORY} &&  \
    echo "source ${VENV_DIRECTORY}/bin/activate" >> "${USER_BASHRC}"
ENV \
    PATH="${VENV_DIRECTORY}/bin:${PATH}" \
    VIRTUAL_ENV=${VENV_DIRECTORY}

# Copy dependency control files and install dependencies
# 1. only the dependency control files for Poetry (equivalent to "requirements.txt"
#    if you have never seen Poetry *.lock and *.toml files before)
#    The unsual square brakets in `poetry.loc[k]` provides GLOB/REGEX syntax which ensures
#    IDEMPOTENCY for the build (i.e. Docker builds correctly the first and every subsequent time).
#    Initially, `poetry.lock` may be absent as it is only created after the first `poetry install`
#    This syntax lets COPY conditionally include it without the build failing if missing.
# 2. .gitignore as it is used by precommit-check in a post build stage of the CICD
WORKDIR $PROJECT_DIRECTORY
COPY --chown="${USER_NAME}:${USER_NAME}" [".gitignore", "pyproject.toml", "poetry.loc[k]", "./"]
RUN \
    poetry install --no-root --no-directory  \
    && chown -R "${USER_NAME}:${USER_NAME}" ${VENV_DIRECTORY}

# Copy the src, test and other code of the project into the container.
# Then install the project itself via Poetry
COPY --chown="${USER_NAME}:${USER_NAME}" . "${PROJECT_DIRECTORY}"
RUN \
    poetry install \
    && chown -R "${USER_NAME}:${USER_NAME}" "${USER_DIRECTORY}" "${PROJECT_DIRECTORY}" "${DATA_DIRECTORY}"

# Use user "admin" to run the build commands below and the server itself.
USER "${USER_NAME}"

###################################
# STAGE 2: development_only_stage #
# - This stage is optional        #
# - It is optimised for VSCode    #
###################################

# To develop from the container we need add some extra directories to work nicely with VSCode
FROM base_stage AS development_only_stage
USER root

# Install hubflow to allow for git flow style development & conditional install
# sudo, giving the user passwordless sudo privileges
WORKDIR "${USER_DIRECTORY}"
ARG HAS_SUDO="${HAS_SUDO:-0}"
RUN git config --global --add safe.directory "${USER_DIRECTORY}/gitflow" \
    && git config --global --add safe.directory "${USER_DIRECTORY}/gitflow/shFlags" \
    && git clone https://github.com/datasift/gitflow \
    && chown -R "${USER_NAME}:${USER_NAME}" gitflow \
    && cd gitflow \
    && ./install.sh \
    && if [ "${HAS_SUDO}" = "1" ]; then \
    apt-get update -y \
    && apt-get install -y sudo \
    && rm -rf /var/lib/apt/lists/* \
    && echo "${USER_NAME:?} ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers; \
    fi

# Return back to the project directory for the entrypoint
WORKDIR "${PROJECT_DIRECTORY}"
USER "${USER_NAME}"

# Prepare the directory for the VSCode server extensions and bash history
RUN mkdir -p "${USER_DIRECTORY}/.vscode-server/extensions" \
    "${USER_DIRECTORY}/.vscode-server-insiders/extensions" \
    && chown -R "${USER_NAME}:${USER_NAME}" \
    "${USER_DIRECTORY}/.vscode-server" \
    "${USER_DIRECTORY}/.vscode-server-insiders" && \
    SNIPPET="export PROMPT_COMMAND='history -a' && export HISTFILE=${USER_DIRECTORY}/.commandhistory/.bash_history" \
    && mkdir "${USER_DIRECTORY}/.commandhistory/" \
    && touch "${USER_DIRECTORY}/.commandhistory/.bash_history" \
    && chown -R "${USER_NAME}:${USER_NAME}" "${USER_DIRECTORY}/.commandhistory/" \
    && echo "$SNIPPET" >> "/home/$USER_NAME/.bashrc"

# Setup the pre-commit hooks so that they are run before each commit
RUN pre-commit install --install-hooks
