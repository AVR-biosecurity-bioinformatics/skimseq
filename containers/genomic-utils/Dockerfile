# syntax=docker/dockerfile:1

ARG UBUNTU_VERSION=24.04

###############################################################################
# Stage 1: Build the software
###############################################################################

FROM ubuntu:${UBUNTU_VERSION} AS builder

ARG HTSLIB_VERSION=1.24
ARG SAMTOOLS_VERSION=1.24
ARG BCFTOOLS_VERSION=1.24
ARG BEDTOOLS_VERSION=2.31.1
ARG SEQKIT_VERSION=2.13.0
ARG MAKE_JOBS=4
ARG TARGETARCH

USER root

# Required for HTTPS access through the DEECA Zscaler proxy.
#
# The certificate must be Base-64/PEM encoded and begin with:
# -----BEGIN CERTIFICATE-----
#
# The certificate is used only in this intermediate build stage.
COPY zscaler-root.cer \
    /usr/local/share/ca-certificates/zscaler-root.crt

# Ubuntu 24.04 stores its default repositories in ubuntu.sources.
#
# Switch to HTTPS and configure APT to tolerate intermittent proxy failures.
RUN sed -i \
        -e 's|http://archive.ubuntu.com|https://archive.ubuntu.com|g' \
        -e 's|http://security.ubuntu.com|https://security.ubuntu.com|g' \
        /etc/apt/sources.list.d/ubuntu.sources && \
    printf '%s\n' \
        'Acquire::Retries "20";' \
        'Acquire::http::Pipeline-Depth "0";' \
        'Acquire::https::Pipeline-Depth "0";' \
        'Acquire::Queue-Mode "access";' \
        'Acquire::https::CaInfo "/usr/local/share/ca-certificates/zscaler-root.crt";' \
        > /etc/apt/apt.conf.d/80-deeca-proxy

# Install build dependencies.
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        ca-certificates \
        curl \
        bzip2 \
        python3 \
        make \
        gcc \
        g++ \
        autoconf \
        automake \
        pkg-config \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libcurl4-openssl-dev \
        libssl-dev \
        libncurses-dev && \
    update-ca-certificates && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /opt/src

###############################################################################
# HTSlib
###############################################################################

RUN curl \
        --fail \
        --show-error \
        --location \
        --proto '=https' \
        --proto-redir '=https' \
        --retry 20 \
        --retry-all-errors \
        --output htslib.tar.bz2 \
        "https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2" && \
    tar -xjf htslib.tar.bz2 && \
    cd "htslib-${HTSLIB_VERSION}" && \
    ./configure --prefix=/usr/local && \
    make -j"${MAKE_JOBS}" && \
    make install && \
    cd .. && \
    rm -rf \
        htslib.tar.bz2 \
        "htslib-${HTSLIB_VERSION}"

###############################################################################
# SAMtools
###############################################################################

RUN curl \
        --fail \
        --show-error \
        --location \
        --proto '=https' \
        --proto-redir '=https' \
        --retry 20 \
        --retry-all-errors \
        --output samtools.tar.bz2 \
        "https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2" && \
    tar -xjf samtools.tar.bz2 && \
    cd "samtools-${SAMTOOLS_VERSION}" && \
    ./configure --prefix=/usr/local && \
    make -j"${MAKE_JOBS}" && \
    make install && \
    cd .. && \
    rm -rf \
        samtools.tar.bz2 \
        "samtools-${SAMTOOLS_VERSION}"

###############################################################################
# BCFtools
###############################################################################

RUN curl \
        --fail \
        --show-error \
        --location \
        --proto '=https' \
        --proto-redir '=https' \
        --retry 20 \
        --retry-all-errors \
        --output bcftools.tar.bz2 \
        "https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2" && \
    tar -xjf bcftools.tar.bz2 && \
    cd "bcftools-${BCFTOOLS_VERSION}" && \
    ./configure --prefix=/usr/local && \
    make -j"${MAKE_JOBS}" && \
    make install && \
    cd .. && \
    rm -rf \
        bcftools.tar.bz2 \
        "bcftools-${BCFTOOLS_VERSION}"

###############################################################################
# BEDTools
###############################################################################

RUN curl \
        --fail \
        --show-error \
        --location \
        --proto '=https' \
        --proto-redir '=https' \
        --retry 20 \
        --retry-all-errors \
        --output bedtools.tar.gz \
        "https://github.com/arq5x/bedtools2/releases/download/v${BEDTOOLS_VERSION}/bedtools-${BEDTOOLS_VERSION}.tar.gz" && \
    tar -xzf bedtools.tar.gz && \
    make -C bedtools2 -j"${MAKE_JOBS}" && \
    cp bedtools2/bin/* /usr/local/bin/ && \
    rm -rf \
        bedtools.tar.gz \
        bedtools2

###############################################################################
# SeqKit
###############################################################################

# Select the release archive matching the target architecture.
RUN case "${TARGETARCH}" in \
        amd64) SEQKIT_ARCH="amd64" ;; \
        arm64) SEQKIT_ARCH="arm64" ;; \
        *) \
            echo "Unsupported architecture: ${TARGETARCH}" >&2; \
            exit 1 \
            ;; \
    esac && \
    curl \
        --fail \
        --show-error \
        --location \
        --proto '=https' \
        --proto-redir '=https' \
        --retry 20 \
        --retry-all-errors \
        --output seqkit.tar.gz \
        "https://github.com/shenwei356/seqkit/releases/download/v${SEQKIT_VERSION}/seqkit_linux_${SEQKIT_ARCH}.tar.gz" && \
    tar -xzf seqkit.tar.gz \
        -C /usr/local/bin \
        seqkit && \
    chmod 0755 /usr/local/bin/seqkit && \
    rm seqkit.tar.gz

###############################################################################
# Validate the builder output
###############################################################################

ENV BCFTOOLS_PLUGINS="/usr/local/libexec/bcftools"

RUN ldconfig && \
    htsfile --version | head -n 1 && \
    samtools --version | head -n 1 && \
    bcftools --version | head -n 1 && \
    bcftools plugin -l | head -n 1 && \
    bedtools --version && \
    seqkit version

# Remove DEECA-specific files before /usr/local is copied to the final image.
RUN rm -f \
        /usr/local/share/ca-certificates/zscaler-root.crt \
        /etc/apt/apt.conf.d/80-deeca-proxy

###############################################################################
# Stage 2: Prepare a minimal runtime filesystem
###############################################################################

FROM ubuntu:${UBUNTU_VERSION} AS runtime-rootfs

USER root

# The certificate is needed temporarily while installing runtime packages.
COPY zscaler-root.cer \
    /usr/local/share/ca-certificates/zscaler-root.crt

RUN sed -i \
        -e 's|http://archive.ubuntu.com|https://archive.ubuntu.com|g' \
        -e 's|http://security.ubuntu.com|https://security.ubuntu.com|g' \
        /etc/apt/sources.list.d/ubuntu.sources && \
    printf '%s\n' \
        'Acquire::Retries "20";' \
        'Acquire::http::Pipeline-Depth "0";' \
        'Acquire::https::Pipeline-Depth "0";' \
        'Acquire::Queue-Mode "access";' \
        'Acquire::https::CaInfo "/usr/local/share/ca-certificates/zscaler-root.crt";' \
        > /etc/apt/apt.conf.d/80-deeca-proxy

# Install runtime dependencies only.
#
# bash and procps are included for Nextflow compatibility.
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        ca-certificates \
        bash \
        procps \
        zlib1g \
        libbz2-1.0 \
        liblzma5 \
        libcurl4t64 \
        libssl3t64 \
        libncursesw6 \
        libtinfo6 \
        libstdc++6 \
        libgcc-s1 && \
    rm -rf /var/lib/apt/lists/*

# Remove the DEECA-specific certificate and rebuild the standard Ubuntu
# certificate store.
#
# The runtime-rootfs stage is subsequently copied into a new scratch-based
# stage, so the earlier certificate layer is not inherited by the final image.
RUN rm -f \
        /usr/local/share/ca-certificates/zscaler-root.crt \
        /etc/apt/apt.conf.d/80-deeca-proxy && \
    update-ca-certificates --fresh

# Create a generic unprivileged runtime user.
#
# The username is descriptive only; permissions are controlled by the numeric
# UID and GID.
RUN groupadd \
        --system \
        --gid 10001 \
        containeruser && \
    useradd \
        --system \
        --uid 10001 \
        --gid containeruser \
        --create-home \
        --shell /bin/bash \
        containeruser && \
    mkdir -p /data && \
    chown containeruser:containeruser /data

###############################################################################
# Stage 3: Assemble the final portable image
###############################################################################

FROM scratch AS runtime

ARG HTSLIB_VERSION=1.24
ARG SAMTOOLS_VERSION=1.24
ARG BCFTOOLS_VERSION=1.24
ARG BEDTOOLS_VERSION=2.31.1
ARG SEQKIT_VERSION=2.13.0
ARG IMAGE_VERSION=1.0.0

LABEL org.opencontainers.image.title="Genomic utilities"
LABEL org.opencontainers.image.description="HTSlib, SAMtools, BCFtools, BEDTools and SeqKit"
LABEL org.opencontainers.image.version="${IMAGE_VERSION}"

# Copy the current, cleaned Ubuntu runtime filesystem.
#
# This copies the filesystem state without inheriting the earlier
# runtime-rootfs layers in which the certificate was present.
COPY --from=runtime-rootfs / /

# Copy only the installed programs and supporting files from the builder.
COPY --from=builder /usr/local/bin/ /usr/local/bin/
COPY --from=builder /usr/local/lib/ /usr/local/lib/
COPY --from=builder /usr/local/libexec/ /usr/local/libexec/
COPY --from=builder /usr/local/share/ /usr/local/share/

ENV PATH="/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"
ENV BCFTOOLS_PLUGINS="/usr/local/libexec/bcftools"

# Validate the installed software, plugin discovery, runtime libraries and
# certificate removal before changing to the non-root user.
RUN ldconfig && \
    htsfile --version | head -n 1 && \
    samtools --version | head -n 1 && \
    bcftools --version | head -n 1 && \
    bcftools plugin -l | head -n 1 && \
    bedtools --version && \
    seqkit version && \
    test ! -e /usr/local/share/ca-certificates/zscaler-root.crt && \
    ! find /etc/ssl/certs \
        -iname '*zscaler*' \
        -print -quit | grep -q . && \
    ! ldd /usr/local/bin/htsfile \
        | grep -q 'not found' && \
    ! ldd /usr/local/bin/samtools \
        | grep -q 'not found' && \
    ! ldd /usr/local/bin/bcftools \
        | grep -q 'not found' && \
    ! ldd /usr/local/bin/bedtools \
        | grep -q 'not found'

WORKDIR /data

USER containeruser

CMD ["/bin/bash"]