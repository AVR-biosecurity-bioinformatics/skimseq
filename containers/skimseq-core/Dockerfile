# syntax=docker/dockerfile:1

ARG UBUNTU_VERSION=24.04
ARG RUST_VERSION=1.94.03

###############################################################################
# Stage 1: Build the software
###############################################################################

FROM ubuntu:${UBUNTU_VERSION} AS builder

ARG HTSLIB_VERSION=1.24
ARG SAMTOOLS_VERSION=1.24
ARG BCFTOOLS_VERSION=1.24
ARG BEDTOOLS_VERSION=2.31.1
ARG BEDOPS_VERSION=2.4.42
ARG SEQKIT_VERSION=2.13.0
ARG BWA_MEM2_VERSION=2.2.1
ARG RUST_VERSION
ARG RIKER_VERSION=0.4.1
ARG SEQTK_VERSION=1.4
ARG PIGZ_VERSION=2.8
ARG DUPBLASTER_VERSION=0.1.0
ARG MINIBWA_VERSION=0.5

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
        make \
        gcc \
        g++ \
        binutils \
        python3 \
        pkg-config \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libcurl4-openssl-dev \
        libssl-dev \
        libncurses-dev \
        libgomp1 && \
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
# BWA-MEM2
###############################################################################

RUN case "${TARGETARCH}" in \
        amd64) BWA_ARCH="x64-linux" ;; \
        *) \
            echo "BWA-MEM2 ${BWA_MEM2_VERSION} binary installation supports amd64 only" >&2; \
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
        --output bwa-mem2.tar.bz2 \
        "https://github.com/bwa-mem2/bwa-mem2/releases/download/v${BWA_MEM2_VERSION}/bwa-mem2-${BWA_MEM2_VERSION}_${BWA_ARCH}.tar.bz2" && \
    tar -xjf bwa-mem2.tar.bz2 && \
    install \
        -m 0755 \
        "bwa-mem2-${BWA_MEM2_VERSION}_${BWA_ARCH}/bwa-mem2" \
        /usr/local/bin/bwa-mem2 && \
    find "bwa-mem2-${BWA_MEM2_VERSION}_${BWA_ARCH}" \
        -maxdepth 1 \
        -type f \
        -name 'bwa-mem2.*' \
        -exec install -m 0755 {} /usr/local/bin/ \; && \
    rm -rf \
        bwa-mem2.tar.bz2 \
        "bwa-mem2-${BWA_MEM2_VERSION}_${BWA_ARCH}"



###############################################################################
# minibwa
###############################################################################

RUN curl \
        --fail \
        --show-error \
        --location \
        --proto '=https' \
        --proto-redir '=https' \
        --retry 20 \
        --retry-all-errors \
        --output minibwa.tar.gz \
        "https://github.com/lh3/minibwa/archive/refs/tags/v${MINIBWA_VERSION}.tar.gz" && \
    tar -xzf minibwa.tar.gz && \
    make \
        -C "minibwa-${MINIBWA_VERSION}" \
        -j"${MAKE_JOBS}" && \
    install \
        -m 0755 \
        "minibwa-${MINIBWA_VERSION}/minibwa" \
        /usr/local/bin/minibwa && \
    rm -rf \
        minibwa.tar.gz \
        "minibwa-${MINIBWA_VERSION}"
        
###############################################################################
# seqtk
###############################################################################

RUN curl \
        --fail \
        --show-error \
        --location \
        --proto '=https' \
        --proto-redir '=https' \
        --retry 20 \
        --retry-all-errors \
        --output seqtk.tar.gz \
        "https://github.com/lh3/seqtk/archive/refs/tags/v${SEQTK_VERSION}.tar.gz" && \
    tar -xzf seqtk.tar.gz && \
    make \
        -C "seqtk-${SEQTK_VERSION}" \
        -j"${MAKE_JOBS}" && \
    install \
        -m 0755 \
        "seqtk-${SEQTK_VERSION}/seqtk" \
        /usr/local/bin/seqtk && \
    rm -rf \
        seqtk.tar.gz \
        "seqtk-${SEQTK_VERSION}"


###############################################################################
# pigz
###############################################################################

RUN curl \
        --fail \
        --show-error \
        --location \
        --proto '=https' \
        --proto-redir '=https' \
        --retry 20 \
        --retry-all-errors \
        --output pigz.tar.gz \
        "https://zlib.net/pigz/pigz-${PIGZ_VERSION}.tar.gz" && \
    tar -xzf pigz.tar.gz && \
    make \
        -C "pigz-${PIGZ_VERSION}" \
        -j"${MAKE_JOBS}" && \
    install \
        -m 0755 \
        "pigz-${PIGZ_VERSION}/pigz" \
        /usr/local/bin/pigz && \
    ln -s pigz /usr/local/bin/unpigz && \
    rm -rf \
        pigz.tar.gz \
        "pigz-${PIGZ_VERSION}"

###############################################################################
# BEDOPS
###############################################################################

RUN case "${TARGETARCH}" in \
        amd64) ;; \
        *) \
            echo "ERROR: BEDOPS binary installation supports amd64 only" >&2; \
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
        --output bedops.tar.bz2 \
        "https://github.com/bedops/bedops/releases/download/v${BEDOPS_VERSION}/bedops_linux_x86_64-v${BEDOPS_VERSION}.tar.bz2" && \
    echo \
        "3a8b4961ff80e201e2a286b91c7a0ee8fff6387ab44c3fc07c3765fdf9aedfef  bedops.tar.bz2" \
        | sha256sum --check - && \
    tar -xf bedops.tar.bz2 && \
    test -d bin && \
    find bin \
        -maxdepth 1 \
        -type f \
        -exec install -m 0755 {} /usr/local/bin/ \; && \
    rm -rf \
        bedops.tar.bz2 \
        bin

###############################################################################
# Rust toolchain - Need latest version for riker and duplblaster
###############################################################################

ENV RUSTUP_HOME="/opt/rustup"
ENV CARGO_HOME="/opt/cargo"
ENV PATH="${CARGO_HOME}/bin:${PATH}"

RUN curl \
        --fail \
        --show-error \
        --location \
        --proto '=https' \
        --proto-redir '=https' \
        --retry 20 \
        --retry-all-errors \
        --output rustup-init \
        "https://static.rust-lang.org/rustup/dist/x86_64-unknown-linux-gnu/rustup-init" && \
    chmod 0755 rustup-init && \
    ./rustup-init \
        -y \
        --no-modify-path \
        --profile minimal \
        --default-toolchain "${RUST_VERSION}" && \
    rm rustup-init && \
    rustc --version && \
    cargo --version

    
###############################################################################
# Riker
###############################################################################

RUN cargo install \
        --locked \
        --version "${RIKER_VERSION}" \
        --root /usr/local \
        riker-ngs && \
    strip /usr/local/bin/riker

###############################################################################
# Dupblaster
###############################################################################

RUN cargo install \
        --locked \
        --version "${DUPBLASTER_VERSION}" \
        --root /usr/local \
        dupblaster && \
    strip /usr/local/bin/dupblaster
    
###############################################################################
# Validate the builder output
###############################################################################

ENV BCFTOOLS_PLUGINS="/usr/local/libexec/bcftools"

# Print versions and exercise basic executable startup.
RUN set -eu; \
    ldconfig; \
    \
    htsfile --version; \
    samtools --version; \
    bcftools --version; \
    bedtools --version; \
    bedops --version; \
    seqkit version; \
    bwa-mem2 version; \
    minibwa --version; \
    riker --version; \
    dupblaster --version; \
    pigz --version; \
    seqtk 2>&1 | grep -m 1 '^Version:'; 

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