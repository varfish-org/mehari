# syntax=docker/dockerfile:1.3

# Note: two-stage building is disabled for now because of the dual nature of
# mehari right now as both CLI and lib; we would have to use features to
# untable the two.

# We don't do cross compilation at the moment but build the dependencies first
# anyway to get the ability to increment.

# ---------------------------------------------------------------------------
# Builder
# ---------------------------------------------------------------------------

# Use ubuntu:noble as the base image
FROM ubuntu:noble AS builder

# Install Rust toolchain and dependencies for compilation of C code (e.g., rocksdb)
RUN apt-get update && \
    apt-get install -y wget curl build-essential clang librocksdb-dev && \
    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y && \
    . $HOME/.cargo/env && \
    rustup component add rustfmt

# Install build dependency `protoc`.
COPY utils/install-protoc.sh /tmp
RUN PREFIX=/usr/local bash /tmp/install-protoc.sh

# Set the working directory.
WORKDIR /usr/src/mehari

# Copy in the rest of the sources.
COPY build.rs Cargo.toml Cargo.lock /usr/src/mehari/
COPY src /usr/src/mehari/src/
COPY protos /usr/src/mehari/protos/
COPY utils/alpine-linker-script.sh /usr/src/mehari/utils/
RUN chmod a+rx /usr/src/mehari/utils/alpine-linker-script.sh
COPY .cargo /usr/src/mehari/.cargo/

# Touch main.rs to prevent cached release build.
RUN touch /usr/src/mehari/src/main.rs

# Build the application
RUN /root/.cargo/bin/cargo build --release

# ---------------------------------------------------------------------------
# Runtime
# ---------------------------------------------------------------------------

FROM ubuntu:noble AS runtime

# Install dependencies (and cleanup afterward)
RUN apt-get update && \
    apt-get install -y libsqlite3-0 && \
    apt-get clean autoclean && \
    apt-get autoremove --yes && \
    rm -rf /var/lib/{apt,dpkg,cache,log}

# Copy application binary from builder image
COPY --from=builder \
    /usr/src/mehari/target/release/mehari \
    /usr/local/bin

# Copy the entrypoint script and make it executable.
COPY utils/docker/entrypoint.sh /
RUN chmod a+rx /entrypoint.sh

# Set the entrypoint.
ENTRYPOINT ["/bin/bash", "/entrypoint.sh"]

# Expose the application port
EXPOSE 8080
