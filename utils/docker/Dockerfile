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
    apt-get install -y unzip wget curl build-essential clang librocksdb-dev libsqlite3-dev && \
    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y && \
    . $HOME/.cargo/env && \
    rustup component add rustfmt

# Install build dependency `protoc`.
COPY utils/install-protoc.sh /tmp
RUN PREFIX=/usr/local bash /tmp/install-protoc.sh

# Set the working directory.
WORKDIR /usr/src/mehari

COPY Cargo.toml Cargo.lock ./
COPY mehari/Cargo.toml mehari/build.rs ./mehari/
COPY mehari/src ./mehari/src/
COPY mehari/protos ./mehari/protos/
COPY mehari-python/Cargo.toml ./mehari-python/
COPY utils/alpine-linker-script.sh ./utils/
RUN chmod a+rx ./utils/alpine-linker-script.sh
COPY .cargo ./.cargo/


# Touch main.rs to prevent cached release build.
RUN touch ./mehari/src/main.rs

# Build the application
RUN /root/.cargo/bin/cargo build --release -p mehari

# ---------------------------------------------------------------------------
# Runtime
# ---------------------------------------------------------------------------

FROM ubuntu:noble AS runtime

# Install dependencies (and cleanup afterward)
RUN apt-get update && \
    apt-get install -y librocksdb8.9 libsqlite3-0 && \
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
