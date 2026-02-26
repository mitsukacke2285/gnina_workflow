############################
# Stage 1 — Build gnina
############################
FROM nvidia/cuda:12.2.2-cudnn8-devel-ubuntu22.04 AS builder

ENV DEBIAN_FRONTEND=noninteractive
WORKDIR /build

# Base dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    git \
    cmake \
    libboost-all-dev \
    libeigen3-dev \
    libopenbabel-dev \
    python3-dev \
    python3-numpy \
    python3-scipy \
    python3-openbabel \
    python3-pytest \
    cython3 \
    libprotobuf-dev \
    protobuf-compiler \
    libjsoncpp-dev \
    ca-certificates \
    gpg \
    wget \
    && rm -rf /var/lib/apt/lists/*

# Install modern CMake (>=3.25) from Kitware
RUN wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc | \
        gpg --dearmor -o /usr/share/keyrings/kitware-archive-keyring.gpg && \
    echo "deb [signed-by=/usr/share/keyrings/kitware-archive-keyring.gpg] https://apt.kitware.com/ubuntu/ jammy main" \
        > /etc/apt/sources.list.d/kitware.list && \
    apt-get update && \
    apt-get install -y cmake && \
    rm -rf /var/lib/apt/lists/*

# Clone gnina
RUN git clone --depth 1 --branch v1.3 https://github.com/gnina/gnina.git

WORKDIR /build/gnina

# Build
RUN mkdir build && cd build && \
    cmake .. \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_POLICY_VERSION_MINIMUM=3.5 && \
    make -j$(nproc)


############################
# Stage 2 — Runtime image
############################
FROM nvidia/cuda:12.2.2-cudnn8-runtime-ubuntu22.04

ENV DEBIAN_FRONTEND=noninteractive
WORKDIR /app

# Runtime dependencies only (no compiler)
RUN apt-get update && apt-get install -y --no-install-recommends \
    libboost-all-dev \
    libopenbabel-dev \
    libprotobuf-dev \
    libjsoncpp-dev \
    python3 \
    python3-numpy \
    python3-scipy \
    python3-openbabel \
    && rm -rf /var/lib/apt/lists/*

# Copy built binary
COPY --from=builder /build/gnina/build /opt/gnina

ENV PATH="/opt/gnina/bin:${PATH}"

WORKDIR /workspace

CMD ["gnina", "--version"]
