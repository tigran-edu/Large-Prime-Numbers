FROM gcc:latest

ARG BOOST_VERSION=1.80.0

WORKDIR /gtest_build

RUN apt update --fix-missing

RUN apt-get update && \
    apt-get install -y \
    libboost-dev libboost-program-options-dev \
    libgtest-dev \
    libssl-dev \
    git \
    curl \
    cmake \
    wget \
    && \
    cmake -DCMAKE_BUILD_TYPE=Release /usr/src/gtest && \
    cmake --build . && \
    mv lib/*.a /usr/lib

RUN cd /tmp && \
    BOOST_VERSION_MOD=$(echo 1.84.0 | tr . _) && \
    wget https://boostorg.jfrog.io/artifactory/main/release/1.84.0/source/boost_${BOOST_VERSION_MOD}.tar.bz2 && \
    tar --bzip2 -xf boost_${BOOST_VERSION_MOD}.tar.bz2 && \
    cd boost_${BOOST_VERSION_MOD} && \
    ./bootstrap.sh --prefix=/usr/local && \
    ./b2 install && \
    rm -rf /tmp/*

WORKDIR /project

COPY ../ ./

CMD cmake . && \
    cmake --build ./ --config RelWithDebInfo && \
    make && \
    CTEST_OUTPUT_ON_FAILURE=TRUE ctest tests --verbose
