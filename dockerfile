# Dockerfile for Willmo3/DiffDomain
# Provides a ready-to-build environment: cmake, gcc/clang, git,
# OpenMP (libgomp/libomp), and a getopt-capable glibc userland
# (Fedora's glibc ships getopt.h + the getopt() function, so nothing
# extra is needed for that beyond the standard build toolchain).

FROM fedora:40

# --- System dependencies ---
RUN dnf install -y \
    gcc \
    gcc-c++ \
    cmake \
    make \
    git \
    ca-certificates \
    libomp-devel \
    && dnf clean all
# gcc/gcc-c++      -> compiler + libgomp (gcc's OpenMP runtime) pulled in automatically
# clang            -> alternative compiler mentioned in the README
# cmake, make      -> build system (repo needs cmake >= 3.5)
# git              -> needed by scripts/update_deps to fetch cereal/eigen
# libomp-devel     -> OpenMP headers/runtime for clang
# getopt() itself lives in glibc, which is always present as a base dependency

# --- Project setup ---
WORKDIR /workspace
COPY . /workspace

RUN cd ./scripts && chmod +x ./update-deps && ./update-deps && cd ..

RUN mkdir -p build && cd build \
    && cmake -DCMAKE_BUILD_TYPE=Release .. \
    && make

CMD ["/bin/bash"]
