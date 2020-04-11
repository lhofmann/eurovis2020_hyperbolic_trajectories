FROM ubuntu:18.04

RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y \
    build-essential \
    cmake \
    libboost-all-dev \
    libeigen3-dev \
    libvtk6-dev \
    python-numpy \
    && rm -rf /var/lib/apt/lists/*

# fix broken vtk package
RUN ln -s /usr/lib/python2.7/dist-packages/vtk/libvtkRenderingPythonTkWidgets.x86_64-linux-gnu.so /usr/lib/x86_64-linux-gnu/libvtkRenderingPythonTkWidgets.so
RUN update-alternatives --install /usr/bin/vtk vtk /usr/bin/vtk6 10
