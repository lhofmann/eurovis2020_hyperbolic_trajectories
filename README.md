# Extraction of Distinguished Hyperbolic Trajectories for 2D Time-Dependent Vector Field Topology 

[![CI](https://github.com/lhofmann/eurovis2020_hyperbolic_trajectories/workflows/CI/badge.svg)](https://github.com/lhofmann/eurovis2020_hyperbolic_trajectories/actions)

Lutz Hofmann and Filip Sadlo, 
Computer Graphics Forum, 
Proceedings of EuroVis 2020.

This source code provides an implementation of the approximate distinguished hyperbolic trajectory method. It is provided as a VTK class `vtkApproximateDHT`, that takes a vector field and a candidate line as input, and outputs the computed approximation.

An example python script (example.py) is provided, that demonstrates the usage of the VTK class. Executing it requires building the VTK module (see below). Optionally, a corresponding ParaView plugin can be built.

## Using docker

A docker environment based on Ubuntu 18.04 is provided. The script `docker.sh` automatically builds the docker image, starts a container, builds the software, and runs the example python script.

## Build requirements

* C++ compiler (tested with GCC v9.2.0)
* CMake (tested with v3.15.4)
* VTK (tested with v8.2.0) or ParaView (tested with v5.6.3 and v5.7.0)
* Eigen3 (tested with v3.3.7)
* Boost (tested with v1.71.0)

## Installation

Run CMake and make sure, that the paths to Eigen3 and Boost are found (header-only libraries suffice). For building a ParaView plugin, enable the option `ENABLE_PARAVIEW`. VTK or ParaView might need their respective paths to be set (`VTK_DIR` or `ParaView_DIR`). 

Running the python script requires the directory `lib/pythonX.Y/site-packages` to be in your PYTHONPATH. The ParaView plugin is built in the directory `bin/plugins` or `lib/plugins`, depending on the version. 

Alternatively, the `install` target copies the required files in locations relative to your VTK or ParaView installation, such that they can be found automatically.
