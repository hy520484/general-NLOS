# General NLOS

This repository contains the code and data for the paper **"Fast Non-line-of-sight Imaging with an Arbitrary Scanning Setup"**.

# Environment Requirements

The code has been tested in the following environments. Other reasonable environments should be OK.

## Hardware Environment

 + OS: Windows 10

 + CPU: Intel i9-9900K

 + GPU: NVIDIA 1080 Ti

## Software Environment

 + Visual Studio 2019

 + CUDA 11.3 Runtime

 + Python 3.7.0

 + pycuda 2022.1

# Quick Start

Please first ensure the software requirements (VS, CUDA, Python, pycuda) are satisfied. Then choose one config file in `reconstruct.py` and run.

# File Description

## Config File

 Most of the parameters of the reconsturction procedure are defined in the config file. The config files corresponding to the results shown in the paper are in the folder `config`. Here we briefly describe the meaning of each parameter.

 + `path_data`: The path to the NLOS measurements.

 + `path_setup`: The path to the NLOS setup, which contains the 3D coordinates of laser points and detection points.

 + `width`: A tunable parameter that control the filters.

 + `N_ld`: A total of `N_ld` * `N_ld` pairs of laser point and detection point.

 + `N_bin`: The number of time bins of the NLOS measurments. The size of the NLOS measurements should be `N_ld` * `N_ld` * `N_bin`.

 + `timeRes`: The time resolution of the NLOS measurements.

 + `N_voxel_s` and `N_voxel_t`: The number of voxels of the hidden space is `N_voxel_s` * `N_voxel_s` * `N_voxel_t`.

 + `x_range0`, ..., `z_range1`: The range of the hidden space.

 + `sigma`: Another tunable parameter that control the filters. If not specified, its default value is 0.23.

 + `th`: A threshold that helps displaying a better reconstructed result (please refer to the function `simple_denoise()` in `utils/show_result.py` for further details). If not specified, its default value is 0.

## Experimental Data

 The experimental data captured by our system are in the folder `experimental_data`. The folder `tuan` contains the data of NLOS imaging with irregular scanning pattern, and the folder `wall` contains the data of NLOS imaging with non-planar relay wall.

## Public Data

 The foler `public_data` contains the data and code of [CC-SOCR](https://doi.org/10.1038/s41467-023-38898-4) [1] and [3D-RSD](https://github.com/ArianaGu/3D-RSD) [2].

# References

1. Liu, Xintong, *et al.* "Non-line-of-sight imaging with arbitrary illumination and detection pattern." *Nature Communications* 14.1 (2023): 3230.

2. Gu, Chaoying, *et al.* "Fast Non-line-of-sight Imaging with Non-planar Relay Surfaces." *2023 IEEE International Conference on Computational Photography (ICCP)*. IEEE, 2023.

