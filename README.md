# Dumpy: A Compact and Adaptive Index for Large Data Series Collections

Dumpy is a data series index that employs a novel, adaptive multi-ary data structure and a node packing algorithm, designed to merge small leaf nodes.
It improves the efficiency and search accuracy of SAX-based indexes in the meanwhile.
Two variants of Dumpy, DumpyFuzzy and Dumpy-Memory that lead to more accurate approximate search results for on-disk and in-memory datasets, respectively.

# License
This archive is free for use for academic and non-profit purposes, but if you use it, please reference it properly.

# Reference
Zeyu Wang, Qitong Wang, Peng Wang, Themis Palpanas, and Wei Wang. 2023. Dumpy: A Compact and Adaptive Index for Large Data Series Collections. Proc. ACM Manag. Data 1, 1, Article 111 (May 2023), 27 pages. https://doi.org/10.1145/3588965

# Disclaimer
The code is provided without warranty of any kind. While we thoroughly tested all code bases on Ubuntu 20.04 LTS (Windows Subsystem of Linux), we do not guarantee that they are exempt from bugs, nor that they will work on other platforms. If you encounter any issues with the code, please feel free to propose them on the ISSUE page of this repo. We will do our best to address your concerns but do not promise to resolve all issues.

# Installation

## Prerequisite

A library is required: boost:serialization.
It is used to serialize the index structure.

## Build

Here is a Cmake project. All the information about compiling is in CMakeList.txt.
We test this project in Cmake 3.20 with C++ language standard 23, however, we SUBJECTIVELY think it can run on relatively lower versions.

*If you can use IDE like Clion*, just open this project use that and build it automatically and then run it.

*Else*, please refer to following instructions.

1. create a "build" directory under the project

2. cd build

3. cmake ..

4. make

5. cd .. && ./bin/Dumpy

## Hints

1. All the configuration is written on the config.ini, including the information about dataset.
To finish your task, please READ and UPDATE config.ini IN DETAIL.

2. A graph skeleton needs to be built beforehand. (index=0)

3. SAX table are strongly recommended to be built before building the whole index. It benefits all types of indexes.

4. SIMD techinques are important for DTW distance. But if your machine doesn't support that or you don't want to use this, just comment the compile command on CMakeList.txt `haswell` and all the codes that use AVX2 commands. Note that all SIMD functions in our repo have a corresponding SISD version.

## Datasets

The datasets we used in the paper (Rand, DNA, ECG, Deep) are now in OneDrive.

See the datasets.txt to find how to access them.

## Issue

If you have any problem, feel free to write your questions in ISSUE.
