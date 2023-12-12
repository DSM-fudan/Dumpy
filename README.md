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

## Reproducibility
After successfully compiling this project, you can reproduce the results in our paper by following procedures.
All the configurations are integrated into `config.ini`, including all the parameters and all the functions of this repo.

0. Build a skeleton graph. Set in `GraphConstruction.cpp` segment number (e.g., 16), and number of bits for fetching the neighbors of a SAX word, (e.g., 3, at most 4).
Set `index = 0` in `config.ini`, and run `./bin/Dumpy`.
After that, update the [other] section in `config.ini`.

1. Build the index. Set `index = 1` if Dumpy, and `index = 2` if Dumpy-Fuzzy. Set `ops = 0` for building index.
`th` is for the leaf node capacity, usually 10,000.
`segmentNum` is the number of segments, also set in `Const.h`, usually 16.
`bitsCardinality` is for the cardinality of bits of SAX words, usually 8.
`fbl_size` is the buffer size, in MB, depending on your machine.
`max_diff`, `fuzzy`, `delta` are for Dumpy-Fuzzy.
`small_perc` is to define a small leaf node for leaf node packing, 1.0 for an exhausted packing.
`max_mask_bit_percentage` is for ensuring the quality of the leaf pack, usually 0.8.
`f_low`, `f_high` are for adaptive split, usually 0.5 and 1.5.
`alpha` is the weighting factor in adaptive split

- [Dataset], the dataset you use should have a section in `config.ini`. Update the paths and the basic info inside.

After confirming these parameters, run `./bin/Dumpy`.

2. Run the queries.
Dumpy supports many types of queries.
- `index = 1 or 5` fast approximate queries, only searching ONE node, ED and DTW, respectively.
- `index = 2 or 8` exact queries, ED and DTW, respectively.
- `index = 3` get the statistics of the index.
- `index = 4 or 6` extended approximate queries. Decide the number of nodes in `nprobe` in `config.ini`, ED and DTW respectively.
- `index = 7` ng-search, to get the accuracy-time curve.
`dtw_window_percent` is to determine the dtw window size.

After confirming these parameters, run `./bin/Dumpy`, then the results will be printed to your screen.

## Datasets

The datasets we used in the paper (Rand, DNA, ECG, Deep) are now in OneDrive.

See the datasets.txt to find how to access them.

## Issue

If you have any problem, feel free to write your questions in ISSUE.
