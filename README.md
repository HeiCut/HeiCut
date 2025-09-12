# Shared-Memory Hypergraph Minimum Cut

## Dependencies

### Mt-KaHyPar
This project depends on the repositories Mt-KaHyPar and kahypar-shared-resources. For more information, see `.gitmodules`. To add these two repositories as submodules, make sure to run the following command:
```
git submodule update --init --recursive
```
Besides, Mt-KaHyPar must be compiled locally, since we do not compile it automatically in this project. The interface `libmtkahypar.so` must be put under `extern/mt-kahypar-library/`.

Since the project depends on Mt-KaHyPar, the following packages must be installed:
```
sudo apt-get install libtbb-dev libhwloc-dev libboost-program-options-dev
```
**Note:** Mt-KaHyPar requires a minimum version of oneTBB, which is sometimes not met by simply installing `libtbb-dev`. In this case, the oneTBB repository must be cloned and built locally using the following commands:
```
git clone https://github.com/oneapi-src/oneTBB.git
cd oneTBB && mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=/usr -DTBB_TEST=OFF ..
sudo cmake --build .
sudo cmake --install .
```

### Gurobi

This projects uses Gurobi as LP solver. To install Gurobi, create a Gurobi account, download [Gurobi for Linux](https://www.gurobi.com/downloads/gurobi-software/) and follow the [installation description](https://support.gurobi.com/hc/en-us/articles/4534161999889-How-do-I-install-Gurobi-Optimizer). Finally download the `gurobi.lic` key and put it under `/opt/gurobi/`.

## Build

You can build the project with CMAKE by using the following command:
```
mkdir build && cd build && cmake -DCMAKE_BUILD_TYPE=Release .. && make
```

## Usage

For each algorithm, you can use `--help` to get a full list of the supported arguments. The default file format of the input hypergraph is expected to be `HMETIS`, but all algorithms also support the `METIS` file format.



### HeiCut

An example usage of HeiCut **without** label propagation and with the tight vertex ordering is:
```
./kernelizer PATH_TO_HYPERGRAPH --ordering_type=tight
```

An example usage of HeiCut **with** one round of label propagation and with the tight vertex ordering is:
```
./kernelizer PATH_TO_HYPERGRAPH --ordering_type=tight --lp_num_iterations=1
```
**Note:** Enable the verbose output via `--verbose` to see the performance of each individual reduction rule.

### BIP

An example usage of the BIP with a timeout of 2 hours is:
```
./ilp PATH_TO_HYPERGRAPH --ilp_timeout=7200
```

### Trimmer

An example usage of Trimmer with the tight vertex ordering is:
```
./trimmer PATH_TO_HYPERGRAPH --ordering_type=tight
```

### Submodular Vertex-Ordering

An example usage of the submodular algorithm with the tight vertex ordering is:
```
./submodular PATH_TO_HYPERGRAPH --ordering_type=tight
```

### k-Core Generator

An example usage the k-core hypergraph generator is:
```
./kcore_generator PATH_TO_HYPERGRAPH PATH_TO_OUTPUT
```