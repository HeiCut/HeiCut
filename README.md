# HeiCut

![C++17](https://img.shields.io/badge/C++-17-blue.svg?style=flat)  
![Linux](https://img.shields.io/badge/OS-Linux-green.svg?style=flat)  
![ALENEX](https://img.shields.io/badge/Conference-ALENEX%2026-orange.svg?style=flat)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17140472.svg)](https://doi.org/10.5281/zenodo.17140472)

---

## Overview

**HeiCut** is a highly efficient, exact solver for the **minimum cut problem in hypergraphs**.  

The algorithm performs repeated rounds of **provably exact reduction rules** that preserve the minimum cut while drastically shrinking instance size. Once no further reductions are possible, an exact solver is applied to compute the minimum cut.

HeiCut is presented in the paper:  
> *Exact Minimum Cuts in Hypergraphs at Scale*, ALENEX 2026.

This repository provides:
- Source code for **HeiCut**  
- Implementations of competing algorithms: **Trimmer**, **Relaxed BIP**, and **ordering-based solvers**  
- Benchmark datasets (`M_{HG}`, `L_{HG}`, and $(k,2)$-core)  
- Scripts to **reproduce all experimental results** from the paper

The source code and benchmark datasets are permanently archived on [Zenodo Software](https://doi.org/10.5281/zenodo.17140472) and [Zenodo Dataset](https://doi.org/10.5281/zenodo.17142170) to ensure long-term availability.

---

## Requirements

HeiCut requires the following:

- A 64-bit **Linux** operating system  
- A modern **C++17** compiler (`g++ ≥ 7` recommended)  
- [CMake][cmake] (≥ 3.16)  
- [Boost.Program_options][Boost.Program_options] (≥ 1.48)  
- [oneTBB][tbb] (≥ 2021.5.0)  
- [hwloc][hwloc]
- [SparseHash](https://github.com/sparsehash/sparsehash) 
- [Gurobi](https://www.gurobi.com/) (used as LP solver)  
- [Mt-KaHyPar](https://github.com/kahypar/mt-kahypar/tree/0ef674ad44c35fb4f601a7eddd3f4f23f0d5d60a), commit `0ef674a`  

Detailed installation instructions are provided below.

---

## Getting Started

Clone the repository **with pre-loaded submodules** (extern/kahypar-shared-resources and extern/mt-kahypar):

```bash
git clone https://github.com/HeiCut/HeiCut.git
```

---

## Dependencies

On Ubuntu, install core dependencies via:

```bash
sudo apt-get install libtbb-dev libhwloc-dev libboost-program-options-dev
```

> **Note:** `libtbb-dev` may be outdated. If so, clone and build [oneTBB](https://github.com/oneapi-src/oneTBB) locally:

```bash
git clone https://github.com/oneapi-src/oneTBB.git
cd oneTBB && mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=/usr -DTBB_TEST=OFF ..
sudo cmake --build .
sudo cmake --install .
```

To install [sparsehash](https://github.com/sparsehash/sparsehash), follow these instructions:
```bash
git clone https://github.com/sparsehash/sparsehash
cd sparsehash
./configure
make install
```

---

### Gurobi

HeiCut uses **Gurobi** as its LP solver.  

1. Create a [Gurobi account](https://www.gurobi.com/).  
2. Download [Gurobi for Linux](https://www.gurobi.com/downloads/gurobi-software/).  
3. Follow the [installation guide](https://support.gurobi.com/hc/en-us/articles/4534161999889-How-do-I-install-Gurobi-Optimizer).
4. Obtain a license to use Gurobi. For instance, you can get a free academic license by following the [official guide](https://www.gurobi.com/academia/academic-program-and-licenses/). 
5. Place your `gurobi.lic` license file in the installation folder (e.g., `/opt/gurobi1203/`).  
6. Add environment variables (update path if needed):  

```bash
export GUROBI_HOME="/opt/gurobi1203/linux64"
export PATH="${PATH}:${GUROBI_HOME}/bin"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${GUROBI_HOME}/lib"
```

---

### Mt-KaHyPar

HeiCut depends on **Mt-KaHyPar** (commit `0ef674a`).  

To install automatically:

```bash
./install_mtkahypar.sh
```

This builds the library and places `libmtkahypar.so` in `HeiCut/extern/mt-kahypar-library/`.

Manual build instructions (if preferred):

```bash
git clone --depth=2 --recursive https://github.com/kahypar/mt-kahypar.git
cd mt-kahypar
git checkout 0ef674a
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make install.mtkahypar   # may need sudo
```
**Note:** If installed locally, the build will exit with an error due to missing permissions. However, the library is still built successfully and is available in the build folder.
Locate `libmtkahypar.so` (usually in `build/lib/`) and copy it into:  
`HeiCut/extern/mt-kahypar-library/`.

#### ⚠️ Known issue (rare): `growt` ref/_mref compilation error

A few users have encountered a **compilation error in Mt-KaHyPar’s `growt` dependency** that looks like:

```
.../external_tools/growt/data-structures/migration_table_iterator.hpp:68:22: error:
‘... migration_table_mapped_reference ...’ has no member named ‘ref’; did you mean ‘_mref’?
68 |                 sref.ref.refresh();
|                      ^~~
|                      _mref
```

If you see this, and the program does not compile, apply the following **manual fix** in the `growt` source and rebuild:

1. Open:
```
<your-mtkahypar-source>/external_tools/growt/data-structures/migration_table_iterator.hpp
````

2. In the class `migration_table_mapped_reference`, **edit the constructor’s initializer list**:

change (possibly line 57) 
```cpp
: _tab(table), _version(ver), _mref(mref)
````
to
```cpp
: _tab(table), _version(ver), _mref(mref), ref(_mref)
```

3. Insert (possibly line 133): 
```cpp
public:
   base_mapped_reference& ref;
````

4. Rebuild Mt-KaHyPar (or rerun your previous build command)
---

## Building HeiCut

After installing all dependencies:

```bash
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

Binaries are placed in `HeiCut/build/`.

---

## Usage

All executables support `--help` to list available arguments.  
Default hypergraph format: **hMETIS** (see [hMETIS manual](https://course.ece.cmu.edu/~ee760/760docs/hMetisManual.pdf)).  
`METIS` format is also supported.

---

### HeiCut

Binary: `kernelizer`

Example (tight ordering, no label propagation):  
```bash
./kernelizer PATH_TO_HYPERGRAPH --ordering_type=tight
```

Example (tight ordering + 1 round of label propagation):  
```bash
./kernelizer PATH_TO_HYPERGRAPH --ordering_type=tight --lp_num_iterations=1
```

> **Tip:** Use `--verbose` to view detailed reduction performance.

---

### Relaxed-BIP

Binary: `ilp`

Example (timeout 2 hours):  
```bash
./ilp PATH_TO_HYPERGRAPH --ilp_timeout=7200
```

---

### Trimmer

Binary: `trimmer`

Example (tight ordering):  
```bash
./trimmer PATH_TO_HYPERGRAPH --ordering_type=tight
```

---

### Vertex-Ordering Solver

Binary: `submodular`

Example (tight ordering):  
```bash
./submodular PATH_TO_HYPERGRAPH --ordering_type=tight
```

---

## Reproducing Paper Results
In the experimental section of our paper, we provide results seperated by dataset.

### Benchmark Datasets

We evaluate on three datasets:  

- **M<sub>HG</sub>** (488 medium instances): Weighted and Unweighted
- **L<sub>HG</sub>** (94 large instances, up to 139M vertices): Weighted and Unweighted  
- **(k,2)-core** (44 synthetic instances with non-trivial cuts): Weighted and Unweighted  

Download all datasets here:  
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17142170.svg)](https://doi.org/10.5281/zenodo.17142170)

Extract into the repository root (`/HeiCut/`):  

- `/HeiCut/med_set/`  
- `/HeiCut/large_set/`  
- `/HeiCut/k,2-core_benchmark/`

**Note**: To save on time and disk space, we only provide weighted versions of the hypergraphs. For reproducing results on the unweighted versions, we pass a flag `--unweighted` to each algorithm to process the weighted hypergraphs as an unweighted (set all edge weights to 1). 

---

### Experimental Scripts

Scripts are located in `/HeiCut/experiments/`, organized per dataset:  

- `medium_weighted`, `medium_unweighted`  
- `large_weighted`, `large_unweighted`  
- `k-core_weighted`, `k-core_unweighted`

Each folder contains a master script that you should run (e.g., `perform_medium_weighted_experiments.sh`). **Important**: If you wish to use our experimental and plotting scripts, avoid using special characters (e.g., spaces, #, %, &) in the HeiCut directory path. Such characters can break file-path parsing in the experimental/plotting scripts and lead to runtime errors.


**Note:** we also provide a global script to run all experiments on all datasets called `perform_all_experiments.sh`. However, we do not recommend using it since it would be simply too time consuming.

> **Dependency:** [GNU parallel](https://www.gnu.org/software/parallel/) and [GNU Time](https://www.gnu.org/software/time/)
> 
> Install on Ubuntu:  
> ```bash
> sudo apt install parallel time
> ```

---

### Setup Notes

- **Time limits:** 2h per algorithm  
- **Memory limits:**  
  - 100GB for `M<sub>HG</sub>`  
  - 300GB for `L<sub>HG</sub>` & `(k,2)-core`  
- For the paper experiments, we use gnu parallel to run 13 instances of experiments on the `M_{HG}` dataset, and 4 instances of experiments on the `L_{HG}` and $(k, 2)$-core datasets concurrently. This was possible since we used a machine with a large amount of memory. As the amount of memory available on the machine you use can vary, by default, in the provided experimetal scripts, we run **1 instance at a time** by default (adjustable).

---

### Results Output

Results are written to:  
`/HeiCut/experiments/<dataset>/generated/all_results/`

Each algorithm produces:  
- **CSV summaries** 
- **Per-instance results**  

Example:  
- `/kernelizer_IT0/all_results.csv` (HeiCut, no LP)  
- `/kernelizer_IT1/all_results.csv` (HeiCut, with LP)
- `/ilp/all_results.csv` (Relaxed-BIP)
- `/trimmer/all_results.csv` (Trimmer)
- `/submodular_tight_single/all_results.csv` (Tight vertex-ordering solver)



In the output `all_results.csv` for each algorithm, each row contains the statistics for an instance. The first three columns correspond to minimum cut, time, and memory respectively. If an algorithm fails on an instance, the minimum cut column is blank. All algorithms except Relaxed-BIP return the exact minimum cut if successful.  

---

### 📊 Plotting Experimental Results

Each experimental dataset folder — for example:
```
experiments/medium_weighted
```
contains a subfolder:
```
experiments/medium_weighted/plot
````
Once you have completed running all experiments and generated the result files, you can produce the **final performance profile plots** exactly as shown in the paper. Note: all algorithms should have attempted to solve each instance and reported their attempt in its all_results.csv, i.e, there may not be missing rows in any CSV. 

#### 🧭 How to plot

1. Navigate to the `plot` directory of the dataset:
```bash
  cd experiments/medium_weighted/plot
```

2. Run the plotting script:
 ```bash
 ./plot_all.sh
 ```
This will generate **performance profile plots** comparing all algorithms on:

* **Memory usage**
* **Runtime**
* **Minimum cut**

The resulting plot corresponds to the dataset associated with the current folder (e.g., `medium_weighted`) and replicates the paper’s figure format.

---

#### 📦 Plotting Dependencies

The plotting pipeline relies on:

* **R** (tested with R 4.3+)
* **LaTeX** with TikZ support
  (ensure `tikzDevice` can compile LaTeX files — we suggest a **full texlive installation** for a smoother experience)

The following R packages are required:
```
ggplot2, plyr, dplyr, RColorBrewer, tikzDevice,
gridExtra, egg, ggpubr, stringr, stringi, ggrepel
```
The R plotting scripts will **automatically attempt to install missing dependencies**.
However, if something goes wrong (e.g., due to system library issues like ICU), simply install the required packages manually:
```r
install.packages(c(
  "ggplot2","plyr","dplyr","RColorBrewer","tikzDevice",
  "gridExtra","egg","ggpubr","stringr","stringi","ggrepel"
), repos = "https://cloud.r-project.org")
```

---

### Custom Experiments (Optional)

Configurable scripts:  

- `generate_experiments.sh` – select algorithms & params  
- `run_experiments.sh` – set time & parallelism  
- `extract_results.sh` – collect results  

Provide hypergraphs in:  
`/HeiCut/experiments/hypergraphs.txt`

---

### k-Core Generator (Optional)

To generate $(k,2)$-core hypergraphs:  

```bash
./kcore_generator PATH_TO_HYPERGRAPH PATH_TO_OUTPUT
```

---

## References

1. C. J. Alpert, *The ISPD98 Circuit Benchmark Suite*, ISPD 1998. [DOI](https://doi.org/10.1145/274535.274546)  
2. N. Viswanathan et al., *The DAC 2012 Routability-Driven Placement Contest*, DAC 2012. [DOI](https://doi.org/10.1145/2228360.2228500)  
3. T. A. Davis and Y. Hu, *The SuiteSparse Matrix Collection*, ACM TOMS 2011. [DOI](https://doi.org/10.1145/2049662.2049663)  
4. A. Belov et al., *The SAT Competition 2014*. [Link](https://satisfiability.org/competition/2014/)  

[cmake]: http://www.cmake.org/ "CMake tool"  
[Boost.Program_options]: http://www.boost.org/doc/libs/1_58_0/doc/html/program_options.html  
[tbb]: https://github.com/oneapi-src/oneTBB  
[hwloc]: https://www.open-mpi.org/projects/hwloc/  
