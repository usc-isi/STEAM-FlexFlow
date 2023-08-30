# USC CARC Discovery

This directory contains job scripts for the [USC Center for Advanced Research Computing (CARC)](https://www.carc.usc.edu/) Discovery cluster.

All commands here are executed on the login nodes.


## Getting Started

First clone this repository, preferably into the `/project` filesystem (in our case, `/project/jpwalter_148/${USER}/`):

```sh
git clone https://github.com/usc-isi/STEAM-FlexFlow.git
cd STEAM-FlexFlow
git submodule init
git submodule update
```

Then apply the `ffsim-opera.diff` in this directory to the `ffsim-opera` submodule in this repository:

```sh
cd ffsim-opera
git apply ../exp/usc-discovery/ffsim-opera.diff
```


## Building

To build STEAM-FlexFlow on Discovery, first purge and load the following software modules:

```sh
module purge
module load usc/8.3.0 cmake/3.16.2 python/3.7.6 cuda/11.2.0 cudnn/8.1.0.77-11.2-cuda hdf5/1.10.6 zlib/1.2.11
```

Then from the top-level directory:

```sh
mkdir build
cd build
../config/config.linux
make -j8
```

Note: do not use more than 8 threads for building, otherwise the login node will kill all of your SSH sessions.


### htsim

First purge and load different software modules:

```sh
module purge
module load usc gcc/11.3
```

Now build the `ffsim-opera` subdirectory, replacing `/path/to/STEAM-FlexFlow` with your own path (in our case, `/project/jpwalter_148/${USER}/STEAM-FlexFlow`).

```sh
cd ffsim-opera/src/clos
FF_HOME=/path/to/STEAM-FlexFlow make -j8
cd datacenter
FF_HOME=/path/to/STEAM-FlexFlow make -j8
```


### htsim for Opera

Despite the directory name above, we use a different htsim for Opera simulations.

First purge and load different software modules:

```sh
module purge
module load usc gcc/11.3
```

Now build the `ffsim-opera-v2` subdirectory, replacing `/path/to/STEAM-FlexFlow` with your own path (in our case, `/project/jpwalter_148/${USER}/STEAM-FlexFlow`).

```sh
cd ffsim-opera-v2/
FF_HOME=/path/to/STEAM-FlexFlow make -j8
cd datacenter
FF_HOME=/path/to/STEAM-FlexFlow make -j8
```


## Creating Opera topologies

Create a new directory (e.g., `opera_topologies`) somewhere convenient to store Opera topology files.
To create these files, first clone a separate repository (tested on revision `c098613`) and load the compatible Matlab software module:

```sh
git clone https://github.com/TritonNetworking/opera-sim.git
cd cd opera-sim/topologies/opera_dynexp_topo_gen
module load matlab/2019a
```

Edit the file `defineTopology_simple_1path.m` to increase the simulated bandwidth: change line 28 from `linkrate=10e9; % bits/s` to `linkrate=100e9; % bits/s`.

Edit the file `pick_topo.m` to fix a typo: change line 8 from `load(sprintf('lifted_Decomp_Nbase=%d_N=%d.mat',Nbase,N),'P');` to `load(sprintf('Lifted_Decomp_Nbase=%d_N=%d.mat',Nbase,N),'P');`

Then edit the file `MAIN.m` to modify the `N` and `k` values:

* For a 16-node topology: `N=4`, `k=8`
* For a 128-node topology: `N=32`, `k=8`
* For a 1024-node topology: `N=128`, `k=16`

After each edit, run the script:

```sh
matlab -nodisplay -nosplash -nodesktop -singleCompThread -sd "$(pwd)" -batch "run('MAIN.m'); exit;"
```

Then move the resulting topology files to your topology directory created above, e.g.:

```sh
mv dynexp_1path_N=4_k=8_G=1.txt dynexp_1path_N=32_k=8_G=1.txt dynexp_1path_N=128_k=16_G=1.txt /path/to/opera_topologies/
```


## Executing Experiments

### Prerequisites

The job script examples in this repository have hardcoded identifiers and paths that you will need to update appropriately.

First, check that the `#SBATCH --account` lines use the account identifier that you've been assigned.

* Update all `FF_HOME` values to point to your STEAM-FlexFlow directory.
* Update any `JSON_DIR` and `TG_DIR` paths to be in your working directory tree where fat-tree simulations produce measurement (.json) and taskgraph (.tg) files, respectively.
* Update `TOPO_DIR` to point to the location of the Opera topology files created above.
* Update `OPERA_HOME` to point to your `ffsim-opera-v2` directory in STEAM-FlexFlow (we had this in a separate `opera-w-ffapp` directory).


### Submitting Jobs

Discovery uses SLURM as its workload manager.
To submit a job, use `sbatch`, e.g.:

```sh
cd exp/usc-discovery/model-a/nodes-16
sbatch model_a_16.job
```

To track submitted/running jobs, check the SLURM queue, e.g.:

```sh
squeue --me
```

Commands that use STEAM-FlexFlow's `dlrmsim` application must run on Discovery's `gpu` partition with the `--gres=gpu:a100:1` SLURM option.
Many jobs also run the htsim applications (i.e., `htsim_tcp_fattree` or `htsim_tcp_flat`) within the same script for convenience.
Some larger simulations and the Opera simulations (i.e., `htsim_ndp_dynexpTopology`) use the `largemem` partition (a few may use the `main` partition instead if they are small enough).
