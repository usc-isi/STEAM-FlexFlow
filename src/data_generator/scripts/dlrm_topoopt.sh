#!/bin/bash

FF_HOME="/home/weiyangw/FlexFlow"

netopt=1
igbw=200
gdbw=256
nlat=1            # network latency 1 us
local_b=128       # local batch size 128
biggpu=4          # Emulate 4 GPU on one machine
d=4               # Node degree 4
b=100             # Bandwidth per link 100Gbps
n=128             # 128 Nodes in the cluster
topo="topoopt"    # TopoOpt cluster
 
globalb=$((n*local_b*biggpu)) # global batch size
mfile="${FF_HOME}/measures/dlrm128.json"
resultdir="a100_dlrm_${topo}_${n}_${b}_${d}_${nlat}_${local_b}_${rid}"
${FF_HOME}/build/Release/examples/cpp/DLRMsim/dlrmsim -ll:gpu 1 -ll:cpu 1 -ll:zsize 20000 -ll:fsize 10000 -ll:util 4 -dm:memoize --embedding-bag-size 100 --arch-sparse-feature-size 128 --arch-embedding-size 10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000 --arch-mlp-bot 2048-2048-2048-2048-2048-2048-2048-2048 --arch-mlp-top 4096-4096-4096-4096-4096-4096-4096-4096-4096-4096-4096-4096-4096-4096-4096-4096-1  --batch-size ${globalb} --interface-bandwidth $b --inter-gpu-bandwidth $igbw --gpu-dram-bandwidth $gdbw --network-latency $nlat --net-opt $netopt --nsimnode $n --search-budget 4000 --mfile $mfile  --enable-propagation --node-degree $d --taskgraph taskgraph_DLRM100G.fbuf --simulator-workspace-size 65536 --big-gpu $biggpu --topology $topo
mkdir $resultdir
mv taskgraph_DLRM100G.fbuf $resultdir/taskgraph.fbuf

bash -c "cd $resultdir && ${FF_HOME}/ffsim-opera/src/clos/datacenter/htsim_tcp_flat -simtime 3600.1 -q 50000 -flowfile ./taskgraph.fbuf -speed $((b*1000)) -ofile nwsim.txt -nodes $n -ssthresh 10000 -rtt 1000"
