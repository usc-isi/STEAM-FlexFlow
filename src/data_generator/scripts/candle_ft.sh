#!/bin/bash

FF_HOME="/home/weiyangw/FlexFlow"

netopt=1
igbw=200
gdbw=256
nlat=1            # network latency 1 us
local_b=256       # local batch size 256
biggpu=4          # Emulate 4 GPU on one machine
d=4               # Node degree 4
b=100             # Bandwidth per link 100Gbps
n=128             # 128 Nodes in the cluster
topo="fattree"    # fattree cluster
 
globalb=$((n*local_b*biggpu)) # global batch size
mfile="${FF_HOME}/measures/candle_128.json"
resultdir="a100_candle_${topo}_${n}_${b}_${d}_${nlat}_${local_b}_${rid}"
${FF_HOME}/build/Release/examples/cpp/candle_unosim/candle_unosim -ll:gpu 1 -ll:cpu 1 -ll:zsize 20000 -ll:fsize 10000 -ll:util 4 -dm:memoize --dense-feature-layers 16384-16384-16384-16384-16384-16384-16384-16384-16384-16384-16384-16384-16384-16384-16384-16384 --dense-layers 16384-16384-16384-16384-16384-16384-16384-16384-1 --batch-size ${globalb} --interface-bandwidth $b --inter-gpu-bandwidth $igbw --gpu-dram-bandwidth $gdbw --network-latency $nlat --net-opt $netopt --nsimnode $n --search-budget 4000 --mfile $mfile  --enable-propagation --node-degree $d --taskgraph taskgraph_CANDLE100GFT.fbuf --simulator-workspace-size 65536 --big-gpu $biggpu --topology $topo

mkdir $resultdir
mv taskgraph_CANDLE100GFT.fbuf $resultdir/taskgraph.fbuf

bash -c "cd $resultdir && ${FF_HOME}/ffsim-opera/src/clos/datacenter/htsim_tcp_fattree -simtime 3600.1 -flowfile ./taskgraph.fbuf -speed $((b*1000)) -ofile nwsim_linkft.txt -nodes $n -ssthresh 10000 -rttrack 1000 -rttnet 1000 -q 10000"

