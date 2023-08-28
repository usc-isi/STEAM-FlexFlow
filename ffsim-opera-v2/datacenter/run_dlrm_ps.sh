#!/bin/sh

#../htsim_ndp_dynexpTopology -cwnd 30 -q 32 -pullrate 1 -thresh 1500 -simtime 10.100001 -utiltime 5.0 -topfile opera_1path_N=108_k=12_G=1.txt -flowfile candle_ps_from_flexflow


./htsim_ndp_dynexpTopology -cutoff 1500 -rlbflow 0 -cwnd 30 -q 8 -simtime 10.0 -pullrate 0.1 -topfile ./dynexp_N=108_k=12_1path.txt -flowfile ../../jsons/tg.fattree.100G.dlrm16.nd4_InOrder.json

#gdb -ex=r --args ../htsim_ndp_dynexpTopology -cwnd 30 -q 8 -simtime .003001 -pullrate 1 -topfile dynexp_N=8_k=8.txt -flowfile flows.txt

