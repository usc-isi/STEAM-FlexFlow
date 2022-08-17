# ideal switch
./htsim_tcp_abssw -simtime 3600.1 -ssthresh 10000 -rtt 1000 -q 10000 -flowfile tg.ideal.1600G.candle16.nd8   -ofile ./2022-07-12_14:18:17/performance.ideal.candle16.speed1600.nd8 -speed 12800000  -nodes 16  
# fat tree switch
./htsim_tcp_fattree -simtime 3600.1 -ssthresh 10000 -rttnet 1000 -rttrack 1000 -q 10000  -flowfile tg.fattree.100G.candle128.nd4   -ofile ./2022-07-12_14:18:17/performance.fattree.candle128.speed100.nd4 -speed 100000  -nodes 128  
# topoopt switch 
./htsim_tcp_flat -simtime 3600.1 -ssthresh 10000 -rtt 1000 -q 50000  -flowfile tg.topoopt.100G.candle128.nd4   -ofile ./2022-07-12_14:18:17/performance.topoopt.candle128.speed100.nd4 -speed 100000  -nodes 128  
