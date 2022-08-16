
appname = ['dlrm', 'candle']
app = ['./DLRMsim/dlrmsim', './candle_unosim/candle_unosim']
prof = ['measures/dlrm', 'measures/candle']
nodes = [16, 32] #, 128, 1024]
topo = ['fattree', 'fattree', 'topoopt']
tglabel = ['ideal', 'fattree', 'topoopt']
degree = [16, 32] # 4, 8, 12]
speed = [50] # , 800, 1600]
# args for measure/
#args = [
#  [	# for dlrm 16, 32, 128
#   " -ll:gpu 1 -ll:cpu 1 -ll:zsize 20000 -ll:fsize 10000 -ll:util 4  -dm:memoize --embedding-bag-size 100 --arch-sparse-feature-size 256 --arch-embedding-size 10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000 --arch-mlp-bot 1024-1024-1024-1024 --arch-mlp-top 2048-2048-2048-2048-2048-2048-2048-2048-1 --search-budget 5000 --inter-gpu-bandwidth 256 --gpu-dram-bandwidth 200 --network-latency 1 --net-opt 1 --enable-propagation --batch-size 16384 --big-gpu 4 --simulator-workspace-size 65536 ",
#" -ll:gpu 1 -ll:cpu 1 -ll:zsize 20000 -ll:fsize 10000 -ll:util 4  -dm:memoize --embedding-bag-size 100 --arch-sparse-feature-size 256 --arch-embedding-size 10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000 --arch-mlp-bot 1024-1024-1024-1024-1024-1024-1024-1024 --arch-mlp-top 2048-2048-2048-2048-2048-2048-2048-2048-2048-2048-2048-2048-2048-2048-2048-2048-1 --search-budget 5000 --inter-gpu-bandwidth 256 --gpu-dram-bandwidth 200 --network-latency 1 --net-opt 1 --enable-propagation --batch-size 32768 --big-gpu 4 --simulator-workspace-size 65536 ",
#" -ll:gpu 1 -ll:cpu 1 -ll:zsize 20000 -ll:fsize 10000 -ll:util 4 -dm:memoize --embedding-bag-size 100 --arch-sparse-feature-size 128 --arch-embedding-size 10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000 --arch-mlp-bot 2048-2048-2048-2048-2048-2048-2048-2048 --arch-mlp-top 4096-4096-4096-4096-4096-4096-4096-4096-4096-4096-4096-4096-4096-4096-4096-4096-1 --batch-size 65536 --inter-gpu-bandwidth 200 --gpu-dram-bandwidth 256 --network-latency 1 --net-opt 1 --search-budget 4000 --enable-propagation --simulator-workspace-size 65536 --big-gpu 4 "
#],
#[	# for candle 16, 32, 128
#" -ll:gpu 1 -ll:cpu 1 -ll:zsize 20000 -ll:fsize 10000 -ll:util 4  -dm:memoize --dense-feature-layers 4096-4096-4096-4096-4096-4096-4096-4096-4096-4096-4096-4096-4096-4096-4096-4096 --dense-layers 4096-4096-4096-4096-4096-4096-4096-4096-1 --search-budget 5000 --inter-gpu-bandwidth 256 --gpu-dram-bandwidth 200 --network-latency 1 --net-opt 1 --enable-propagation --batch-size 16384 --big-gpu 4 --simulator-workspace-size 65536 ",
#" -ll:gpu 1 -ll:cpu 1 -ll:zsize 20000 -ll:fsize 10000 -ll:util 4  -dm:memoize --dense-feature-layers 8192-8192-8192-8192-8192-8192-8192-8192-8192-8192-8192-8192-8192-8192-8192-8192-8192-8192-8192-8192-8192-8192-8192-8192-8192-8192-8192-8192-8192-8192-8192-8192 --dense-layers 8192-8192-8192-8192-8192-8192-8192-8192-8192-8192-8192-8192-8192-8192-8192-8192-1 --search-budget 5000 --inter-gpu-bandwidth 256 --gpu-dram-bandwidth 200 --network-latency 1 --net-opt 1 --enable-propagation --batch-size 32768 --big-gpu 4 --simulator-workspace-size 65536 " ,
#" -ll:gpu 1 -ll:cpu 1 -ll:zsize 20000 -ll:fsize 10000 -ll:util 4 -dm:memoize --dense-feature-layers 16384-16384-16384-16384-16384-16384-16384-16384-16384-16384-16384-16384-16384-16384-16384-16384 --dense-layers 16384-16384-16384-16384-16384-16384-16384-16384-1 --batch-size 131072 --inter-gpu-bandwidth 200 --gpu-dram-bandwidth 256 --network-latency 1 --net-opt 1 --search-budget 4000 --enable-propagation --simulator-workspace-size 65536 --big-gpu 4 "
#]
#]

# args for candle16_large.json, candle_32_large.json, dlrm16_large.json, dlrm16_large.json
args = [
 [ # for dlrm 16, 32
   " -ll:gpu 1 -ll:cpu 1 -ll:zsize 20000 -ll:fsize 10000 -ll:util 4  -dm:memoize --embedding-bag-size 100 --arch-sparse-feature-size 128 --arch-embedding-size 10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000 --arch-mlp-bot 4096-4096-4096-4096-4096-4096-4096-4096 --arch-mlp-top 8192-8192-8192-8192-8192-8192-8192-8192-8192-8192-8192-8192-8192-8192-8192-8192-1 --search-budget 5000 --inter-gpu-bandwidth 256 --gpu-dram-bandwidth 200 --network-latency 1 --net-opt 1 --enable-propagation --batch-size 8192 --big-gpu 4 --simulator-workspace-size 65536 ",
" -ll:gpu 1 -ll:cpu 1 -ll:zsize 20000 -ll:fsize 10000 -ll:util 4  -dm:memoize --embedding-bag-size 100 --arch-sparse-feature-size 128 --arch-embedding-size 10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000-10000000 --arch-mlp-bot 4096-4096-4096-4096-4096-4096-4096-4096 --arch-mlp-top 8192-8192-8192-8192-8192-8192-8192-8192-8192-8192-8192-8192-8192-8192-8192-8192-1 --search-budget 5000 --inter-gpu-bandwidth 256 --gpu-dram-bandwidth 200 --network-latency 1 --net-opt 1 --enable-propagation --batch-size 16384 --big-gpu 4 --simulator-workspace-size 65536 ",
 ],
 [ # for candle 16, 32
" -ll:gpu 1 -ll:cpu 1 -ll:zsize 20000 -ll:fsize 10000 -ll:util 4  -dm:memoize --dense-feature-layers 16384-16384-16384-16384-16384-16384-16384-16384-16384-16384-16384-16384-16384-16384-16384-16384 --dense-layers 16384-16384-16384-16384-16384-16384-16384-16384-1 --search-budget 5000 --inter-gpu-bandwidth 256 --gpu-dram-bandwidth 200 --network-latency 1 --net-opt 1 --enable-propagation --batch-size 16384 --big-gpu 4 --simulator-workspace-size 65536 ",
" -ll:gpu 1 -ll:cpu 1 -ll:zsize 20000 -ll:fsize 10000 -ll:util 4  -dm:memoize --dense-feature-layers 16384-16384-16384-16384-16384-16384-16384-16384-16384-16384-16384-16384-16384-16384-16384-16384 --dense-layers 16384-16384-16384-16384-16384-16384-16384-16384-1 --search-budget 5000 --inter-gpu-bandwidth 256 --gpu-dram-bandwidth 200 --network-latency 1 --net-opt 1 --enable-propagation --batch-size 32768 --big-gpu 4 --simulator-workspace-size 65536 " ,
 ]
]

for ni, n in enumerate(nodes):
    print("# nodes (%d) " % n)
    nodestr = ' --nsimnode ' + str(n)
    for ai, a in enumerate(app):
        print("# app (%s) " % a)
        for ti, t in enumerate(topo):
            topostr = ' --topology ' + t
            for s in speed:
                speedstr = ' --interface-bandwidth ' + str(s)
                for d in degree:
                    degreestr = ' --node-degree ' + str(d) 
                    if (tglabel[ti] == 'ideal'):
                        speedstr = ' --interface-bandwidth ' + str(s*d)
                    profstr = ' --mfile ' + prof[ai] + str(n) + '_large.json'
                    tgraphstr = ' --taskgraph tg.' + tglabel[ti] + '.' + str(s) + 'G.' + appname[ai] + str(nodes[ni]) + ".nd" + str(d)
                    cmd = app[ai] + args[ai][ni] + profstr + tgraphstr + degreestr + speedstr + nodestr + topostr + " --no-gpu "
                    print(cmd)


