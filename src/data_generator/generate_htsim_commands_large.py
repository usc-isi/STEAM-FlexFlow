
from datetime import datetime
import os

apps = ["dlrm", "candle"]
#nodes = [8, 16, 32, 128]
bandwidth = [400000]
#links = [4, 8, 16]
links = [4] # , 32]
links_fattree = [1, 2, 3] 
topologies = ["ideal", "fattree", "topoopt", "sip-ml", "intrepid8", "intrepid16"]
topo_dict = { "ideal":0, "fattree": 1, "topoopt": 2, "sip-ml": 3, "intrepid8": 4, "intrepid16": 5}
bandwidth_link_weight = [links, links_fattree, [], [], [], [], []]	# ideal fat tree: d * b
#fattree_nodes = [16, 16, 32, 128]
#topoopt_nodes = [8, 16, 32, 128]
#sip_ml_nodes = [8, 16, 32, 128]
#intrepid8_nodes = [8]
#intrepid16_nodes = [16]
ideal_nodes = [16, 32, 128, 1024]
fattree_nodes = [16, 32, 128, 1024]
topoopt_nodes = [16, 32, 128, 1024]
sip_ml_nodes = []
intrepid8_nodes = []
intrepid16_nodes = []
#fattree_nodes_allocated = [16, 16, 54, 128]
ideal_nodes_allocated = [16, 32, 128, 1024]
fattree_nodes_allocated = [16, 54, 128, 1024]
topoopt_nodes_allocated = topoopt_nodes
sip_ml_nodes_allocated = sip_ml_nodes
intrepid8_nodes_allocated = intrepid8_nodes
intrepid16_nodes_allocated = intrepid16_nodes
topo_nodes = [ideal_nodes, fattree_nodes, topoopt_nodes, sip_ml_nodes, intrepid8_nodes, intrepid16_nodes]
topo_nodes_allocated = [ideal_nodes_allocated, fattree_nodes_allocated, topoopt_nodes_allocated, sip_ml_nodes_allocated, intrepid8_nodes_allocated, intrepid16_nodes_allocated]

ideal_topologies = ["ideal"]
fattree_topologies = ["fattree"]
topoopt_topologies = ["topoopt"]
sip_ml_topologies = ["sip-ml"]
intrepid8_topologies = ["ring1", "ring2", "ring3"]
intrepid16_topologies = ["ring11"]
topo_types = [ideal_topologies, fattree_topologies, topoopt_topologies, sip_ml_topologies, intrepid8_topologies, intrepid16_topologies]

topo_file_prefix = ["ideal", "fattree", "topoopt", "fattree", "fattree", "fattree"]
topo_file_prefix_opt = [[], [], [], [], intrepid8_topologies, intrepid16_topologies]

topo_htsim = ['./htsim_tcp_abssw', './htsim_tcp_fattree', './htsim_tcp_flat', './htsim_tcp_dyn_flat', './htsim_tcp_flat', './htsim_tcp_flat']
sim_args_ideal = " -simtime 3600.1 -ssthresh 10000 -rtt 1000 -q 10000"
sim_args_fattree = " -simtime 3600.1 -ssthresh 10000 -rttnet 1000 -rttrack 1000 -q 10000 "
sim_args_flat = " -simtime 3600.1 -ssthresh 10000 -rtt 1000 -q 50000 "
sim_args_dyn_flat = " -simtime 3600.1 -rdelay 25 -omethod dheu -ssthresh 10000 -rtt 1000 -q 50000 "
topo_sim_args = [sim_args_ideal, sim_args_fattree, sim_args_flat, sim_args_dyn_flat, sim_args_flat, sim_args_flat]
option_link_deg = [[], [], [], links, [], []]
option_link_deg_ifile = [[], [], links, [], [], []]


# datetime object containing current date and time
now = datetime.now()
# dd/mm/YY H:M:S
dt_string = now.strftime("%Y-%m-%d_%H:%M:%S")
dt_string = "./" + dt_string + "/"
os.mkdir(dt_string)
#print("directory created!")

for a_index, a in enumerate(apps):
    for t in topologies:
        t_index = topo_dict[t]
        nodes = topo_nodes[t_index]
        nodes_allocated = topo_nodes_allocated[t_index]
        topo = topo_types[t_index]
        for tt in topo:
            for j, n in enumerate(nodes):
                node_str = " -nodes " + str(nodes_allocated[j]) + " "
                for b in bandwidth:
                    for i, d in enumerate(links):
                        l = bandwidth_link_weight[t_index]
                        ifile = 'tg.' + topo_file_prefix[t_index] + "." + str(int(b/1000)) + 'G.' + a + str(n) + '_large.nd' + str(d)
#                        print("ifile(%s), 'tg.', topo_file_prefix(%s), a(%s); tt(%s)" % (ifile, topo_file_prefix[t_index], a, tt))
                        option_ifile_deg = " "
#                        if (option_link_deg_ifile[t_index] != []):
#                            option_ifile_deg = ".nd" + str(option_link_deg_ifile[t_index][i])
                        ifile = " -flowfile " + ifile + option_ifile_deg + " "
                        ofile = " -ofile " + dt_string + "performance." + tt + "." + a + str(n) + ".speed" + str(int(b/1000)) + ".nd" + str(d)
                        actual_bandwidth = b
                        if (l != []):
                            actual_bandwidth = b * l[i]
                        bandwidth_str = " -speed " + str(int(actual_bandwidth)) + " "
                        option_deg = " "
                        if (option_link_deg[t_index] != []):
                            option_deg = " -deg " + str(option_link_deg[t_index][i]) + " "
                        full_cmd = topo_htsim[t_index] + topo_sim_args[t_index] + ifile + ofile + bandwidth_str + node_str + option_deg
                        print("%s" % full_cmd) 
