# python code to read a task graph in flatbuffers format and replace its topology with rings for INTREPID network architecture
# namespace: FlatBufTaskGraph
import flatbuffers as fb
import FlatBufTaskGraph.TaskGraph as tg
import FlatBufTaskGraph.Connection as conn
import FlatBufTaskGraph.Device as device
import FlatBufTaskGraph.Operator as operator
import FlatBufTaskGraph.Route as route
import FlatBufTaskGraph.Path as path
import FlatBufTaskGraph.Task as task
import FlatBufTaskGraph.Topology as tp
import FlatBufTaskGraph.Rings as rings
import FlatBufTaskGraph.RingDescriptor as ringdescriptor

import copy
import sys
import argparse

builder = fb.Builder(0)


g_task = []
g_conn = []
g_route = []
g_device = []

def make_vector(_builder, elm_l):
    for e in elm_l:
        _builder.PrependUOffsetTRelative(e)
    return builder.EndVector()

def add_conns(_builder, _conn):
    conn_l = []
    for c in _conn:
        conn.Start(_builder)
        conn.AddFromnode(_builder, c[0])
        conn.AddTonode(_builder, c[1])
        conn.AddNconn(_builder, c[2])
        conn_i = conn.End(_builder)
        conn_l.append(conn_i)
    tg.StartConnVector(_builder, len(conn_l))
    conn_v = make_vector(_builder, conn_l)
    return conn_v

def add_hopnode(_builder, _hopnode):
#    _hopnode.reverse()
    path.StartHopnodeVector(_builder, len(_hopnode))
    for p in _hopnode:
        _builder.PrependUint32(p)
    hopnode = _builder.EndVector()
    return hopnode 

def add_paths(_builder, _path):
    path_l = []
    #print("_path = %s" % str(_path))
    for p in _path:  # path = [hopnode, chance]
        #print("p[0] = %s" % str(p[0]))
        hopnode_v = add_hopnode(_builder, p[0])
        # add chance
        path.Start(_builder)
        path.AddHopnode(_builder, hopnode_v)
        path.AddChance(_builder, p[1])
        path_i = path.End(_builder)
        path_l.append(path_i)
    route.StartPathsVector(_builder, len(path_l))
    path_v = make_vector(_builder, path_l)
    return path_v

def add_routes(_builder, _route):
    route_l = []
    for r in _route:
        #print(r)
        path_v = add_paths(_builder, r[2])
        route.Start(_builder)
        route.AddFromnode(_builder, r[0])
        route.AddTonode(_builder, r[1])
        route.AddPaths(_builder, path_v)
        route_i = route.End(_builder)
        route_l.append(route_i)
    tg.StartRoutesVector(_builder, len(route_l))
    route_v = make_vector(_builder, route_l)
    return route_v

def pr2dlist(title, plist):
    print("%s", title)
    print(plist)
    for p in plist:
        print(p)

def get_send_dest_id(src, dest):

    devproperties = int(src * Ntotalnodes + dest)
    return devproperties
    
def create_compute_task(taskid, nodeid, compute, nexttasks):
    task_i = []
    task_i.append(0)    # type: TASK_FORWARD
    task_i.append(taskid)    # task_id:
    task_i.append(nodeid)    # device_id: gpu id 
    task_i.append(compute)  # runtime: 0.1 sec
    task_i.append(1)  # xfersize : just a nonzero value
    task_i.append(nexttasks)
    return task_i

def create_communication_task(taskid, src, dest, dev_id, communication, nexttasks):
    task_i = []
    task_i.append(5)    # type: TASK_NOMINALCOMM  # TASK_COMM is not supported
    task_i.append(taskid)    # task_id:
    task_i.append(dev_id)
    task_i.append(0)  # runtime: 0.1 sec
    task_i.append(communication)  # runtime: 0.1 sec
    task_i.append(nexttasks)    # nexttasks
    return task_i

# Rule for numbering devices and tasks
#   switches: 0 - Nswitches
#   gpu/cpu : Nswitchs - (Nswitches + Nnodes - 1)
# 8 tasks, task graph having ring comm. pattern
def generate_ring1_task():
    global g_task
    global g_device
    global Nnodes

    g_task = []
    dev_l = []
    dev_id = Ntotalnodes

    task_id = 0
    for i in range(Nnodes):
        task_i = []
        j = i + 1
        if (j == Nnodes): j = 0
        src = i
        dest = j
        nexttasks = []
        nexttasks.append(i*2+1)
        task_i = create_compute_task(i*2, src, 0.1, nexttasks)
        g_task.append(task_i)
        dev_id = dev_id + 1
        nexttasks = []
        # nexttasks.append(999)   # empty list is not acceptable. trial. Not working!
        #nexttasks.append(14)   # empty list is not acceptable. the last task. 
        nexttasks.append(i*2+2)   # empty list is not acceptable. next compute task
        dev_l.append([dev_id, src, dest])
        task_i = create_communication_task(i*2+1, src, dest, dev_id, 1000000000+1, nexttasks)
        g_task.append(task_i)
        task_id = i*2+2
    nexttasks = [task_id]
    task_i = create_compute_task(task_id, 0, 0.1, nexttasks)
    g_task.append(task_i)

    #print("dev_l = %s" % str(dev_l))
    #print("g_task = %s" % str(g_task))

    for i in range(Nnodes): # GPU
        node_id = i 
        dev_i = [0, node_id, node_id,  node_id, 0]
        g_device.append(dev_i)

    for d in dev_l: # COMM
        device_properties = get_send_dest_id(d[1], d[2])
        dev_i = [13, d[0], device_properties, device_properties, 1000]
        g_device.append(dev_i)
    #print("g_device = %s" % str(g_device))

def normalize_nodeid(ref_node, offset, in_group):
    new_node = ref_node + offset
    if (in_group):
        ref_gid = int(ref_node / Nodespergroup)
        if (offset < 0):
            if (new_node < ref_gid * Nodespergroup): new_node = new_node + Nodespergroup
        else:
            if (new_node >= (ref_gid + 1) * Nodespergroup): new_node = new_node - Nodespergroup
        assert(new_node >= (ref_gid) * Nodespergroup and new_node < (ref_gid +1) * Nodespergroup)
    else:
        if (new_node < 0): new_node = new_node + Nnodes
        if (new_node >= Nnodes): new_node = new_node - Nnodes
        assert(new_node >= 0 and new_node < Nnodes)

    print("normalize_nodeid: ref_node(%d), offset(%d), in_group(%s), return(%d)" % (ref_node, offset, str(in_group), new_node))

    return new_node

def normalized_total (node_id):
    if (node_id < 0): node_id = node_id + Nnodes
    if (node_id >= Nnodes): node_id = node_id - Nnodes
    return node_id

def normalized_diff(diff, in_group):
    if diff >= 0: return (False, diff)
    if (in_group): diff = diff + Nodespergroup
    else: diff = diff + Nnodes
    return (True, diff)

def same_group(n1, n2):
    if (int(n1 / Nodespergroup) == int(n2 / Nodespergroup)): return True
    return False

def same_ring(n1, n2):
    print("same_ring: %d offset(%d), %d offset(%d)" % (n1, n1 % Nodespergroup, n2, n2 % Nodespergroup))
    if (n1 % Nodespergroup == n2 % Nodespergroup): return True
    return False

def get_dist(src, dest, in_group):
    nnodes = Nnodes
    if (in_group): nnodes = Nodespergroup
    half_ring = int(nnodes / 2)

    raw_diff = dest - src
    p_diff = raw_diff
    r_diff = nnodes - raw_diff
    reverse_direction = False
    if (raw_diff <= 0): 
        reverse_direction = True
        p_diff = raw_diff + nnodes
        r_diff = - raw_diff
    return (raw_diff, p_diff, r_diff)

def get_route_ring(src, dest, d_list, in_group):
    nnodes = Nnodes
    if (in_group): nnodes = Nodespergroup
    half_ring = int(nnodes / 2)

    (raw_diff, p_diff, r_diff) = get_dist(src, dest, in_group)

    print("get_route_ring: src(%d), dest(%d), in_group(%s), d_list = %s, raw_diff(%d), p_diff(%d), r_diff(%d)" % (src, dest, str(in_group), str(d_list), raw_diff, p_diff, r_diff))

    t_list = d_list.copy()
    gh = []
    loop_count = len(t_list)
    while (loop_count != 0):
        t_src = src
        h = [src]
        t_src = src
        t_dest = dest
        print("start %d " % (t_src))
        for s in t_list:
            o_s = s
            (raw_diff, p_diff, r_diff) = get_dist(t_src, t_dest, in_group)
            if s > half_ring:	# reverse link
                diff = r_diff
                s = nnodes - s
            else: diff = p_diff
            print("inside for: s(%d), diff(%d), src(%d), dest(%d), in_group(%s), d_list = %s, raw_diff(%d), p_diff(%d), r_diff(%d)" % (s, diff, src, dest, str(in_group), str(d_list), raw_diff, p_diff, r_diff))
            while diff >= s:
                if (in_group):
                    t_src = normalize_nodeid(t_src, o_s, True)
                else:
                    t_src = normalize_nodeid(t_src, o_s, False)
                diff = diff - s
                print("next %d, diff (%d)" % (t_src, diff))
                h.append(t_src)
        assert(t_src == dest)
        loop_count = loop_count - 1
        gh.append(h)
        del t_list[0]
        t_list.append(s)

    # get the shortest one
    min_len = 99999999
    for i, l in enumerate(gh):
        print("candidate %d: %s" % (i, str(l)))
        if (len(l) < min_len):
            min_len = len(l)
            best_id = i
    
    print("best : %s" % str(gh[best_id]))
    return gh[best_id]

def get_droute(src, dest):

    print("dist_l = %s" % str(dist_l)) 
    print("dist_g = %s" % str(dist_g)) 
    diff = normalized_total(dest - src)
    t_id = src
    if (Nodespergroup == 1): # single ring, only arg_topo[1] matters
        h = get_route_ring(src, dest, dist_l, True)
        return h
    else:	# multiple level ring, assuming that each ring connects all components
        src_g_id = int(src / Nodespergroup)
        dest_g_id = int(dest / Nodespergroup)
        src_offset = src % Nodespergroup	# outside ring-id
        dest_offset = dest % Nodespergroup	# outside ring-id
        diff_offset = dest_offset - src_offset	
        if (same_group(src, dest)): # on the same group
            h = get_route_ring(src, dest, dist_l, True)
            print("h: same group: %d --> %d: %s" % (src, dest, str(h)))
            return h
        elif (same_ring(src, dest)): # check 2 different paths
            src_t_id = normalize_nodeid(src, 1, True)
            dest_t_id = normalize_nodeid(dest, 1, True)
            other_offset = src_t_id % Nodespergroup
            # single ring route
            h1 = get_route_ring(src, dest, dist_g[src_offset], False)
            # double ring route 
            h2 = get_route_ring(src_t_id, dest_t_id, dist_g[other_offset], False)
            h2.append(dest)
            h2.insert(0, src)
            if (len(h2) < len(h1)): return h2
            return h1
        else:
            src_t_id = normalize_nodeid(src, diff_offset, True)
            dest_t_id = normalize_nodeid(dest, - diff_offset, True)
            print("src(%d), dest(%d), src_t_id(%d), dest_t_id(%d)" % (src, dest, src_t_id, dest_t_id))
            assert(same_group(src, src_t_id))
            assert(same_group(dest_t_id, dest))
            assert(same_ring(src, dest_t_id))
            assert(same_ring(src_t_id, dest))
            # 1. route src --> src_t_id --> dest
            h1 = get_route_ring(src_t_id, dest, dist_g[dest_offset], False)
            if (src != src_t_id): h1.insert(0, src)
            print("h1: %d --> %d: %s" % (src_t_id, dest, str(h1)))
            # 2. route src --> dest_t_id --> dest
            h2 = get_route_ring(src, dest_t_id, dist_g[src_offset], False)
            if (dest != dest_t_id): h2.append(dest)
            print("h2: %d --> %d: %s" % (src, dest_t_id, str(h2)))
            if (len(h1) < len(h2)):
                h = h1
            else:
                h = h2
            print("return: %d --> %d: %s" % (src, dest, str(h)))
            return h

def get_route(src, dest, distance_l):
    if (src == dest): return []
    (src_nodeid, src_switchid) = get_rel_ids(src)
    (dest_nodeid, dest_switchid) = get_rel_ids(dest)
    h = [src]

    # node is part of the switch (HPPN)
    if (src_switchid < 0 or dest_switchid < 0):
        assert(dest_switchid < 0 and src_switchid < 0)
        src_switchid = src
        dest_switchid = dest

    # node and the switch is separated
    if (src_switchid == dest_switchid):
        h = [src, src_switchid + Nnodes, dest]
        return h

    nnodes = Nswitches
    if (Nswitches == 0): nnodes = Nnodes

    distance_l.sort(reverse = True)
    diff = dest_switchid - src_switchid
    if (Nswitches > 0): diff = diff + Nswitches
    elif (diff < 0): diff = diff + Nnodes 

    g_switchid = get_switch_gid(src_switchid)
    t_switchid = src_switchid
    if (src != g_switchid):
        h.append(g_switchid)
    while diff > 0:
        for d in distance_l:
            while diff >= d:
                t_switchid = (t_switchid + d) % nnodes
                #print("diff = %d, d = %d, t_switchid = %d" % (diff, d, t_switchid))
                h.append(get_switch_gid(t_switchid))
                diff = diff - d
    #print("src %d, dest %d, src_switchid %d, dest_switchid %d, t_switchid %d, dest_switchid %d: h = %s" % (src, dest, src_switchid, dest_switchid, t_switchid, dest_switchid, str(h)))
    assert (t_switchid == dest_switchid)
    if (dest != dest_switchid):
        h.append(dest)
    return h

# return global (nodeid, switchid) of gpu node
def get_ids(i):
    assert(i < Nnodes)
    nodeid = i
    switchid = int(i/NodeperSwitch) + Nnodes
    if (Nswitches == 0): switchid = -1
    return (nodeid, switchid)

# return (nodeid, switchid) of gpu node
def get_rel_ids(i):
    assert(i < Nnodes)
    nodeid = i
    switchid = int(i/NodeperSwitch)
    if (Nswitches == 0): switchid = -1
    return (nodeid, switchid)

# get global id of gpu node 
def get_node_gid(i):
    assert(i < Nnodes)
    return i
    
# get global id of switch from its relative id
def get_switch_gid(i):
    if (Nswitches == 0): return i
    assert(i < Nswitches)
    return i + Nnodes


def generate_dring(arg_topo):
    global g_conn
    global g_route
    global dist_g
    global dist_l

    g_conn = []
    g_route = []
    dist_l = []
    dist_g = []

    print("arg_topo[1] = %s", str(arg_topo[1]))
    print("arg_topo[2] = %s", str(arg_topo[2]))
    # assume that all the nodes in a group is connected to inter-group ring
    assert(len(arg_topo[1]) == Nodespergroup)	
    assert(Ngroups * Nodespergroup == Nnodes)
 
    # node_id numbering convention
    # i-th node in the j-th group = j * Nodespergroup + i
    # group internal connection
    intra_group_conn = arg_topo[2]
    for t in intra_group_conn:
        if (t[0] < 0): t[0] = t[0] + Nodespergroup

    # inter-group connection
    inter_group_conn = arg_topo[1]
    for g in inter_group_conn:	# for each group
        for t in g:	# actual connection info
            if (t[0] < 0): t[0] = t[0] + Ngroups

    for l in inter_group_conn:
        l.sort(key=lambda x: x[0], reverse=True) # sort according to the distance
        dist_g.append([i[0] * Nodespergroup for i in l])

    intra_group_conn.sort(key=lambda x: x[0], reverse=True) # sort according to the distance
    dist_l = [i[0] for i in intra_group_conn]

    # intra-group
    for l in intra_group_conn:
        for g in range(Ngroups):
            for t in range(Nodespergroup):
                srcid = g * Nodespergroup + l[0] + t
                diff = l[0]
                assert(l[0] > 0)
                destid = normalize_nodeid(srcid, diff, True)
                g_conn.append([srcid, destid, l[1]])

    # inter-group
    for out_ring_id, lg in enumerate(inter_group_conn): # for each outside ring
        s_gid_start = out_ring_id
        print("-- connection of ring %d --" % (out_ring_id))
        for l in lg:	# for each connection within the ring
            for j in range(Ngroups):
                s_gid = s_gid_start + j * Nodespergroup
                d_gid = normalized_total(s_gid + l[0]*Nodespergroup)
                g_conn.append([s_gid, d_gid, l[1]])
                print(g_conn[-1])

    # route
    for i in range(Nnodes):
        for j in range(Nnodes):
            if (i == j): continue
            h = get_droute(i, j)
            pv = [[h, 1.0]]   # single route is assumed for now
            r = [i, j, pv] # from_node, to_node, path
            g_route.append(r)
    #pr2dlist("route matrix", g_route)

# 8 switches, ring 1
# 1 node/switch
# total numbers of devices = 16
# encoding of node ID:
#   switch ID: (switch id) * (node/switch)
#   node ID :  (switch id) + (node id)
# connection from a node to a switch:
# connection from a switch to a switch:
# arg_topo: list of [distance, no of connections]
def generate_ring(arg_topo):
    global g_conn
    global g_route

    g_conn = []
    g_route = []

    for t in arg_topo:
        assert(t[0] > 0)
    # ids of gpu nodes are always from 0 to (ngpus - 1)
    # then, let's design like this
    # gpu i is connected to switch i/(gpus per node) + (ngpus)

    # switch to a node (independent of arg_topo)
    if (Nswitches != 0):
        for i in range(Nnodes):
            (node_id, switch_id) = get_ids(i)
            g_conn.append([switch_id, node_id, 10])
            g_conn.append([node_id, switch_id, 10])

    arg_topo.sort(key=lambda x: x[0]) # sort according to the distance
    distances = [i[0] for i in arg_topo]
    #print("distances = %s" % str(distances))
    # switch to switch
    if (Nswitches == 0): nnodes = Nnodes;
    for _t in arg_topo:
        first_switch = Nnodes
        if (Nswitches == 0): first_switch = 0
        distance = _t[0]
        nconn = _t[1]
        for i in range(nnodes):
            j = i + distance
            if (j < 0): j = j + nnodes 
            if (j >= nnodes): j = j - nnodes
            g_conn.append([i+first_switch, j+first_switch, nconn])
        #pr2dlist("connection matrix", g_conn)

    # route
    for i in range(Nnodes):
        for j in range(Nnodes):
#            if (i == j): 
#                continue
            h = get_route(i, j, distances)
#            h.reverse()
            pv = [[h, 1.0]]   # single route is assumed for now
            r = [i, j, pv] # from_node, to_node, path
            g_route.append(r)
    #pr2dlist("route matrix", g_route)

def add_ringdescriptorjumps(_builder, _rd):
    if (_rd == []): return None
    rd_l = []
#    _rd.reverse()
    ringdescriptor.StartJumpsVector(_builder, len(_rd))
    print("_rd = %s" % str(_rd))
    for r in _rd:
        _builder.PrependInt32(r)
    jumps = _builder.EndVector()
    return jumps

def add_ringpaths(_builder, _rd):
    if (_rd == []): return None
    rd_l = []
    for r in _rd:
        rd_v = add_ringdescriptorjumps(_builder, r)
        rd_l.append(rd_v)
    rings.StartRingpathsVector(_builder, len(_rd))
    paths_v = make_vector(_builder, rd_l)
    return paths_v

def add_rings(_builder, _rings):
    print("rings = %s" % str(_rings))
    if (_rings == []): return None
    rings_l = []
    for r in _rings:
        rd_v = add_ringpaths(_builder, r[1])
        rings.Start(_builder)
        rings.AddRingsz(_builder, r[0])
        rings.AddRingpaths(_builder, rd_v)
        r_i = rings.End(_builder)
        rings_l.append(r_i)
    tg.StartRingsVector(_builder, len(rings_l))
    rings_v = make_vector(_builder, rings_l)
    return rings_v

def add_ops(_builder, _ops):
    if (_ops == []): return None
    o_l = []
    for o in _ops:
        print("operator = %s" % o)
        str_name = _builder.CreateString(o[2])
        operator.Start(_builder)
        operator.AddOpid(_builder, o[0])
        operator.AddType(_builder, o[1])
        operator.AddName(_builder, str_name)
        o_i = operator.End(_builder)
        o_l.append(o_i)
    tg.StartOpsVector(_builder, len(o_l))
    o_v = make_vector(_builder, o_l)
    return o_v

def add_devices(_builder, _devices):
    #print("---devices = %s" % str(_devices))
    device_l = []
    for d in _devices:
        device.Start(_builder)
        device.AddType(_builder, d[0])
        device.AddDeviceid(_builder, d[1])
        device.AddNodeid(_builder, d[2])
        device.AddDeviceproperty(_builder, d[3])
        device.AddBandwidth(_builder, d[4])
        device_i = device.End(_builder)
        device_l.append(device_i)
    tg.StartDevicesVector(_builder, len(device_l))
    device_v = make_vector(_builder, device_l)
    return device_v

def add_nexttasks(_builder, _nt):
    #print("nexttasks = %s" % str(_nt))
    task.StartNexttasksVector(_builder, len(_nt))
    for n in _nt:
        #print("add %s" % str(n))
        _builder.PrependUint64(n)
    nt_v = _builder.EndVector()
    return nt_v

def add_tasks(_builder, _tasks):
    task_l = []
    for t in _tasks:
        #print("t = %s" % str(t))
        nexttask_v = []
        if (len(t[5]) > 0):
            nexttask_v = add_nexttasks(_builder, t[5])
        task.Start(_builder)
        task.AddType(_builder, t[0])
        task.AddTaskid(_builder, t[1])
        task.AddDeviceid(_builder, t[2])
        task.AddRuntime(_builder, t[3])
        task.AddXfersize(_builder, t[4])
        if (len(t[5]) > 0):
            task.AddNexttasks(_builder, nexttask_v)
        task_i = task.End(_builder)
        task_l.append(task_i)
    tg.StartTasksVector(_builder, len(task_l))
    task_v = make_vector(_builder, task_l)
    return task_v

#build taskgraph from vectors
def build_taskgraph(builder, ngpupernode, Nnodes, Nswitches, gpubw, drambw, netbw, conn_v, ops_v, tasks_v, devices_v, routes_v, rings_v, fname):
    tg.Start(builder)
    tg.AddNgpupernode(builder, ngpupernode)
    tg.AddNnode(builder, Nnodes)
    tg.AddNswitch(builder, Nswitches)
    tg.AddIntergpubw(builder, gpubw) 
    tg.AddDrambw(builder, drambw)
    tg.AddNetbw(builder, netbw)
    tg.AddConn(builder, conn_v)
    tg.AddOps(builder, ops_v)
    tg.AddTasks(builder, tasks_v)
    tg.AddDevices(builder, devices_v)
    tg.AddRoutes(builder, routes_v)
    if (rings_v != None):
        tg.AddRings(builder, rings_v)
    tg1 = tg.End(builder)
    builder.Finish(tg1)
    buffer = builder.Output()
    #print("file write")
    fw = open(fname, 'wb')
    fw.write(buffer)
    fw.close()

# build task graph from a list
def write_flatbuffer_tg(builder, tg_l, fname):
    
    for i in range(len(tg_l)):
        if (isinstance(tg_l[i], list)):
            tg_l[i].reverse()
            if (i == 8): # Tasks
                for t in tg_l[i]:
                    t[5].reverse()
            elif (i == 10): # Routes
                for r in tg_l[i]: # route
                    for p in r[2]: # path
                        p[0].reverse()
                    r[2].reverse()
                
    conn_v = add_conns(builder, tg_l[6])
    ops_v = add_ops(builder, tg_l[7])
    tasks_v = add_tasks(builder, tg_l[8])
    devices_v = add_devices(builder, tg_l[9])
    routes_v = add_routes(builder, tg_l[10])
    rings_v = add_rings(builder, tg_l[11])
    build_taskgraph(builder, tg_l[0], tg_l[1], tg_l[2], tg_l[3], tg_l[4], tg_l[5], conn_v, ops_v, tasks_v, devices_v, routes_v, rings_v, fname)

#build topology
def build_topology(builder, conn_v, route_v, fname):
    tp.Start(builder)
    tp.AddConn(builder, conn_v)
    tp.AddRoutes(builder, route_v)
    tp1 = tp.End(builder)
    builder.Finish(tp1)
    buffer = builder.Output()
    #print("file write")
    fw = open(fname, 'wb')
    fw.write(buffer)
    fw.close()

def read_conn(mytg):
    myconn_length = mytg.ConnLength()
#    print("myconn_length = %d" % myconn_length)
    conn_l = []
    for i in range(myconn_length):
        myconn = mytg.Conn(i)
        fnode = myconn.Fromnode()
        tonode = myconn.Tonode()
        nconn = myconn.Nconn()
        print("fnode %d tonode %d nconn %d" % (fnode, tonode, nconn))
        conn_l.append([fnode, tonode, nconn])
    return conn_l

def read_routes(mytg):
    myroutes_length = mytg.RoutesLength()
#    print("myroutes_length = %d" % myroutes_length)
    routes_l = []
    for i in range(myroutes_length):
        myroutes = mytg.Routes(i)
        fnode = myroutes.Fromnode()
        tonode = myroutes.Tonode()
        mypath_length = myroutes.PathsLength()
#        print("from(%d) to (%d)" % (fnode, tonode))
        p_l = []
        for j in range(mypath_length):
            mypath = myroutes.Paths(j)
            hopnode = mypath.HopnodeAsNumpy()
            h_l = []
            if (not mypath.HopnodeIsNone()):
                h_l = hopnode.tolist()
            chance = mypath.Chance()
#            print("\thopnode = %s : prob %f" % (str(hopnode), chance))
            p_l.append([h_l, chance])
        routes_l.append([fnode, tonode, p_l])
    return routes_l

def read_ops(mytg):
    myops_length = mytg.OpsLength()
    l = []
    for i in range(myops_length):
        myops = mytg.Ops(i)
        opid = myops.Opid()
        mytype = myops.Type()
        name = myops.Name()
        l.append([opid, mytype, name])
#        print("\tops = %s " % (str(l)))
    return l

def read_tasks(mytg):
    mytasks_length = mytg.TasksLength()
    tasks_l = []
    for i in range(mytasks_length):
        l = []
        mytask = mytg.Tasks(i)
        ttype = mytask.Type()
        taskid = mytask.Taskid()
        deviceid = mytask.Deviceid()
        runtime = mytask.Runtime()
        xfersize = mytask.Xfersize()
        nexttasks = mytask.NexttasksAsNumpy()
        nexttasks_l = []
        if (not mytask.NexttasksIsNone()): 
            nexttasks_l = nexttasks.tolist()
        tasks_l.append([ttype, taskid, deviceid, runtime, xfersize, nexttasks_l])
#        print("type %d, taskid %d, deviceid %d, runtime %d, xfersize %d" % (ttype, taskid, deviceid, runtime, xfersize))
#        print("next tasks = %s" % str(nexttasks))
    return tasks_l

def read_devices(mytg):
    mydevices_length = mytg.DevicesLength()
    devices_l = []
    for i in range(mydevices_length):
        mydevice = mytg.Devices(i)
        dtype = mydevice.Type()
        deviceid = mydevice.Deviceid()
        nodeid = mydevice.Nodeid()
        deviceproperty = mydevice.Deviceproperty()
        bandwidth = mydevice.Bandwidth()
        devices_l.append([dtype, deviceid, nodeid, deviceproperty, bandwidth])
#        print("type %d, deviceid %d, node id %d, device property %d, bandwitdh %d" % (dtype, deviceid, nodeid, deviceproperty, bandwidth))
    return devices_l

def read_rings(mytg):
    myrings_length = mytg.RingsLength()
    rings_l = []
    for i in range(myrings_length):	# for each ring
        myring = mytg.Rings(i)
        ringsz = myring.Ringsz()
        myringpaths_len = myring.RingpathsLength()
        l = []
        for j in range(myringpaths_len): # for each ringdescriptor
            myringdescriptor = myring.Ringpaths(j)
            if (myringdescriptor.JumpsIsNone()): 
                print("jumpis None")
                continue
            print("jumps length = %d" % myringdescriptor.JumpsLength())
            jumps = myringdescriptor.JumpsAsNumpy()
            print("jumps = %s" % str(jumps))
            l.append(jumps.tolist())
        rings_l.append([ringsz, l])
#    print("Rings = %s" % str(rings_l))
    return rings_l

def print_list(name, l):
    print("%s = " % name)
    for t in l:
        print(t)

def do_diff(tg1, tg2):
    for i in range(len(tg1)):
        if tg1[i] != tg2[i]:
            print("Different at %d-th element" % i)
            print_list("original", tg1[i])
            print_list("copy", tg2[i])

# read taskgraph and returns the list of the taskgraph
def read_taskgraph(fname):
    buf = open(fname, 'rb').read()
    buf = bytearray(buf)
    mytg = tg.TaskGraph.GetRootAs(buf, 0)
    ngpupernode = mytg.Ngpupernode()
    nnode = mytg.Nnode()
    nswitch = mytg.Nswitch()
    intergpubw = mytg.Intergpubw()
    drambw = mytg.Drambw()
    netbw = mytg.Netbw()
    print("ngpupernode(%d), nnode(%d), nswitch(%d), ingergpubw(%f), drambw(%f), netbw(%f)" % (ngpupernode, nnode, nswitch, intergpubw, drambw, netbw))
    conn_l = read_conn(mytg)
    print_list("connection " , conn_l)
    ops_l = read_ops(mytg)
    print_list("ops " , ops_l)
    tasks_l = read_tasks(mytg)
    print_list("tasks " , tasks_l)
    devices_l = read_devices(mytg)
    print_list("devices " , devices_l)
    routes_l = read_routes(mytg)
    print_list("routes " , routes_l)
    rings_l = read_rings(mytg)
    print_list("rings " , rings_l)
    tg_l = [ngpupernode, nnode, nswitch, intergpubw, drambw, netbw, conn_l, ops_l, tasks_l, devices_l, routes_l, rings_l]
    return tg_l

def usage():
    print("Usage: python3 tgconverter.py <input flatbuffer task graph> <output flatbuffer task graph> <Number of nodes> <Number of swtiches> <Ring type {1, 2, 3, 11}>")
    print("     : Ring type 1 -3 is for single ring, 11 or higher is for double ring")
    sys.exit()

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('-i',  '--input', type=str, help='input taskgraph file in flabuffers', required=True)
parser.add_argument('-o',  '--output', type=str, help='output taskgraph file in flabuffers', required=True)
parser.add_argument('-N',  '--Nodes', type=int, help='number of nodes in the taskgraph file in flabuffers', required=True)
parser.add_argument('-S',  '--Switches', type=int, help='number of switches in the taskgraph file in flabuffers', required=True)
parser.add_argument('-T',  '--Type', type=int, help='Ring type of the new taskgraph file in flabuffers', required=True)

# describe topology
# num_nodes, connection
Ngpupernode = 2
Nnodes = 16 
Nodespergroup = 1
Ngroups = 1
Nswitches = 8
Gpubw = 0.0
Drambw = 0.0
Netbw = 1000.0
output_fname = 'ringtopo1.bin'

args = vars(parser.parse_args())
print(args)
Nnodes = args['Nodes']
Nswitches = args['Switches']
NodeperSwitch = 1
if (Nswitches != 0):	# Node itself is a switch
    NodeperSwitch = int(Nnodes / Nswitches)
Ntotalnodes = Nnodes + Nswitches
input_fname = args['input']
output_fname = args['output']
ring_type = args['Type']

assert(Ntotalnodes == Nnodes + Nswitches)

#generate_ring1_task()
ring_topo1 = [[1, 1], [2, 1], [3, 1], [4, 1], [5, 1], [6, 1], [7, 1]] # Ring 1
ring_topo2 = [[1, 8]] # Ring 2
ring_topo3 = [[1, 7], [2,1]] # Ring 3
ring_topo = [ring_topo1, ring_topo2, ring_topo3]

# double ring topology description
# [number of nodes/group, [connection between the group from different nodes in the group]], [connection insdie of the group], 
dring_topo1 = [1, [ring_topo1], []]
dring_topo2 = [1, [ring_topo2], []]
dring_topo3 = [1, [ring_topo3], []]
dring_topo11 = [2, [[[1, 8]], [[-1, 8]]], [[1, 8]]]
num_single_topologies = 3
dring_topo = [dring_topo1, dring_topo2, dring_topo3, dring_topo11]
dist_g = []
dist_l = []
if (ring_type > 10): # dring
     topo = dring_topo[ring_type - 10 -1 + num_single_topologies]
     Nodespergroup = topo[0]
     Ngroups = int(Nnodes / Nodespergroup)
     generate_dring(topo) 
else:
     generate_ring(ring_topo[ring_type - 1])
conn_v = add_conns(builder, g_conn)
route_v = add_routes(builder, g_route)
#task_v = add_tasks(builder, g_task)
#device_v = add_devices(builder, g_device)

#build_taskgraph(builder, Ngpupernode, Nnodes, Nswitches, Gpubw, Drambw, Netbw, conn_v, route_v, task_v, device_v, 'test-tg.bin')
#build_topology(builder, conn_v, route_v, 'test-tp.bin')

tg_l = read_taskgraph(input_fname)

#sys.exit()

tg_l_save = copy.deepcopy(tg_l)
#do_diff(tg_l, tg_l_save)
#print_list("old conn", tg_l[6])
#print_list("old route", tg_l[10])
tg_l[6] = g_conn
tg_l[10] = g_route

write_flatbuffer_tg(builder, tg_l, output_fname)
tg_l2 = read_taskgraph(output_fname)
print_list("tg_l after write_flatbuffer_tg", tg_l[6])
print_list("tg_l_save again", tg_l_save[6])
do_diff(tg_l_save, tg_l2)

#print_list("new conn", tg_l[6])
#print_list("new route", tg_l[10])
