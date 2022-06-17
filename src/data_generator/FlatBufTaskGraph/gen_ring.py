
import flatbuffers as fb
import TaskGraph as tg
import Connection as conn
import Device as device
import Route as route
import Path as path
import Task as task
import Topology as tp

builder = fb.Builder(0)

# describe topology
# num_nodes, connection
Nnodes = 8
Nswitches = 8
NodeperSwitch = int(Nnodes / Nswitches)
NodesperCluster = (NodeperSwitch + 1)
Ntotalnodes = Nnodes + Nswitches
assert(Ntotalnodes == Nnodes + Nswitches)

g_conn = []
g_route = []
g_device = []

# conn
# two dimensional array conn[i][j] = n, if there are 'n' links

#conn.Start(builder)
#conn.AddFromnode(builder, 0)
#conn.AddTonode(builder, 1)
#conn.AddNconn(builder, 8)
#conn1 = conn.End(builder)

#tg.StartConnVector(builder,1)
#builder.PrependUOffsetTRelative(conn1)
#conn_v1 = builder.EndVector(1)

def make_vector(_builder, elm_l):
    for e in elm_l:
        _builder.PrependUOffsetTRelative(e)
    return builder.EndVector(len(elm_l))

# route
#path.StartHopnodeVector(builder, 2)
#builder.PrependUint32(1)
#builder.PrependUint32(2)
#hopnode1 = builder.EndVector(2)

#path.Start(builder)
#path.AddHopnode(builder, hopnode1)
#path.AddChance(builder, 1.0)
#path1 = path.End(builder)

#route.StartPathsVector(builder, 1)
#builder.PrependUOffsetTRelative(path1)
#route_path = builder.EndVector(1)

#route.Start(builder)
#route.AddFromnode(builder, 0)
#route.AddTonode(builder, 2)
#route.AddPaths(builder, route_path)
#route1 = route.End(builder)

#tg.StartRoutesVector(builder, 1)
#builder.PrependUOffsetTRelative(route1)
#route_vector = builder.EndVector(1)




#tg.AddConn(builder, conn_v1)
#tg.AddRoutes(builder, route_vector)
#tg1 = tg.End(builder)

# ops

# tasks

# devices

# routes

# rings


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
    path.StartHopnodeVector(_builder, len(_hopnode))
    for p in _hopnode:
        _builder.PrependUint32(p)
    hopnode = _builder.EndVector(len(_hopnode))
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
    global Nnodes
    global Nswitches
    global NodeperSwitch
    global NodesperCluster

    devproperties = int(src * Ntotalnodes + dest)
    return devproperties
    
def get_nodeid(i):
    global Nnodes
    global Nswitches
    global NodeperSwitch
    global NodesperCluster

    node_offset = int(i % NodeperSwitch)
    switch_id = int(i / NodeperSwitch)
    node_id = int(switch_id * NodesperCluster + node_offset)
    return node_id

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
        switch_id = i * NodesperCluster
        node_id = switch_id + 1
        j = i + 1
        if (j == Nnodes): j = 0
        src = i
        dest = j
        #src_id = get_nodeid(i)
        #dest_id = get_nodeid(i)
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

def get_route(src, dest, distance_l):
    (src_nodeid, src_switchid) = get_rel_ids(src)
    (dest_nodeid, dest_switchid) = get_rel_ids(dest)
    h = [src]
   
    distance_l.sort(reverse = True)
    diff = dest_switchid - src_switchid
    if diff < 0:
        diff = diff + Nswitches 
    g_switchid = get_switch_gid(src_switchid)
    t_switchid = src_switchid
    h.append(g_switchid)
    while diff > 0:
        for d in distance_l:
            while diff >= d:
                t_switchid = (t_switchid + d) % Nswitches
                #print("diff = %d, d = %d, t_switchid = %d" % (diff, d, t_switchid))
                h.append(get_switch_gid(t_switchid))
                diff = diff - d
    #print("src %d, dest %d, src_switchid %d, dest_switchid %d, t_switchid %d, dest_switchid %d: h = %s" % (src, dest, src_switchid, dest_switchid, t_switchid, dest_switchid, str(h)))
    assert (t_switchid == dest_switchid)
    h.append(dest)
    return h

# return global (nodeid, switchid) of gpu node
def get_ids(i):
    assert(i < Nnodes)
    nodeid = i
    switchid = int(i/NodeperSwitch) + Nnodes
    return (nodeid, switchid)

# return (nodeid, switchid) of gpu node
def get_rel_ids(i):
    assert(i < Nnodes)
    nodeid = i
    switchid = int(i/NodeperSwitch)
    return (nodeid, switchid)

# get global id of gpu node 
def get_node_gid(i):
    assert(i < Nnodes)
    return i
    
# get global id of switch from its relative id
def get_switch_gid(i):
    assert(i < Nswitches)
    return i + Nnodes

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
    global g_device
    global Ntotalnodes
    global NodesperCluster
    global NodeperSwitch

    g_conn = []
    g_route = []

    for t in arg_topo:
        assert(t[0] > 0)
    # ids of gpu nodes are always from 0 to (ngpus - 1)
    # then, let's design like this
    # gpu i is connected to switch i/(gpus per node) + (ngpus)

    # switch to a node (independent of arg_topo)
    for i in range(Nnodes):
        (node_id, switch_id) = get_ids(i)
        g_conn.append([switch_id, node_id, 1])
        g_conn.append([node_id, switch_id, 1])

    arg_topo.sort(key=lambda x: x[0]) # sort according to the distance
    distances = [i[0] for i in arg_topo]
    #print("distances = %s" % str(distances))
    # switch to switch
    for _t in arg_topo:
        first_switch = Nnodes
        distance = _t[0]
        nconn = _t[1]
        for i in range(Nswitches):
            j = i + distance
            if (j < 0): j = j + Nswitches
            if (j >= Nswitches): j = j - Nswitches
            g_conn.append([i+first_switch, j+first_switch, nconn])
        #pr2dlist("connection matrix", g_conn)

    # route
    for i in range(Nnodes):
        for j in range(Nnodes):
            if (i == j): 
                continue
            h = get_route(i, j, distances)
            h.reverse()
            pv = [[h, 1.0]]   # single route is assumed for now
            r = [i, j, pv] # from_node, to_node, path
            g_route.append(r)
    #pr2dlist("route matrix", g_route)

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
    tg.StartTasksVector(_builder, len(device_l))
    device_v = make_vector(_builder, device_l)
    return device_v

def add_nexttasks(_builder, _nt):
    #print("nexttasks = %s" % str(_nt))
    task.StartNexttasksVector(_builder, len(_nt))
    for n in _nt:
        #print("add %s" % str(n))
        _builder.PrependUint64(n)
    nt_v = _builder.EndVector(len(_nt))
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

generate_ring1_task()
arg_topo1 = [[1, 1], [2, 1], [3, 1], [4, 1], [5, 1], [6, 1], [7, 1]] # Ring 1
arg_topo2 = [[1, 8]] # Ring 2
arg_topo3 = [[1, 7], [2,1]] # Ring 3
generate_ring(arg_topo3)
conn_v = add_conns(builder, g_conn)
route_v = add_routes(builder, g_route)
task_v = add_tasks(builder, g_task)
#print("g_device again = %s" % str(g_device))
device_v = add_devices(builder, g_device)

#build taskgraph
tg.Start(builder)
tg.AddNgpupernode(builder, 1)
tg.AddNnode(builder, Nnodes)
tg.AddNswitch(builder, Nswitches)
tg.AddIntergpubw(builder, 0.0)
tg.AddDrambw(builder, 0.0)
tg.AddNetbw(builder, 1000.0)
tg.AddConn(builder, conn_v)
tg.AddRoutes(builder, route_v)
tg.AddTasks(builder, task_v)
tg.AddDevices(builder, device_v)
tg1 = tg.End(builder)
builder.Finish(tg1)
buffer = builder.Output()
#print("file write")
fw = open('test-tg.bin', 'wb')
fw.write(buffer)
fw.close()

#build topology
tp.Start(builder)
tp.AddConn(builder, conn_v)
tp.AddRoutes(builder, route_v)
tp1 = tp.End(builder)
builder.Finish(tp1)
buffer = builder.Output()
#print("file write")
fw = open('test-tp.bin', 'wb')
fw.write(buffer)
fw.close()

