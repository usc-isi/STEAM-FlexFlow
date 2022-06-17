
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
import copy

builder = fb.Builder(0)

# describe topology
# num_nodes, connection
Ngpupernode = 2
Nnodes = 16 
Nswitches = 8
NodeperSwitch = int(Nnodes / Nswitches)
Ntotalnodes = Nnodes + Nswitches
Gpubw = 0.0
Drambw = 0.0
Netbw = 1000.0

assert(Ntotalnodes == Nnodes + Nswitches)

g_conn = []
g_route = []
g_device = []

def make_vector(_builder, elm_l):
    for e in elm_l:
        _builder.PrependUOffsetTRelative(e)
    return builder.EndVector(len(elm_l))

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

def get_route(src, dest, distance_l):
    if (src == dest): return []
    (src_nodeid, src_switchid) = get_rel_ids(src)
    (dest_nodeid, dest_switchid) = get_rel_ids(dest)
    h = [src]

    if (src_switchid == dest_switchid):
        h = [src, src_switchid + Nnodes, dest]
        return h

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
    for r in _rd:
        _builder.PrependInt32(r)
    jumps = _builder.EndVector(len(_rd))
    return jumps

def add_rings(_builder, _rings):
    print("rings = %s" % str(_rings))
    if (_rings == []): return None
    rings_l = []
    for r in _rings:
        rd_v = add_ringdescriptorjumps(_builder, r[1])
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
            chance = mypath.Chance()
#            print("\thopnode = %s : prob %f" % (str(hopnode), chance))
            p_l.append([hopnode.tolist(), chance])
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
    for i in range(myrings_length):
        myring = mytg.Rings(i)
        ringsz = mytg.Ringsz
        myringpaths_len = myring.RingpathsLength()
        l = []
        for j in range(myringpaths_len):
            myringpath = myring.Ringpaths(j)
            jumps = myringpath.JumpsAsNumpy()
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

output_fname = 'ringtopo1.bin'
generate_ring1_task()
arg_topo1 = [[1, 1], [2, 1], [3, 1], [4, 1], [5, 1], [6, 1], [7, 1]] # Ring 1
arg_topo2 = [[1, 8]] # Ring 2
arg_topo3 = [[1, 7], [2,1]] # Ring 3
generate_ring(arg_topo1)
conn_v = add_conns(builder, g_conn)
route_v = add_routes(builder, g_route)
task_v = add_tasks(builder, g_task)
device_v = add_devices(builder, g_device)

#build_taskgraph(builder, Ngpupernode, Nnodes, Nswitches, Gpubw, Drambw, Netbw, conn_v, route_v, task_v, device_v, 'test-tg.bin')
#build_topology(builder, conn_v, route_v, 'test-tp.bin')

tg_l = read_taskgraph('output-tg.fattree')
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
