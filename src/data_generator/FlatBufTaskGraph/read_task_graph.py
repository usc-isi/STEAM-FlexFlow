
# namespace: FlatBufTaskGraph
import flatbuffers 
import TaskGraph as tg
#import Connection as conn
import Device as device
import Route as route
import Path as path
import Task as task
import Topology as tp


buf = open('output-tg.fattree', 'rb').read()
buf = bytearray(buf)
mytg = tg.TaskGraph.GetRootAs(buf, 0)
ngpupernode = mytg.Ngpupernode()
print("ngpupernode = %d" % ngpupernode)
myconn_length = mytg.ConnLength()
print("myconn_length = %d" % myconn_length)
for i in range(myconn_length):
    myconn = mytg.Conn(i)
#    print("c = %s" % str(c))
#tasks = mytg.Tasks
#myconn = conn.Connection.GetRootAsConnection(conn,0)

#for c in conn:
#   print("tasks = %s" % str(c))
#print("conn = %s" % str(conn))


