# automatically generated by the FlatBuffers compiler, do not modify

# namespace: FlatBufTaskGraph

import flatbuffers
from flatbuffers.compat import import_numpy
np = import_numpy()

class Task(object):
    __slots__ = ['_tab']

    @classmethod
    def GetRootAs(cls, buf, offset=0):
        n = flatbuffers.encode.Get(flatbuffers.packer.uoffset, buf, offset)
        x = Task()
        x.Init(buf, n + offset)
        return x

    @classmethod
    def GetRootAsTask(cls, buf, offset=0):
        """This method is deprecated. Please switch to GetRootAs."""
        return cls.GetRootAs(buf, offset)
    # Task
    def Init(self, buf, pos):
        self._tab = flatbuffers.table.Table(buf, pos)

    # Task
    def Type(self):
        o = flatbuffers.number_types.UOffsetTFlags.py_type(self._tab.Offset(4))
        if o != 0:
            return self._tab.Get(flatbuffers.number_types.Int16Flags, o + self._tab.Pos)
        return 0

    # Task
    def Taskid(self):
        o = flatbuffers.number_types.UOffsetTFlags.py_type(self._tab.Offset(6))
        if o != 0:
            return self._tab.Get(flatbuffers.number_types.Uint64Flags, o + self._tab.Pos)
        return 0

    # Task
    def Deviceid(self):
        o = flatbuffers.number_types.UOffsetTFlags.py_type(self._tab.Offset(8))
        if o != 0:
            return self._tab.Get(flatbuffers.number_types.Uint64Flags, o + self._tab.Pos)
        return 0

    # Task
    def Runtime(self):
        o = flatbuffers.number_types.UOffsetTFlags.py_type(self._tab.Offset(10))
        if o != 0:
            return self._tab.Get(flatbuffers.number_types.Float32Flags, o + self._tab.Pos)
        return 0.0

    # Task
    def Xfersize(self):
        o = flatbuffers.number_types.UOffsetTFlags.py_type(self._tab.Offset(12))
        if o != 0:
            return self._tab.Get(flatbuffers.number_types.Uint64Flags, o + self._tab.Pos)
        return 0

    # Task
    def Nexttasks(self, j):
        o = flatbuffers.number_types.UOffsetTFlags.py_type(self._tab.Offset(14))
        if o != 0:
            a = self._tab.Vector(o)
            return self._tab.Get(flatbuffers.number_types.Uint64Flags, a + flatbuffers.number_types.UOffsetTFlags.py_type(j * 8))
        return 0

    # Task
    def NexttasksAsNumpy(self):
        o = flatbuffers.number_types.UOffsetTFlags.py_type(self._tab.Offset(14))
        if o != 0:
            return self._tab.GetVectorAsNumpy(flatbuffers.number_types.Uint64Flags, o)
        return 0

    # Task
    def NexttasksLength(self):
        o = flatbuffers.number_types.UOffsetTFlags.py_type(self._tab.Offset(14))
        if o != 0:
            return self._tab.VectorLen(o)
        return 0

    # Task
    def NexttasksIsNone(self):
        o = flatbuffers.number_types.UOffsetTFlags.py_type(self._tab.Offset(14))
        return o == 0

def Start(builder): builder.StartObject(6)
def TaskStart(builder):
    """This method is deprecated. Please switch to Start."""
    return Start(builder)
def AddType(builder, type): builder.PrependInt16Slot(0, type, 0)
def TaskAddType(builder, type):
    """This method is deprecated. Please switch to AddType."""
    return AddType(builder, type)
def AddTaskid(builder, taskid): builder.PrependUint64Slot(1, taskid, 0)
def TaskAddTaskid(builder, taskid):
    """This method is deprecated. Please switch to AddTaskid."""
    return AddTaskid(builder, taskid)
def AddDeviceid(builder, deviceid): builder.PrependUint64Slot(2, deviceid, 0)
def TaskAddDeviceid(builder, deviceid):
    """This method is deprecated. Please switch to AddDeviceid."""
    return AddDeviceid(builder, deviceid)
def AddRuntime(builder, runtime): builder.PrependFloat32Slot(3, runtime, 0.0)
def TaskAddRuntime(builder, runtime):
    """This method is deprecated. Please switch to AddRuntime."""
    return AddRuntime(builder, runtime)
def AddXfersize(builder, xfersize): builder.PrependUint64Slot(4, xfersize, 0)
def TaskAddXfersize(builder, xfersize):
    """This method is deprecated. Please switch to AddXfersize."""
    return AddXfersize(builder, xfersize)
def AddNexttasks(builder, nexttasks): builder.PrependUOffsetTRelativeSlot(5, flatbuffers.number_types.UOffsetTFlags.py_type(nexttasks), 0)
def TaskAddNexttasks(builder, nexttasks):
    """This method is deprecated. Please switch to AddNexttasks."""
    return AddNexttasks(builder, nexttasks)
def StartNexttasksVector(builder, numElems): return builder.StartVector(8, numElems, 8)
def TaskStartNexttasksVector(builder, numElems):
    """This method is deprecated. Please switch to Start."""
    return StartNexttasksVector(builder, numElems)
def End(builder): return builder.EndObject()
def TaskEnd(builder):
    """This method is deprecated. Please switch to End."""
    return End(builder)