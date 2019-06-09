import numpy as np
import os
import readfiles as readf

debug = True

def LaterBetter(_data, _col):
    shape = _data.shape
    _t = _data[:, _col]
    now = 0
    Prev = shape[0] - 1
    Next = shape[0] - 2
    data = np.zeros(shape)
    data[now, :] = _data[Prev, :]
    while Next >= 0:
        if _t[Next] < _t[Prev]:
            Prev = Next
            now = now + 1
            data[now, :] = _data[Prev, :]
        Next = Next - 1
    data = data[0:now+1, :]
    data = data[::-1, :]
    return data


def FindSteps(_route, _caseName, _maxStep):
    Steps = []
    for step in range(_maxStep):
        filename = "%s/%s.horiz_avg.0.%d" % (_route, _caseName, step)
        if os.path.isfile(filename):
            Steps.append(step)
    Steps = np.array(Steps)
    return Steps

def CheckDir(_dirName):
    if os.path.isdir(_dirName) is False:
        os.mkdir(_dirName)

def IbcLayer(_paraDict):
    Range = []
    IbcInterfaceIndex = 1
    blur = 0.8
    InterF = readf.get_variable(_paraDict, 'z_interface', 'float_list')
    zLith = readf.get_variable(_paraDict, 'z_lith', 'float')
    ro = readf.get_variable(_paraDict, 'radius_outer', 'float')
    upper = ro - zLith
    lower = InterF[IbcInterfaceIndex]
    thick = upper - lower
    Range.append(lower-blur*thick)
    Range.append(upper+blur*thick)
    return Range

def LowerMantle(_paraDict):
    tiny = 1e-6
    Range = []
    ro = readf.get_variable(_paraDict, 'radius_outer', 'float')
    ri = readf.get_variable(_paraDict, 'radius_inner', 'float')
    Range.append(ri-tiny)
    Range.append((ri+ro)/2.0+tiny)
    return Range
