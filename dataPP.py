import os

#------------------------------lhy 20180530---------------------
#derive index in paraview saved images for a input steps
#this is useful when saving interval changes in CitcomS output
def map_steps_to_index_pv(route,caseName,maxStep):
    d = dict()
    index = 0
    for step in range(0,maxStep+1):
        filename = "%s/%s.horiz_avg.0.%d" % (route,caseName,step)
        if os.access(filename,os.F_OK):
            d[step] = index
            index+=1
    return d

