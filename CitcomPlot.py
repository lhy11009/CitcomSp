import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as Mplt
import readfiles as Mread
import plotfigures as Mpf
import post_processing as pp
from HorizPlot import HorizAvg
from HorizPlot import PlotHaSum
from TimeFile import TimeDat
from Manipulate import FindSteps
from Manipulate import IbcLayer
from Manipulate import LowerMantle
from Files import fFileColD
from Files import HzFileColD

Route = "/home/p300/Desktop/SpaImpactOverturn/I420deg60ibc50v1e-3R65c192"
CaseName = "a"
InputFile = "%s/input_solidus_i/in_65moon" % (Route)

ppFile = "%s/post_process" % (Route)
plogfile = '%s/plog' % (Route)
MaxStep = 100001
oRoute = './test'
# read input #
inD = Mread.read_input(InputFile)
# get step list #
Steps = FindSteps(Route, CaseName, MaxStep)
TD = TimeDat(inD, fFileColD, Steps)
ToF = TD.Check(Route)
TD.Read(Route)
timeArray = TD.GetTime()
TD.PlotMachineTime(oRoute)
# horiz_avg file #
IbcCol = HzFileColD['chemical'][2]
VrCol = HzFileColD['vr']
VthCol = HzFileColD['vth']
VCol = [VrCol, VthCol]
rIbc = IbcLayer(inD)
rLM = LowerMantle(inD)
HA = HorizAvg(inD, HzFileColD, Steps, timeArray)
HA.Read(Route, 0)
rr = HA.GetR()
i = 0
size = Steps.size
Vrms = np.zeros(size)
Ibc = np.zeros(size)
IbcInLayer = np.zeros(size)
IbcLM = np.zeros(size)
for step in Steps:
    HA.CheckStep(Route, step)
    HA.Read(Route, step)
    HA.Plot(oRoute, None)
    Vrms[i] = HA.SumInRange(VCol, None, 'sqrt')
    Ibc[i] = HA.SumInRange(IbcCol, None, 'normal')
    IbcInLayer[i] = HA.SumInRange(IbcCol, rIbc, 'normal')
    IbcLM[i] = HA.SumInRange(IbcCol, rLM, 'normal')
    i = i+1
IbcOri = Ibc[0]
PlotHaSum(oRoute, timeArray, Steps, Ibc/IbcOri, 'Total_IBC')
PlotHaSum(oRoute, timeArray, Steps, IbcInLayer/Ibc, 'Retained_IBC', 'normal')
PlotHaSum(oRoute, timeArray, Steps, IbcLM/Ibc, 'Overturned_IBC', 'normal')
input()
#---------------qb.dat and qs.dat----------------------#
#pp.heat_flux(inD,time_array,Steps,Route,oroute,CaseName,cols=5,colT=None)
#------------------mf file------------------------------#
settings={'col':10}
settings['column']={'step':0,'tracer':[3,6,7,8,9],'eular':4} #first in melting is total melting
settings['color']={'tracer':'k','eular':'y','melting':['k','r','b','g','c']}
pp.plot_mf_melting(inD,time_array,Route,oroute,CaseName,settings)
cdict={'tracer':[0,2,3,4,5],'eular':1}
print("mf done")
print("continue?")
input()
#-------------------MF file-----------------------------#
nprocx = int(inD['nprocx'])
nprocy = int(inD['nprocy'])
nprocz = int(inD['nprocz'])
nodex = int(inD['nodex'])
nodez = int(inD['nodez'])
nox = (nodex+1)//nprocx
noz = (nodez+1)//nprocz
nproc = 12*nprocx*nprocy*nprocz
plot_MF_surf=pp.plot_MF_surf(cdict,nproc,nprocz,nox,noz)
plot_MF_surf(Route,CaseName,7400,oroute)
input()
#--------------------remove files----------------------------#
settings={'col':10}
settings['column']={'step':0,'temp':1,'melting':[5,6,7,8,9]} #first in melting is total melting
settings['color']={'temp':'r','melting':['k','r','b','g','c']}
#settings={'col':6}
#settings['column']={'step':0,'temp':1,'melting':[5]} #first in melting is total melting
#settings['color']={'temp':'r','melting':['k']}
#pp.plot_volume_1(inD,time_array,Route,oroute,CaseName,settings)
p_dict=Mread.read_input(ppFile)
e_steps=Mread.get_variable(p_dict,'episode_steps','int_list')
e_ends = e_steps[len(e_steps)-1]
#try:
#    steps_remain=Mread.get_variable(p_dict,'steps_remain','int_list')
#except KeyError:
#    pp.remove_extra_file(Route,CaseName,inD,time_array,Steps,ppFile,plogfile)
#else:
#    steps_remain_end = steps_remain[len(steps_remain)-1]
#    if steps_remain_end != Steps[len(Steps)-1]:
#        pp.remove_extra_file1(Route,CaseName,inD,time_array,Steps,steps_remain_end,ppFile,plogfile)
#pp.sph_plot(inD,rr,Route,oroute,CaseName,ppFile,6000,'T')
#--------------------------sph file-------------------------#
for step in e_steps:
#for step in steps_remain:
    pp.sph_plot(inD,rr,Route,oroute,CaseName,ppFile,step,'T')
#step = steps_remain[-1]
#pp.sph_plot(inD,rr,Route,oroute,CaseName,ppFile,step,'T')

#-----------------------surf file---------------------------#

#print('sph ploting done')
#input()
#temp0 = [0.0,0.0]
#for step in steps_remain:
#    temp = pp.surf_plot(inD,Route,oroute,CaseName,step,vtype='total',bound=[-25000,25000,0.0],minor=[5000,5000,1000],check_bound=False,rot=[0,90],angle=[-45,-15])
#,cpt='colorS.cpt')
#rot=[0,90],angle=[-45,-15],
#    temp0[0] = min(temp[0],temp0[0])
#    temp0[1] = max(temp[1],temp0[1])
#print(temp0)
#dfile = pp.hoz_combine(inD,Route,CaseName,time_array,Steps,HzFileColD)
#tval = [0.0,500.0,50.0,10.0]
#tval = [0.0,200.0,20.0,4.0]
#tval = [0.0,100.0,10.0,2.0]
#tval = [0.0,50.0,5.0,1.0]
#tval = [0.0,20.0,2.0,0.4]
#pp.plot_combine(inD,oRoute,CaseName,dfile,'temp',cdict,tval,cpt='colorT.cpt')
#pp.plot_combine(inD,oRoute,CaseName,dfile,'rho',cdict,tval,rhodminor=1.0,rhocoff=10.0)

