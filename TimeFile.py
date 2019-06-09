import os
import numpy as np
import RaiseError as rerr
from Manipulate import LaterBetter
from Manipulate import CheckDir
from matplotlib import pyplot as plt
import readfiles as readf

hour = 3600
year = 365*24*3600
debug = True


class TimeDat():

    def __init__(self, _paraDict, _colsDict, _stepsTuple):
        self.cD = _colsDict
        self.sT = _stepsTuple
        self.dataFile = _paraDict['datafile']
        Ro = readf.get_variable(_paraDict, 'radius', 'float')
        Kappa = readf.get_variable(_paraDict, 'thermdiff', 'float')
        self.tScale = Ro**2.0/Kappa/(1e6*year)

    def Check(self, _route):
        dataFile = self.dataFile
        filename = os.path.join(_route, "%s.time" % (dataFile))
        return os.path.isfile(filename)

    def Read(self, _route):
        cD = self.cD
        StepCol = cD['step']
        dataFile = self.dataFile
        filename = os.path.join(_route, "%s.time" % (dataFile))
        if os.path.isfile(filename) is False:
            raise rerr.NoFileError(filename)
        Mdata = np.genfromtxt(filename)
        Mdata = LaterBetter(Mdata, StepCol)
        self.Mdata = Mdata

    def GetTime(self):
        cD = self.cD
        tCol = cD['time']
        return self.Mdata[:, tCol]

    def PlotMachineTime(self, _route):
        CheckDir(_route)
        cD = self.cD
        tCol = cD['time']
        Mdata = self.Mdata
        t = Mdata[:, tCol]*self.tScale
        MtCol = cD['machine_time']
        Mt = Mdata[:, MtCol]/hour
        if debug:
            np.savetxt('./debugMt', Mt)
        StepCol = cD['step']
        Step = Mdata[:, StepCol]
        Color = ['c', 'k']
        LineType = '-'
        fig, ax1 = plt.subplots()   # plot
        ax1.plot(Step, t, color=Color[0], linestyle=LineType, label='Model Time')
        ax1.set(xlabel='Step', ylabel='Model Time [Ma]')
        ax1.set_ylabel('Model Time', color=Color[0])
        ax1.tick_params(axis='y', labelcolor=Color[0])
        LineType = '-.'
        ax2 = ax1.twinx()
        ax2.plot(Step, Mt, color=Color[1], linestyle=LineType, label='Machine Time')
        ax2.set_ylabel('Machine Time [h]', color=Color[1])
        ax2.tick_params(axis='y', labelcolor=Color[1])
        fig.tight_layout()
        filename = os.path.join(_route, "time.eps")
        fig.savefig(filename)
