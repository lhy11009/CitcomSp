import os
import numpy as np
import RaiseError as rerr
import readfiles as readf
from matplotlib import pyplot as plt
from Manipulate import CheckDir

year = 365*24*3600
debug = True
quiet = False


class HorizAvg():

    def __init__(self, _paraDict, _colsDict, _stepsTuple, _tArray):
        self.cD = _colsDict
        self.sT = _stepsTuple
        self.tA = _tArray
        self.dataFile = _paraDict['datafile']
        self.nProcZ = readf.get_variable(_paraDict, 'nprocz', 'int')
        Alpha = readf.get_variable(_paraDict, 'thermexp', 'float')
        Tref = readf.get_variable(_paraDict, 'reftemperature', 'float')
        self.Rho = readf.get_variable(_paraDict, 'density', 'float')
        self.RhoScale = Alpha*Tref*self.Rho
        self.bA = readf.get_variable(_paraDict, 'buoyancy_ratio', 'float_list')
        Ro = readf.get_variable(_paraDict, 'radius', 'float')
        Kappa = readf.get_variable(_paraDict, 'thermdiff', 'float')
        self.tScale = Ro**2.0/Kappa/(1e6*year)
        self.vScale = Kappa/Ro*year*1e2

    def CheckStep(self, _route, _step):
        nProcZ = self.nProcZ
        dataFile = self.dataFile
        for i in range(nProcZ):
            filename = os.path.join(_route, '%s.horiz_avg.%d.%d' % (dataFile, i, _step))
            if debug:
                print("checking", filename)
            if os.path.isfile(filename) is False:
                return False
        return True

    def Read(self, _route, _step):
        self.step = _step
        dataFile = self.dataFile
        filename = os.path.join(_route, '%s.horiz_avg.0.%d' % (dataFile, _step))
        if os.path.isfile(filename) is False:
            raise rerr.NoFileError(filename)
        if debug:
            print('reading Horiz_avg: ', filename)   # debug
        Mdata = np.genfromtxt(filename)
        nProcZ = self.nProcZ
        for i in range(1, nProcZ):
            filename = os.path.join(_route, '%s.horiz_avg.%d.%d' % (dataFile, i, _step))
            if os.path.isfile(filename) is False:
                raise rerr.NoFileError(filename)
            if debug:
                print('reading Horiz_avg: ', filename)   # debug
            Mdata1 = np.genfromtxt(filename)
            end = Mdata1.shape[0]
            Mdata = np.vstack((Mdata, Mdata1[1:end, :]))
        self.Mdata = Mdata

    def Plot(self, _route, _type):
        if quiet is False:
            print("plot Horiz_avg file, step %d" % (self.step))
        CheckDir(_route)
        step = self.step
        tScale = self.tScale
        time = self.tA[step] * tScale
        self.TPlot('temp', time, _route, 'r')
        self.TPlot('visc', time, _route, 'b', 'logx')
        self.RhoPlot(time, _route, 'k')
        self.VPlot(time, _route, ['c', 'g'])

    def PlotValue(self, _rad, _T, _Name, _time, _route, _type, _Color):
        _step = self.step
        LineType = '-'
        MarkType = '.'
        rUnit = '1.0'
        TUnit = '1.0'   # get data and max min values#
        rMin = min(_rad)
        rMax = max(_rad)
        TMin = min(_T)
        TMax = max(_T)
        fig, ax = plt.subplots()   # plot
        if _type is 'normal':
            ax.plot(_T, _rad, color=_Color, linestyle=LineType,
                    label=_Name, marker=MarkType)
        elif _type is 'logx':
            ax.semilogx(_T, _rad, color=_Color, linestyle=LineType,
                        label=_Name, marker=MarkType)
        rRange = rMax - rMin    # settings
        TRange = TMax - TMin
        if _type is 'n':
            plt.xlim(TMin-0.05*TRange, TMax+0.05*TRange)
        plt.ylim(rMin, 1.0)
        ax.set(xlabel='%s [%s]' % (_Name, TUnit),
               ylabel='Radius [%s]' % (rUnit), title='step %s' % (_step))
        ax.legend()
        if _type is 'normal':
            plt.text(TMin+0.1*TRange, rMin+0.1*rRange, 'Time = %.2f Myr' % (_time))
        if _type is 'logx':
            plt.text(TMin*(TMax/TMin)**0.1, rMin+0.1*rRange, 'Time = %.2f Myr' % (_time))
        filename = os.path.join(_route, "%s_%06d.eps" % (_Name, _step))
        fig.savefig(filename)
        plt.close(fig)

    def TPlot(self, _Name, _time, _route, _Color, _type='normal'):
        cD = self.cD
        rad = self.Mdata[:, cD['radius']]
        T = self.Mdata[:, cD[_Name]]
        oRoute = os.path.join(_route, _Name)
        CheckDir(oRoute)
        self.PlotValue(rad, T, _Name, _time, oRoute, _type, _Color)

    def GetR(self):
        col = self.cD['radius']
        Radius = self.Mdata[:, col]
        return Radius

    def GetRho(self):
        cL = self.cD['chemical']
        Rho = self.Rho
        RhoScale = self.RhoScale
        bA = self.bA
        shape = self.Mdata.shape
        rho = np.zeros(shape[0])
        for i in range(len(cL)):
            col = cL[i]
            rho = rho + self.Mdata[:, col] * bA[i]
        rho = rho * RhoScale + Rho
        return rho

    def RhoPlot(self, _time, _route, _Color):
        cD = self.cD
        rad = self.Mdata[:, cD['radius']]
        rho = self.GetRho()
        PlotType = 'normal'
        Name = 'rho'
        oRoute = os.path.join(_route, Name)
        CheckDir(oRoute)
        self.PlotValue(rad, rho, Name, _time, oRoute, PlotType, _Color)

    def VPlot(self, _time, _route, _Color):
        vScale = self.vScale
        cD = self.cD
        rad = self.Mdata[:, cD['radius']]
        vr = self.Mdata[:, cD['vr']] * vScale
        vth = self.Mdata[:, cD['vth']] * vScale
        step = self.step
        Name = 'velo'
        LineType = '-'
        MarkType = '.'
        rUnit = '1.0'
        TUnit = 'cm/yr'   # get data and max min values#
        oRoute = os.path.join(_route, Name)
        CheckDir(oRoute)
        rMin = min(rad)
        rMax = max(rad)
        TMin = min(min(vr), min(vth))
        TMax = max(max(vr), max(vth))
        TRange = TMax - TMin
        fig, ax = plt.subplots()   # plot
        ax.plot(vr, rad, color=_Color[0], linestyle=LineType,
                label='vr', marker=MarkType)
        ax.plot(vth, rad, color=_Color[1], linestyle=LineType,
                label='vth', marker=MarkType)
        rRange = rMax - rMin    # settings
        plt.xlim(TMin-0.05*TRange, TMax+0.05*TRange)
        plt.ylim(rMin, 1.0)
        ax.set(xlabel='%s [%s]' % (Name, TUnit),
               ylabel='Radius [%s]' % (rUnit), title='step %s' % (step))
        ax.legend()
        plt.text(TMin+0.1*TRange, rMin+0.1*rRange, 'Time = %.2f Myr' % (_time))
        filename = os.path.join(oRoute, "%s_%06d.eps" % (Name, step))
        fig.savefig(filename)
        plt.close(fig)

    def SumInRange(self, _col, _range, _type):
        Mdata = self.Mdata
        shape = Mdata.shape
        cD = self.cD
        rad = self.Mdata[:, cD['radius']]
        Sum = 0
        for i in range(shape[0]-1):
            rb = rad[i]
            ru = rad[i+1]
            rr = (rb+ru)/2.0
            if _range is None:
                Cond = True
            elif (rr > _range[0] and rr < _range[1]):
                Cond = True
            else:
                Cond = False
            if Cond:
                weight = 4.0/3.0 * np.pi * (ru**3.0 - rb**3.0)
                if _type is 'normal':
                    vv = (Mdata[i+1, _col] + Mdata[i, _col])/2.0
                elif _type is 'abs':
                    vMean = (Mdata[i+1, _col] + Mdata[i, _col])/2.0
                    vv = abs(vMean)/2.0
                elif _type is 'sqrt':
                    vv = 0.0
                    for col in _col:
                        vMean = (Mdata[i+1, col] + Mdata[i, col])/2.0
                        vv = vv + vMean**2.0
                    vv = vv**0.5
                else:
                    ErrMessage = "SumInRange type cannot be %s" % (_type)
                    raise rerr.FuncInputError(ErrMessage)
                Sum = Sum + weight * vv
        if debug:
            print('Sum =', Sum)
        return Sum


def PlotHaSum(_route, _time, _steps, _sum, _Name, _type=None):
    if quiet is False:
        print("plot sum from Horiz_avg file:", _Name)
    oRoute = os.path.join(_route, 'ha_sum')
    CheckDir(oRoute)
    time = _time[_steps]
    fig, ax = plt.subplots()
    ax.plot(time, _sum)
    ax.set(xlabel='Time [Ma]', ylabel=_Name)
    if _type is 'normal':
        plt.ylim((-0.05, 1.05))
    filename = os.path.join(oRoute, _Name)
    fig.tight_layout()
    fig.savefig(filename)
