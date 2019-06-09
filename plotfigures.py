import matplotlib.pyplot as Mplt
import os
import sys
sys.path.append('../public')

#----plot horizontal average result at a single step----#
def plot_horiz_avg(xdata,ydata,Cstep,Ctime,Csettings):
    def _conf(name,default):
        try:
            return Csettings[name]
        except KeyError:
            return default
    name=Csettings['name']
    color=Csettings['color']
    x_unit=Csettings['x_unit']
    y_unit=Csettings['y_unit']
    o_dir=Csettings['o_dir']
    plot_type=Csettings['plot_type']
    line_type=_conf('line','-')
    mark_type=_conf('mark',None)
    try:
        col = xdata.shape[1] 
    except IndexError:
        col = 1
    #----get data and max min values----#
    x_mm = []
    y_mm = []
    if col==1:
        x_mm.append(min(xdata))
        x_mm.append(max(xdata))
    else:
        x_mm.append(min(xdata[:,0]))
        x_mm.append(max(xdata[:,0]))
        for c in range(1,col):
            temp=min(xdata[:,c])
            if temp<x_mm[0]:
                x_mm[0]=temp
            temp=max(xdata[:,c])
            if temp>x_mm[1]:
                x_mm[1]=temp
    y_mm.append(min(ydata))
    y_mm.append(max(ydata))
    #----plot----#
    fig, ax = Mplt.subplots()
    if col==1:
        if plot_type is 'n':
            ax.plot(xdata, ydata,color='%s'%(color),linestyle=line_type,label=name,marker=mark_type)
        else:
            ax.semilogx(xdata, ydata,color='%s'%(color),linestyle=line_type,label=name,marker=mark_type)
    else:
        for c in range(col):
            if plot_type is 'n':
                ax.plot(xdata[:,c], ydata,color='%s'%(color[c]),linestyle=line_type,label=name[c+1],marker=mark_type)
            else:
                ax.semilogx(xdata[:,c], ydata,color='%s'%(color[c]),linestyle=line_type,label=name[c+1],marker=mark_type)
    #----settings----#
    x_range=x_mm[1]-x_mm[0]
    y_range=y_mm[1]-y_mm[0]
    if plot_type is 'n':
        Mplt.xlim(x_mm[0]-0.05*x_range, x_mm[1]+0.05*x_range)
    Mplt.ylim(y_mm[0], 1.0)
    ax.set(xlabel='[%s]'%(x_unit), ylabel='Radius [%s]'%(y_unit),
       title='step %s'%(Cstep))
    ax.legend()
    Mplt.text(x_mm[0]+0.1*x_range,y_mm[0]+0.1*y_range,'Time = %.2f Myr'%(Ctime))
    #print("%s/%s_%06d.png"%(o_dir,name,Cstep)) #debug
    if col==1:
        fig.savefig("%s/%s_%06d.eps"%(o_dir,name,Cstep))
    else:
        fig.savefig("%s/%s_%06d.eps"%(o_dir,name[0],Cstep))
    Mplt.close(fig)
    #Mplt.show()




#----creat dirname if not exit----#
def assign_output_dir(*args):
    i = 0
    for mem in args:
        if i == 0:
            dirname = mem
        else:
            dirname = os.path.join(dirname,mem)
        i = i+1
    if os.path.isdir(dirname) is False:
        os.mkdir(dirname)
    return dirname

def data_plot(xdata,ydata,settings):
    def _set(key,default):
        try:
            return settings[key]
        except KeyError:
            return default
    try:
        col = xdata.shape[1] 
    except IndexError:
        col = 1
    name=_set('name','Data') #bug here when default is not equal to col
    color=_set('color','r')
    plot_type=_set('plot_type','n')
    text=_set('text',None)
    x_mm = []
    y_mm = []
    x_mm.append(xdata.min())
    x_mm.append(xdata.max())
    y_mm.append(ydata.min())
    y_mm.append(ydata.max())
    fig, ax = Mplt.subplots()
    if col==1:
        if plot_type is 'n':
            ax.plot(xdata, ydata,color='%s'%(color),linestyle=line_type,label=name,marker=mark_type)
        else:
            ax.semilogx(xdata, ydata,color='%s'%(color),linestyle=line_type,label=name,marker=mark_type)
    else:
        for c in range(col):
            if plot_type is 'n':
                ax.plot(xdata[:,c], ydata,color='%s'%(color[c]),linestyle=line_type,label=name[c+1],marker=mark_type)
            else:
                ax.semilogx(xdata[:,c], ydata,color='%s'%(color[c]),linestyle=line_type,label=name[c+1],marker=mark_type)
    #----settings----#
    x_range=x_mm[1]-x_mm[0]
    y_range=y_mm[1]-y_mm[0]
    if plot_type is 'n':
        Mplt.xlim(x_mm[0]-0.05*x_range, x_mm[1]+0.05*x_range)
    Mplt.ylim(y_mm[0], 1.0)
    ax.set(xlabel='[%s]'%(x_unit), ylabel='Radius [%s]'%(y_unit),
       title='step %s'%(Cstep))
    ax.legend()
    if text is not None:
        Mplt.text(x_mm[0]+0.1*x_range,y_mm[0]+0.1*y_range,text)
    return fig
