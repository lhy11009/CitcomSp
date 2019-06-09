import sys
sys.path.append('../public')
import numpy as np
import matplotlib.pyplot as plt
import math
import plotfigures as Mpf
import readfiles as readf
import globalp as Gp
import os
import surf_combine as surfc

def horiz_avg(in_dict,time_array,step_tuple,column_dict,route,o_dir,case_name):
#plot horizontal_averagye result
    nprocz=readf.get_variable(in_dict,'nprocz','int')
    noz=readf.get_variable(in_dict,'nodez','int')
    Ra=readf.get_variable(in_dict,'rayleigh','float')
    Buoy=readf.get_variable(in_dict,'buoyancy_ratio','float_list')
    Alpha=readf.get_variable(in_dict,'thermexp','float')
    Tref=readf.get_variable(in_dict,'reftemperature','float')
    Tsurf=readf.get_variable(in_dict,'surftemperature','float')
    Rho=readf.get_variable(in_dict,'density','float')
    Eta=readf.get_variable(in_dict,'refvisc','float')
    Ro=readf.get_variable(in_dict,'radius','float')
    Kappa=readf.get_variable(in_dict,'thermdiff','float')
    lnoz=int((noz-1)/nprocz)+1
    cols=column_dict['cols']
    n_c=len(column_dict['chemical'])
    #----derive scalings----#
    year=365*24*3600
    r_scaling=Ro/1e3
    time_scaling=pow(Ro,2.0)/Kappa/(1e6*year)
    rho_scaling=Alpha*Tref*Rho
    velo_scaling=Kappa/Ro*year*1e2
    #---output informations---#
    case_dir = Mpf.assign_output_dir("%s/%s"%(o_dir,case_name))
    c_dir = Mpf.assign_output_dir("%s/C"%(case_dir))
    temp_dir = Mpf.assign_output_dir("%s/temp"%(case_dir))
    velo_dir = Mpf.assign_output_dir("%s/velo"%(case_dir))
    visc_dir = Mpf.assign_output_dir("%s/visc"%(case_dir))
    rho_dir = Mpf.assign_output_dir("%s/rho"%(case_dir))
    for step in step_tuple:
        print(step,end='\r')
        data_in = readf.read_horig_output(route,case_name,step,nprocz,lnoz,cols)
        time=time_array[step]*time_scaling
        ##----plot figures----#
        col_r=column_dict['radius']
        ydata=data_in[:,col_r]
        settings={'name':'temp','x_unit':'1','y_unit':'1','o_dir':temp_dir,'plot_type':'n','color':'r','mark':'.'}
        col_t=column_dict['temp']
        xdata=data_in[:,col_t]
        Mpf.plot_horiz_avg(xdata,ydata,step,time,settings)
        temp=['C']
        #----plot chemicals----#
        for n in range(n_c):
            temp.append('chemical %d'%(n))
        settings={'name':temp,'x_unit':'1','y_unit':'1','o_dir':c_dir,'plot_type':'n','color':['r','b','g','c'],'line':':'}
        xdata=np.zeros((noz,n_c))
        n=0
        for col_c in column_dict['chemical']:
            xdata[:,n]=data_in[:,col_c]
            n=n+1
        Mpf.plot_horiz_avg(xdata,ydata,step,time,settings)
        settings={'name':'rho','x_unit':'kg','y_unit':'1','o_dir':rho_dir,'plot_type':'n','color':'k','mark':'.'}
        i=0
        xdata=np.zeros(noz)
        for col_c in column_dict['chemical']:
            xdata=xdata+data_in[:,col_c]*Buoy[i]
            i=i+1
        xdata=xdata*rho_scaling+Rho
        Mpf.plot_horiz_avg(xdata,ydata,step,time,settings)
        settings={'name':'visc','x_unit':'1','y_unit':'1','o_dir':visc_dir,'plot_type':'l','color':'b','mark':'.'}
        col_v=column_dict['visc']
        xdata=data_in[:,col_v]
        Mpf.plot_horiz_avg(xdata,ydata,step,time,settings)
        settings={'name':['velo','vr','vth'],'x_unit':'cm/yr','y_unit':'1.0','o_dir':velo_dir,'plot_type':'n','color':('c','g'),'mark':'.'}
        xdata=np.zeros((noz,2))
        col_vr=column_dict['vr']
        col_vth=column_dict['vth']
        xdata[:,0]=data_in[:,col_vr]
        xdata[:,1]=data_in[:,col_vth]
        Mpf.plot_horiz_avg(xdata,ydata,step,time,settings)

def horiz_avg_1(in_dict,time_array,step_tuple,column_dict,route,o_dir,case_name):
#plot horizontal_averagye result, no chemical
    nprocz=readf.get_variable(in_dict,'nprocz','int')
    noz=readf.get_variable(in_dict,'nodez','int')
    Ra=readf.get_variable(in_dict,'rayleigh','float')
    Buoy=readf.get_variable(in_dict,'buoyancy_ratio','float_list')
    Alpha=readf.get_variable(in_dict,'thermexp','float')
    Tref=readf.get_variable(in_dict,'reftemperature','float')
    Tsurf=readf.get_variable(in_dict,'surftemperature','float')
    Rho=readf.get_variable(in_dict,'density','float')
    Eta=readf.get_variable(in_dict,'refvisc','float')
    Ro=readf.get_variable(in_dict,'radius','float')
    Kappa=readf.get_variable(in_dict,'thermdiff','float')
    lnoz=int((noz-1)/nprocz)+1
    cols=column_dict['cols']
    n_c=len(column_dict['chemical'])
    #----derive scalings----#
    year=365*24*3600
    r_scaling=Ro/1e3
    time_scaling=pow(Ro,2.0)/Kappa/(1e6*year)
    velo_scaling=Kappa/Ro*year*1e2
    #---output informations---#
    case_dir = Mpf.assign_output_dir("%s/%s"%(o_dir,case_name))
    temp_dir = Mpf.assign_output_dir("%s/temp"%(case_dir))
    velo_dir = Mpf.assign_output_dir("%s/velo"%(case_dir))
    visc_dir = Mpf.assign_output_dir("%s/visc"%(case_dir))
    for step in step_tuple:
        print(step,end='\r')
        data_in = readf.read_horig_output(route,case_name,step,nprocz,lnoz,cols)
        time=time_array[step]*time_scaling
        ##----plot figures----#
        col_r=column_dict['radius']
        ydata=data_in[:,col_r]
        settings={'name':'temp','x_unit':'1','y_unit':'1','o_dir':temp_dir,'plot_type':'n','color':'r','mark':'.'}
        col_t=column_dict['temp']
        xdata=data_in[:,col_t]
        Mpf.plot_horiz_avg(xdata,ydata,step,time,settings)
        settings={'name':'visc','x_unit':'1','y_unit':'1','o_dir':visc_dir,'plot_type':'l','color':'b','mark':'.'}
        col_v=column_dict['visc']
        xdata=data_in[:,col_v]
        Mpf.plot_horiz_avg(xdata,ydata,step,time,settings)
        settings={'name':['velo','vr','vth'],'x_unit':'cm/yr','y_unit':'1.0','o_dir':velo_dir,'plot_type':'n','color':('c','g'),'mark':'.'}
        xdata=np.zeros((noz,2))
        col_vr=column_dict['vr']
        col_vth=column_dict['vth']
        xdata[:,0]=data_in[:,col_vr]
        xdata[:,1]=data_in[:,col_vth]
        Mpf.plot_horiz_avg(xdata,ydata,step,time,settings)


def chemical(in_dict,time_array,step_tuple,column_dict,route,o_dir,case_name,pp_file):
    get_end=1
    split=6
    init_lower_percent=0.01
    pi=np.pi
    tiny=1e-8
    year=365*24*3600
    case_dir = Mpf.assign_output_dir("%s/%s"%(o_dir,case_name))
    o_dir = Mpf.assign_output_dir("%s/chemical"%(case_dir))
    if os.path.isfile(pp_file):
        pp_dict = readf.read_input(pp_file)
        try:
            pp_dict['end_time']
        except KeyError:
            get_end = 1
        else:
            get_end = 0
    else:
        get_end = 1
    #----read parameters----#
    cols=column_dict['cols']
    nprocz=readf.get_variable(in_dict,'nprocz','int')
    noz=readf.get_variable(in_dict,'nodez','int')
    Ro=readf.get_variable(in_dict,'radius','float')
    interface=readf.get_variable(in_dict,'z_interface','float_list')
    r_inter=interface[-1]
    r_o=readf.get_variable(in_dict,'radius_outer','float')
    r_in=readf.get_variable(in_dict,'radius_inner','float')
    r_lith=r_o-readf.get_variable(in_dict,'z_lith','float')
    Kappa=readf.get_variable(in_dict,'thermdiff','float')
    r_half=(r_o+r_in)/2.0
    lnoz=int((noz-1)/nprocz)+1
    time_scaling=pow(Ro,2.0)/Kappa/(1e6*year)
    total_ibc=4*pi/3*(pow(r_lith,3.0)-pow(r_lith-20e3/Ro,3.0))
    #----prepare chemical data----#
    chemical=np.zeros((len(step_tuple),3))
    time1 = np.zeros(len(step_tuple))
    n = 0
    for step in step_tuple:
        time1[n] = time_array[step]
        print(step,end='\r')
        data_in = readf.read_horig_output(route,case_name,step,nprocz,lnoz,cols)
        rr = data_in[:,column_dict['radius']]
        C = data_in[:,column_dict['ic']]
        for i in range(noz-1):
            volume = 4.0/3.0*pi*(pow(rr[i+1],3.0)-pow(rr[i],3.0))
            if rr[i]>r_inter-tiny and rr[i]<r_lith+tiny:
                chemical[n,0] = chemical[n,0]+(C[i]+C[i+1])/2.0*volume
            if rr[i]>r_half-tiny:
                chemical[n,1] = chemical[n,1]+(C[i]+C[i+1])/2.0*volume
            if rr[i]<r_half+tiny:
                chemical[n,2] = chemical[n,2]+(C[i]+C[i+1])/2.0*volume
        n = n+1
    total_chemical = chemical[0,1]+chemical[0,2]
    chemical_low = chemical[:,2] #chemical_low is lower chemical percent
    #print(chemical[0,1]+chemical[0,2])
    #print(chemical[:,0]) #debug
    #print(chemical[:,1])
    #print(chemical[:,2])
    #----plot figures----#
    #----figure 1 : retained chemical----#
    fig, ax = plt.subplots()
    line1, = ax.plot(time1*time_scaling, chemical[:,0]/chemical[0,0],'r-')
    ax.set(xlabel='Time [ma]', ylabel='Chemical Percent[1]', title='Retained Chemicals in Initial Layer')
    fig.savefig("%s/retained_ibc.eps"%(o_dir))
    #plt.show()
    plt.close(fig)
    #----figure 2 : overturned chemical----#
    time2 = time1*time_scaling
    fig, ax = plt.subplots()
    line1, = ax.plot(time2, chemical_low/total_chemical,'b-')
    ax.set(xlabel='Time [ma]', ylabel='Chemical Percent[1]', title='Overturned Chemical')
    fig.savefig("%s/overturned_chemical.eps"%(o_dir))
    #plt.show()
    plt.close(fig)
    #----figure 3 : overturned ibc volume----#
    fig, ax = plt.subplots()
    line1, = ax.plot(time2, chemical_low/total_chemical*total_ibc,'c-')
    ax.set(xlabel='Time [ma]', ylabel='IBC Volume[1]', title='Overturned IBC')
    fig.savefig("%s/overturned_IBC.eps"%(o_dir))
    if get_end is 1:
        temp = input("input scheme to continue getting end for overturn: 'a(auto)/m(manual):")
        if temp is 'a':
            cid = fig.canvas.mpl_connect('button_press_event', on_press)
            plt.show()
            end_time = get_end_time(Gmc[0],Gmc[1],Gmc[2])
            diff_t = abs(time2-end_time)
            end_n = np.argmin(diff_t)
            end_c = chemical_low[end_n]
            episode_c = np.linspace(0.0,1.0,split)*end_c
            episode_c[0] = init_lower_percent*end_c
            episode_steps = []
            for cc in episode_c:
                diff_c = abs(chemical_low-cc)
                tempn = np.argmin(diff_c)
                episode_steps.append(step_tuple[tempn])
            end_c_percent = end_c/total_chemical
            episode_completness = episode_c/end_c
            fig.canvas.mpl_disconnect(cid)
        elif temp is 'm':
            plt.show()
            temp = input('endstep for episode:')
            end_step = int(temp)
            end_time = time_array[end_step]*time_scaling
            diff_t = abs(time2-end_time)
            end_n = np.argmin(diff_t)
            episode_t = np.linspace(0.0,1.0,split)*end_time
            end_c = chemical_low[end_n]
            end_c_percent = end_c/total_chemical
            episode_completness = episode_t/end_time
            episode_steps = []
            for tt in episode_t:
                diff_t = abs(time2-tt)
                tempn = np.argmin(diff_t)
                episode_steps.append(step_tuple[tempn])
        readf.write_to_file(pp_file,'end_time',end_time,'float')
        readf.write_to_file(pp_file,'end_c_percent',end_c_percent,'float')
        readf.write_to_file(pp_file,'episode_completness',episode_completness,'float','array')
        readf.write_to_file(pp_file,'episode_steps',episode_steps,'int','array')
    else:
        plt.show()
    plt.close(fig)
#----plot heat_flux----#
def heat_flux(in_dict,time_array,step_tuple,route,o_dir,case_name,**kwargs):

    def _config(name,default):
        try:
            return kwargs[name]
        except KeyError:
            return default

    year=365*24*3600
    write_q = readf.get_variable(in_dict,'write_q_files','int')
    Kappa = readf.get_variable(in_dict,'thermdiff','float')
    Cp = readf.get_variable(in_dict,'cp','float')
    Rho = readf.get_variable(in_dict,'density','float')
    Tref = readf.get_variable(in_dict,'reftemperature','float')
    Tsurf = readf.get_variable(in_dict,'surftemperature','float')
    Ro=readf.get_variable(in_dict,'radius','float')
    case_dir = Mpf.assign_output_dir("%s/%s"%(o_dir,case_name))
    o_dir = Mpf.assign_output_dir("%s/heat_flux"%(case_dir))
    Therm_k = Kappa*Cp*Rho
    time_scaling=pow(Ro,2.0)/Kappa/(1e6*year)
    hq_scaling = Therm_k*Tref/Ro*1e3
    #----read data----#
    q_file = "%s/%s.qb.dat"%(route,case_name)
    cols = _config('cols',5)
    q_data,q_steps = readf.read_q_file(q_file,time_array,cols)
    #for i in q_data[:,0]:
    #    print(i)#debug
    #print(q_steps) #debug
    colt = _config('colt',0)
    colq = _config('colq',1)
    colT = _config('colT',3)
    print(q_data.size)
    q_time = q_data[:,0]*time_scaling
    #print(q_data[:,0]) #debug
    q_hq = q_data[:,1]*hq_scaling
    if colT is not None:
        Tc = q_data[:,3]*Tref+Tsurf
    else:
        Tc = Tref*np.ones(q_time.size)+Tsurf
    B,qad = magnetism_intensiy(in_dict,q_hq/1e3,Tc)
    qad = qad*1e3
    B = B*1e6
    #----plot figures----#
    fig, ax = plt.subplots()
    ax = plt.subplot(2,1,1)
    #fig, ax = plt.subplots()
    plt.plot(q_time, q_hq,'r-')
    ax.set(ylabel='Heat Flux[mw/m3]')
    ax = plt.subplot(2,1,2)
    #fig, ax = plt.subplots(212)
    plt.plot(q_time, B,'b-')
    ax.set(xlabel='Time [ma]', ylabel='Magnetism [/muT]')
    #plt.show()
    fig.savefig("%s/cmb_heatflux.eps"%(o_dir))
    plt.close(fig)
#----plot sph expansion----#
def sph_expansion(in_dict,time_array,o_dir,case_name,pp_file):
    case_dir = Mpf.assign_output_dir("%s/%s"%(o_dir,case_name))
    o_dir = Mpf.assign_output_dir("%s/field"%(case_dir))
    ep_steps
    for step in ep_steps:
        sph_plot(in_dict,rr,route,o_dir,case_name,pp_file,step)
#----plot a single sph expansion----#
def sph_plot(in_dict,rr,route,o_dir,case_name,pp_file,step,vtype='T'):
    print('plot sph expansion figures, step %d\n'%step)
    color_mini = 0.05
    r_minor = 4.0
    r_mini = 2.0
    gmt_dir = Mpf.assign_output_dir("./gmt")
    file_out=('%s/%s.field.%d.dat'%(gmt_dir,case_name,step))
    case_dir = Mpf.assign_output_dir("%s/%s"%(o_dir,case_name))
    field_dir=Mpf.assign_output_dir('%s/field'%(case_dir))
    radius=readf.get_variable(in_dict,'radius','float')
    rin=readf.get_variable(in_dict,'radius_inner','float')
    nprocz=readf.get_variable(in_dict,'nprocz','int')
    noz=readf.get_variable(in_dict,'nodez','int')
    ll_max=readf.get_variable(in_dict,'output_ll_max','int')

    ofile = '%s/%s_field_%s_%d.ps'%(field_dir,case_name,vtype,step)
    ofile1 = '%s/%s_field_%s_%dl.eps'%(field_dir,case_name,vtype,step)
    lnoz=int((noz-1)/nprocz)+1
    nn_max = int((ll_max+1)*(ll_max+2)/2)
    boundary = np.zeros(noz)
    boundary[0] = 1
    boundary[noz-1] = 1
    #----read in data from .field output----#
    data = readf.read_field(route,case_name,nprocz,noz,ll_max,step)
    #----output file for gmt plot----#
    sph_max = np.zeros(noz)
    fout = open(file_out,'w')
    for nz in range(noz):
        temp = data[nz,1:ll_max+1]
        ll0 = np.argmax(temp)+1
        if boundary[nz]:
            sph_max[nz] = 0.0
            temp = temp*0.0
        else:
            sph_max[nz] = temp[ll0-1]
            temp = temp/sph_max[nz]
        for ll in range(1,ll_max+1):
            fout.write('%.4e %d %.4e\n'%(rr[nz],ll,temp[ll-1]))
    fout.close()
    #----plot spectrum output with gmt----#
    ro = radius/1e5
    ri = radius*rin/1e5
    fig,ax = plt.subplots()
    plt.plot(sph_max,rr*ro,'b-')
    ax.set(xlabel='Amplitude [1.0]', ylabel='Radius [100km]')
    plt.xlim(0.0,)
    plt.ylim(ri,ro)
    #plt.show()
    fig.savefig(ofile1)
    plt.close(fig)
    bashcommand = []
    bashcommand.append('gmt nearneighbor -i1,0[s%.1f],2 -R1/%d/%.1f/%.1f -I1.0/0.3 -S1 -N1 -G%s/result.grd -V %s'%(ro,ll_max,ri,ro,gmt_dir,file_out))
    bashcommand.append('gmt makecpt -Cno_green -T0/1/%.3f > %s/color.cpt'%(color_mini,gmt_dir))
    bashcommand.append('gmt grdimage %s/result.grd -R1/%d/%.1f/%.1f/ -C%s/color.cpt -JX6i/3i -X1.0i -Y2.0i -Ba5f1:Degree:/a%.1ff%.1f:Radius/100km:WSen -P -K > %s'%(gmt_dir,ll_max,ri,ro,gmt_dir,r_minor,r_mini,ofile))
    bashcommand.append('gmt psscale -D3.5i/-0.75i/15.0c/0.5ch -C%s/color.cpt -P -O -U -V -B:Power:0.1f0.1WSne >> %s'%(gmt_dir,ofile))
    bashcommand.append('gmt psconvert -A -Tj %s'%ofile)
#    bashcommand.append('gmt psconvert -A -Tf %s'%ofile)
    ofile2 = '%s/%s_field_%s_%d.jpg'%(field_dir,case_name,vtype,step)
    ofile3 = '%s/%s_field_%s_%dl.png'%(field_dir,case_name,vtype,step)
    bashcommand.append('convert %s -resize x1368 %s' % (ofile1,ofile3))
    bashcommand.append('convert %s %s +append %s' % (ofile2,ofile3,ofile2))
    #print(bashcommand) #debug
    #input()
    readf.run_bash(bashcommand)
    #os.remove(ofile1)
    os.remove(ofile3)
#----plot surf result----#
def surf_plot(in_dict,route,o_dir,case_name,step,**kwargs):
    def _config(name,default):
        try:
            return kwargs[name]
        except KeyError:
            return default
    gmt_dir = _config('gmt_dir','./gmt')
    rot = _config('rot',None)
    angle = _config('angle',[0.0,0.0])
    vtype = _config('vtype','total')
    cpt = _config('cpt',None)
    check_bound = _config('check_bound',False)
    cc_file = './cc'
    gravacc=readf.get_variable(in_dict,'gravacc','float')
    Rho=readf.get_variable(in_dict,'density','float')
    Ro=readf.get_variable(in_dict,'radius','float')
    Eta=readf.get_variable(in_dict,'refvisc','float')
    Kappa=readf.get_variable(in_dict,'thermdiff','float')
    if vtype is 'total':
        deg = -1
    case_dir = Mpf.assign_output_dir("%s/%s"%(o_dir,case_name))
    surf_dir = Mpf.assign_output_dir("%s/surf"%(case_dir))
    file0 = '%s/%s.surf0.%d'%(cc_file,case_name,step)
    if rot is None:
        ofile = '%s/topo_%d_%d.ps'%(surf_dir,step,deg)
        file1 = '%s/%s_%d_surf_combine'%(cc_file,case_name,step)
    else:
        ofile = '%s/topo_%d_%drt%d_%d_%d_%d.ps'%(surf_dir,step,deg,rot[0],rot[1],angle[0],angle[1])
        file1 = '%s/%s_%d_surf_combine_rt%d_%d_%d_%d'%(cc_file,case_name,step,rot[0],rot[1],angle[0],angle[1])
    scl_ang=180.0/np.pi
    if deg is -1:
        col = 2
    else:
        col = deg+2
    scl_topo = Eta*Kappa/Ro**2.0/(Rho*gravacc)

    topo_minor = 200
    sminor = 500
    smini = 100

    #----combine surf output----#
    if not os.path.isfile(file1):
        readf.write_surf_runfile(in_dict,route,case_name,step,vtype='total')
        bashcommand = []
        bashcommand.append('%s/images.x %s/runfile'%(cc_file,cc_file))
        if rot is not None:
            bashcommand.append('./lib/eula_polar %s %s %s %s %s %s'%(file0,file1,rot[0],rot[1],angle[0],angle[1]))
        readf.run_bash(bashcommand)
        if rot is None:
            os.rename(file0,file1)
    #----read head value---#
    with open(file1) as fin:
        sdata = fin.readline().split()
    data = []
    for s in sdata:
        data.append(float(s)*scl_topo)
    try:
        bound = kwargs['bound'] #enlist bound in kwargs to creat a liste
        topomm = [bound[0],bound[1]] #that stands for boundary value for
        off_topo = bound[2] #dynamic topography
    except KeyError: #else take from input
        idat = input('(min,max,offset),reference:(%.4e %.4e %.4e)\n'%(data[0],data[1],data[2]));
        idat = idat.split(',')
        topomm = []
        topomm.append(float(idat[0]))
        topomm.append(float(idat[1]))
        off_topo = float(idat[2])
    if data[0]<topomm[0] or data[1]>topomm[1]:
        print("at step %s, default value for bound is (%.4e %.4e), given is %.4e %.4e"
        %(step,data[0],data[1],topomm[0],topomm[1]))
    if check_bound:
        return data
    try:
        minor = kwargs['minor'] #enlist minor in kwargs to creat a list that
        topo_minor = minor[0] #stand for interval in color bar, interval of
        sminor = minor[1] #scale
        smini = minor[2]
    except KeyError:
        isdo = input('change (topo_minor,sminor,smini) value,default is\
            (%.4e %.4e,%.4e),y/n\n'%(topo_minor,sminor,smini))
        if isdo is 'y':
            nval = input('new value:')
            nval = nval.split(',')
            topo_minor = float(nval[0])
            sminor = float(nval[1])
            smini = float(nval[2])


    #----plot figure using gmt----#
    bashcommand = []
    bashcommand.append('gmt nearneighbor -i1[s%.2f],0[s%.2f][o90],%d[s%.2f][o%.2f] -R0/360/-90/90 -I1/1 -S10 -N1 -G%s/topo.grd -V -h1 %s'%(scl_ang,-scl_ang,col,scl_topo,off_topo,gmt_dir,file1))
    if cpt is None:
        cpt = '%s/color.cpt'%(gmt_dir)
        bashcommand.append('gmt makecpt -Cno_green -T%.2f/%.2f/%.2f > %s/color.cpt'%(topomm[0],topomm[1],topo_minor,gmt_dir))
    else:
        cpt = '%s/%s'%(gmt_dir,cpt)
    bashcommand.append('gmt grdimage %s/topo.grd -R0/360/-90/90 -C%s -JW90/6i -X0.35i -Y6.0i -B120g30/60g30 -K -P > %s'%(gmt_dir,cpt,ofile))
    bashcommand.append('gmt psscale -D3.5i/-0.5i/15.0c/0.5ch -C%s -P -O -U -V -B:Topography/m:%.2ff%.2fWSne >> %s'%(cpt,sminor,smini,ofile))
    bashcommand.append('gmt psconvert -A -Tj %s'%ofile)
    readf.run_bash(bashcommand)


#----record press position----#
def on_press(event):
    print('you pressed', event.button, event.xdata, event.ydata)
    global Gmc
    try:
        Gmc
    except NameError:
        Gmc=[]
    Gmc.append([event.xdata,event.ydata])

#----calculate overturn end----#
def get_end_time(x1,x2,x3):
    A=np.zeros((3,3))
    A[0,:] = np.array([x1[0]*x1[0],x1[0],1])
    A[1,:] = np.array([x2[0]*x2[0],x2[0],1])
    A[2,:] = np.array([x3[0]*x3[0],x3[0],1])
    b = np.array([x1[1],x2[1],x3[1]])
    x = np.linalg.solve(A,b.T)
    return -x[1]/2.0/x[0]

def magnetism_intensiy(in_dict,q_data,Tc):
    pi = np.pi
    G = 6.67e-11
    alphac = 6e-5
    rhoc = 7400
    kc = 40
    Ro = readf.get_variable(in_dict,'radius','float')
    rc = readf.get_variable(in_dict,'radius_inner','float')
    Rc = Ro*rc
    cpc = 800
    f = 1.0/7
    fohm = 1.0
    c = 0.63
    mu0 = 4*pi*1e-7
    qad = 4*pi*G*alphac*rhoc*Tc*kc*Rc/(3*cpc)
    qdiff = (q_data-qad)*(q_data-qad>0.0)
    B = f*(Rc/Ro)**(3.0)*(2*mu0*fohm*c*rhoc**(1.0/3))**(0.5)*(qdiff*qad*Rc/(kc*Tc))**(1.0/3)
    return B,qad

def remove_extra_file(route,case_name,in_dict,time_array,step_tuple,p_file,plogfile):
    caps = 12
    steps_inter = 4
    steps_after = 10
    nprocz=readf.get_variable(in_dict,'nprocz','int')
    nprocx=readf.get_variable(in_dict,'nprocx','int')
    nprocy=readf.get_variable(in_dict,'nprocy','int')
    nproc = caps*nprocx*nprocy*nprocz

    p_dict = readf.read_input(p_file)
    e_steps=readf.get_variable(p_dict,'episode_steps','int_list')
    length = len(step_tuple)
    time_at_steps = [time_array[step] for step in step_tuple]
    time0 = time_array[e_steps[0]]
    time1 = time_array[e_steps[len(e_steps)-1]]
    plog = open(plogfile,'a')
    #----get time of steps that are worth keeping----#
    time_query = []
    for n in range(1,len(e_steps)):
        stime = time_array[e_steps[n-1]]
        ftime = time_array[e_steps[n]]
        for time in np.linspace(stime,ftime,steps_inter):
            time_query.append(time)
    for time in np.linspace(time1,time_array[len(time_array)-1],steps_after):
        time_query.append(time)
    steps_in = []
    #----get steps that are worth keeping----#
    for step in step_tuple:
        if time_array[step]<time0:
            steps_in.append(step)
    for time in time_query:
        n = np.argmin(abs(time_at_steps-time))
        step = step_tuple[n]
        if step not in steps_in:
            steps_in.append(step)
    is_see = input('deleting steps: see steps to keep?(y/n)')
    if is_see:
        print(steps_in)
    input()
    #----delete extra files----#
    for n in range(len(step_tuple)):
        step = step_tuple[n]

#        bash_command = []
#        bash_command.append('rm %s/%s.proc0.%d.vts'%(route,case_name,step))
#        bash_command.append('rm %s/%s.%d.vtm'%(route,case_name,step))
        if step not in steps_in:
            for proc in range(nproc):
                filename = '%s/%s.proc%d.%d.vts'%(route,case_name,proc,step)
                if os.path.isfile(filename):
                    os.remove(filename)
            filename = '%s/%s.%d.vtm'%(route,case_name,step)
            if os.path.isfile(filename):
                os.remove(filename)
#            readf.run_bash(bash_command)
            plog.write('remove file at step %d\n'%(step))
    plog.close()
    #----write remain steps to p_file----#
    readf.write_to_file1(p_file,'steps_remain',steps_in,vtype='int',vvtype='array')

def remove_extra_file1(route,case_name,in_dict,time_array,step_tuple,step_start,p_file,plogfile):
    caps = 12
    steps_inter = 30
    nprocz=readf.get_variable(in_dict,'nprocz','int')
    nprocx=readf.get_variable(in_dict,'nprocx','int')
    nprocy=readf.get_variable(in_dict,'nprocy','int')
    nproc = caps*nprocx*nprocy*nprocz

    p_dict = readf.read_input(p_file)
    e_steps = readf.get_variable(p_dict,'episode_steps','int_list')
    steps_in = readf.get_variable(p_dict,'steps_remain','int_list')
    length = len(step_tuple)
    time_at_steps = [time_array[step] for step in step_tuple]
    step_end = step_tuple[len(step_tuple)-1]
    stime = time_array[step_start]
    ftime = time_array[step_end]
    plog = open(plogfile,'a')
    #----get time of steps that are worth keeping----#
    time_query = []
    for time in np.linspace(stime,ftime,steps_inter):
        time_query.append(time)
    #----get steps that are worth keeping----#
    for time in time_query:
        n = np.argmin(abs(time_at_steps-time))
        step = step_tuple[n]
        if step not in steps_in:
            steps_in.append(step)
    is_see = input('deleting steps: see steps to keep?(y/n)')
    if is_see:
        print(steps_in)
    input()
    #----delete extra files----#
    for n in range(len(step_tuple)):
        step = step_tuple[n]
        if step < step_start:
            continue

#        bash_command = []
#        bash_command.append('rm %s/%s.proc0.%d.vts'%(route,case_name,step))
#        bash_command.append('rm %s/%s.%d.vtm'%(route,case_name,step))
        if step not in steps_in:
            for proc in range(nproc):
                filename = '%s/%s.proc%d.%d.vts'%(route,case_name,proc,step)
                if os.path.isfile(filename):
                    os.remove(filename)
                if proc%nprocz is nprocz-1:
                    filename = '%s/%s.surf.%d.%d'%(route,case_name,proc,step)
                    if os.path.isfile(filename):
                        os.remove(filename)
                    filename = '%s/%s.surf_ori.%d.%d'%(route,case_name,proc,step)
                    if os.path.isfile(filename):
                        os.remove(filename)
            filename = '%s/%s.surf_sph.%d'%(route,case_name,step)
            if os.path.isfile(filename):
                os.remove(filename)
            filename = '%s/%s.%d.vtm'%(route,case_name,step)
            if os.path.isfile(filename):
                os.remove(filename)
#            readf.run_bash(bash_command)
            plog.write('remove file at step %d\n'%(step))
    plog.close()
    #----write remain steps to p_file----#
    readf.overwrite_to_file(p_file,'steps_remain',steps_in,vtype='int',vvtype='array')
    input()

#----plot volume results, T and melting----#
def plot_volume(in_dict,time_array,route,o_dir,case_name):
    #----parameters----#
    Ro = readf.get_variable(in_dict,'radius','float')
    Kappa = readf.get_variable(in_dict,'thermdiff','float')
    year=365*24*3600;
    time_scaling=pow(Ro,2.0)/Kappa/(1e6*year)
    melt_scaling=pow(Ro,3.0)/1e9
    v_file='%s/%s.volume.dat'%(route,case_name)
    if os.path.isfile(v_file) is False:
        print('no file %s'%(v_file))
        return
    #----read data----#
    v_data,v_steps = readf.read_q_file(v_file,time_array,6)
    v_time = v_data[:,0]*time_scaling
    v_T = v_data[:,1]
    v_melt = v_data[:,5]*melt_scaling
    #print(v_time) #debug
    #print(v_T) #debug
    #----plot figures----#
    case_dir = Mpf.assign_output_dir("%s/%s"%(o_dir,case_name))
    o_dir = Mpf.assign_output_dir("%s/vol"%(case_dir))
    fig, ax = plt.subplots()
    ax = plt.subplot(2,1,1)
    plt.plot(v_time, v_T,'r-')
    ax.set(ylabel='Temperature[1.0]')
    ax = plt.subplot(2,1,2)
    plt.plot(v_time, v_melt,'c-')
    ax.set(xlabel='Time [ma]', ylabel='Total melt volume [km^3]')
    fig.savefig("%s/volume_value.eps"%(o_dir))
    plt.close(fig)
#----another method to plot volume result----#
def plot_volume_1(in_dict,time_array,route,o_dir,case_name,settings):
    #----parameters----#
    column_dict=settings['column']
    color_dict=settings['color']
    Ro = readf.get_variable(in_dict,'radius','float')
    Kappa = readf.get_variable(in_dict,'thermdiff','float')
    year=365*24*3600;
    time_scaling=pow(Ro,2.0)/Kappa/(1e6*year)
    melt_scaling=pow(Ro,3.0)/1e9
    v_file='%s/%s.volume.dat'%(route,case_name)
    if os.path.isfile(v_file) is False:
        print('no file %s'%(v_file))
        return
    data = readf.read_data_1(v_file)
    temp = column_dict['step']
    step_array = data[:,temp]
    xdata = np.zeros(step_array.shape)
    for n in range(step_array.size):
        step = int(step_array[n])
        xdata[n] = time_array[step]*time_scaling
    #----plot figures----#
    case_dir = Mpf.assign_output_dir("%s/%s"%(o_dir,case_name))
    o_dir = Mpf.assign_output_dir("%s/vol"%(case_dir))
    fig, ax = plt.subplots()
    #----temperature----#
    ax = plt.subplot(2,1,1)
    temp = column_dict['temp']
    ydata = data[:,temp]
    color_type = color_dict['temp']
    plt.plot(xdata, ydata,color=color_type,linestyle='-',label='average T')
    ax.legend()
    ax.set(ylabel='Temperature[1.0]')
    #----melting----#
    ax = plt.subplot(2,1,2)
    temp = column_dict['melting'][0]
    color_type = color_dict['melting'][0]
    n_c = len(column_dict['melting'])-1
    ydata = data[:,temp]*melt_scaling
    plt.plot(xdata, ydata,color=color_type,linestyle='-',label='total')
    for n in range(n_c):
        name = 'chemical %d'%(n)
        temp = column_dict['melting'][n+1]
        ydata = data[:,temp]*melt_scaling
        color_type = color_dict['melting'][n+1]
        plt.plot(xdata, ydata,color=color_type,linestyle='--',label=name)
    ax.legend()
    ax.set(xlabel='Time [ma]', ylabel='Total melt volume [km^3]')
    fig.savefig("%s/volume_value.eps"%(o_dir))
    plt.close(fig)

#----mf_melting----#
def plot_mf_melting(in_dict,time_array,route,o_dir,case_name,settings):
    #----parameters----#
    column_dict=settings['column']
    color_dict=settings['color']
    Ro = readf.get_variable(in_dict,'radius','float')
    Kappa = readf.get_variable(in_dict,'thermdiff','float')
    year=365*24*3600;
    time_scaling=pow(Ro,2.0)/Kappa/year
    filename='%s/%s.mf.dat'%(route,case_name)
    if os.path.isfile(filename) is False:
        print('no file %s'%(filename))
        return
    data = np.genfromtxt(filename)
    temp = column_dict['step']
    step_array = data[:,temp]
    temp = column_dict['eular']
    eu_rate = data[:,temp]
    temp = column_dict['tracer']
    csize = len(temp)-1
    tr_rate = data[:,temp[0]]
    tr_rate_comp = np.zeros((tr_rate.size,csize))
    for i in range(csize):
        tr_rate_comp[:,i] = data[:,temp[i+1]]
    eu = np.zeros(step_array.size)
    tr = np.zeros(step_array.size)
    tr_comp = np.zeros((step_array.size,csize))
    time = np.zeros(step_array.size)
    dtime = np.zeros(step_array.size)
    for n in range(step_array.size):
        step = int(step_array[n])
        time[n] = time_array[step]*time_scaling
        if n>0:
            dtime = time[n]-time[n-1]
            eu[n] = eu[n-1]+eu_rate[n-1]*dtime
            tr[n] = tr[n-1]+tr_rate[n-1]*dtime
            tr_comp[n,:] = tr_comp[n-1,:]+tr_rate_comp[n-1,:]*dtime
    #----plot figures----#
    case_dir = Mpf.assign_output_dir("%s/%s"%(o_dir,case_name))
    o_dir = Mpf.assign_output_dir("%s/melting"%(case_dir))
    fig, ax = plt.subplots()
    #----temperature----#
    ax = plt.subplot(2,1,1)
    color_type = color_dict['tracer']
    plt.plot(time/1e6, tr_rate, color=color_type, linestyle='--',label='Tracer Method')
    color_type = color_dict['eular']
    plt.plot(time/1e6, eu_rate, color=color_type, linestyle='--',label='Eular Method')
    color_type = color_dict['melting']
    for i in range(csize):
        plt.plot(time/1e6, tr_rate_comp[:,i], color=color_type[i+1], linestyle='-',label='comp%d'%(i))
    ax.set(ylabel='Melting Producting Rate [km^3/yr]')

    ax = plt.subplot(2,1,2)
    color_type = color_dict['tracer']
    plt.plot(time/1e6, tr, color=color_type, linestyle='--',label='Tracer Method')
    color_type = color_dict['eular']
    plt.plot(time/1e6, eu, color=color_type, linestyle='--',label='Eular Method')
    color_type = color_dict['melting']
    for i in range(csize):
        plt.plot(time/1e6, tr_comp[:,i], color=color_type[i+1], linestyle='-',label='comp%d'%(i))
    ax.legend()
    ax.set(ylabel='Melting Volume[km^3/yr]')
    ax.set(xlabel='time [ma]')
    fig.tight_layout()
    fig.savefig("%s/mml_melting.eps"%(o_dir))

#----plot mf melting on surfuce-----#
class  plot_MF_surf():
    def __init__(self,cdict,nproc,nprocz,nox,noz):
        self.cdict = cdict
        self.nproc = nproc
        self.nprocz = nprocz
        self.scl_MF = 1e6/(4*math.pi*1740**2.0/(12*(nox-1)**2.0))
        print(self.scl_MF)
        self.combine_surf_file = surfc.combine_surf_file("MF",cdict,nproc,nprocz,nox,noz)
        pass
    def __call__(self,route,case_name,step,oroute):
        self.combine_surf_file(route,case_name,step,self.nprocz-1)
        self.gmt_plot(route,oroute,case_name,step)
        pass
    def gmt_plot(self,route,oroute,case_name,step,cpt=None):
    #----plot figure using gmt----#
        if os.path.isdir(oroute) is False:
            os.path.mkdir(oroute)
        mfroute = os.path.join(oroute,case_name,"melting")
        if os.path.isdir(mfroute) is False:
            os.path.mkdir(mfroute)
        filename = os.path.join(route,"%s.MF0.%d"%(case_name,step))
        grdname = os.path.join(mfroute,"MF_%d.grd"%(step))
        outname = os.path.join(mfroute,"MF_%d"%(step))
        scl_ang=180.0/np.pi
        col = 2
        mm = np.zeros(2)
        with open(filename,'r') as f:
            inputs = f.readline().split()
        mm[0] = float(inputs[0])*self.scl_MF
        mm[1] = float(inputs[1])*self.scl_MF
        minor = (mm[1]-mm[0])/10.0
        sminor = minor*2.5
        smini = sminor/5.0
        bashcommand = []
        bashcommand.append('gmt nearneighbor -i1[s%.2f][o180],0[s%.2f][o90],%d[s%.2f] -R0/360/-90/90 -I1/1 -S10 -N1 -G%s -V -h1 %s'%(scl_ang,-scl_ang,col,self.scl_MF,grdname,filename))
        if cpt is None:
            cpt = os.path.join(mfroute,"color.cpt")
            bashcommand.append('gmt makecpt -Cno_green -T%.2f/%.2f/%.2f > %s/color.cpt'%(mm[0],mm[1],minor,mfroute))
        else:
            cpt = os.path.join(mfroute,cpt)
        bashcommand.append('gmt grdimage %s -R0/360/-90/90 -C%s -JW0/6i -X0.35i -Y6.0i -B120g30/60g30 -K -P > %s'%(grdname,cpt,outname))
        bashcommand.append('gmt psscale -D3.5i/-0.5i/15.0c/0.5ch -C%s -P -O -U -V -B:Melting_Rate/km_per_Myr:%.2ff%.2fWSne >> %s'%(cpt,sminor,smini,outname))
        bashcommand.append('gmt psconvert -A -Tj %s'%outname)
        print(bashcommand) #debug
        readf.run_bash(bashcommand)
        pass

#----plot sph expansion----#

def hoz_combine(in_dict,route,cname,ta,steps,hdict):
    fname='%s_hoz_combine'%(cname)
    dir0 = './gmt'
    temp = hdict['temp']
    rho = hdict['rho']
    col = hdict['col']
    nprocz = readf.get_variable(in_dict,'nprocz','int')
    noz=readf.get_variable(in_dict,'nodez','int')
    lnoz=int((noz-1)/nprocz)+1
    foname = '%s/%s'%(dir0,fname)
    if os.path.isfile(foname):
        os.remove(foname)
    for i in range(len(steps)):
        step = steps[i]
        t = ta[step]
        for pz in range(nprocz):
            filename = '%s/%s.horiz_avg.%d.%d'%(route,cname,pz,step)
            data = readf.read_data(filename)
            with open(foname,'a') as fo:
                for lnz in range(lnoz):
                    if pz>=1 and lnz==0:
                        continue
                    nn = lnz*col
                    fo.write('%.4e %.4e %.4e %.4e\n'%(t,data[nn],data[nn+temp],data[nn+rho]))
    return fname

def plot_combine(in_dict,oroute,cname,dfile,tp,cdict,tval,**kwarg):

    def _config(name,default):
        try:
            return kwarg[name]
        except KeyError:
            return default

    dcwd = './gmt'
    Ro=readf.get_variable(in_dict,'radius','float')
    Kappa=readf.get_variable(in_dict,'thermdiff','float')
    ri=readf.get_variable(in_dict,'radius_inner','float')
    Buoy = readf.get_variable(in_dict,'buoyancy_ratio','float')
    Alpha = readf.get_variable(in_dict,'thermexp','float')
    year = 365*24*3600;
    Time = Ro**2.0/Kappa/(1e6*year)
    R = Ro/1e5
    Temp=readf.get_variable(in_dict,'reftemperature','float')
    Toff=readf.get_variable(in_dict,'surftemperature','float')
    Tab = 273.15
    Toff = Toff-Tab
    Rho = readf.get_variable(in_dict,'density','float')
    Drho = Buoy*Alpha*Temp*Rho
    tmin = tval[0]
    tmax = tval[1]
    tminor = tval[2]
    tmini = tval[3]
    rmin = R*ri
    rmax = R
    odir0 = Mpf.assign_output_dir('%s/%s'%(oroute,cname))
    odir = Mpf.assign_output_dir('%s/combine'%(odir0))
    ofile = '%s/combine_%s_%.0f_%.0f.ps'%(odir,tp,tmin,tmax)
    if tp is 'rho':
        Dat = Drho
        Doff = Rho
        dmin = Rho
        dmax = Rho+Drho
        dminor = _config('rhodminor',5.0)
        coff = _config('rhocoff',50.0)
    elif tp is 'temp':
        Dat = Temp
        Doff = Toff
        dmin = Toff
        dmax = Toff+Temp
        dminor = 50.0
        coff = 500.0
    t = cdict['time']
    r = cdict['radius']
    dat = cdict[tp]
    abash = []
    abash.append('gmt nearneighbor -i%d[s%.4f],%d[s%.4f],%d[s%.4fo%.4f] -R%.4f/%.4f/%.4f/%.4f -I0.1/0.1 -S30 -N1 -Gresult.grd -V %s'%(t,Time,r,R,dat,Dat,Doff,tmin,tmax,rmin,rmax,dfile))
    cpt = _config('cpt',None)
    if cpt is None:
        abash.append('gmt makecpt -Cno_green -T%.4f/%.4f/%.4f > color.cpt'%(dmin,dmax,dminor))
        cfile = 'color.cpt'
    else:
        cfile = cpt
    abash.append('gmt grdimage result.grd -R%.4f/%.4f/%.4f/%.4f -C%s -JX4i/4i -Ba%.4ff%.4f:Time/Ma:/a4.0f2.0:Radius/100km:WSen -K -P > %s'%(tmin,tmax,rmin,rmax,cfile,tminor,tmini,ofile))
    abash.append('gmt psscale -C%s -D4.4i/2.0i/4.0i/0.2i -P -O -V -B%.4f >> %s'%(cfile,coff,ofile))
    abash.append('gmt psconvert -A -Tj %s'%ofile)
    print(abash) #debug
    readf.run_bash(abash,dcwd)

def make_color(infile,ofile,total,arr):
    col = 4
    fout = open(ofile,'w')
    with open(infile,'r') as fin:
        i = 0
        while 1:
            line = fin.readline()
            if line is '':
                break
            dat = line.split()
            if i <= total-2:
                fout.write('%d %s %d %s\n'%(arr[i],dat[1],arr[i+1],dat[3]))
                i = i+1
            else:
                fout.write(line)
    fout.close()



