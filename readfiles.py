import subprocess 
import numpy as Mnp
import os
import sys

#----read in radius----#
def read_radius(in_dict,route,casename):
    nprocz = get_variable(in_dict,'nprocz','int')
    noz = get_variable(in_dict,'nodez','int')
    lnoz=int((noz-1)/nprocz)+1
    data = read_horig_output(route,casename,0,nprocz,lnoz,1)
    return data[:,0]

#----read in horizontal average result---#
def read_horig_output(Croute,Ccasename,Cstep,Cnprocz,Clnoz,cols):
    noz = (Clnoz-1)*Cnprocz+1
    Rdata = Mnp.zeros((noz,cols))
    for pz in range(Cnprocz):
        filename = "%s/%s.horiz_avg.%d.%d"%(Croute,Ccasename,pz,Cstep)
        fin = open(filename,'r')
        offset = (Clnoz-1)*pz
        nl = 0
        for line in fin:
            temp = line.split(' ')
            for nc in range(cols):
                Rdata[nl+offset][nc] = float(temp[nc])
            nl = nl+1
        fin.close()
    return Rdata

#----read in time----#
def read_time(Croute,Ccasename):
    filename = "%s/%s.time"%(Croute,Ccasename)
    fin = open(filename,'r')
    step = list()
    time = list()
    for line in fin:
        temp = line.split(' ')
        step.append(int(temp[0]))
        time.append(float(temp[1]))
    #----drop duplicate data----#
    fstep = step[len(step)-1]
    Rdata = Mnp.zeros(fstep+1)
    temp = fstep+1
    for n in range(len(step)-1,-1,-1):
        n_s = step[n]
        if n_s<temp:
            Rdata[n_s] = time[n]
            temp = n_s
    fin.close();
    return Rdata
#----read in q_file----#
def read_q_file(filename,time_array,column,colt=0):
    tiny = 1e-4
    ttiny = 1e-9
    #----read data from qfile----#
    fin = open(filename,'r')
    data = read_data(filename)
    lines = len(data)//column
    #----eliminate overlapping data----#
    step = list()
    ilist = list()
    time_after = 1.1*data[column*(lines-1)+colt]
    for i in range(lines-1,-1,-1):
        time = data[column*i+colt]
        if time<time_after:
            ilist.append(i)
            time_after = time
            if time<ttiny:
                step.append(0)
                continue
            for n in range(len(time_array)):
                if abs(time-time_array[n])/time<tiny:
                    step.append(n)
                    break
    #----generate data matrix for return----#
    step = step[::-1]
    Rdata = Mnp.zeros((len(ilist),column))
    for n in range(len(ilist)-1,-1,-1):
        i = ilist[n]
        for j in range(column):
            Rdata[n,j] = data[column*i+j]
    return Rdata,step
#----read input_file----#
def read_input(Cfilename):
    Rdict = dict()
    fin = open(Cfilename,'r')
    for line in fin:
        #----elimate space at head----#
        for i in range(len(line)):
            if line[i]!=' ' and line[i]!='\t':
                break
        line=line[i:len(line)] 
        #----skip some lines----#
        if (line is '\n') or (line[0] is '#'):
            continue
        #----elimate content after '#'----#
        for i in range(len(line)):
            if line[i] is '#':
                break
        line=line[0:i] 
        temp = line.split('=')
        if temp[1][0] is '"':
            temp[1] = temp[1][1:len(temp[1])-1]
        #print(line) #debug
        Rdict[temp[0]]=temp[1]
    fin.close()
    return Rdict

#----get value of a variable from an existing dictionary----#
def get_variable(Cdict,Cname,Ctype):
    try:
        temp=Cdict[Cname]
    except NameError:
        print('no name %s in dictionary\n',Cname)
        os.exit()
    #print(temp) #debug
    if Ctype is 'string':
        Rvalue=temp
    elif Ctype is 'int':
        Rvalue=int(temp)
    elif Ctype is 'float':
        Rvalue=float(temp)
    elif Ctype is 'int_list':
        Rvalue=temp.split(',')
        for i in range(len(Rvalue)):
            Rvalue[i] = int(Rvalue[i])
    elif Ctype is 'float_list':
        Rvalue=temp.split(',')
        for i in range(len(Rvalue)):
            Rvalue[i] = float(Rvalue[i])
    return Rvalue

def read_field(route,case_name,nprocz,noz,ll_max,step):
    header = 5
    inter = 2
    columns = 6
    pcolumns = 2
    Dname = ('T')
    f_dict = {'T':[2,3],'Vr':[4,5]} #columns, starts from 0
    lnoz=int((noz-1)/nprocz)+1
    nn_max = int((ll_max+1)*(ll_max+2)/2)
    Rdata = Mnp.zeros((noz,ll_max+1))
    for procz in range(nprocz):
        f_file = "%s/%s.field.%d.%d"%(route,case_name,procz,step)
        data = read_data(f_file)
        for i in range(lnoz):
            nz = i+(lnoz-1)*procz
            Atemp = Mnp.zeros(ll_max+1)
            for ll in range(ll_max+1):
                for mm in range(ll+1):
                    pp = int(ll*(ll+1)/2+mm)
                    for name in Dname:
                        idx0 = header+(i*nn_max+pp)*columns+inter*i+f_dict[name][0]
                        idx1 = header+(i*nn_max+pp)*columns+inter*i+f_dict[name][1]
                        Atemp[ll] = Atemp[ll]+data[idx0]**2.0+data[idx1]**2.0
                Atemp[ll] = Atemp[ll]**0.5
            Rdata[nz,:] = Atemp
    return Rdata

def write_surf_runfile(in_dict,route,case_name,step,vtype='total'):
    #----patameters----#
    ofile = './cc/runfile'
    caps = 12 
    nox = get_variable(in_dict,'nodex','int')
    noy = get_variable(in_dict,'nodey','int')
    noz = get_variable(in_dict,'nodez','int')
    nprocx = get_variable(in_dict,'nprocx','int')
    nprocy = get_variable(in_dict,'nprocy','int')
    nprocz = get_variable(in_dict,'nprocz','int')
    inputf0="%s/%s"%(route,case_name) #
    outputf0="./cc/%s"%(case_name) # output result in the same as input
    cpu_xy = caps*nprocx*nprocz
    cpu_z = nprocz
    file_type=1
    get_deltaT=3
    comp_temp=1
    get_slab_center=40
    slab_t=-0.1
    special_value=1.0
    #----output to runfile----#
    with open(ofile,'w') as fout:
        fout.write('%s\n'%(inputf0))
        fout.write('%s\n'%(inputf0))
        fout.write('%s\n'%(outputf0))
        fout.write('%d %d %d\n'%(nox,noy,noz))
        fout.write('%d %d\n'%(cpu_xy,cpu_z))
        fout.write('%d %d %d 0 %d\n'%(file_type,get_deltaT,step,comp_temp))
        fout.write('%d %.4f 0.0 %.4f\n'%(get_slab_center,slab_t,special_value))
        fout.write('-1')


#----read from a pure data file----#
def read_data(filename):
    data=[]
    with open(filename) as f:
        for line in f:
            line.strip()
            for part in line.split():
                if part is not ' ':
                    data.append(float(part))
    return data
#----read data and return a matrix----#
def read_data_1(filename):
    with open(filename,'r') as fin:
        line=fin.readline()
        s_list=line.split(' ')
    col=0
    for s in s_list:
        if s == ' ' or s == '\n' or s == '':
            pass
        else:
            col=col+1
    in_data=read_data(filename)
    row=len(in_data)//col
    data=Mnp.array(in_data)
    data.resize((row,col))
    return data
#----get steps by finding horiz file----#

def get_steps(Croute,Ccase_name,Cmax_step):
    Rtuple = ()
    for step in range(Cmax_step):
        filename = "%s/%s.horiz_avg.0.%d"%(Croute,Ccase_name,step)
        if os.path.isfile(filename):
            Rtuple=Rtuple+(step,)
    return Rtuple

#----write to file when file or variable doesn't exit----#
def write_to_file1(pp_file,name,value,vtype='int',vvtype='simple'):
    if os.path.isfile(pp_file): 
        pp_dict = read_input(pp_file)
        try:
            pp_dict[name]
        except KeyError:
            write = 1
        else:
            write = 0
    else:
        write = 1
    if write:
        write_to_file(pp_file,name,value,vtype,vvtype)

def overwrite_to_file(pp_file,name,value,vtype='int',vvtype='simple'):
    temp_file = '%s_temp'%(pp_file)
    if os.path.isfile(pp_file): 
        pp_dict = read_input(pp_file)
    else:
        pp_dict = []
    for key in pp_dict:
        if key == name:
            write_to_file(temp_file,name,value,vtype,vvtype)
        else:
            try:
                fp = open(temp_file,'a')
            except IOError:
                print('cannot open %s for writing',pp_file)
            else:
                fp.write('%s=%s\n'%(key,pp_dict[key]))
                fp.close()
    os.rename(temp_file,pp_file)
                
#----write to post-process file----#
def write_to_file(pp_file,name,value,vtype='int',vvtype='simple',vvvtype='value'):
    with open(pp_file,'a') as fout:
        if vtype is 'int':
            s='%d'
        elif vtype is 'float':
            s='%.4e'
        elif vtype is 'str':
            s='%s'
        if vvtype is 'array':
            string = ''
            for i in range(len(value)):
                string = string+s%(value[i])
                if i is not len(value)-1:
                    string = string+','
        else:
            string = s%value
        if vvvtype is 'quote':
            string = '\''+string+'\''
        fout.write("%s="%(name)+string+'\n')

#----run bash file----#
def run_bash(bashcommand,dcwd=None,shll=False):
    #----output output to a specific file----#
    for command in bashcommand:
        if shll is True:
            process = subprocess.Popen(command, stdout=subprocess.PIPE, cwd=dcwd,shell=True)
        else:
            process = subprocess.Popen(command.split(), stdout=subprocess.PIPE,cwd=dcwd)
        output, error = process.communicate()
    return

#----read data from vtk files, return as data and coordinate as numpy arrays----#
def read_vtk_file(filename):
    def _data_type(line):
        d_type={}
        d_list=line.split(' ')
        name=None
        comp=1
        for string in d_list:
            index=string.find('=')
            if index > -1:
                head=string[:index]
                value=string[index+2:len(string)-1]
                if head == 'Name':
                    name=value
                elif head == 'NumberOfComponents':
                    comp=int(value)
        return name,comp
    d_dict={}
    with open(filename,'r') as fin:
        for i in range(5):
            line=fin.readline()
        while 1:
            line=fin.readline()
            name,comp=_data_type(line)
            if name==None:
                break
            data=[]
            while 1:
                line=fin.readline()
                if line.find('<') > -1:
                    break
                d_list=line.split(' ')
                length=len(d_list)
                for j in range(length):
                    data.append(float(d_list[j]))
            data=Mnp.array(data)
            row=len(data)//comp
            data.resize((row,comp))
            d_dict[name]=data
        #read coordinate#
        data=[]
        comp=3
        for i in range(4):
            line=fin.readline()
        while 1:
            line=fin.readline()
            if line.find('<') > -1:
                break
            d_list=line.split(' ')
            length=len(d_list)
            for j in range(length):
                data.append(float(d_list[j]))
        data=Mnp.array(data)
        row=len(data)//comp
        data.resize((row,comp))
        coord=data
    return d_dict,coord
#----write vtk file from data and coord, should be numpy arrays----#
def write_vtk_file(filename,o_data,o_coord,step,proc,nx,ny,nz):
    with open(filename,'w') as fin:
        #vtk head
        fin.write('<?xml version=\"1.0\"?>\n')
        fin.write('<VTKFile type=\"StructuredGrid\" version=\"0.1\" \
compressor=\"vtkZLibDataCompressor\" byte_order=\"LittleEndian\">\n')
        fin.write('  <StructuredGrid WholeExtent="%d %d %d %d %d %d">\n'\
                %(nz[0],nz[1],nx[0],nx[1],ny[0],ny[1]))
        fin.write('    <Piece Extent="%d %d %d %d %d %d">\n'\
                %(nz[0],nz[1],nx[0],nx[1],ny[0],ny[1]))
        fin.write('      <PointData Scalars=\"temperature\" Vectors=\"velocity\">\n')
        #vtk data
        for key in o_data:
            data=o_data[key]
            if data.shape[1]>1:
                comp=data.shape[1]
                fin.write('        <DataArray type=\"Float32\" Name=\"%s\"\
 NumberOfComponents=\"%d\" format=\"ascii\">\n'%(key,comp))
            else:
                fin.write('        <DataArray type=\"Float32\" Name=\"%s\" format=\"ascii\">\n'%(key))
            for i in range(data.shape[0]):
                for j in range(data.shape[1]):
                    if j == data.shape[1]-1:
                        fin.write('%.4e\n'%(data[i,j]))
                    else:
                        fin.write('%.4e '%(data[i,j]))
            fin.write('        </DataArray>\n')
        fin.write('      </PointData>\n')
        fin.write('      <CellData>\n')
        fin.write('      </CellData>\n')
        fin.write('      <Points>\n')
        #vtk coordinate
        fin.write('        <DataArray type=\"Float32\" Name=\"coordinate\"\
 NumberOfComponents=\"3\" format=\"ascii\">\n') 
        for i in range(o_coord.shape[0]):
            for j in range(o_coord.shape[1]):
                if j == o_coord.shape[1]-1:
                    fin.write('%.4e\n'%(o_coord[i,j]))
                else:
                    fin.write('%.4e '%(o_coord[i,j]))
        fin.write('        </DataArray>\n')
        fin.write('     </Points>\n')
        fin.write('    </Piece>\n')
        fin.write('  </StructuredGrid>\n')
        fin.write('</VTKFile>\n')

    
