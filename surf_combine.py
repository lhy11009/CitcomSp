import sys
import os
sys.path.append('../public')
import numpy as np
import readfiles as readf
import mathf as mathf


def bound_points(coord0,noz):
    for row in range(coord0.shape[0]):
        if row%noz == 0:
            nxy = row//noz
            coord1[nxy,:] = coord0[row,:]
    return coord1

def bound_elements(coord0,nox,noz):
    coord1 = np.zeros(((nox-1)**2,3))
    n0 = np.zeros(4)
    for i in range(nox-1):
        for j in range(nox-1):
            n0[0] = i*noz*nox+j*noz
            n0[1] = (i+1)*noz*nox+j*noz
            n0[2] = i*noz*nox+(j+1)*noz
            n0[3] = (i+1)*noz*nox+(j+1)*noz
            n = i*(nox-1)+j
            for k in n0:
                coord1[n,:] = coord1[n,:]+coord0[int(k),:]
            coord1[n,:] = coord1[n,:]/4.0
    return coord1

def line_prepender(filename, line):
    with open(filename, 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        f.write(line.rstrip('\r\n') + '\n' + content)

class combine_surf_file():
    
    def __init__(self,fname,cdict,nproc,nprocz,nox,noz):
        self.fname=fname
        self.cdict=cdict
        self.nproc=nproc
        self.nprocz=nprocz
        self.nox=nox
        self.noz=noz
        pass

    def __call__(self,route,case_name,step,procz):
        nprocxy = self.nproc//self.nprocz
        out_filename=os.path.join(route,"%s.MF0.%d"%(case_name,step))
        mm0 = np.zeros(2)
        if os.path.isfile(out_filename):
            os.remove(out_filename)
        for i in range(nprocxy):
            MF0,coord0 = self.read(route,case_name,step,i)
            coord1 = bound_elements(coord0,self.nox,self.noz)
            #print(coord1.shape[0])
            #print(coord1)
            #input()
            mm1 = self.maxmin(MF0,0)
            mm0[0] = min(mm0[0],mm1[0])
            mm0[1] = max(mm0[1],mm1[1])
            #print(coord1) #debug
            coord1 = mathf.cart2sph(coord1)
            #print(coord1)
            #input()
            self.output(out_filename,MF0,coord1)
            #print(MF0.size) #debug
            #print(MF0)
#            input()
        #print(mm0)
        line_prepender(out_filename,"%s %s"%(mm0[0],mm0[1]))
        pass 

    def read(self,route,case_name,step,procxy,flag="add"):
    #read in MF file, flag is unused#
        proc = procxy*self.nprocz
        filename=os.path.join(route,"%s.MF.%d.%d"%(case_name,proc,step))
        MF0 = np.genfromtxt(filename)
        for j in range(1,self.nprocz):
            proc = procxy*self.nprocz+j
            filename=os.path.join(route,"%s.proc%d.%d.vts"%(case_name,proc,step))
            ddict,coord=readf.read_vtk_file(filename)
            filename=os.path.join(route,"%s.MF.%d.%d"%(case_name,proc,step))
            MF1 = np.genfromtxt(filename)
            MF0 = MF0 + MF1
        return MF0,coord

    def output_header(selt,filename,mm):
        with open(filename,'w') as f:
            f.write("%.4e %.4e\n"%(mm[0],mm[1]))

    def output(self,filename,MF0,coord):
    #output to a MF0 file for gmt#
        cshape = coord.shape
        Mshape = MF0.shape
        with open(filename,'a') as f:
            for i in range(Mshape[0]):
                #print(coord[i,:]) #debug
                #input()
                f.write("%.4e"%coord[i,1])
                for j in range(2,cshape[1]):
                    f.write(" %.4e"%coord[i,j])
                for j in range(Mshape[1]):
                    f.write(" %.4e"%MF0[i,j])
                f.write("\n")

    def maxmin(self,MF,column):
        mm = np.zeros(2)
        mm[0] = np.min(MF[:,column])
        mm[1] = np.max(MF[:,column])
        return mm
