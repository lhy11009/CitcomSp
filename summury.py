import os
import sys
sys.path.append('../public')
import imgPP #self writen module for dealing images
import dataPP #self writen module dealing with CitcomS output files
import readfiles as Mread

route = "/home/p300/Documents/Research/work/CitcomS2/MgH30V5e20E100L"
#dataRoute = "/home/p300/Desktop/lunar_project/m1t_1"
dataRoute = "/media/p300/Seagate0/mg_project/MgH30V5e20E100L"
caseName = "MgH30V5e20E100L"
pp_file = '%s/post_process'%(dataRoute)
use_eps = 1
pp_dict = Mread.read_input(pp_file)
PVarray = []
MLarray = []
SGarray = []
steparray = []
#image types from paraview goes here
#PVarray.append('volume')
#PVarray.append('slice')
#PVarray.append('composition')
#PVarray.append('cslice0')
PVarray.append('cslice2')
PVarray.append('melting')
#image types from matlab goes here
MLarray.append('temp')
MLarray.append('rho')
MLarray.append('visc')
MLarray.append('velo')
include_field = 1
#image that is singular (not base on step) goes here
SGarray.append('chemical/overturned_IBC')
SGarray.append('chemical/retained_ibc')
SGarray.append('heat_flux/cmb_heatflux')
SGarray.append('vol/volume_value')
#steps of interest
steparray.append(0)
e_steps = Mread.get_variable(pp_dict,'episode_steps','int_list')
for step in e_steps:
    if step > 0:
        steparray.append(step)
try:
    rsteps = Mread.get_variable(pp_dict,'steps_remain','int_list')
except KeyError:
    map_to_index = 1
else:
    map_to_index = 0
    steparray.append(rsteps[-1])

picN = 0
fileO = '%s/summury.jpg' % (route)
maxStepPv = 100000 # max step for searching exitsting output file from CitcomS
#----map step to index of 3d figures----#
if map_to_index:
    stepIndexPv = dataPP.map_steps_to_index_pv(dataRoute,caseName,maxStepPv) #get index for paraview output of a specified step
else:
    stepIndexPv = {}
    for n in range(len(rsteps)):
        stepIndexPv[rsteps[n]] = n
#print(stepIndexPv)
#input()
for step in steparray:
#convert images within specified steps horizontally and 
#convert results from different steps vertically 
    picNN = 0
    pvIndex = stepIndexPv[step]
    fileOO = '%06d.jpg' % (step)
    for typename in PVarray:
    #convert images form paraview
        filein = '%s/%s/%s.%04d.jpg' % (route,typename,typename,pvIndex) #careful for paraview outputs, that's why I use this Index for filenames
        print("Reading image: %s" %(filein))
        imgPP.convert_image(filein,fileOO,picNN,'h')
        picNN += 1
    for typename in MLarray:
    #convert images form matlab
        if use_eps:
            filein = '%s/%s/%s_%06d.eps' % (route,typename,typename,step)
        else:
            filein = '%s/%s/%s_%06d.png' % (route,typename,typename,step)
        print("Reading image: %s" %(filein))
        fileinin = '%s_%06d.png' % (typename,step)
        imgPP.convert_image(filein,fileinin,0,'h',1148,872)
        imgPP.convert_image(fileinin,fileOO,picNN,'h')
        os.remove(fileinin)
        picNN += 1   
    if include_field and step is not 0:
        filein = '%s/field/%s_field_T_%d.jpg' % (route,caseName,step)
        print("Reading image: %s" %(filein))
        fileinin = 'field_%d.png' % (step)
        imgPP.convert_image(filein,fileinin,0,'h',2600,872)
        imgPP.convert_image(fileinin,fileOO,picNN,'h')
        os.remove(fileinin)
        picNN += 1   
    imgPP.convert_image(fileOO,fileO,picN,'v')
    os.remove(fileOO)
    picN += 1
#convert singular file
picNN = 0
fileOO = 'Simage.jpg'
for fileN in SGarray:
    if use_eps:
        filein = '%s/%s.eps' % (route,fileN)
    else:
        filein = '%s/%s.png' % (route,fileN)
    print("Reading image: %s" %(filein))
    fileinin = 'Simage%d.jpg' % (picNN)
    imgPP.convert_image(filein,fileinin,0,'h',1148,872)
    imgPP.convert_image(fileinin,fileOO,picNN,'h')
    os.remove(fileinin)
    picNN += 1
imgPP.convert_image(fileOO,fileO,picN,'v')
os.remove(fileOO)
picN += 1

    

#bashCommand = 'convert %s %s +append ./test.jpg' % (filename,filename1)
#process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
#output, error = process.communicate()


