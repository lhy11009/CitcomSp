import numpy as np

def cart2sph(xyz):
    pts = np.zeros((xyz.shape[0],3))
    xy = xyz[:,0]**2 + xyz[:,1]**2
    pts[:,0] = np.sqrt(xy + xyz[:,2]**2)
    pts[:,1] = np.arctan2(np.sqrt(xy), xyz[:,2]) # for elevation angle defined from Z-axis down
    pts[:,2] = np.arctan2(xyz[:,1], xyz[:,0])
    return pts

def eula_trans(xyz,pts):
    xyz1 = 0
    return xyz1
    pass
