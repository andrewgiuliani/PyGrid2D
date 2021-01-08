import numpy as np
import ipdb

def output_ply(vertices, cells, cells_ij, domain, Nx, Ny):
    f = open(domain.name+"_"+str(Nx)+"_"+str(Ny)+".ply", 'w')
    f.write("ply\n")
    f.write("format ascii 1.0\n")
    f.write("comment this file is a " +str(Nx)+"x"+str(Ny)+" " + domain.name + "\n")
    f.write("element vertex " + str(vertices.shape[0]) + "\n")
    f.write("property float x\n")
    f.write("property float y\n")
    tot_count = 0
    for c in cells:
        tot_count = tot_count + c.shape[0]
        
    f.write("element face " + str(tot_count) +"\n")
    f.write("end_header\n")
    np.savetxt(f,vertices)
    for c in cells:
        np.savetxt(f,c,fmt='%i')

    f.close()

