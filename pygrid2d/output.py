import numpy as np

def output_nor(vertices, mesh_data, face_data, domain, Nx, Ny, q):
    assert q == 1

    f = open(domain.name+"_"+str(Nx)+"_"+str(Ny)+"_q" + str(q)+".nor", 'w')
    num_faces = len(face_data)
    nor = np.zeros((num_faces,3))
    count = 0
    for face in face_data:
        v = face['vertices']
        x = np.zeros((2,))
        y = np.zeros((2,))
        for i in range(2):
            x[i] = vertices[v[i],0]
            y[i] = vertices[v[i],1]
        dx = x[1]-x[0]
        dy = y[1]-y[0]
        nn = np.sqrt(dx**2 + dy**2)
        nor[count,0] = dy/nn
        nor[count,1] = -dx/nn
        count+=1
    np.savetxt(f,nor)
    f.close()
    
def output_ply(vertices, cells, domain, Nx, Ny, q):
    f = open(domain.name+"_"+str(Nx)+"_"+str(Ny)+"_q" + str(q)+".ply", 'w')
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

def output_unstructured(vertices, mesh_data, face_data, domain, Nx, Ny, q):
    tot_elem = len(mesh_data)
    tot_face = len(face_data)
    
    f = open(domain.name+"_"+str(Nx)+"_"+str(Ny)+"_q" + str(q)+".unstr", 'w')
    f.write("vertices " + str(vertices.shape[0]) + "\n")
    np.savetxt(f,vertices)
    
    f.write("cells " + str(tot_elem) + "\n")
    for c in mesh_data:
        cdata = np.concatenate( ([c['cut']], c['fidx']) )
        np.savetxt(f, cdata.reshape((1,-1)), fmt='%i')

    f.write("faces " + str(tot_face) + "\n")
    for face in face_data:
        fdata = np.concatenate( ([face['cut']], face['vertices']) )
        np.savetxt(f, fdata.reshape((1,-1)), fmt='%i')

    f.write("faces left right " + str(tot_face) + "\n")
    for face in face_data:
        np.savetxt(f, np.array(face['lr']).reshape((1,-1)), fmt='%i')



def output_ag(vertices, cells, cells_ij, nv, ncf, domain, Nx, Ny, q, vertex_idx, vert_in):
    
    irr = -np.ones( (Nx, Ny) ).astype(int)
    irr[cells_ij[0]] = 1
    
    irr_num = 2
    for c in range(1, len(cells_ij)):
        if len(cells_ij[c]) > 0:
            irr[cells_ij[c]] = np.arange( cells_ij[c][0].size ) + irr_num
            irr_num = irr_num + cells_ij[c][0].size

    f = open(domain.name + "_" + str(Nx) + "x" + str(Ny)+"_q"+str(q)+".dat","w")
    f.write(str(Nx) + " " + str(Ny) + "\n")
    f.write("{:.16E} {:.16E} {:.16E} {:.16E}\n".format(domain.left, domain.bottom, domain.right, domain.top)) 
    f.write("{:.16E} {:.16E}\n".format(  (domain.right - domain.left) / Nx, (domain.top - domain.bottom) / Ny  ) )
    np.savetxt(f, irr, fmt = '%i')
    f.write("{}\n{}\n{}\n".format(4,0,q) )
    f.write("{}\n".format( 0 )  )
    f.write("{} {}\n".format( -1,-1 ))
    dx = (domain.right - domain.left) / Nx
    dy = (domain.top - domain.bottom) / Ny
    np.savetxt(f, np.array([ [0,0], [0,dy], [dx,dy],[dx,0] ]) , fmt = '%.16E')  
    
    # compute side data
    hb_list = [None]*len(cells)
    ht_list = [None]*len(cells)
    vl_list = [None]*len(cells)
    vr_list = [None]*len(cells)
    tot_cell_count = 0
    for cid1 in range( 1, len(cells) ):
        c = cells[cid1]
        v = nv[cid1]
        tot_cell_count = tot_cell_count + c.shape[0]
        if c.shape[0] == 0:
            continue

        ij = cells_ij[cid1]
        
        v1 = ij
        v2 = (v1[0]+1, v1[1]  )
        v3 = (v1[0]+1, v1[1]+1)
        v4 = (v1[0]  , v1[1]+1)

        
        where_v1 = np.where( np.logical_and(vertex_idx[v1][:,None] == c, vert_in[v1][:,None]) )
        where_v2 = np.where( np.logical_and(vertex_idx[v2][:,None] == c, vert_in[v2][:,None]) )
        where_v3 = np.where( np.logical_and(vertex_idx[v3][:,None] == c, vert_in[v3][:,None]) )
        where_v4 = np.where( np.logical_and(vertex_idx[v4][:,None] == c, vert_in[v4][:,None]) )

        idx_v1 = -np.ones(c.shape[0]).astype(int)  
        idx_v2 = -np.ones(c.shape[0]).astype(int)  
        idx_v3 = -np.ones(c.shape[0]).astype(int)  
        idx_v4 = -np.ones(c.shape[0]).astype(int)  

        idx_v1[where_v1[0]] = where_v1[1] 
        idx_v2[where_v2[0]] = where_v2[1] 
        idx_v3[where_v3[0]] = where_v3[1] 
        idx_v4[where_v4[0]] = where_v4[1] 

        hb = np.where( np.logical_and(idx_v1 > -1, idx_v2 > -1)  , idx_v1,   -1)  
        hb = np.where( np.logical_and(idx_v1 > -1, idx_v2 == -1) , idx_v1,   hb)  
        hb = np.where( np.logical_and(idx_v1 == -1, idx_v2 > -1) , idx_v2-1, hb) 
        
        ht = np.where( np.logical_and(idx_v3 > -1,  idx_v4 > -1)  , idx_v3,   -1)  
        ht = np.where( np.logical_and(idx_v3 > -1,  idx_v4 == -1) , idx_v3,   ht)  
        ht = np.where( np.logical_and(idx_v3 == -1, idx_v4 > -1)  , idx_v4-1, ht) 
 

        vl = np.where( np.logical_and(idx_v4 > -1,  idx_v1 > -1)  , idx_v4,   -1)  
        vl = np.where( np.logical_and(idx_v4 > -1,  idx_v1 == -1) , idx_v4,   vl)  
        vl = np.where( np.logical_and(idx_v4 == -1, idx_v1 > -1)  , idx_v1-1, vl) 
        
        vr = np.where( np.logical_and(idx_v2 > -1,  idx_v3 > -1)  , idx_v2,   -1)  
        vr = np.where( np.logical_and(idx_v2 > -1,  idx_v3 == -1) , idx_v2,   vr)  
        vr = np.where( np.logical_and(idx_v2 == -1, idx_v3 > -1)  , idx_v3-1, vr) 
 
#        hb_list[cid1] = hb
#        ht_list[cid1] = ht
#        vl_list[cid1] = vl
#        vr_list[cid1] = vr
        # need to flip the sides
        num = v - (ncf[cid1] * (q-1) + ncf[cid1] - 1) - 2
        hb_list[cid1] = np.where(hb>-1, num-hb, hb)
        ht_list[cid1] = np.where(ht>-1, num-ht, ht)
        vl_list[cid1] = np.where(vl>-1, num-vl, vl)
        vr_list[cid1] = np.where(vr>-1, num-vr, vr)



        for cid2 in range(c.shape[0]):           
            X1 = vertices[c[cid2,0],0]
            Y1 = vertices[c[cid2,0],1]
            if ncf[cid1] == 1:
                X2 = vertices[c[cid2, -ncf[cid1] *(q-1) - 1 ],0]
                Y2 = vertices[c[cid2, -ncf[cid1] *(q-1) - 1 ],1]
            else:
                X2 = vertices[c[cid2, -ncf[cid1] *(q-1) - 1 -1 ],0]
                Y2 = vertices[c[cid2, -ncf[cid1] *(q-1) - 1 -1 ],1]

            bc1 = domain.bc(X1,Y1)
            bc2 = domain.bc(X2,Y2)
            bid1 = domain.bc_id(bc1)
            bid2 = domain.bc_id(bc2)
            f.write("{}\n{}\n{}\n".format(v - (ncf[cid1] * (q-1) + ncf[cid1] - 1), ncf[cid1] * (q-1) + ncf[cid1] - 1 , q))
            f.write("{}\n".format( ncf[cid1] )  )
            f.write("{} {}\n".format( bid1, bid2 ))
            cl = c[cid2,:]
            num = v - (ncf[cid1] * (q-1) + ncf[cid1] - 1)
            cl[:num] = np.flip(cl[:num])
            cl[num:] = np.flip(cl[num:])
            np.savetxt(f, vertices[cl,:], fmt = '%.16E')
             
    np.savetxt( f, np.array( [ 3, 2, 1, 0] ).reshape( (1,-1) ).astype(int) , fmt = "%i")            
    for cid1 in range(1,len(cells)):
        c = cells[cid1]
        if c.shape[0] == 0:
            continue
        np.savetxt(f, np.hstack( (hb_list[cid1][:,None], 
                                  vr_list[cid1][:,None], 
                                  ht_list[cid1][:,None],  
                                  vl_list[cid1][:,None] ) ) , fmt = '%i')
        
    f.close()

        
#    from matplotlib import collections  as mc
#    from matplotlib.collections import LineCollection
#    import matplotlib.pyplot as plt
#    for cid1 in range(1, len(cells)):
#        c = cells[cid1]
#        v = nv[cid1]
#        if c.shape[0] == 0:
#            continue
#
#        for cid2 in range(c.shape[0]):           
#            XX = vertices[c[cid2,:],0]
#            YY = vertices[c[cid2,:],1]
#            XX = np.hstack( (XX, XX[0] ) )
#            YY = np.hstack( (YY, YY[0] ) )
#
#            plt.plot(XX,YY)
#            if hb_list[cid1][cid2] > -1:
#                plt.scatter(np.mean(XX[hb_list[cid1][cid2]:hb_list[cid1][cid2]+2]),
#                            np.mean(YY[hb_list[cid1][cid2]:hb_list[cid1][cid2]+2])+0.01, c = 'black',s=2)
#            if ht_list[cid1][cid2] > -1:
#                plt.scatter(np.mean(XX[ht_list[cid1][cid2]:ht_list[cid1][cid2]+2]),
#                            np.mean(YY[ht_list[cid1][cid2]:ht_list[cid1][cid2]+2])-0.01, c = 'green', s=2)
#            if vl_list[cid1][cid2] > -1:
#                plt.scatter(np.mean(XX[vl_list[cid1][cid2]:vl_list[cid1][cid2]+2])+0.01,
#                            np.mean(YY[vl_list[cid1][cid2]:vl_list[cid1][cid2]+2]), c = 'blue',s=2)
#            if vr_list[cid1][cid2] > -1:
#                plt.scatter(np.mean(XX[vr_list[cid1][cid2]:vr_list[cid1][cid2]+2])-0.01,
#                            np.mean(YY[vr_list[cid1][cid2]:vr_list[cid1][cid2]+2]), c = 'purple', s=2)    
#        plt.show()


   
