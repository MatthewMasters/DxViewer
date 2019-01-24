import numpy as np
import math
import sys

# python writeDx.py 3bkl_protein.pdb
# returns 3bkl_C.dx, 3bkl_O.dx, 3bkl_N.dx, 3bkl_S.dx

def writeDx(filename, grid, origin=[0,0,0]):
    delta = [1.0,1.0,1.0]

    out_dx = open(filename, 'w')
    out_dx.write("object 1 class gridpositions counts %i %i %i\n" % (grid.shape[0],grid.shape[1],grid.shape[2]))
    out_dx.write("origin %.3f %.3f %.3f\n" % (origin[0],origin[1],origin[2]))
    out_dx.write("delta %.5f   0   0\n" % delta[0])
    out_dx.write("delta 0   %.5f   0\n" % delta[1])
    out_dx.write("delta 0   0   %.5f\n" % delta[2])
    out_dx.write("object 2 class gridconnections counts %i %i %i\n" % (grid.shape[0],grid.shape[1],grid.shape[2]))
    out_dx.write("object 3 class array type double rank 0 items %i data follows\n" % grid.size)

    count=0
    #print(grid)
    for data in grid.flat:
        if count==3: out_dx.write('\n'); count=0
        out_dx.write('%8.5f\t'%(float(data)))
        count+=1
    out_dx.write('\n')
    out_dx.close()

origin = [999,999,999]
fo = open(sys.argv[1])
for line in fo:
    if line[:4] == 'ATOM':
        atm_type = str(line[76:78]).strip()
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        if x < origin[0]:
            origin[0] = x
        if y < origin[1]:
            origin[1] = y
        if z < origin[2]:
            origin[2] = z
fo.close()

gsize = 96 #48

data = np.zeros((4,gsize,gsize,gsize))

fo = open(sys.argv[1])
for line in fo:
    if line[:4] == 'ATOM':
        atm_type = str(line[76:78]).strip()
        x = float(line[30:38]) - origin[0]
        y = float(line[38:46]) - origin[1]
        z = float(line[46:54]) - origin[2]

        x_b = math.floor(x)
        y_b = math.floor(y)
        z_b = math.floor(z)

        if 0 <= x <= gsize-1 and 0 <= y <= gsize-1 and 0 <= z <= gsize-1:
            x_space = [m for m in range(x_b-1,x_b+3)]
            y_space = [m for m in range(y_b-1,y_b+3)]
            z_space = [m for m in range(z_b-1,z_b+3)]
            for x_s in x_space:
                if x_s < 0 or x_s >= gsize:
                    continue
                for y_s in y_space:
                    if y_s < 0 or y_s >= gsize:
                        continue
                    for z_s in z_space:
                        if z_s < 0 or z_s >= gsize:
                            continue
                        dis = math.sqrt((x_s - x) ** 2 + (y_s - y) ** 2 + (z_s - z) ** 2)
                        density = math.e**(-0.6*dis**2)
                        if atm_type == 'C':
                            data[0,x_s,y_s,z_s] = data[0,x_s,y_s,z_s] + density
                        elif atm_type == 'O':
                            data[1,x_s,y_s,z_s] = data[1,x_s,y_s,z_s] + density
                        elif atm_type == 'N':
                            data[2,x_s,y_s,z_s] = data[2,x_s,y_s,z_s] + density
                        elif atm_type == 'S':
                            data[3,x_s,y_s,z_s] = data[3,x_s,y_s,z_s] + density
fo.close()

writeDx('3bkl_C.dx', data[0], origin)
writeDx('3bkl_O.dx', data[1], origin)
writeDx('3bkl_N.dx', data[2], origin)
writeDx('3bkl_S.dx', data[3], origin)



