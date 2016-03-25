import numpy as np

def xdatcar_2_xyz(path_to_xdatcar_folder):
    xdatcar = open(path_to_xdatcar_folder+"/XDATCAR", 'r')
    xyz = open(path_to_xdatcar_folder+"/XDATCAR.xyz", 'w')
    xyz_fract = open(path_to_xdatcar_folder+"/XDATCAR_fract.xyz", 'w')

    system = xdatcar.readline()
    scale = float(xdatcar.readline().rstrip('\n'))

    #get lattice vectors
    a1 = np.array([ float(s)*scale for s in xdatcar.readline().rstrip('\n').split() ])
    a2 = np.array([ float(s)*scale for s in xdatcar.readline().rstrip('\n').split() ])
    a3 = np.array([ float(s)*scale for s in xdatcar.readline().rstrip('\n').split() ])

    #Save scaled lattice vectors
    lat_rec = open(path_to_xdatcar_folder+"/lattice.vectors", 'w')
    lat_rec.write(str(a1[0])+' '+str(a1[1])+' '+str(a1[2])+'\n')
    lat_rec.write(str(a2[0])+' '+str(a2[1])+' '+str(a2[2])+'\n')
    lat_rec.write(str(a3[0])+' '+str(a3[1])+' '+str(a3[2]))
    lat_rec.close()


    #Read xdatcar
    element_names = xdatcar.readline().rstrip('\n').split()

    element_dict = {}
    element_numbers = xdatcar.readline().rstrip('\n').split()

    i = 0
    N = 0
    for t in range(len(element_names)):
        element_dict[element_names[t]] = int(element_numbers[i])
        N += int(element_numbers[i])
        i += 1

    while True:
        line = xdatcar.readline()
        if len(line) == 0:
            break
        xyz.write(str(N) + "\ncomment\n")
        xyz_fract.write(str(N)+"\ncomment\n")
        for el in element_names:
            for i in range(element_dict[el]):
                p = xdatcar.readline().rstrip('\n').split()
                coords = np.array([ float(s) for s in p ])
                cartesian_coords = coords[0]*a1+coords[1]*a2+coords[2]*a3
                xyz.write(el+ " " + str(cartesian_coords[0])+ " " + str(cartesian_coords[1]) + " " + str(cartesian_coords[2]) +"\n")
                xyz_fract.write(el+ " " + str(coords[0])+ " " + str(coords[1]) + " " + str(coords[2]) +"\n")
    xdatcar.close()
    xyz.close()
    xyz_fract.close()