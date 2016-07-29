import numpy as np
from copy import deepcopy
import json

# This function reads an XYZ file and a list of lattice vectors L = [x,y,z] and gives MSD + unwrapped coordinates

# This function should take an xdatcar object! not xyz file...

def MSD(xyz_file_folder,L=None):
    if L==None:
        L = read_lat_vec(xyz_file_folder)
    a = []; l = [];
    a.append(L[0]); a.append(L[1]); a.append(L[2]); #basis vectors in cartesian coords
    l.append(np.sqrt(np.dot(a[0],a[0]))); l.append(np.sqrt(np.dot(a[1],a[1]))); l.append(np.sqrt(np.dot(a[2],a[2]))); #basis vector lengths

    file = open(xyz_file_folder+"/XDATCAR_fract.xyz", 'r')
    recorder = open(xyz_file_folder+"/msd.out", 'w')
    coord_rec = open(xyz_file_folder+"/unwrapped.xyz", 'w')

    origin_list = [] # Stores the origin as [element,[coords]]
    prev_list = [] # Stores the wrapped previous step
    unwrapped_list = [] # Stores the instantenous unwrapped

    msd = [] #Stores atom-wise MSD  Stores msd as [msd]
    msd_dict ={} #Stores element-wise MSD
    msd_lattice = []
    msd_dict_lattice ={}

    element_list = [] # element list
    element_dict = {} # number of elements stored

    content = file.readline()
    N = int(content)

    for i in range(N):
        msd.append(np.float64('0.0'))

        msd_lattice.append([0.0, 0.0, 0.0 ])

    file.readline()
    step = 0

    while True:
        step += 1
# Get and store the origin coordinates in origin_dict at first step
        if step == 1:
            for i in xrange(N):
                t = file.readline().rstrip('\n').split()
                element = t[0]
                if element not in element_list:
                    element_list.append(element)
                if element not in element_dict:
                    element_dict[element] = 1.0
                else:
                    element_dict[element] += 1.0
                coords = np.array( [ float(s) for s in t[1:] ] )
                origin_list.append([element,coords])
      # Copy the first set of coordinates as prev_dict and unwrapped
            unwrapped_list = deepcopy(origin_list)
            prev_list = deepcopy(origin_list)
            recorder.write("step ")
            for element in element_list:
                recorder.write(element+" ")
                for dir in ['_x','_y','_z']:
                    recorder.write(element+dir+" ")
            recorder.write("\n")

# Read wrapped coordinates into wrapped_dict
        content = file.readline()
        if len(content) == 0:
            # print "\n---End of file---\n"
            break
        N = int(content)
        file.readline()
        wrapped_list = [] # Erease the previous set of coordinates
        for i in xrange(N):
            t = file.readline().rstrip('\n').split()
            element = t[0]
            coords = np.array( [ float(s) for s in t[1:] ] )
            wrapped_list.append([element,coords])

        coord_rec.write(str(N)+ "\ncomment\n")

# Unwrap coodinates and get MSD

        for atom in range(N):

            msd[atom] = 0.0

            coord_rec.write(wrapped_list[atom][0])

            # decompose wrapped atom coordinates to onto lattice vectors:
            w1 = wrapped_list[atom][1][0]
            w2 = wrapped_list[atom][1][1]
            w3 = wrapped_list[atom][1][2]

            # decompose prev atom coordinates to onto lattice vectors:
            p1 = prev_list[atom][1][0]
            p2 = prev_list[atom][1][1]
            p3 = prev_list[atom][1][2]

            #get distance between periodic images and use the smallest one
            if np.fabs(w1 - p1) > 0.5:
                 u1 = w1 - p1 - np.sign(w1 - p1)
            else:
                 u1 = w1 - p1

            if np.fabs(w2 - p2) > 0.5:
                 u2 = w2 - p2 - np.sign(w2 - p2)
            else:
                 u2 = w2 - p2

            if np.fabs(w3 - p3) > 0.5:
                 u3 = w3 - p3 - np.sign(w3 - p3)
            else:
                 u3 = w3 - p3

            #add unwrapped displacements to unwrapped coords

            unwrapped_list[atom][1][0] += u1
            unwrapped_list[atom][1][1] += u2
            unwrapped_list[atom][1][2] += u3

            uw = unwrapped_list[atom][1][0]*a[0] + unwrapped_list[atom][1][1]*a[1] +unwrapped_list[atom][1][2]*a[2]
            ol = origin_list[atom][1][0]*a[0] + origin_list[atom][1][1]*a[1] + origin_list[atom][1][2]*a[2]

            msd[atom] = np.linalg.norm(uw-ol)**2
            msd_lattice[atom] = [np.linalg.norm(uw[0]-ol[0])**2,np.linalg.norm(uw[1]-ol[1])**2,np.linalg.norm(uw[2]-ol[2])**2]
#            print msd[atom]
#            print msd_lattice[atom]
            coord_rec.write(" " + np.array_str(uw).replace("[","").replace("]",""))
            coord_rec.write("\n")

        prev_list = [] # Store current wrapped coordinates for the next step
        prev_list = deepcopy(wrapped_list)
# record msd
        recorder.write(str(step) + " ")

        for el in element_list:
             msd_dict[el] = 0.0
             msd_dict_lattice[el]=[0.,0.,0.]

        for atom in range(len(msd)):
            msd_dict[wrapped_list[atom][0]] += msd[atom]/element_dict[wrapped_list[atom][0]]
            for i in range(3):
                msd_dict_lattice[wrapped_list[atom][0]][i] += msd_lattice[atom][i]/element_dict[wrapped_list[atom][0]]

        for el in element_list:
            recorder.write(str(msd_dict[el])+ " " + str(msd_dict_lattice[el][0])+ " " + str(msd_dict_lattice[el][1])+ " " + str(msd_dict_lattice[el][2])+ " ")

        recorder.write("\n")
        # if step % 10 == 0:
        #    print step
    recorder.close()
    file.close()
    coord_rec.close()
    with open(xyz_file_folder+"/composition.json","w") as f:
        json.dump(element_dict,f)

def read_lat_vec(path):
    lat_file = open(path+'/lattice.vectors','r')
    line = []
    for i in xrange(3):
        line.append([float(x) for x in lat_file.readline().rstrip('\n').split()])
        # print line[i]
    lattice = np.array([line[0],line[1],line[2]])
    return lattice

#You can read the lattice vectors from lattice.vector file
#lattice = read_lat_vec()

#Or define the lattice vector manually as in
#lattice =np.array([[-12.181156,-4.306689,7.459404],[0.000000,-12.920067,7.459404],[0.000000,0.000000,14.918808]])

#Run the MSD calculator with XDATCAR_fract.xyz and the lattice vector defined above
#MSD("XDATCAR_fract.xyz",lattice)