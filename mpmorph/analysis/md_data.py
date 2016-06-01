import numpy as np
import re
import os

def get_MD_data(outcar_path, search_keys=None, search_data_column=None):
    '''
    Extracts the pressure, kinetic energy and total energy data from
    VASP MD OUTCAR.

    Args:
          outcar_path:
          search_keys:
          search_keys:
          search_data_column:
        - outcar_path = path to OUTCAR to be parsed
    Returns:
        - A nested list of MD steps where each search key value is
          listed.
    '''
    # Initial map of keywords to serach for and data to map out from that line in OUTCAR
    # search_keys = ['external', 'kinetic energy EKIN', 'ETOTAL']
    # index of stripped column of data in that line, starts from 0
    # search_data_column = [3, 4, 4]
    if search_data_column is None:
        search_data_column = [3, 4, 4, 4]
    if search_keys is None:
        search_keys = ['external', 'kinetic energy EKIN', '% ion-electron', 'ETOTAL']
    outcar = open(outcar_path)
    print "OUTCAR opened"
    data_list = []
    md_step = 0
    for line in outcar:
        line = line.rstrip()
        for key_index in range(len(search_keys)):
            if re.search(search_keys[key_index],line):
                if key_index == 0:
                    data_list.append([float(line.split()[search_data_column[key_index]])])
                else:
                    try:
                        data_list[md_step].append(float(line.split()[search_data_column[key_index]]))
                    except IndexError:
                        break
                if key_index == len(search_keys)-1:
                    md_step +=1
    print "Requested information parsed."
    outcar.close()
    return data_list

def autocorrelation(data_list):
    """
    TODO
    Args:
        data_list:

    Returns:

    """

def get_MD_stats(data_list):
    """
    Args: data_list is the list of MD data returned by get_MD_data
    Returns: means and standard deviations
    """
    data_list = np.array(data_list)
    stats = []
    for col in range(data_list.shape[1]):
        data_col = data_list[:,col]
        stats.append( ( np.mean(data_col), np.std(data_col) ) )
    return stats

def plot_md_data(data_list):
    '''
    Args:
        data_list:

    Returns:
        matplotlib plt object

    '''

def parse_pressure(path, averaging_fraction=0.5):
    os.system("grep external " + path + "/OUTCAR o| awk '{print $4}' > "+path +"/pres")
    os.system("grep volume/ion " + path + "/OUTCAR | awk '{print $5}' > "+path +"/vol")
    if os.path.isfile(path+"/OUTCAR"):
        with open(path+"/pres") as f:
            p = [float(line.rstrip()) for line in f]
        with open(path+"/vol") as f:
            vol = [float(line.rstrip()) for line in f][0]
        pressure = np.array(p)
        avg_pres = np.mean( pressure[int(averaging_fraction*(len(pressure)-1)):] )
    else:
        raise ValueError("No OUTCAR found.")
    return avg_pres, vol, pressure