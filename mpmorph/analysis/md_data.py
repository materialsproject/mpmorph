import numpy as np
import re

def get_MD_data(outcar_path, search_keys = ['external', 'kinetic energy EKIN', '% ion-electron', 'ETOTAL'], search_data_column = [3, 4, 4, 4] ):
    '''
    Extracts the pressure, kinetic energy and total energy data from
    VASP MD OUTCAR.

    Args:
        - outcar_path = path to OUTCAR to be parsed
    Returns:
        - A nested list of MD steps where each search key value is
          listed.
    '''
    # Initial map of keywords to serach for and data to map out from that line in OUTCAR
    # search_keys = ['external', 'kinetic energy EKIN', 'ETOTAL']
    # index of stripped column of data in that line, starts from 0
    # search_data_column = [3, 4, 4]
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