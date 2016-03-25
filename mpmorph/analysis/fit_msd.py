import numpy as np
import matplotlib.pyplot as plt
import json

class DiffusionCoeffs:
    def __init__(self, path_to_msd, t_step=2.0):
        '''
        :param path_to_msd: path to MSD output file (see msd module)
        :param t_step: (time step in femtoseconds)
        :return: None
        '''
        msd_file = open(path_to_msd+"/msd.out", "r")
        self.header = msd_file.readline().rstrip("\n").split()
        self.msd_data = np.loadtxt(msd_file)
        self.msd_data[:, 1:] *= 1e-16 # convert from Angs**2 to cm**2
        self.msd_data[:, 0] *= t_step * 1e-15  # convert to seconds.
        self.n_steps = self.msd_data.shape[0]
        self.n_species = (self.msd_data.shape[1]-1)//4
        self.skip_first = self.n_steps//3
        self.diff_coeffs = {}
        with open(path_to_msd+"/composition.json","r") as f:
            self.composition = json.load(f)
        print self.composition


    def compute_diff_coeff(self):
        '''
        Calculates diffusion coefficients
        '''
        self.diff_coeffs = {} # Clear any existing coeffs
        x = self.msd_data[:, 0]

        for col in range(1, len(self.header)):

            if "_" in self.header[col]:
                d = 2.0 # for 1 d diffision axes
            else:
                d = 6.0 # for 3 d diffusion
            y = self.msd_data[:, col]
            self.diff_coeffs[self.header[col]] = np.polyfit(x[self.skip_first:]*d, y[self.skip_first:], 1)

        print "element D(cm^2/s) intercept"
        for col in self.header[1:]:
            print col, self.diff_coeffs[col][0], self.diff_coeffs[col][1]

        return self.diff_coeffs

    def plot_MSD(self, temp=1000, plot_directions=False):

        if not len(self.diff_coeffs):
            print("MSD fitting must be performed prior to plotting.")
            self.compute_diff_coeff()

        plt.figure()

        x = self.msd_data[:, 0]*1e12  #convert to picoseconds for plotting
        legend_str = []
        for col in range(1, self.n_species*4+1):
            if "_" in self.header[col]:
                if plot_directions == False:
                    continue
            legend_str.append(self.header[col])
            y = self.msd_data[:, col]*1e16
            plt.plot(x,y)
            plt.xlabel("time (picoseconds)",fontsize=14)
            plt.ylabel('MSD ($\\times 10^{-16}$ cm$^2$/s)', fontsize=14)

        formula = '' #supercell composition
        for k,v in self.composition.items():
            formula += ''.join((k,"$_{",str(int(v)),"}$"))

        plt.legend(legend_str, bbox_to_anchor=(1.03, 1.0), loc=2, borderaxespad=0.1,fontsize=14)
        plt.title(formula+" T = " + str(temp) + " K",fontsize=14)
        plt.tick_params(labelsize=13)
        plt.show()
        return plt

class Activation:
    def __init__(self,dt_dict):
        '''
        dt_dict

        '''
        self.dt_dict = {}
        for el in dt_dict:
            self.dt_dict[el] = np.array(dt_dict[el])
        # check lenghts
        self.activation_barriers = {}
        self.d_rt = {}
        self.cov = {}

    def compute_activation(self):
        for el in self.dt_dict:
            self.dt_dict[el] = np.array([ [1/T, np.log(D)] for T,D in self.dt_dict[el] ] )
            d_fit = np.polyfit(self.dt_dict[el][:,0],self.dt_dict[el][:,1],1,cov=True)
            print "d_ft", d_fit
            d_RT = d_fit[0][0]/298.0+d_fit[0][1]
            print "-Q/k ln(D0) ln[D(RT)] D(RT)[cm2/s]"
            print d_fit[0], d_RT, np.exp(d_RT)
            print d_fit
            self.activation_barriers[el] = d_fit[0]
            self.d_rt[el] = d_RT
            self.cov[el] = d_fit[1]
        return self.activation_barriers


    def plot_activation(self, colors, title="Self-diffusion"):
        self.title = title
        plt.figure()
        legend_list=[]

        for el in self.dt_dict:
            plt.scatter(self.dt_dict[el][:,0]*1000, self.dt_dict[el][:,1],s=100,c=colors[el])
            plt.scatter(1.0/298.0*1000,self.d_rt[el],s=100,c=colors[el],marker='s')
            x_low, x_high = min(self.dt_dict[el][:,0]-0.0005), 1.0/298.0+0.0005
            y_low = self.activation_barriers[el][0]*x_low + self.activation_barriers[el][1]
            y_high = self.activation_barriers[el][0]*x_high + self.activation_barriers[el][1]
            plt.plot([x_low*1000, x_high*1000], [y_low, y_high],
                     color=colors[el], linestyle='-', linewidth=1)
            plt.xlim(x_low*1000,x_high*1000)
            legend_list.append(el
                               +" ($Q=$"+str('{0:.1f}'.format(-self.activation_barriers[el][0]))
                               +"$k_B$)")
        plt.legend(legend_list,bbox_to_anchor=(0.975, 0.975), loc=1, borderaxespad=0.1,fontsize=12)
        plt.title(self.title,fontsize=14)
        plt.ylabel("ln[$D$ (cm$^2$/s)]")
        plt.xlabel(r'$\frac{\mathrm{1000 K}}{T}$', fontsize=25)
        plt.grid(True)
        plt.show()


from collections import OrderedDict
from itertools import cycle

def plot_D_comp(D_list,plot_elements=None, title="", comp_axis_label="x"):
    if plot_elements == None:
        l = []
        for comp in D_list:
            for el in comp['comp'].keys() :
                l.append(el)
        plot_elements = list(set(l))

    plot_dict = {}
    for i in plot_elements:
        plot_dict[i]=[]
    for comp in D_list:
        for key in plot_dict:
            if key in comp['comp']:
                plot_dict[key].append((comp['comp'].get('Cr',0.0)/16.0, comp['coeffs'][key][0]))
    for key in plot_dict:
        arr = np.array(plot_dict[key])
        plot_dict[key]= arr[arr[:,0].argsort()]

    plot_dict = OrderedDict(plot_dict)
    markers = cycle(['o','s','^','+','.'])

    plt.figure()
    for key in plot_dict:
        plt.grid(True)
        plt.semilogy(plot_dict[key][:,0],plot_dict[key][:,1], linestyle='-',
                 marker=next(markers), markersize = 10, alpha=0.8)
        plt.xlabel(comp_axis_label,fontsize=16)
        plt.ylabel('$D_i$ cm$^2$/s',fontsize=16)
        plt.title(title)
        plt.xlim([-0.1,2.1])
        plt.tick_params(labelsize=15)

    plt.legend(plot_dict.keys(),bbox_to_anchor=(1.025, 0.975), loc=2, borderaxespad=0.1,fontsize=14)
    plt.show()