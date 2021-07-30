# KG 7/11/2021

import numpy as np
import copy
import matplotlib.pyplot as plt
from code2_ev_and_freq import phononBC
from code4_modal_calc import En_hist_n_center

class contributing_phonons_n_PBC:
    def __init__(self,val1, val2, val3, val4, val5, val6):
        self.w = val1 # Phonon frequencies
        self.v = val2 # Phonon eigenvectors
        self.En = val3 # modal contributions
        self.atomsp = val4 # atom species (string)
        self.atomNs = val5 # atom species numbers (int)
        self.exclude_hopping_ion = val6

    def plot_En_PBC(self):
        freq = self.w.copy()
        self.En /= np.sum(self.En) # Normalizing En
        En = self.En.copy()

        modalEn_instance = En_hist_n_center(freq, En, 50)
        [w_En, counter_En, w_En_center] = modalEn_instance.hist_for_En_n_center()
        phononBC_instance = phononBC(50, freq, self.v, self.atomNs, self.atomsp)
        [w_TotalDOS, counter_TotalDOS, TotalPBC] = phononBC_instance.TotalDOS() # Calculate Total DOS histogram
        [w_partialDOS, counter_partialDOS, partialPBC] = phononBC_instance.partialDOS('O') # Calculate partial DOS histogram
        counter_partialDOS /= np.sum(counter_TotalDOS) # Normalizing partial DOS first
        counter_TotalDOS /= np.sum(counter_TotalDOS) # Then, normalizing total DOS

        # plot series 1
        plt.figure('Phonon contributions')
        plt.rcParams.update({'font.size': 13})
        plt.fill_between(w_En, 0, self.smooth(counter_En,3), color="black", alpha=0.5, label='Phonon contributions')
        plt.fill_between(w_TotalDOS, 0, self.smooth(counter_TotalDOS,3), color="blue", alpha=0.5, label='Total phonon DOS')
        plt.fill_between(w_partialDOS, 0, self.smooth(counter_partialDOS,3), color="red", alpha=0.5, label='Oxygen pDOS')
        plt.plot(freq,En, linestyle='none', marker='o', markerfacecolor='none', markersize=7, markeredgecolor='k', label='Phonon contributions')

        minx = np.floor(np.min(freq))
        maxx = np.ceil(np.max(freq))
        miny = 0
        maxy = np.ceil(np.max(En)*25.)/25.

        #plt.plot([w_En_center,w_En_center],[miny,0.6*maxy],'k--')
        #plt.plot([TotalPBC,TotalPBC],[miny,0.6*maxy],'b--')
        #plt.plot([partialPBC,partialPBC],[miny,0.6*maxy],'r--')

        #plt.text(w_En_center, 0.63*maxy, 'En_center', fontsize=10, rotation = 'vertical', color = 'k')
        #plt.text(TotalPBC, 0.63*maxy, 'tPBC', fontsize=10, rotation = 'vertical', color = 'b')
        #plt.text(partialPBC, 0.63*maxy, 'pPBC', fontsize=10, rotation = 'vertical', color = 'r')
        
        plt.xlabel('Frequency (THz)')
        plt.ylabel('Normlized excitation energy')
        plt.xlim(minx, maxx)
        plt.ylim(miny, maxy)

        plt.legend(fontsize=11)
        filename = 'phonon_exclude_' + str(self.exclude_hopping_ion) + '.png'
        plt.savefig(filename)
        #plt.show()
        plt.close()

    def smooth(self, y, box_pts):
        box = np.ones(box_pts)/box_pts
        y_smooth = np.convolve(y, box, mode='same')
        return y_smooth

