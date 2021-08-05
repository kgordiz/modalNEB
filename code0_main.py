# KG 7/29/2021
##
# The main code on the modalNEB package
##
import os
import numpy as np
import pandas as pd
import glob

from code1_POSCARs_info import POSCARs_info
from code2_ev_and_freq import ev_freq_info
from code2_ev_and_freq import phononBC
from code3_PBC_correct import PBC_correct
from code4_modal_calc import modal_calc
from code4_modal_calc import En_hist_n_center
import code5_energetics_n_MB_info as Eng_n_MB
from code100_plotting import contributing_phonons_n_PBC

########################### Some input info related to folder addresses ###########################
# If you need to activate a specific Python environment, make sure you do it before tyou run Python
Phonopy_path = "" # If phonopy is not accessible systemwide, give the path to it here

image1 = 0 # The image containing the hop origin
image2 = 2 # The image containing the after the hop image
image3 = 3 # This will correspond to the TS image number

exclude_hopping_ion = 0 # exclude the hopping ion from modal analysis? [1: yes] [0: no]

# start the dataframe
data = {
  "MB (eV)": [],
  "En Max freq (THz)": [],
  "En center freq (THz)": [],
  "Total PBC (THz)": [],
  "O PBC (THz)": []
}
df = pd.DataFrame(data)

########################### folder addresses
DFPT_vasprun_path = 'DFPT_00'
NEB_path = "./" 
###########################
# run Phonopy to get the Hessian matrix
os.system(Phonopy_path + "phonopy --fc " + DFPT_vasprun_path + "/vasprun.xml > out.txt") # FORCE_CONSTANTS files should be created in the current folder

# Read the structure info (atomic masses, xyz of atoms, etc.)
struct_info = POSCARs_info(image1, image2, image3, NEB_path) # The passed arguments are the image numbers (folder numbers), TS image number and the address to the NEB folder
struct_info.read_masses_and_xyzs()

LD_info = ev_freq_info(struct_info.mass) # The passed argument is the atomic mass vector
LD_info.perform_LD()

# calculate diff and apply PBC
dx_direct = struct_info.x1 - struct_info.x2
diff_correcting = PBC_correct(dx_direct, struct_info.L, struct_info.coeff, struct_info.Natoms) # The passed argument are the diff coords and the L vector & maybe the coeff
diff_correcting.perform_PBC_correct()

# calculate the modal contributions
modal_info = modal_calc(diff_correcting.dx_cartesian, LD_info.w, LD_info.v, struct_info.mass, struct_info.Natoms, exclude_hopping_ion) # The passed argument are ...
modal_info.perform_modal()

# Calculate the MB
MB = Eng_n_MB.MB_determine(NEB_path)

### Extract some other numbers
# Related to modal contributions
modalEn_instance = En_hist_n_center(LD_info.w, modal_info.En, 50)
[w_En, counter_En, w_En_center] = modalEn_instance.hist_for_En_n_center()
w_En_max = modalEn_instance.freq_max_En()
# Related to Phonon BCs
phononBC_instance = phononBC(50, LD_info.w, LD_info.v, struct_info.atomNs_passed, struct_info.atomsp_passed)
[w_TotalDOS, counter_TotalDOS, TotalPBC] = phononBC_instance.TotalDOS()
[w_partialDOS, counter_partialDOS, partialPBC] = phononBC_instance.partialDOS('O')

# add the info to the dataframe
df.loc[len(df.index)] = [MB,
                         w_En_max, 
                         w_En_center,
                         TotalPBC,
                         partialPBC]

# plot modal contributions and phonon band centers
plot1 = contributing_phonons_n_PBC(LD_info.w, LD_info.v, modal_info.En, struct_info.atomsp_passed, struct_info.atomNs_passed, exclude_hopping_ion)
plot1.plot_En_PBC()

filename = 'Results_' + str(image1) + '_' + str(image2) + '.csv'
df.index += 1 
df.to_csv(filename)

