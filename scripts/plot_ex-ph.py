## plot color map of 
import numpy as np 
from netCDF4 import Dataset
from matplotlib import pyplot as plt 
from io_exph import *
import scienceplots
import sys
import yaml
# ---- Plotting options -----
plt.style.use(['science','nature'])
plt.rcParams.update({"font.family": "Helvetica", 'font.size': 12,
    'axes.labelsize': 8,'axes.titlesize':10,
    'legend.fontsize': 8, 'xtick.labelsize': 8, 'ytick.labelsize': 8})
width = 3.375; height = 2.2
fig = plt.figure(); fig.set_size_inches(width, height)

"""
# Example input file
calc_folder: "./"
SAVE_dir: "gw_bse/SAVE"
elph_file: "elph/ndb.elph"
Gamma_states: [0,4]
Q_states: [0,4]
exph_file= "./Ex-ph.npy"
gauge_invar = true
"""
# --- Read input file ---
if len(sys.argv) < 2:
    print("Usage: python3 xx.py inputfile.yaml")
    sys.exit(1)

input_file = sys.argv[1]

with open(input_file, "r") as f:
    params = yaml.safe_load(f)

# Read inputs
calc_folder = params.get("calc_folder", ".")
SAVE_dir = calc_folder + params.get("SAVE_dir", "SAVE")
elph_file = calc_folder + params.get("elph_file", "ndb.elph")
Gamma_states = params.get("Gamma_states", [0,4])
Q_states = params.get("Q_states", [0,4])
exph_file= params.get("exph_file","Ex-ph.npy")
gauge_invar = params.get("gauge_invar", True)
#
## read the lattice data
print('*' * 30, ' Program started ', '*' * 30)
#
print("\n===== Input Parameters =====")
print(f"calc_folder : {calc_folder}")
print(f"SAVE_dir    : {SAVE_dir}")
print(f"elph_file   : {elph_file}")
print(f"Gamma states: {Gamma_states}")
print(f"Q states: {Gamma_states}")
print(f"Ex-ph file: {exph_file}")
print(f"gauge_invar : {gauge_invar}")
print("============================\n")
#
print('Reading Lattice data')
lattice, nibz, symm_mats, ele_time_rev, _ = get_SAVE_Data(save_folder=SAVE_dir)
blat = 2*np.pi*np.linalg.inv(lattice.T)
kpt = Dataset(elph_file)['qpoints'][...].data.astype(numpy_float)
kpt_cart = np.einsum('ij,nj->ni',blat,kpt,optimize=True)
exph = np.load(exph_file)

if gauge_invar:
    exph_abs = np.abs(exph)**2
    exph_abs = np.sum(exph_abs[:,:,Gamma_states[0]:Gamma_states[1],Q_states[0]:Q_states[1]],axis=(1,2,3))
    print(exph_abs.shape)
    exph_abs = np.sqrt(exph_abs)*27.2114
else :
    exph_abs = np.abs(exph)**2
    exph_abs = np.sum(exph_abs[:,:,Gamma_states[0]:Gamma_states[1],Q_states[0]:Q_states[1]],axis=(2,3))
    print(exph_abs.shape)
    exph_abs = np.sum(np.sqrt(exph_abs)*27.2114,axis=-1)
#print(exph.shape)


for i in range(-1,0):
    for j in range(-1,0):
        kp1 = kpt.copy()
        kp1[:,0] = kp1[:,0]+i
        kp1[:,1] = kp1[:,1]+j
        kpt_cart1 = np.einsum('ij,nj->ni',blat,kp1,optimize=True)
        plt.scatter(kpt_cart1[:,0],kpt_cart1[:,1],c=exph_abs,s=25,linewidth=0.5,vmin=0)

plt.axis('equal')
plt.xlabel('$q_x$ ($\\AA^{-1}$)')
plt.ylabel('$q_y$ ($\\AA^{-1}$)')
plt.title('$\\mathcal{G} = \\Sigma_{S^\prime}^{4}\\Sigma_{S=1}^{2}\\Sigma_{\\nu} \\big | \\langle \\mathbf{Q}=\mathbf{q}, S^\prime | \\partial_{\\mathbf{q}}^\\nu V_{\\text{scf}} | \\mathbf{Q}=\\mathbf{0}, S\\rangle \\big |^2$',fontsize=7)
cb = plt.colorbar()
cb.ax.set_title('$\\mathcal{G}$ (eV)',fontsize=8)
fig.savefig("plot_ex_ph.pdf",bbox_inches='tight',pad_inches = 0.03)
plt.show()

