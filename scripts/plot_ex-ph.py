## plot color map of 
import numpy as np 
from matplotlib import pyplot as plt 
import scienceplots
plt.style.use(['science','nature'])
plt.rcParams.update({"font.family": "Helvetica", 'font.size': 12,
    'axes.labelsize': 8,'axes.titlesize':10,
    'legend.fontsize': 8, 'xtick.labelsize': 8, 'ytick.labelsize': 8})
width = 3.375; height = 2.2
fig = plt.figure(); fig.set_size_inches(width, height)



Gamma_states=[0,4]
Q_states = [0,4]
gauge_invar = True
kpt = np.loadtxt('kpts_48_48')
lattice = np.array([[6.277054, -3.138527, 0], [0, 5.436089, 0], [0, 0, 32.1253] ])






lattice = lattice*0.529177 ## bohr to angstrom
blat = 2*np.pi*np.linalg.inv(lattice.T)

kpt_cart = np.einsum('ij,nj->ni',blat,kpt,optimize=True)
exph = np.load('Ex-ph.npy')

if gauge_invar:
    exph_abs = np.abs(exph)**2
    exph_abs = np.sum(exph_abs[:,:,Gamma_states[0]:Gamma_states[1],Q_states[0]:Q_states[1]],axis=(1,2,3))
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

