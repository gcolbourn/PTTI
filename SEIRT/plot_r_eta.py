
import numpy as np
np.random.seed(0)
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from cpyment import CModel
from SEIR_x_UD import seirxud_abm_gill, SEIRxUD

def plot(gridsize, output, **params):
    fontP = FontProperties()
    fontP.set_size('small')

    fig, (ax1, ax2) = plt.subplots(2,1, figsize=(8,8))
    ax1.set_xlabel("t")
    ax1.set_ylabel("R(t)")
    ax2.set_xlabel("t")
    ax2.set_ylabel("E(t) + I(t)")
   
    beta, c, gamma = params["beta"], params["c"], params["gamma"]
 
    rr = []
    for i in range(gridsize+1):
        eta = float(i)/gridsize
        p = params.copy()
        p["eta"] = eta
        model = SEIRxUD(**p)
        traj = model.run_cmodel(etadamp=1, has_memory_states=True)["y"]
    
        R0 = beta*c/gamma
        print("eta = {0} R0 = {1}".format(eta, R0))
    
        EI = traj[:,2] + traj[:,3] + traj[:,4] + traj[:,5]
    
        r = model.rseries(traj)
    
        ax1.plot(model.t, r, lw=0.75, label="eta = {}".format(eta))
        ax2.plot(model.t, EI, lw=0.75, label="eta = {}".format(eta))
   
    tmax = params["tmax"]
 
    ax1.axhline(1, c=(0,0,0), lw=0.5, ls='--')
    ax1.set_xlim(0,tmax)
    ax1.set_ylim(0,3.5)
    ax2.set_xlim(0,tmax)
    ax1.legend(prop=fontP, loc="upper right")
    ax2.legend(prop=fontP, loc="upper right")
    plt.subplots_adjust(right=0.85)
    plt.savefig("{0}-r_eta.png".format(output))
