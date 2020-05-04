#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from SEIR_x_UD import SEIRxUD, INDEX_SU, INDEX_SD, INDEX_EU, INDEX_ED, INDEX_IU, INDEX_ID, INDEX_RU, INDEX_RD
import sys

def _recompute(gridsize, output, **params):
    etarange = np.linspace(0.0, 1.0, gridsize)
    thetarange = np.linspace(0.0, 0.5, gridsize)

    rows = []
    for eta in etarange:
        for theta in thetarange:
            p = params.copy()
            p["eta"] = eta
            p["theta"] = theta
            print(p)
    
            model = SEIRxUD(**p)
    
            traj = model.run_cmodel(has_memory_states=True)['y']
    
            S = traj[:,INDEX_SU] + traj[:,INDEX_SD]
            E = traj[:,INDEX_EU] + traj[:,INDEX_ED]
            I = traj[:,INDEX_IU] + traj[:,INDEX_ID]
            R = traj[:,INDEX_RU] + traj[:,INDEX_RD]
    
            Rs = model.rseries(traj)
            peak = np.argmax(E+I)
    
            row = np.array([eta, theta, np.log10(np.min(S)), np.log10(np.max(E+I)), np.log10(np.max(R)), model.t[peak], Rs[60], Rs[-1]])
            print("\t".join(map(str, row)))
            rows.append(row)
    
    field = np.array(rows)
    np.savetxt("{0}-scan_eta_theta.tsv".format(output), field, delimiter="\t")
    return field

def plot(gridsize, output, **params):
    try:
        field = np.loadtxt("{0}-scan_eta_theta.tsv".format(output), delimiter="\t")
    except:
        field = _recompute(gridsize, output, **params)

    extent = [np.min(field[:,0]), np.max(field[:,0]), np.min(field[:,1]), np.max(field[:,1])]
    
    print(extent)
    labels = {
      "S": "log(uninfected)",
      "I": "log(peak infections)",
      "R": "log(total infections)",
      "T": "days to peak",
      "P": "R(60 days)",
    }
    
    for i, j in enumerate("SIRTP"):
        print(i,j)
        fig, ax = plt.subplots()
        if j in "P":
            extra = { "norm": mpl.colors.Normalize(vmin=-1, vmax=3) }
        else:
            extra = {}
        im = ax.imshow(np.reshape(field[::, 2+i], (gridsize, gridsize)).T[::-1,:], extent=extent, cmap=plt.get_cmap("inferno"), **extra)
        ax.set_xlabel("eta - tracing success")
        ax.set_ylabel("theta - testing rate")
        bar = fig.colorbar(im, shrink=0.55)
        bar.set_label(labels[j])
        fig.savefig("{0}-{1}_phase.png".format(output, j))
