#!/usr/bin/env python

import numpy as np
from SEIR_x_UD import SEIRxUD

N = 1000
I0 = 0.01
c = 13
beta = 0.018
alpha = 0.2
gamma = 0.0714
theta = 0.15
kappa = 0.03
eta = 0.1
chi = 0.1

tmax = 600

# Scan ranges?
chirange = np.linspace(0.1, 1.0, 10)
etarange = np.linspace(0.1, 1.0, 10)

f = open('scan_results_theta{0}.dat'.format(theta), 'w')

for chi in chirange:
    for eta in etarange:
        params = {
            'c': c,
            'beta': beta,
            'alpha': alpha,
            'gamma': gamma,
            'theta': theta,
            'kappa': kappa,
            'eta': eta,
            'chi': chi
        }

        model = SEIRxUD(N=N, tmax=tmax, I0=I0, **params)

        trajsABM = model.run_abm(samples=10)

        trajABMAvg = np.average(trajsABM, axis=0)
        trajABMStd = np.std(trajsABM, axis=0)
        trajODE = model.run_cmodel(has_memory_states=True)['y']

        # Measure
        finalDiff = trajODE[-1,:8]-trajABMAvg[:,-1]
        finalRelDiff = finalDiff/trajABMStd[:,-1]
        peakiODE = np.argmax(trajODE[:,:8], axis=0)
        peakiABM = np.argmax(trajABMAvg, axis=1)
        peakDiff = trajODE[peakiODE, range(8)]-trajABMAvg[range(8), peakiABM]
        peakRelDiff = peakDiff/trajABMStd[range(8), peakiABM]
        peakTDiff = model.t[peakiODE]-model.t[peakiABM]

        f.write('{0} {1} '.format(chi, eta))
        f.write(' '.join(map(str, finalDiff)) + ' ')
        f.write(' '.join(map(str, finalRelDiff)) + ' ')
        f.write(' '.join(map(str, peakDiff)) + ' ')
        f.write(' '.join(map(str, peakRelDiff)) + ' ')
        f.write(' '.join(map(str, peakTDiff)) + ' ')

        f.write('\n')

    f.write('\n')
        



