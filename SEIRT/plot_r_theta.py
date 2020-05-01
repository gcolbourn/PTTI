
import numpy as np
np.random.seed(0)
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from cpyment import CModel
from SEIR_x_UD import seirxud_abm_gill, SEIRxUD

N = 67000000
I0 = 1.0/67000000
c = 13
beta = 0.033
alpha = 0.2
gamma = 1.0/7
kappa = 1.0/14
eta = 0.0
chi = 0.0

tmax = 600

params = {
    'c': c,
    'beta': beta,
    'alpha': alpha,
    'gamma': gamma,
    'kappa': kappa,
    'eta': eta,
    'chi': chi
}

fontP = FontProperties()
fontP.set_size('small')

fig, (ax1, ax2) = plt.subplots(2,1, figsize=(8,8))
ax1.set_xlabel("t")
ax1.set_ylabel("R(t)")
ax2.set_xlabel("t")
ax2.set_ylabel("E(t) + I(t)")

rr = []
for i in range(10):
    theta = float(i)/10
    params["theta"] = theta
    model = SEIRxUD(N=N, tmax=tmax, tsteps=30000, I0=I0, **params)
    traj = model.run_cmodel(etadamp=1, has_memory_states=True)["y"]

    R0 = beta*c/gamma
    print("theta = {0} R0 = {1}".format(theta, R0))

    EI = traj[:,2] + traj[:,3] + traj[:,4] + traj[:,5]

    r = model.rseries(traj)

    ax1.plot(model.t, r, lw=0.75, label="theta = {}".format(theta))
    ax2.plot(model.t, EI, lw=0.75, label="theta = {}".format(theta))

ax1.set_xlim(0,400)
ax1.set_ylim(0,3.5)
ax2.set_xlim(0,400)
ax1.legend(prop=fontP, loc="upper right")
ax2.legend(prop=fontP, loc="upper right")
plt.subplots_adjust(right=0.85)
plt.savefig("r_theta.png")
