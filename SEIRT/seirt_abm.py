import numpy as np
from numba import jit

PROG_S = 0
PROG_E = 1
PROG_I = 2
PROG_R = 3
PROG_T = 4


@jit(nopython=True)
def _gillespie(tmax, agents, beta, gamma, delta, theta, th_delay, etaCT,
               tblock=100):

    times = np.zeros(tblock)
    traj = np.zeros((tblock, len(agents), 3))
    traj[0] = agents

    i = 0
    t = 0
    while t < tmax:
        t1 = traj[i].copy()
        S = np.sum(t1[:, 0] == PROG_S)
        E = np.sum(t1[:, 0] == PROG_E)
        I = np.sum(t1[:, 0] == PROG_I)
        R = np.sum(t1[:, 0] == PROG_R)
        T = np.sum(t1[:, 0] == PROG_T)

        wSE = beta*S*I
        wEI = gamma*E
        wIR = delta*I
        wIT = theta*I
        wTR = delta*T

        W = wSE + wEI + wIR + wIT + wTR

        if W <= 0:
            break

        dt = -np.log(np.random.random())/W

        # What is the next upcoming CT?
        dCT = np.inf
        iCT = -1
        for a_i in range(len(agents)):
            dd = th_delay-t+t1[a_i, 2]
            if (t1[int(t1[a_i, 1]), 0] == PROG_T and
                    dd < dCT and dd >= 0):
                dCT = dd
                iCT = a_i

        if iCT >= 0 and dCT < dt:

            t += dCT

            if np.random.random() >= etaCT:
                # Failed tracing
                t1[iCT, 2] = np.inf
            else:
                # Traced!
                t1[iCT, 0] = PROG_T
                t1[iCT, 2] = np.inf
        else:

            t += dt

            # Pick an event type
            p = np.array([wSE, wEI, wIR, wIT, wTR])/W
            p = np.cumsum(p)
            echoice = np.random.random()

            if echoice < p[0]:
                # S to E
                Sis = np.where(t1[:, 0] == PROG_S)[0]
                Si = np.random.choice(Sis)
                Iis = np.where(t1[:, 0] == PROG_I)[0]
                Ii = np.random.choice(Iis)
                t1[Si, 0] = PROG_E
                t1[Si, 1] = Ii
                t1[Si, 2] = t
            elif echoice < p[1]:
                # E to I
                Eis = np.where(t1[:, 0] == PROG_E)[0]
                Ei = np.random.choice(Eis)
                t1[Ei, 0] = PROG_I
            elif echoice < p[2]:
                # I to R
                Iis = np.where(t1[:, 0] == PROG_I)[0]
                Ii = np.random.choice(Iis)
                t1[Ii, 0] = PROG_R
            elif echoice < p[3]:
                # I to T
                Iis = np.where(t1[:, 0] == PROG_I)[0]
                Ii = np.random.choice(Iis)
                t1[Ii, 0] = PROG_T
            elif echoice < p[4]:
                Tis = np.where(t1[:, 0] == PROG_T)[0]
                Ti = np.random.choice(Tis)
                t1[Ti, 0] = PROG_R

        i += 1
        if i >= len(traj):
            traj = np.concatenate((traj, np.zeros((tblock, len(agents), 3))),
                                  axis=0)
            times = np.concatenate((times, np.zeros(tblock)), axis=0)
        traj[i] = t1
        times[i] = t

    return i, times, traj


class ABMSEIRT(object):
    """
    beta:   S->E
    gamma:  E->I
    delta:  I->R, T->NT
    theta:  I->T
    """

    def __init__(self, beta, gamma, delta, theta, th_delay, etaCT=1,
                 N=1000):

        # Initialise a population, everyone
        # The columns are:
        # - state
        # - infector
        # - time at which infector has been detected

        self.beta = beta
        self.gamma = gamma
        self.delta = delta
        self.theta = theta
        self.th_delay = th_delay
        self.etaCT = etaCT
        self.N = N
        self.agents = np.zeros((N, 3))
        self.agents[:, 0] = PROG_S
        self.agents[:, 2] = np.inf

        # No known infector for anyone
        self.agents[:, 1] = -1

    def reset(self, N_I):

        self.agents[:, 0] = PROG_S
        # Pick N_I individuals at random
        ii = np.random.choice(range(self.N), size=N_I, replace=False)
        self.agents[ii, 0] = PROG_I

    def gillespie(self, tmax, samples=100):

        trajs = []
        for s_i in range(samples):
            trajs .append(_gillespie(tmax, self.agents, self.beta, self.gamma, self.delta,
                                     self.theta, self.th_delay, self.etaCT))

        return trajs
