import numpy as np
from numba import jit

# States
STATE_S = 0
STATE_E = 1
STATE_I = 2
STATE_R = 3
STATE_T = 4


@jit(nopython=True)
def count_states(states):
    S = np.sum((states == STATE_S))
    E = np.sum((states == STATE_E))
    I = np.sum((states == STATE_I))
    R = np.sum((states == STATE_R))
    T = np.sum((states == STATE_T))
    return S, E, I, R, T


@jit(nopython=True)
def random_agent_i(states, tstate):
    return np.random.choice(np.where(states == tstate)[0])


@jit(nopython=True)
def seirt_abm_gill(tmax, N, I0, c, beta, alpha, gamma, theta, kappa, eta, chi):

    # Generate states
    states = np.zeros(N)
    # Generate infector record
    infectors = -np.ones(N)
    # Generate traceable status
    traceable = np.zeros(N)

    # Infect I0 patients
    istart = np.random.choice(np.arange(N), size=I0, replace=False)
    for i in istart:
        states[i] = STATE_I

    times = []
    Ss = []
    Es = []
    Is = []
    Rs = []
    Ts = []

    t = 0
    while t < tmax:

        S, E, I, R, T = count_states(states)
        Ss.append(S)
        Es.append(E)
        Is.append(I)
        Rs.append(R)
        Ts.append(T)
        times.append(t)

        # Traceable agents?
        CT = np.sum(traceable)

        # Probability of various transitions
        wSE = beta*c*S*I/N
        wEI = alpha*E
        wIR = gamma*I
        wIT = theta*I
        wTR = kappa*T
        wCT = chi*CT

        W = wSE + wEI + wIR + wIT + wTR + wCT
        if W <= 0:
            break

        p = np.cumsum(np.array([wSE, wEI, wIR, wIT, wTR, wCT]))/W

        dt = -np.log(np.random.random())/W

        rn = np.random.random()
        if rn < p[0]:
            # Infect a random S
            si = random_agent_i(states, STATE_S)
            ii = random_agent_i(states, STATE_I)
            states[si] = STATE_E
            infectors[si] = ii
        elif rn < p[1]:
            # E becomes I
            ei = random_agent_i(states, STATE_E)
            states[ei] = STATE_I
        elif rn < p[2]:
            # I becomes R
            ii = random_agent_i(states, STATE_I)
            states[ii] = STATE_R
            traceable[ii] = False
        elif rn < p[3]:
            # I becomes T
            ii = random_agent_i(states, STATE_I)
            states[ii] = STATE_T
            traceable[ii] = False
            # Also set all those who have it as an infector as traceable
            ctis = np.where(infectors == ii)[0]
            traceable[ctis] = (np.random.random(len(ctis)) < eta)
        elif rn < p[4]:
            # T becomes R
            ti = random_agent_i(states, STATE_T)
            traceable[ti] = False
            states[ti] = STATE_R
        elif rn < p[5]:
            # Random traceable?
            cti = np.random.choice(np.where(traceable)[0])
            traceable[cti] = False
            states[cti] = STATE_T

        t += dt

    S, E, I, R, T = count_states(states)
    Ss.append(S)
    Es.append(E)
    Is.append(I)
    Rs.append(R)
    Ts.append(T)
    times.append(t)

    return times, Ss, Es, Is, Rs, Ts
