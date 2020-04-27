import numpy as np
from numba import jit
from scipy.interpolate import interp1d
from cpyment import CModel

# States
STATE_S = 0
STATE_E = 1
STATE_I = 2
STATE_R = 3

INDEX_SU = 0
INDEX_EU = 2
INDEX_IU = 4
INDEX_RU = 6

INDEX_SD = 1
INDEX_ED = 3
INDEX_ID = 5
INDEX_RD = 7


@jit(nopython=True)
def count_states(states, diagnosed):
    SU = np.sum((states == STATE_S)*(1-diagnosed))
    SD = np.sum((states == STATE_S)*diagnosed)
    EU = np.sum((states == STATE_E)*(1-diagnosed))
    ED = np.sum((states == STATE_E)*diagnosed)
    IU = np.sum((states == STATE_I)*(1-diagnosed))
    ID = np.sum((states == STATE_I)*diagnosed)
    RU = np.sum((states == STATE_R)*(1-diagnosed))
    RD = np.sum((states == STATE_R)*diagnosed)
    return SU, SD, EU, ED, IU, ID, RU, RD


@jit(nopython=True)
def random_agent_i(states, diagnosed, tstate, tdiag=None):
    if tdiag is None:
        return np.random.choice(np.where(states == tstate)[0])
    else:
        return np.random.choice(np.where((states == tstate) *
                                         (diagnosed == tdiag))[0])


@jit(nopython=True)
def seirxud_abm_gill(tmax=10,
                     N=1000,
                     I0=10,
                     c=5,
                     beta=0.05,
                     alpha=0.2,
                     gamma=0.1,
                     theta=0.0,
                     kappa=0.05,
                     eta=0,
                     chi=0):

    # Generate states
    states = np.zeros(N)
    # Generate diagnosed state
    diagnosed = np.zeros(N)
    # Generate traceable status
    traceable = np.zeros(N)
    # Contact matrix
    contactM = np.zeros((N, N))

    # Infect I0 patients
    istart = np.random.choice(np.arange(N), size=I0, replace=False)
    for i in istart:
        states[i] = STATE_I

    times = []
    traj = []

    t = 0
    while t < tmax:

        counts = count_states(states, diagnosed)
        traj.append(counts)
        times.append(t)

        # Traceable agents?
        CT = np.sum(traceable)

        E = counts[INDEX_EU] + counts[INDEX_ED]
        I = counts[INDEX_IU] + counts[INDEX_ID]

        # Possible contacts
        wSIc = c*counts[INDEX_SU]*counts[INDEX_IU]/N
        # E becomes I
        wEI = alpha*E
        # I becomes R
        wIR = gamma*I
        # I is diagnosed
        wIUID = theta*counts[INDEX_IU]
        # Diagnosed S is released
        wSDSU = kappa*counts[INDEX_SD]
        # Diagnosed R is released
        wRDRU = kappa*counts[INDEX_RD]
        # Someone who's traceable gets quarantined
        wCT = chi*CT

        Wtot = wSIc + wEI + wIR + wIUID + wSDSU + wRDRU + wCT
        if Wtot <= 0:
            break

        wp = np.array([wSIc, wEI, wIR, wIUID, wSDSU, wRDRU, wCT])
        wp = np.cumsum(wp)/Wtot

        dt = -np.log(np.random.random())/Wtot

        rn = np.random.random()

        if rn < wp[0]:
            # Contact between a random SU and a random IU
            si = random_agent_i(states, diagnosed, STATE_S, False)
            ii = random_agent_i(states, diagnosed, STATE_I, False)
            contactM[si][ii] = True
            if np.random.random() <= beta:
                states[si] = STATE_E
        elif rn < wp[1]:
            # E becomes I
            ei = random_agent_i(states, diagnosed, STATE_E)
            states[ei] = STATE_I
        elif rn < wp[2]:
            # I becomes R
            ii = random_agent_i(states, diagnosed, STATE_I)
            states[ii] = STATE_R
        elif rn < wp[3]:
            # Diagnosis
            ii = random_agent_i(states, diagnosed, STATE_I, False)
            diagnosed[ii] = True
            traceable[ii] = False
            # Also set all those who have it as an infector as traceable
            ctis = np.where((contactM[:, ii]*(diagnosed == 0)))[0]
            traceable[ctis] = np.logical_or(traceable[ctis],
                                            np.random.random(len(ctis)) < eta)
        elif rn < wp[4]:
            si = random_agent_i(states, diagnosed, STATE_S, True)
            diagnosed[si] = False
        elif rn < wp[5]:
            ri = random_agent_i(states, diagnosed, STATE_R, True)
            diagnosed[ri] = False
        else:
            # Contact tracing
            # Random traceable?
            cti = np.random.choice(np.where(traceable)[0])
            diagnosed[cti] = True
            traceable[cti] = False

        t += dt

    counts = count_states(states, diagnosed)
    traj.append(counts)
    times.append(t)

    return times, traj


class SEIRxUD(object):

    default_params = {
        'c': 4,
        'beta': 0.3,
        'alpha': 0.1,
        'gamma': 0.3,
        'theta': 0.0,
        'kappa': 0.03,
        'eta': 0.0,
        'chi': 0.7
    }

    def __init__(self, N=1000, I0=0.01, tmax=100, tsteps=1000, params={}):

        self.t = np.linspace(0, tmax, tsteps)
        self.N = N
        self.I0 = I0
        self.params = dict(self.default_params)
        self.params.update(params)

    def run_abm(self, samples=10):
        trajs = []

        for i in range(samples):
            print('{0}/{1}'.format(i+1, samples))
            traj = seirxud_abm_gill(tmax=self.t[-1],
                                    N=self.N,
                                    I0=int(self.N*self.I0),
                                    **self.params)
            trajs.append((traj[0], np.array(traj[1]).T))

        # Do interpolation and averaging
        intptrajs = np.array([interp1d(tr[0], tr[1], kind='previous',
                                       bounds_error=False,
                                       fill_value=(tr[1][:, 0], tr[1][:, -1])
                                       )(self.t)
                              for tr in trajs])

        return intptrajs

    def make_cmodel(self, etadamp=1):

        cm = CModel(['SU', 'SD', 'EU', 'ED', 'IU', 'ID', 'RU', 'RD'])

        N = self.N
        beta = self.params['beta']
        c = self.params['c']
        chi = self.params['chi']
        eta = self.params['eta']*etadamp
        alpha = self.params['alpha']
        gamma = self.params['gamma']
        theta = self.params['theta']
        kappa = self.params['kappa']

        cm.set_coupling_rate('SU*IU:SU=>EU', beta*c/N)
        cm.set_coupling_rate('SU*IU:SU=>SD', chi*eta*c/gamma*(1-beta)*theta/N)
        cm.set_coupling_rate('SD:SD=>SU', kappa)

        cm.set_coupling_rate('EU:EU=>IU', alpha)
        cm.set_coupling_rate('ED:ED=>ID', alpha)
        cm.set_coupling_rate('EU:EU=>ED', eta*chi*theta)

        cm.set_coupling_rate('IU:IU=>RU', gamma)
        cm.set_coupling_rate('ID:ID=>RD', gamma)
        cm.set_coupling_rate('IU:IU=>ID', theta*(1+eta*chi))

        cm.set_coupling_rate('RU*IU:RU=>RD', chi*eta*c/gamma*theta/N)
        cm.set_coupling_rate('RD:RD=>RU', kappa)

        self.cm = cm

    def run_cmodel(self, t0=0, y0=None,  etadamp=1):

        t0i = np.where(self.t >= t0)[0][0]

        self.make_cmodel(etadamp)

        if y0 is None:
            y0 = np.zeros(8)
            y0[INDEX_SU] = self.N*(1-self.I0)
            y0[INDEX_IU] = self.N*self.I0

        return self.cm.integrate(self.t[t0i:], y0)
