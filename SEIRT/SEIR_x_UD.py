import numpy as np
from numba import jit
from math import sqrt
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from cpyment import CModel
from glob import glob
import argparse
import yaml
import time
import logging as log
import sys

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

log.basicConfig(stream=sys.stdout, level=log.INFO,
                format='%(asctime)s - %(levelname)s - %(message)s')


@jit(nopython=True)
def count_states(states, diagnosed):
    SU = np.sum((states == STATE_S)*(diagnosed == 0))
    SD = np.sum((states == STATE_S)*diagnosed)
    EU = np.sum((states == STATE_E)*(diagnosed == 0))
    ED = np.sum((states == STATE_E)*diagnosed)
    IU = np.sum((states == STATE_I)*(diagnosed == 0))
    ID = np.sum((states == STATE_I)*diagnosed)
    RU = np.sum((states == STATE_R)*(diagnosed == 0))
    RD = np.sum((states == STATE_R)*diagnosed)
    return SU, SD, EU, ED, IU, ID, RU, RD


@jit(nopython=True)
def count_pcis(states, diagnosed, contactM):
    SUi = np.where((states == STATE_S)*(1-diagnosed))[0]
    SDi = np.where((states == STATE_S)*diagnosed)[0]
    EUi = np.where((states == STATE_E)*(1-diagnosed))[0]
    EDi = np.where((states == STATE_E)*diagnosed)[0]
    IUi = np.where((states == STATE_I)*(1-diagnosed))[0]
    IDi = np.where((states == STATE_I)*diagnosed)[0]
    RUi = np.where((states == STATE_R)*(1-diagnosed))[0]
    RDi = np.where((states == STATE_R)*diagnosed)[0]

    N = states.shape[0]
    SUpci = np.sum(contactM[SUi])/(len(SUi)*N)
    SDpci = np.sum(contactM[SDi])/(len(SDi)*N)
    EUpci = np.sum(contactM[EUi])/(len(EUi)*N)
    EDpci = np.sum(contactM[EDi])/(len(EDi)*N)
    IUpci = np.sum(contactM[IUi])/(len(IUi)*N)
    IDpci = np.sum(contactM[IDi])/(len(IDi)*N)
    RUpci = np.sum(contactM[RUi])/(len(RUi)*N)
    RDpci = np.sum(contactM[RDi])/(len(RDi)*N)

    return SUpci, SDpci, EUpci, EDpci, IUpci, IDpci, RUpci, RDpci


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
                     chi=0,
                     return_pcis=False):

    # Generate states
    states = np.zeros(N, dtype=np.int8)
    # Generate diagnosed state
    diagnosed = np.zeros(N, dtype=np.bool8)
    # Generate traceable status
    traceable = np.zeros(N, dtype=np.bool8)
    # Contact matrix
    contactM = np.zeros((N, N), dtype=np.bool8)

    # Infect I0 patients
    istart = np.random.choice(np.arange(N), size=I0, replace=False)
    for i in istart:
        states[i] = STATE_I

    times = []
    traj = []
    pcis = []

    t = 0
    while t < tmax:

        counts = count_states(states, diagnosed)
        traj.append(counts)
        times.append(t)
        if return_pcis:
            pcis.append(np.sum(np.sum(contactM, axis=0) > 0)/N)
        else:
            pcis.append(0)

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
            contactM[si, ii] = True
            if np.random.random() <= beta:
                states[si] = STATE_E
        elif rn < wp[1]:
            # E becomes I
            ei = random_agent_i(states, diagnosed, STATE_E)
            states[ei] = STATE_I
        elif rn < wp[2]:
            # I becomes R
            ii = random_agent_i(states, diagnosed, STATE_I)
            contactM[:, ii] = False
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
    if return_pcis:
        pcis.append(np.sum(np.sum(contactM, axis=0) > 0)/N)
    else:
        pcis.append(0)

    return times, traj, pcis


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

    def __init__(self, N=1000, I0=0.01, tmax=100, tsteps=1000, **params):

        self.t = np.linspace(0, tmax, tsteps)
        self.N = N
        self.I0 = I0
        self.params = dict(self.default_params)
        self.params.update((k, v) for (k, v) in params.items()
                           if k in self.default_params)

    def run_abm(self, samples=10, return_pcis=False):
        trajs = []

        for i in range(samples):
            log.info('Sampling trajectory {0}/{1}'.format(i+1, samples))
            traj = seirxud_abm_gill(tmax=self.t[-1],
                                    N=self.N,
                                    I0=int(self.N*self.I0),
                                    return_pcis=return_pcis,
                                    **self.params)
            trajs.append((traj[0], np.array(traj[1]).T, traj[2]))

        # Do interpolation and averaging
        intptrajs = np.array([interp1d(tr[0], tr[1], kind='previous',
                                       bounds_error=False,
                                       fill_value=(tr[1][:, 0], tr[1][:, -1])
                                       )(self.t)
                              for tr in trajs])
        if return_pcis:
            intppcis = np.array([interp1d(tr[0], tr[2], kind='previous',
                                          bounds_error=False,
                                          fill_value=(tr[2][0], tr[2][-1])
                                          )(self.t)
                                 for tr in trajs])
            return intptrajs, intppcis
        else:
            return intptrajs

    def make_cmodel(self, etadamp=1, has_memory_states=False):

        states = ['SU', 'SD', 'EU', 'ED', 'IU', 'ID', 'RU', 'RD']

        if has_memory_states:
            states += ['CIS', 'CIE', 'CII', 'CIR']

        cm = CModel(states)

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
        cm.set_coupling_rate('SD:SD=>SU', kappa)

        cm.set_coupling_rate('EU:EU=>IU', alpha)
        cm.set_coupling_rate('ED:ED=>ID', alpha)

        cm.set_coupling_rate('IU:IU=>RU', gamma)
        cm.set_coupling_rate('ID:ID=>RD', gamma)

        cm.set_coupling_rate('RD:RD=>RU', kappa)

        cm.set_coupling_rate('EU:EU=>ED', eta*chi*theta)
        cm.set_coupling_rate('IU:IU=>ID', theta*(1+eta*chi))
        # Now the stuff that depends on memory
        if has_memory_states:

            cm.set_coupling_rate('IU*SU:=>CIS', c*(1-beta)/N)
            cm.set_coupling_rate('IU*CIS:CIS=>CIE', c*beta/N)
            cm.set_coupling_rate('CIS:CIS=>', gamma+theta*(1+eta*chi))

            cm.set_coupling_rate('IU*SU:=>CIE', c*beta/N)
            cm.set_coupling_rate('CIE:CIE=>CII', alpha)
            cm.set_coupling_rate('CIE:CIE=>', gamma+theta*(1+eta*chi))

            cm.set_coupling_rate('IU*IU:=>CII', c/N)
            cm.set_coupling_rate('CII:CII=>CIR', gamma)
            cm.set_coupling_rate('CII:CII=>', gamma+theta*(1+eta*chi))

            cm.set_coupling_rate('IU*RU:=>CIR', c/N)
            cm.set_coupling_rate('CIR:CIR=>', gamma+theta*(1+eta*chi))

            cm.set_coupling_rate('CIS:SU=>SD', chi*eta*theta)
            cm.set_coupling_rate('CIS*CIS:SU=>SD', chi*(1-(1-eta)**2)*theta/N)
            cm.set_coupling_rate('CIR:RU=>RD', chi*eta*theta)
            cm.set_coupling_rate('CIR*CIR:RU=>RD', chi*(1-(1-eta)**2)*theta/N)
        else:
            cm.set_coupling_rate('SU*IU:SU=>SD', chi*eta *
                                 c/(gamma+theta*(1+eta*chi))*(1-beta)*theta/N)
            cm.set_coupling_rate('RU*IU:RU=>RD', chi*eta *
                                 c/(gamma+theta*(1+eta*chi))*theta/N)

        self.cm = cm

    def run_cmodel(self, t0=0, y0=None, etadamp=1, has_memory_states=False):

        t0i = np.where(self.t >= t0)[0][0]

        self.make_cmodel(etadamp, has_memory_states)

        L = 8 + 4*has_memory_states
        if y0 is None:
            y0 = np.zeros(L)
            y0[INDEX_SU] = self.N*(1-self.I0)
            y0[INDEX_IU] = self.N*self.I0
        else:
            if L != len(y0):
                raise ValueError(('Invalid size of starting state! Use {0}'
                                  ' for simulations with{1} memory states').format(L,
                                                                                   '' if has_memory_states else 'out'))

        return self.cm.integrate(self.t[t0i:], y0)


if __name__ == '__main__':
    parser = argparse.ArgumentParser("SEIRT Agent-Based Simulator")
    parser.add_argument("-N", type=int, default=1000, help="Population size")
    parser.add_argument("-I", type=int, default=10, help="Infected at start")
    parser.add_argument("-c", type=float, default=13.0, help="Contact rate")
    parser.add_argument("-b", "--beta", type=float,
                        default=0.025, help="Infection per contact")
    parser.add_argument("-a", "--alpha", type=float,
                        default=0.2, help="Progression rate E -> I")
    parser.add_argument("-g", "--gamma", type=float,
                        default=0.1, help="Progression I->R")
    parser.add_argument("-t", "--theta", type=float,
                        default=0.1, help="Testing rate")
    parser.add_argument("-e", "--eta", type=float,
                        default=0.75, help="Tracing efficiency")
    parser.add_argument("-x", "--chi", type=float,
                        default=0.25, help="Testing rate")
    parser.add_argument("--tmax", type=float, default=100.0,
                        help="Simulation end time")
    parser.add_argument("--steps", type=int, default=1000,
                        help="Time steps (ODEs)")
    parser.add_argument("--ode", default=False,
                        action="store_true", help="Run ODE simulation")
    parser.add_argument("--abm", default=False, action="store_true",
                        help="Run mechanistic ABM simulation")
    parser.add_argument("--samples", type=int, default=10,
                        help="Number of samples for ABM")
    parser.add_argument("--yaml", type=str, default=None,
                        help="Read parameters from file")
    parser.add_argument("--seed", type=int,
                        default=time.time(), help="Random seed")
    parser.add_argument("-o", "--output", type=str,
                        default="simdata", help="Output file")
    parser.add_argument("--plot", default=False,
                        action="store_true", help="Generate plots")
    parser.add_argument("--memory", default=False,
                        action="store_true", help="Use memory states in ODE model")
    
    args = parser.parse_args()

    params = {
        "N": args.N, "I0": float(args.I)/args.N,
        "c": args.c, "beta": args.beta, "alpha": args.alpha, "gamma": args.gamma,
        "theta": args.theta, "eta": args.eta, "chi": args.chi,
        "tmax": args.tmax, "tsteps": args.steps
    }

    if args.yaml is not None:
        with open(args.yaml) as fp:
            ydata = yaml.load(fp)
            params.update(ydata.get("sim", {}))
            args.output = ydata.get("meta", {}).get("output", args.output)
            args.seed = ydata.get("meta", {}).get("seed", time.time())

    R0 = params["beta"]*params["c"]/params["gamma"]
    N = params["N"]

    log.info("Params: {0}".format(params))
    log.info("R0={0}".format(R0))

    np.random.seed(int(args.seed))

    if args.ode:
        sim = SEIRxUD(**params)
        log.info("Running deterministic model with contact tracing")
        traj = sim.run_cmodel(has_memory_states=args.memory)
        sim.t.dump("{0}.t".format(args.output))
        traj["y"].dump("{0}.y".format(args.output))
        log.info("Running deterministic model without contact tracing")
        trajNoCT = sim.run_cmodel(etadamp=0, has_memory_states=args.memory)
        trajNoCT["y"].dump(args.output + ".n")
    if args.abm:
        sim = SEIRxUD(**params)
        log.info("Running mechanistic agent-based model")
        trajs = sim.run_abm()  # args.samples)
        sim.t.dump("{0}.t".format(args.output))
        trajs.dump("{0}.trajs".format(args.output))
    if args.plot:
        log.info("Generating plots")
        t = np.load(args.output + ".t", allow_pickle=True)
        traj = {"y": np.load(args.output + ".y", allow_pickle=True)}
        trajNoCT = {"y": np.load(args.output + ".n", allow_pickle=True)}

        trajs = [np.load(f, allow_pickle=True)
                 for f in glob(args.output + "*.trajs")]
        trajsGill = np.vstack(trajs)
        trajsGillAvg = np.average(trajsGill, axis=0)
        trajsGillStd = np.std(trajsGill, axis=0)

        ##
        # plot trajectories and envelope
        ##
        fig, ax = plt.subplots()
        colors = {
            'S': np.array((0.0, 0.5, 1.0)),
            'E': np.array((1.0, 0.6, 0.0)),
            'I': np.array((1.0, 0.2, 0.1)),
            'R': np.array((0.0, 0.7, 0.4)),
        }

        stdn = 1
        t0i = 0

        for i, s in enumerate('SEIR'):
            for j, d in enumerate('UD'):
                col = colors[s]*(1 if j == 0 else 0.5)
                ax.plot(t[t0i:], traj['y'][:, 2*i+j], c=col, label=s+d, lw=2.0)
                ax.plot(t, trajNoCT['y'][:, 2*i+j], c=col, lw=2.0, ls='--')
                ax.plot(t, trajsGillAvg[2*i+j], c=col, lw=0.7)
                ax.fill_between(t,
                                trajsGillAvg[2*i+j]+stdn*trajsGillStd[2*i+j],
                                trajsGillAvg[2*i+j]-stdn*trajsGillStd[2*i+j],
                                color=list(col) + [0.3])

        ax.axhline((1-1/R0)*N, c=(0, 0, 0), lw=0.5, ls='--')
        ax.legend()
        plt.savefig(args.output + ".png")

        ##
        # plot diagnosed comparison
        ##
        fig, ax = plt.subplots()
        ax.set_title('Percentage of diagnosed people')
        ax.set_xlabel('t')
        ax.set_ylabel('*D (%)')
        ax.plot(t[t0i:], np.sum(traj['y'][:, 1::2], axis=1)*100/N, label='ODE')
        ax.plot(t, np.sum(trajsGillAvg[1::2], axis=0)*100/N, label='ABM')
        ax.legend()
        plt.savefig(args.output + "-diagnosed.png")

    log.info("Done")
