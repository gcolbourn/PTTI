import numpy as np
from cpyment import CModel

# A constructor to turn our sensible parameters into rates for the CModel


def build_model(N=1000,             # Population
                beta=0.06,          # Contact term between S and various I
                c=5,                # Average number of daily contacts
                t_E=3.7,            # Time spent in E
                t_Ip=1.5,           # Time spent in Ip
                t_I=2.3,            # Time spent in normal infective phases
                f_Ia=0.5,           # Fraction of asymptomatic patients
                f_H=0.068,          # Probability of being hospitalised
                f_ICU=0.022,        # Probability of ending in the ICU
                t_H=19.11,          # Time spent hospitalised
                f_HR=0.89,          # Fraction hospitalised who recovers
                t_ICU=22.63,        # Time spent in ICU
                f_ICUR=0.8,         # Fraction in ICU who recovers
                theta0=0.0,         # Fraction of people tested daily regularly
                thetaI=0.1,         # Fraction of infected, symptomatic people who are tested
                # Efficiency of contact tracing (contacts traced per non-contact positive)
                etaCT=0.0,
                r=1,                # Recall of testing protocol
                s=1,                # Specificity of testing protocol
                t_Q=14.0,           # Time spent in quarantine unless symptomatic
                # Efficiency of Quarantine-Quarantine isolation (if 1 => no infection)
                eta_QQ=1.0,
                # Efficiency of Quarantine-NonQuarantine isolation (if 1 => no infection)
                eta_QNQ=1.0
                ):

    # States
    states = ['S', 'E', 'Ip', 'Ia', 'Is', 'R',
              'TS', 'TE', 'TIp', 'TIa', 'TIs', 'TR',
              'H', 'ICU', 'D']

    cm = CModel(states)

    # INFECTION PROCESSES

    # Outside of tracing
    cm.set_coupling_rate('S*Ip:S=>E', beta*c/N)
    cm.set_coupling_rate('S*Ia:S=>E', beta*c/N)
    cm.set_coupling_rate('S*Is:S=>E', beta*c/N)

    # Between traced/quarantined cases
    cm.set_coupling_rate('TS*TIp:TS=>TE', (1-eta_QQ)*beta*c/N)
    cm.set_coupling_rate('TS*TIa:TS=>TE', (1-eta_QQ)*beta*c/N)
    cm.set_coupling_rate('TS*TIs:TS=>TE', (1-eta_QQ)*beta*c/N)

    # Cross-infection rates
    cm.set_coupling_rate('TS*Ip:TS=>TE', (1-eta_QNQ)*beta*c/N)
    cm.set_coupling_rate('TS*Ia:TS=>TE', (1-eta_QNQ)*beta*c/N)
    cm.set_coupling_rate('TS*Is:TS=>TE', (1-eta_QNQ)*beta*c/N)
    cm.set_coupling_rate('S*TIp:TS=>TE', (1-eta_QNQ)*beta*c/N)
    cm.set_coupling_rate('S*TIa:TS=>TE', (1-eta_QNQ)*beta*c/N)
    cm.set_coupling_rate('S*TIs:TS=>TE', (1-eta_QNQ)*beta*c/N)

    # PROGRESSION

    # E->Ip
    cm.set_coupling_rate('E:E=>Ip', 1.0/t_E)
    cm.set_coupling_rate('TE:TE=>TIp', 1.0/t_E)

    # Ip->Ia/Is
    cm.set_coupling_rate('Ip:Ip=>Ia', f_Ia/t_Ip)
    cm.set_coupling_rate('Ip:Ip=>Is', (1.0-f_Ia)/t_Ip)
    cm.set_coupling_rate('TIp:TIp=>TIa', f_Ia/t_Ip)
    cm.set_coupling_rate('TIp:TIp=>TIs', (1.0-f_Ia)/t_Ip)

    # Recovery from I states
    cm.set_coupling_rate('Ia:Ia=>R', 1.0/t_I)
    cm.set_coupling_rate('Is:Is=>R', (1.0-f_H-f_ICU)/t_I)
    cm.set_coupling_rate('TIa:TIa=>TR', 1.0/t_I)
    cm.set_coupling_rate('TIs:TIs=>TR', (1.0-f_H-f_ICU)/t_I)

    # Hospitalisation and ICU
    cm.set_coupling_rate('Is:Is=>H', f_H/t_I)
    cm.set_coupling_rate('TIs:TIs=>H', f_H/t_I)
    cm.set_coupling_rate('Is:Is=>ICU', f_ICU/t_I)
    cm.set_coupling_rate('TIs:TIs=>ICU', f_ICU/t_I)

    # Outcomes from hospital
    cm.set_coupling_rate('H:H=>R', f_HR/t_H)
    cm.set_coupling_rate('ICU:ICU=>R', f_ICUR/t_ICU)
    cm.set_coupling_rate('H:H=>D', (1-f_HR)/t_H)
    cm.set_coupling_rate('ICU:ICU=>D', (1-f_ICUR)/t_ICU)

    # TESTING AND TRACING

    n = c*(t_Ip+t_I)
    R0 = beta*n

    # S=>TS
    # False positives
    cm.set_coupling_rate('S:S=>TS', (1-s)*theta0)

    # Contact tracing, true positives
    cm.set_coupling_rate('S*Ip:S=>TS', (1-beta)*n*etaCT*thetaI*r/N)
    cm.set_coupling_rate('S*Ia:S=>TS', (1-beta)*n*etaCT*thetaI*r/N)
    cm.set_coupling_rate('S*Is:S=>TS', (1-beta)*n*etaCT*thetaI*r/N)

    # E=>TE
    # False positives + contact tracing
    cm.set_coupling_rate('E:E=>TE', (1-s)*theta0+n*etaCT*thetaI*r)

    # I=>TI
    for s1 in ('Ip', 'Ia', 'Is'):
        # True positives + contact tracing
        cm.set_coupling_rate('{0}:{0}=>T{0}'.format(s1), r*thetaI*(1+n*etaCT))

    # R=>TR
    # False positives
    cm.set_coupling_rate('R:R=>TR', (1-s)*theta0)

    # Contact tracing, true positives
    cm.set_coupling_rate('R*Ip:R=>TR', n*etaCT*thetaI*r/N)
    cm.set_coupling_rate('R*Ia:R=>TR', n*etaCT*thetaI*r/N)
    cm.set_coupling_rate('R*Is:R=>TR', n*etaCT*thetaI*r/N)

    # False positives for everyone
    for s1 in ('S', 'E', 'Ip', 'Ia', 'Is', 'R'):
        for s2 in ('S', 'E', 'R'):
            cm.set_coupling_rate('{0}*{1}:{0}=>T{0}',
                                 (1-s)*theta0*n*etaCT/N)

    # Quarantine expiration
    cm.set_coupling_rate('TS:TS=>S', 1.0/t_Q)
    cm.set_coupling_rate('TR:TR=>R', 1.0/t_Q)

    return cm


if __name__ == '__main__':
    m = build_model()
    m.pprint()
