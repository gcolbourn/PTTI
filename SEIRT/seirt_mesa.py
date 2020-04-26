import numpy as np
from mesa import Agent, Model
from mesa.time import RandomActivation
from mesa.datacollection import DataCollector

# States
STATE_S = 0
STATE_E = 1
STATE_I = 2
STATE_R = 3
STATE_T = 4


class SEIRTAgent(Agent):
    def __init__(self, unique_id, model, state=STATE_S):
        super().__init__(unique_id, model)

        self.state = state
        self.traced = False
        self.infected = []

    def step(self):

        # Now this depends on state
        if self.state == STATE_E:
            pi = 1-np.exp(-self.model.gamma*self.model.dt)
            if self.model.random.random() < pi:
                self.model.counts[STATE_E] -= 1
                self.model.counts[STATE_I] += 1
                self.state = STATE_I

            if self.traced:
                pt = 1-np.exp(-self.model.eta*self.model.dt)
                if self.model.random.random() < pi:
                    self.state = STATE_T
                    self.model.counts[STATE_E] -= 1
                    self.model.counts[STATE_T] += 1

        elif self.state == STATE_I:
            # First, try to infect
            S = self.model.counts[STATE_S]
            pi = 1-np.exp(-self.model.beta*self.model.dt *
                          S/self.model.num_agents)
            if self.model.random.random() < pi:
                # Pick one to infect!
                si = self.model.random.choice([i for i, a in
                                               enumerate(
                                                   self.model.schedule.agents)
                                               if a.state == STATE_S])
                self.model.schedule.agents[si].state = STATE_E
                self.model.counts[STATE_S] -= 1
                self.model.counts[STATE_E] += 1
                self.infected.append(si)

            # Then try to recover, or be traced
            pr = 1-np.exp(-self.model.delta*self.model.dt)
            if self.traced:
                pt = 1-np.exp(-(self.model.theta+self.model.eta)*self.model.dt)
            else:
                pt = 1-np.exp(-self.model.theta*self.model.dt)
            if self.model.random.random() < pr:
                self.state = STATE_R
                self.model.counts[STATE_I] -= 1
                self.model.counts[STATE_R] += 1
            elif self.model.random.random() < pt:
                self.state = STATE_T
                self.model.counts[STATE_I] -= 1
                self.model.counts[STATE_T] += 1
                # Label all infected individuals as traced!
                for i in self.infected:
                    self.model.schedule.agents[i].traced = True
        elif self.state == STATE_T:
            pr = 1-np.exp(-self.model.kappa*self.model.dt)
            if self.model.random.random() < pr:
                self.state = STATE_R
                self.model.counts[STATE_T] -= 1
                self.model.counts[STATE_R] += 1


class SEIRTModel(Model):

    def __init__(self, N,
                 I0=0.1,
                 beta=0.3,
                 gamma=0.1,
                 delta=0.1,
                 theta=0.1,
                 kappa=0.03,
                 eta=0.2,
                 dt=0.001):

        self.num_agents = N
        self.beta = beta
        self.gamma = gamma
        self.delta = delta
        self.theta = theta
        self.kappa = kappa
        self.eta = eta
        self.dt = dt
        self.t = 0

        self.schedule = RandomActivation(self)

        # Starting counts
        self.counts = [0]*5

        # Create agents
        for i in range(self.num_agents):
            state = STATE_I if self.random.random() < I0 else STATE_S
            a = SEIRTAgent(i, self, state)
            self.counts[state] += 1
            self.schedule.add(a)

        # Now, data collection
        self.datacollector = DataCollector(
            model_reporters={'t': 't',
                             'S': lambda m: m.counts[STATE_S],
                             'E': lambda m: m.counts[STATE_E],
                             'I': lambda m: m.counts[STATE_I],
                             'R': lambda m: m.counts[STATE_R],
                             'T': lambda m: m.counts[STATE_T]
                             }
        )

    def step(self):
        self.t += self.dt
        self.datacollector.collect(self)
        self.schedule.step()
