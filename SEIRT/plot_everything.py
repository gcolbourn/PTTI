import argparse
import yaml

import plot_r_eta
import plot_r_theta
import plot_eta_theta

if __name__ == '__main__':
    parser = argparse.ArgumentParser("SEIRT Agent-Based Simulator")
    parser.add_argument("-N", type=int, default=1000, help="Population size")
    parser.add_argument("-I", type=int, default=10, help="Infected at start")
    parser.add_argument("-c", type=float, default=13.0, help="Contact rate")
    parser.add_argument("-b", "--beta", type=float,
                        default=0.033, help="Infection per contact")
    parser.add_argument("-a", "--alpha", type=float,
                        default=0.2, help="Progression rate E -> I")
    parser.add_argument("-g", "--gamma", type=float,
                        default=1.0/7, help="Progression I->R")
    parser.add_argument("-t", "--theta", type=float,
                        default=1.0/7, help="Testing rate")
    parser.add_argument("-e", "--eta", type=float,
                        default=0.25, help="Tracing efficiency")
    parser.add_argument("-x", "--chi", type=float,
                        default=1.0, help="Testing rate")
    parser.add_argument("-k", "--kappa", type=float,
                        default=1.0/14, help="Isolation time")
    parser.add_argument("--tmax", type=float, default=600.0,
                        help="Simulation end time")
    parser.add_argument("--steps", type=int, default=1000,
                        help="Time steps (ODEs)")
    parser.add_argument("--yaml", type=str, default=None,
                        help="Read parameters from file")
    parser.add_argument("-o", "--output", type=str,
                        default="simdata", help="Output file")
    parser.add_argument("--gridsize", type=int,
                        default=10, help="Steps between parameter extremes")

    args = parser.parse_args()

    params = {
        "N": args.N, "I0": float(args.I)/args.N,
        "c": args.c, "beta": args.beta, "alpha": args.alpha, "gamma": args.gamma,
        "theta": args.theta, "eta": args.eta, "chi": args.chi,
        "tmax": args.tmax, "tsteps": args.steps, "gridsize": args.gridsize,
        "output": args.output
    }

    if args.yaml is not None:
        with open(args.yaml) as fp:
            ydata = yaml.load(fp)
            params.update(ydata.get("sim", {}))
            params["output"] = ydata.get("meta", {}).get("output", args.output)

    plot_r_eta.plot(**params)
    plot_r_theta.plot(**params)
    plot_eta_theta.plot(**params)
