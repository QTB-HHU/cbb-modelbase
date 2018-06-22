import modelbase
import model
import parameters
import reactions
import numpy as np
import matplotlib.pyplot as plt


def plotSteadyState(t, y, m, groups):
    legend = m.cpdNames
    nrows = int(np.ceil((len(groups))/2))
    fig, ax = plt.subplots(nrows, 2, figsize=[16, 3*nrows])
    ax = ax.ravel()
    for plot, g in enumerate(groups):
        for i in g:
            ax[plot].plot(t, y[:, legend.index(i)], label=i)
        if plot % 2 == 0:
            ax[plot].legend(bbox_to_anchor=[-0.15, 1], loc="upper right", borderaxespad=0)
        else:
            ax[plot].legend(bbox_to_anchor=[1.15, 1], loc="upper left", borderaxespad=0)
            ax[plot].yaxis.tick_right()
    plt.tight_layout()


if __name__ == '__main__':
    # Poolman model
    r = reactions.Reactions()
    p = parameters.Parameters
    m = model.Poolman(p, r)

    init = {"PGA": 0.6437280277346407,
            "BPGA": 0.001360476366780556,
            "GAP": 0.011274125311289358,
            "DHAP": 0.24803073890728228,
            "FBP": 0.019853938009873073,
            "F6P": 1.0950701164493861,
            "G6P": 2.5186612678035734,
            "G1P": 0.14608235353185037,
            "SBP": 0.09193353265673603,
            "S7P": 0.23124426886012006,
            "E4P": 0.028511831060903877,
            "X5P": 0.036372985623662736,
            "R5P": 0.06092475016463224,
            "RUBP": 0.24993009253928708,
            "RU5P": 0.02436989993734177,
            "ATP": 0.43604115800259613}
    y0 = np.array([init[i] for i in m.cpdNames])

    s = modelbase.Simulator(m)
    s.integrator.verbosity = 50
    s.timeCourse(np.linspace(0, 200, 1000), y0)

    t, y = s.results[0]["t"], s.results[0]["y"]

    groups = [["G6P", "PGA", "F6P", "S7P"],
              ["RU5P", "X5P", "SBP", "G1P"],
              ["ATP", "DHAP", "RUBP", "R5P"],
              ["BPGA", "E4P", "FBP", "GAP"]]

    plotSteadyState(t, y, m, groups)

    # Poolman model with dynamic NADPH
    r = reactions.ReactionsNADPH()
    p = parameters.ParametersNADPH
    m = model.PoolmanNADPH(p, r)

    init = {"PGA": 0.599645270373,
            "BPGA": 0.000907499521924,
            "GAP": 0.011839616887,
            "DHAP": 0.260471552645,
            "FBP": 0.021895569623,
            "F6P": 1.2456290719,
            "G6P": 2.86494686535,
            "G1P": 0.166166918189,
            "SBP": 0.1120019621,
            "S7P": 0.233467059202,
            "E4P": 0.0330766864679,
            "X5P": 0.0374527459593,
            "R5P": 0.0627333486958,
            "RUBP": 0.261466058509,
            "RU5P": 0.0250933393445,
            "ATP": 0.414993685612,
            "NADPH": 0.281543418344}

    y0 = np.array([init[i] for i in m.cpdNames])

    s = modelbase.Simulator(m)
    s.set_initial_value(y0)
    s.integrator.verbosity = 50
    s.timeCourse(np.linspace(0, 200, 1000), y0)
    t, y = s.results[0]["t"], s.results[0]["y"]

    groups = [["G6P", "PGA", "F6P", "S7P"],
              ["RU5P", "X5P", "SBP", "G1P"],
              ["ATP", "DHAP", "RUBP", "R5P", "NADPH"],
              ["BPGA", "E4P", "FBP", "GAP"]]

    plotSteadyState(t, y, m, groups)
    plt.show()
