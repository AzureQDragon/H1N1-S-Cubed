# Quarantine Model, Daniel Moon and Cy Pabis, 4/9/20

# Importing numpy and matplotlib
import numpy as np
import matplotlib.pyplot as plt 

# Euler Function for Graphing
def Euler(S0, I0, Q0, R0, d, dt, a):
    global s, i, q, r

    # Setting the 
    b = .1

    # Setup of the t array
    t = np.arange(0, 100, dt)

    # Create arrays filled with zeros
    suscept = np.zeros(len(t))
    infect = np.zeros(len(t))
    quar = np.zeros(len(t))
    recover = np.zeros(len(t))

    # Set up the initial populations
    suscept[0] = S0
    infect[0] = I0
    quar[0] = Q0
    recover[0] = R0

    # Setting up relations for systems of equations and graphing
    for i in range(1, len(t)):

        suscept[i] =  (suscept[i-1]) - dt*(a * suscept[i-1] * infect[i-1])
        infect[i] = infect[i-1]+ dt*(a*suscept[i-1] * infect[i-1] - b*infect[i-1] - d*infect[i-1])
        quar[i] = quar[i-1] + dt*(b*infect[i-1] - (1/6) * quar[i-1])
        recover[i] = recover[i-1] + dt*(d*infect[i-1] + (1/6) *quar[i-1])

    # Graph the individual systems 
    s, = plt.plot(t, suscept)
    i, = plt.plot(t, infect)
    q, = plt.plot(t, quar)
    r, = plt.plot(t, recover)

    # Labels for graph legend and labels
    plt.xlabel("Time (Days)")
    plt.ylabel("Population")
    plt.legend(("Susceptible Pop", "Infected Pop", "Quarantined Pop", "Recovered Pop"))
    plt.title("H1N1 Infection")
    plt.grid(True)
    plt.show()
    return


# Setting starting values for suscept, infect, quar, recover, and the recovery rate
S0 = 18223
I0 = 11
Q0 = 0
R0 = 0
d = 1/6
a = 25 * (10**-5)
# Running the Euler function
Euler(S0, I0, Q0, R0, d, .01, a)