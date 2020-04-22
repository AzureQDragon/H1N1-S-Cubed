# Quarantine Model, Daniel Moon and Cy Pabis, 4/9/20
# Importing numpy and matplotlib
import numpy as np
import matplotlib.pyplot as plt 

# Euler Function for Graphing
def Euler(S0, I0, Q0, R0, dt):
    global s, i, q, r

    # Setting the constant variables
    b = .1
    a_0 = 25 * (10**-5)
    d_0 = 1/6
    d = d_0
    a = a_0
    g = (1/6) * (1/.8)
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
    P = .01147*I0
    # Setting up relations for systems of equations and graphing
    for i in range(1, len(t)):

        a = .05741*.00025*.5025 + .00025*.4975
        print(a)
        #a =  .05741*(.00025)*(.01147)*(infect[i - 1]) + .00025*(.98853*infect[i - 1]) + .05741*(.00025)*(quar[i - 1])
        #d = d_prev * (1-P) + g * P
        #d_prev = d 
        suscept[i] =  (suscept[i-1]) - dt*(a * suscept[i-1] * infect[i-1])
        infect[i] = infect[i-1]+ dt*(a*suscept[i-1] * infect[i-1] - b*infect[i-1] - d*infect[i-1])
        quar[i] = quar[i-1] + dt*(b*infect[i-1] - g * quar[i-1])
        recover[i] = recover[i-1] + dt*(d*infect[i-1] + g*quar[i-1])

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


# Running the Euler function
Euler(S0, I0, Q0, R0, .01)
