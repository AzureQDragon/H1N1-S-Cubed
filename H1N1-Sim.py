# Quarantine Model, Daniel Moon and Cy Pabis, 4/9/20
# Importing numpy and matplotlib
import numpy as np
import matplotlib.pyplot as plt 

# Euler Function for Graphing
def Euler(S0, I0, Q0, R0, dt):
    global s, i, q, r

    # Setting the constant variables
    b = .1
    a =  (1.8/31.5) *.00025*.5025 + .00025*.4975
    print(a)
    d = .5025*1.25*(1/6) + .4975*(1/6)
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

        suscept[i] =  (suscept[i-1]) - dt*(a * suscept[i-1] * infect[i-1])
        infect[i] = infect[i-1]+ dt*(a*suscept[i-1] * infect[i-1] - b*infect[i-1] - d*infect[i-1])
        quar[i] = quar[i-1] + dt*(b*infect[i-1] - g * quar[i-1])
        recover[i] = recover[i-1] + dt*(d*infect[i-1] + g*quar[i-1])
    
    # Find the maximums of infected population and quarantine population and print
    print(infect.max())
    print(quar.max())
    print(len(quar))
    # Find the index of 
    for i in range(5, 10000):
        if infect[i] <= 1:
            output = i
            break
    print(output)
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
