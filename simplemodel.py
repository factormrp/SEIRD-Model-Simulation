from scipy.integrate import odeint
import numpy as np
import sys

def get_hyperparams():
    """
    Asks the user to enter in custom hyperparameter values
    """
    incub = input("What is the incubation rate? ")
    infec = input("What is the infection rate? ")
    R_not = input("What is the R_0 factor? ")
    N = input("What is the total population size? ")
    return incub,infec,R_not,N

def get_consts(incub=None,infec=None,R_0=None):
    """
    Takes some of the hyperparameters to compute the constants used in the system of diff eqs
    """
    # raise exception on improper use of the function
    if incub is None or infec is None or R_0 is None:
        raise Exception("Failed to pass in necessary hyperparameter in constants calculation")

    # compute constants
    alpha = 1/incub
    gamma = 1/infec
    beta = R_0*gamma
    return alpha,beta,gamma

def get_initials():
    """
    Asks the user to enter in the initial values to be used in the IVP
    """
    s_0 = input("Enter the inital population of susceptibles: ")
    e_0 = input("Enter the inital population of exposed: ")
    i_0 = input("Enter the inital population of infected: ")
    r_0 = input("Enter the inital population of recovered: ")
    return [s_0,e_0,i_0,r_0]

def plot_results(s=None,e=None,i=None,r=None,t=None):
    """
    Plots the relative population sizes over time
    """
    # import libraries and specify the backend
    import matplotlib
    matplotlib.use("TkAgg")
    import matplotlib.pyplot as plt 

    # raise exception on improper use of the function
    if s is None or e is None or i is None or r is None or t is None:
        raise Exception("Failed to pass in necessary data values in plotting")

    # plot the results
    plt.figure(figsize=(16,9))
    plt.subplot(2,1,1)
    plt.plot(t,s,color="blue",lw=2,label="susceptible")
    plt.plot(t,r,color="green",lw=2,label="recovered")
    plt.ylabel("Fraction of Pop")
    plt.legend()
    plt.subplot(2,1,2)
    plt.plot(t,e,color="orange",lw=2,label="exposed")
    plt.plot(t,i,color="red",lw=2,label="infected")
    plt.ylabel("Fraction of Pop")
    plt.xlabel("Time in Days")
    plt.legend()
    plt.show()


if __name__ == "__main__":
    
    # initialize hyperparameters
    incubation,infection,R_not,N = None,None,None,None

    # if user specified, then ask user for hyperparameters, otherwise use hardcoded values 
    if "--get-args" in sys.argv or "-g" in sys.argv:
        incubation,infection,R_not,N,u = get_hyperparams()
    else:
        incubation,infection,R_not,N,u = 5.1,3.3,2.4,33517,0.1 # hard-coded values

    # create derived constants for the system of diff eqs using the hyperparameters
    alpha,beta,gamma = get_consts(incub=incubation,infec=infection,R_0=R_not)

    # get initial values
    x_not = get_initials()
    
    # define the time range over which simulation is to be run
    t = np.linspace(0,200,101)

    # define system
    def sim(x,t):
        """
        Runs the SEIR model simulation over time
        """
        # raise expection on improper use of the function
        if alpha is None or beta is None or gamma is None or u is None:
            raise Exception("Failed to run simulation")

        # get populations
        s,e,i,r = x

        # initialize the system of diff eqs and define the equations
        dx = np.zeros(4)
        dx[0] = -(1-u)*beta*s*i
        dx[1] = (1-u)*beta*s*i - alpha*e
        dx[2] = alpha*e - gamma*i
        dx[3] = gamma*i
        return dx

    # run the SEIR model over the time range
    x = odeint(sim,x_not,t)
    s = x[:,0]
    e = x[:,1]
    i = x[:,2]
    r = x[:,3]

    # plot the results
    plot_results(s,e,i,r,t)
