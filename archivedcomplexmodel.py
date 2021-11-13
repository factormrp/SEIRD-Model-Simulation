def estimateParams(N,E,I,R,D):
    """
    Defines a whole set of parameters to search best fit from
    """
    # import dependencies
    import numpy as np
    from lmfit import Parameters

    # get values for susceptible population group
    S = [N - (e+i+r+d) for e,i,r,d in zip(E,I,R,D)]

    # initialize lists
    dE,dI,dR,dD = [],[],[],[]
    alphas,betas,gammas,deltas = [],[],[],[]
    
    # create lists of all possible parameter values
    for i in range(len(S)-1):
        dE.append(E[i+1]-E[i])
        dI.append(I[i+1]-I[i])
        dR.append(R[i+1]-R[i])
        dD.append(D[i+1]-D[i])
        alphas.append((dI[-1]+dR[-1]+dD[-1])/E[i])
        betas.append((N*(dE[-1]+dI[-1]+dR[-1]+dD[-1]))/(S[i]*I[i]))
        gammas.append(dR[-1]/I[i])
        deltas.append(dD[-1]/I[i])

    # define initial parameters
    params = Parameters()
    params.add("alpha",np.average(alphas),min=min(alphas),max=max(alphas))
    params.add("beta",np.average(betas),min=min(betas),max=max(betas))
    params.add("gamma",np.average(gammas),min=min(gammas),max=max(gammas))
    if min(deltas) != max(deltas):
        params.add("delta",np.average(deltas),min=min(deltas),max=max(deltas))
    else:
        params.add("delta",np.average(deltas),vary=False)        
    return params

def error(params,inits,timespan,data):
    """
    Calculates the residual values of estimates
    """
    # calls the solver to grab solution
    soln = solver(timespan,inits,params)

    return (soln[:,1:]-data).ravel()

def solver(t,vals,params):
    """
    Defines the system of ODEs and uses passed values to update the model
    """
    # import dependency        
    from scipy.integrate import odeint

    # grab values from passed parameters
    n,e,i,r,d = vals

    # calculate s population
    s = n - (e+i+r+d)

    # define params
    alpha,beta = params["alpha"].value,params["beta"].value
    gamma,delta = params["gamma"].value,params["delta"].value

    return odeint(SEIRDModel,[s,e,i,r,d],t,args=(alpha,beta,gamma,delta))
    
def SEIRDModel(z,t,alpha,beta,gamma,delta):
    """
    Defines the system of ODEs
    """
    # grab values
    s,e,i,r,d = z
    N = s+e+i+r+d

    # define the differential system
    dS = -(beta*s*e)/N
    dE = (beta*s*e)/N - alpha*e
    dI = alpha*e - gamma*i - delta*i
    dR = gamma*i
    dD = delta*i

    return [dS,dE,dI,dR,dD]

def loadData(path,county,year=2021,months=[10],records=14):
    """
    Loads in source data from json file and returns distinct population lists
    """
    # import dependencies
    from time import strptime
    import json

    # initialize lists of distinct population groups
    tested,confirmed,recovered,deceased = [],[],[],[]

    # open the file, iterate over each nested dictionary, and append appropriate data
    with open(path) as f:
        data = json.load(f)

        for day,data in data[county]["dates"].items():
            day = strptime(day,"%Y-%m-%d")
            is_from_last_month = (day.tm_year == year and day.tm_mon in months) 
            if data.get("total") and is_from_last_month:
                tested.append(data["total"].get("tested") if data["total"].get("tested") else 0)
                confirmed.append(data["total"].get("confirmed") if data["total"].get("confirmed") else 0)
                recovered.append(data["total"].get("recovered") if data["total"].get("recovered") else 0)
                deceased.append(data["total"].get("deceased") if data["total"].get("deceased") else 0)
    return tested[-records:],confirmed[-records:],recovered[-records:],deceased[-records:]

def plotCompare(actual,predicted,t,title="Population over Time"):
    """
    Plots simulated values against source data values
    """
    # import dependencies
    import matplotlib
    matplotlib.use("TkAgg")
    import matplotlib.pyplot as plt

    # plot the comparison
    plt.figure(figsize=(16,9))
    plt.plot(t,actual,color="blue",lw=2,label="actual")
    plt.plot(t,predicted,color="orange",lw=2,label="predicted")
    plt.ylabel("Proportion of Population")
    plt.xlabel("Time in days")
    plt.title(title)
    plt.legend()
    plt.show()

#   #   #   #   #   #

if __name__ == "__main__":
    # import dependencies
    import os
    import numpy as np
    from lmfit import minimize

    # define filters
    YEAR = 2021
    MONTHS = [10]
    RECORDS = 14
    COUNTIES = ["AP","CH","DL","MH"]
    names = {"AP":"Andhra Pradesh","CH":"Chandigarh","DL":"Delhi","MH":"Maharashtra"}
    TRAINDAYS = 7

    # define population constant
    N = 1398705031

    # define path to data
    path = os.path.join("data","india-data.json")

    # Run simulation for each desired county
    for county in COUNTIES:
        # load data into lists
        E,I,R,D = loadData(path,county,year=YEAR,months=MONTHS,records=RECORDS)

        # get parameter estimation ranges based on tdays of data
        params = estimateParams(N,E[:TRAINDAYS],I[:TRAINDAYS],R[:TRAINDAYS],D[:TRAINDAYS])

        # convert input data to array to feed into optimizer for residuals calculations
        data = np.concatenate([np.array(i) for i in [E,I,R,D]]).reshape((4,14)).T

        # define the initial conditions and period of simulation
        x0 = [N,E[0],I[0],R[0],D[0]]
        period = list(range(len(E)))

        # Run the model through optimizer to get best parameters
        result = minimize(error,params,args=(x0,period,data),method="leastsq")

        # Run the model using the best parameters
        soln = solver(period,x0,result.params)
        S_p,E_p,I_p,R_p,D_p = soln[:,0],soln[:,1],soln[:,2],soln[:,3],soln[:,4]

        plotCompare(I,I_p,period,title=f"{names[county]} Infected Pop over Time")
        plotCompare(R,R_p,period,title=f"{names[county]} Recovered Pop over Time")
        plotCompare(D,D_p,period,title=f"{names[county]} Deceased Pop over Time")
