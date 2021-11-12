
class SEIRDModel():
    """
    Encapsulates the system of equations with helper functions and class variables for ease of
    time-dependent calculations
    """

    def __init__(self,population,susceptible,tested,confirmed,recovered,deceased):
        """
        Initializes the model with source data
        """
        self.N = population
        self.S = susceptible
        self.E = tested
        self.I = confirmed
        self.R = recovered
        self.D = deceased
        self.calc_deltas()

    def calc_deltas(self):
        """
        Calculates all daily differences in source data populations
        """
        self.S_diff = [0]+[self.S[i]-self.S[i-1] for i in range(1,len(self.S))]
        self.E_diff = [0]+[self.E[i]-self.E[i-1] for i in range(1,len(self.E))]
        self.I_diff = [0]+[self.I[i]-self.I[i-1] for i in range(1,len(self.I))]
        self.R_diff = [0]+[self.R[i]-self.R[i-1] for i in range(1,len(self.R))]
        self.D_diff = [0]+[self.D[i]-self.D[i-1] for i in range(1,len(self.D))]

    def differentials(self,vals,t):
        """
        Defines the system of ODEs and uses passed values to update the model
        """
        # grab values from passed parameters
        s,e,i,r,d = vals
        t = int(t) if t<len(self.E_diff) else len(self.E_diff)-1

        # grab the differences between previous paremeters and current ones
        de,di,dr,dd = self.E_diff[t],self.I_diff[t],self.R_diff[t],self.D_diff[t]

        # calculate the time dependent constants
        alpha = (di+dr+dd)/e
        beta = (de+di+dr+dd)*self.N/(s*i)
        gamma = dr/i
        delta = dd/i

        # define the differential system
        dS = -(beta*s*i)/self.N
        dE = (beta*s*i)/self.N - alpha*e
        dI = alpha*e - gamma*i - delta*i
        dR = gamma*i
        dD = delta*i

        return [dS,dE,dI,dR,dD]
        
    def run(self):
        """
        Leverages the scipy.integrate.odeint function to run stepwise updates
        """
        # import dependency        
        from scipy.integrate import odeint

        # define the time range, initial values, and run it
        t = list(range(len(self.S)))
        x0 = [self.S[0],self.E[0],self.I[0],self.R[0],self.D[0]]
        x = odeint(self.differentials,x0,t)
        return x

#   #   #   #   #   # 

def loadData(path,N,year=2021,months=[10]):
    """
    Loads in source data from json file and returns distinct population lists
    """
    # import dependency
    from time import strptime
    import json

    # initialize lists of distinct population groups
    susceptible,tested,confirmed,recovered,deceased = [],[],[],[],[]

    # open the file, iterate over each nested dictionary, and append appropriate data
    with open(path) as f:
        data = json.load(f)

        for day,data in data["MH"]["dates"].items():
            day = strptime(day,"%Y-%m-%d")
            is_from_last_month = (day.tm_year == year and day.tm_mon in months) 
            if data.get("total") and is_from_last_month:
                tested.append(data["total"].get("tested") if data["total"].get("tested") else 0)
                confirmed.append(data["total"].get("confirmed") if data["total"].get("confirmed") else 0)
                recovered.append(data["total"].get("recovered") if data["total"].get("recovered") else 0)
                deceased.append(data["total"].get("deceased") if data["total"].get("deceased") else 0)
                susceptible.append(N - (tested[-1]+confirmed[-1]+recovered[-1]+deceased[-1]))
    return susceptible,tested,confirmed,recovered,deceased

def plot_compare(actual,predicted,t,title="Population over Time"):
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
    # import dependency
    import os

    # define date filters
    YEAR = 2021
    MONTHS = [8,9,10]

    # get data to run the model
    N = 1398705031
    path = os.path.join("data","india-data.json")
    susceptible,tested,confirmed,recovered,deceased = loadData(path,N,year=YEAR,months=MONTHS)

    # create the model with given data and run the ODE system
    model = SEIRDModel(N,susceptible,tested,confirmed,recovered,deceased)
    X = model.run()

    # collect results
    s,e,i,r,d = X[:,0],X[:,1],X[:,2],X[:,3],X[:,4]
    t = list(range(len(s)))

    # plot results for each group of individuals in Maharashtra District
    plot_compare(tested,e,t,title="Exposed Population over Time for Maharashtra District")
    plot_compare(confirmed,i,t,title="Infected Population over Time for Maharashtra District")
    plot_compare(recovered,r,t,title="Recovered Population over Time for Maharashtra District")
    plot_compare(deceased,d,t,title="Deceased Population over Time for Maharashtra District")
