def plot_compare(actual, predicted, t, title="Population over Time"):
    """
    Plot simulated values against source data values
    """
    # import dependencies
    import matplotlib
    matplotlib.use("TkAgg")
    import matplotlib.pyplot as plt

    # plotting
    plt.figure(figsize=(16, 9))
    # source data plot
    plt.plot(t, actual, color="blue", lw=2, label="actual")
    # simulated/predict plot
    plt.plot(t, predicted, color="orange", lw=2, label="predicted")
    # add labels, title, and legend
    plt.ylabel("Proportion of Population")
    plt.xlabel("Time in days")
    plt.title(title)
    plt.legend()
    plt.show()

if __name__ == "__main__":
    # import dependency
    from tdseird import SEIRD
    import os

    # define the list of files for which to run the simulation on
    states = ["andaman","andhra","delhi","manipur"]
    input_files = ["india-"+state+".csv" for state in states]

    # initialize population
    N = 1398705031
    # get data to run the model

    for state,filename in zip(states,input_files):
        path = os.path.join("data",filename)
        outpath = os.path.join("data","output")

        # create model with data
        model = SEIRD(path,N,outpath)

        # get initials
        start,T,I,R,D = model.feed_data()

        # get results of model
        df = model.final_run()

        # define range for which to display values
        t = list(range(len(I[-21:])))

        # plot results
        #plot_compare(I[-21:],list(df["Active"]),t,title=f"Infected Population over Simulation Period for {state}")
        #plot_compare(R[-21:],list(df["Recovered"]),t,title=f"Recovered Population over Simulation Period for {state}")
        #plot_compare(D[-21:],list(df["Deceased"]),t,title=f"Deceased Population over Simulation Period for {state}")

        plot_compare(I[-21:],R[-21:],D[-21:],t,title=f"Simulation over Period for {state}")
