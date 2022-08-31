## Script to plot the accuracy of each method
import pandas as pd
import matplotlib.pyplot as plt
from os import listdir
from os.path import isfile, join
import configparser
from functools import reduce

def load_data():
    """
    All simulation files should be in the "output" folder. This function loads in all files in that folder.
    """
    folder = "output"
    results_files = [f for f in listdir(folder) if isfile(join(folder, f))]
    ## Ensure we only have relevant files
    results_files = [f for f in results_files if "output" in f]
    ## Combine the replicates
    no_reps = len(results_files)

    ## Quick check: is this value the same as in the parameter file?
    config = configparser.ConfigParser()
    config.read('parameters.cfg')
    no_reps_params = int(config.get("Simulation parameters","no_reps"))
    if no_reps != no_reps_params:
        print("Error! Number of files does not equal the number of replicates")
        SystemExit()
    
    ## Load in the files
    df_list = []
    for file in results_files:
        df = pd.read_csv("output/" + file)
        df_list.append(df)

    ## Sum all values
    df_summed = reduce(lambda x, y: x.add(y, fill_value=0), df_list)

    ## Create the TP, TN, FP, FN
     ## Add in the columns for the results
    df_summed["TP_bf"] = 0
    df_summed["TN_bf"] = 0
    df_summed["FP_bf"] = 0
    df_summed["FN_bf"] = 0

    df_summed["TP_Method_A"] = 0
    df_summed["TN_Method_A"] = 0
    df_summed["FP_Method_A"] = 0
    df_summed["FN_Method_A"] = 0

    df_summed["TP_Method_B"] = 0
    df_summed["TN_Method_B"] = 0
    df_summed["FP_Method_B"] = 0
    df_summed["FN_Method_B"] = 0

    ## Bayes factor: summed across replicates

    df_summed.loc[(df_summed.mobile > 0) & (df_summed.log10BF>=1),"TP_bf"] = 1
    df_summed.loc[(df_summed.mobile == 0) & (df_summed.log10BF<1),"TN_bf"] = 1
    df_summed.loc[(df_summed.mobile == 0) & (df_summed.log10BF>=1),"FP_bf"] = 1
    df_summed.loc[(df_summed.mobile > 0) & (df_summed.log10BF<1),"FN_bf"] = 1

    ## Define how many replicates need to be positive (default = all)

    df_summed.loc[(df_summed.mobile > 0) & (df_summed.Method_A==no_reps),"TP_Method_A"] = 1
    df_summed.loc[(df_summed.mobile == 0) & (df_summed.Method_A<no_reps),"TN_Method_A"] = 1
    df_summed.loc[(df_summed.mobile == 0) & (df_summed.Method_A==no_reps),"FP_Method_A"] = 1
    df_summed.loc[(df_summed.mobile > 0) & (df_summed.Method_A<no_reps),"FN_Method_A"] = 1

    df_summed.loc[(df_summed.mobile > 0) & (df_summed.Method_B == no_reps),"TP_Method_B"] = 1
    df_summed.loc[(df_summed.mobile == 0) & (df_summed.Method_B<no_reps),"TN_Method_B"] = 1
    df_summed.loc[(df_summed.mobile == 0) & (df_summed.Method_B == no_reps),"FP_Method_B"] = 1
    df_summed.loc[(df_summed.mobile > 0) & (df_summed.Method_B<no_reps),"FN_Method_B"] = 1

    ## Check that all values have been assigned
    if sum(df_summed[["TP_bf","TN_bf","FP_bf","FN_bf"]].sum()) != len(df_summed):
        print("Error! BF values not assigned")
    if sum(df_summed[["TP_Method_A","TN_Method_A","FP_Method_A","FN_Method_A"]].sum()) != len(df_summed):
        print("Error! Method A values not assigned")
    if sum(df_summed[["TP_Method_B","TN_Method_B","FP_Method_B","FN_Method_B"]].sum()) != len(df_summed):
        print("Error! Method B values not assigned!")

    return df_summed

def plot_data_all(df, func_parameter):
    """ This function plots the accuracy of each method as a function of the passed parameter """
    func_parameter_list = set(df[func_parameter].to_list())
    val_list = []
    bf_list = []
    methoda_list = []
    methodb_list = []

    for val in func_parameter_list:
        val_list.append(val)
        plot_df = df[df[func_parameter]==val]
        tp = plot_df["TP_bf"].sum()
        fp = plot_df["FP_bf"].sum()
        tn = plot_df["TN_bf"].sum()
        fn = plot_df["FN_bf"].sum()
        bf_accuracy = (tp + tn) / (tp + tn + fp + fn)
        bf_list.append(bf_accuracy)
        plt.plot(val, bf_accuracy,"kx")


        tp = plot_df["TP_Method_A"].sum()
        fp = plot_df["FP_Method_A"].sum()
        tn = plot_df["TN_Method_A"].sum()
        fn = plot_df["FN_Method_A"].sum()
        abs_accuracy = (tp + tn) / (tp + tn + fp + fn)
        methoda_list.append(abs_accuracy)
        plt.plot(val, abs_accuracy,"rx")

        tp = plot_df["TP_Method_B"].sum()
        fp = plot_df["FP_Method_B"].sum()
        tn = plot_df["TN_Method_B"].sum()
        fn = plot_df["FN_Method_B"].sum()
        uni_accuracy = (tp + tn) / (tp + tn + fp + fn)
        methodb_list.append(uni_accuracy)
        plt.plot(val, uni_accuracy,"bx")

    plt.plot(val, bf_accuracy, "kx", label = "BF")
    plt.plot(val, abs_accuracy, "rx", label = "Method A")
    plt.plot(val, uni_accuracy, "bx",label = "Method B")

    df = pd.DataFrame([val_list, bf_list, methoda_list, methodb_list])
    df = df.T
    df.columns = [func_parameter, "bf","methoda","methodb"]
    df.to_csv("output/plot_results" + func_parameter + "_.csv", index=None)

    plt.xlabel(func_parameter)
    plt.ylabel("Accuracy")
    plt.legend()
    plt.savefig("output/results.png",dpi=300)
    plt.show()

def plot_data_bf(df, func_parameter):
    """ This function plots the accuracy of Bayes factors as a function of the passed parameter """
    func_parameter_list = set(df[func_parameter].to_list())
    val_list = []
    bf_list = []
    methoda_list = []
    methodb_list = []

    for val in func_parameter_list:
        val_list.append(val)
        plot_df = df[df[func_parameter]==val]
        tp = plot_df["TP_bf"].sum()
        fp = plot_df["FP_bf"].sum()
        tn = plot_df["TN_bf"].sum()
        fn = plot_df["FN_bf"].sum()
        bf_accuracy = (tp + tn) / (tp + tn + fp + fn)
        bf_list.append(bf_accuracy)
        plt.plot(val, bf_accuracy,"kx")

    plt.plot(val, bf_accuracy, "kx", label = "BF")

    df = pd.DataFrame([val_list, bf_list])
    df = df.T
    df.columns = [func_parameter, "bf"]
    df.to_csv("output/plot_results" + func_parameter + "_.csv", index=None)

    plt.xlabel(func_parameter)
    plt.ylabel("Accuracy")
    plt.savefig("output/results.png",dpi=300)
    plt.show()