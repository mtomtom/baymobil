## Script to plot the accuracy of each method
import pandas as pd
import matplotlib.pyplot as plt

def plot_data(df, func_parameter):
    """ This function plots the accuracy of each method as a function of the passed parameter """
    func_parameter_list = set(df[func_parameter].to_list())
    for val in func_parameter_list:
        plot_df = df[df[func_parameter]==val]
        tp = plot_df["TP_bf"].sum()
        fp = plot_df["FP_bf"].sum()
        tn = plot_df["TN_bf"].sum()
        fn = plot_df["FN_bf"].sum()
        bf_accuracy = (tp + tn) / (tp + tn + fp + fn)
        plt.plot(val, bf_accuracy,"kx")

        tp = plot_df["TP_Method_A"].sum()
        fp = plot_df["FP_Method_A"].sum()
        tn = plot_df["TN_Method_A"].sum()
        fn = plot_df["FN_Method_A"].sum()
        abs_accuracy = (tp + tn) / (tp + tn + fp + fn)
        plt.plot(val, abs_accuracy,"rx")

        tp = plot_df["TP_Method_B"].sum()
        fp = plot_df["FP_Method_B"].sum()
        tn = plot_df["TN_Method_B"].sum()
        fn = plot_df["FN_Method_B"].sum()
        uni_accuracy = (tp + tn) / (tp + tn + fp + fn)
        plt.plot(val, uni_accuracy,"bx")

    plt.plot(val, bf_accuracy, "kx", label = "BF")
    plt.plot(val, abs_accuracy, "rx", label = "Method A")
    plt.plot(val, uni_accuracy, "bx",label = "Method B")

    plt.xlabel("N")
    plt.ylabel("Accuracy")
    plt.legend()
    #plt.savefig("N_comp_p01q.png",dpi=300)
    plt.show()