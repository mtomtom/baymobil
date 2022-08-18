## Functions to create simulated datasets
import random
import pandas as pd
import numpy as np
import itertools
import configparser
import json
from distutils.util import strtobool
import baymobil as baymob

def load_parameters():
    config = configparser.ConfigParser()
    config.read('parameters.cfg')
    N_values = json.loads(config.get("Simulation parameters","N_values"))
    q_values = json.loads(config.get("Simulation parameters","q_values"))
    N2_values = json.loads(config.get("Simulation parameters","N2_values")) 

    ## Check if the homografts should have separate read depths. Else make them the same as the heterografts
    constant_Nhom = bool(strtobool((config.get("Simulation parameters","constant_Nhom"))))
    constant_Nhom_value = int(config.get("Simulation parameters","constant_Nhom_value"))

    ## Set the random flags
    random_N = bool(strtobool((config.get("Simulation parameters","random_N"))))
    random_Nhom = bool(strtobool((config.get("Simulation parameters","random_Nhom"))))
    random_q = bool(strtobool((config.get("Simulation parameters","random_q"))))

    ## These should be integers
    no_transcripts = int(config.get("Simulation parameters","no_transcripts"))
    return N_values, q_values, N2_values, no_transcripts, constant_Nhom, constant_Nhom_value, random_N, random_Nhom, random_q

def setup_df_parameters(N_values, q_values, N2_values, no_transcripts, constant_Nhom, constant_Nhom_value,random_N, random_Nhom, random_q):
    ## Create iterations of all parameter values
    settings_list = [N_values, q_values, N2_values]
    data = (list(itertools.product(*settings_list)))
    df = pd.DataFrame(data, columns = ['N','q',"N2_func"])
    df = pd.concat([df]*no_transcripts)

    if constant_Nhom:
        df["Nhom1"] = constant_Nhom_value
        df["Nhom2"] = constant_Nhom_value
    else:
        df["Nhom1"] = df["N"]
        df["Nhom2"] = df["N"]

    ## Define our mobile transcripts - each transcript should have one unique definition which is kept the same for all parameter values
    mobile_def = random.choices([True, False], weights=[0.5, 0.5], k=no_transcripts)
    mobile_def = np.repeat(mobile_def, len(data))
    df["mobile"] = mobile_def

    ## Apply random flags (if appropriate)
    if random_N:
        df["N"] = df["N"].apply(lambda x: np.random.randint(low = 10, high = x))
    if random_Nhom:
        df["Nhom1"] = df["Nhom1"].apply(lambda x: np.random.randint(low = 10, high = x))
        df["Nhom2"] = df["Nhom2"].apply(lambda x: np.random.randint(low = 10, high = x))
    if random_q:
        df["q"] = df["q"].apply(lambda x: np.random.uniform(low = 0, high = x))
    return df

def create_errors(row, q):
        rand_numbs = np.random.rand(int(row))
        nhom1 = (rand_numbs < q).sum()
        return int(nhom1)

def generate_errors(df):
    ## Create the sequencing errors
    df["nhom1"] = df.apply(lambda x: create_errors(x.Nhom1, x.q),axis=1)
    df["nhom2"] = df.apply(lambda x: create_errors(x.Nhom2, x.q),axis = 1)
    df["n"] = df.apply(lambda x: create_errors(x.N, x.q), axis = 1)
    df["N2"] = df["N2_func"] * df["mobile"] * df["q"] * df["N"]

    df["n"] = df["n"] + df["N2"]
    return df

def run_bayes_analysis(df):
    ## Run Bayes analysis
    df["log10BF"] = df.apply(lambda x: baymob.fasterpostN2(x.Nhom1,x.nhom1,x.Nhom2,x.nhom2,x.N,x.n,10)[2], axis=1)

    ## Add in the columns for the results
    df["TP_bf"] = 0
    df["TN_bf"] = 0
    df["FP_bf"] = 0
    df["FN_bf"] = 0

    df.loc[(df.mobile == True) & (df.log10BF>=1),"TP_bf"] = 1
    df.loc[(df.mobile == False) & (df.log10BF<1),"TN_bf"] = 1
    df.loc[(df.mobile == False) & (df.log10BF>=1),"FP_bf"] = 1
    df.loc[(df.mobile == True) & (df.log10BF<1),"FN_bf"] = 1

    return df

## Run the analysis for Method A and Method B
def run_method_a(df):
    df["TP_Method_A"] = 0
    df["TN_Method_A"] = 0
    df["FP_Method_A"] = 0
    df["FN_Method_A"] = 0

    df["Method_A"] = df["n"]
    df.loc[df.Method_A>0, "Method_A"] = 1

    df.loc[(df.mobile == True) & (df.Method_A==1),"TP_Method_A"] = 1
    df.loc[(df.mobile == False) & (df.Method_A==0),"TN_Method_A"] = 1
    df.loc[(df.mobile == False) & (df.Method_A==1),"FP_Method_A"] = 1
    df.loc[(df.mobile == True) & (df.Method_A==0),"FN_Method_A"] = 1

    return df

def run_method_b(df):
    df["TP_Method_B"] = 0
    df["TN_Method_B"] = 0
    df["FP_Method_B"] = 0
    df["FN_Method_B"] = 0

    ## Universal pipeline
    df["Method_B"] = df["n"]
    df.loc[df.Method_B>0, "Method_B"] = 1
    df.loc[(df.nhom1>0)|(df.nhom2>0),"Method_B"] = 0

    df.loc[(df.mobile == True) & (df.Method_B==1),"TP_Method_B"] = 1
    df.loc[(df.mobile == False) & (df.Method_B==0),"TN_Method_B"] = 1
    df.loc[(df.mobile == False) & (df.Method_B==1),"FP_Method_B"] = 1
    df.loc[(df.mobile == True) & (df.Method_B==0),"FN_Method_B"] = 1

    return df

def create_simulated_data():
    N_values, q_values, N2_values, no_transcripts, constant_Nhom, constant_Nhom_value, random_N, random_Nhom, random_q = load_parameters()
    df = setup_df_parameters(N_values, q_values, N2_values, no_transcripts, constant_Nhom, constant_Nhom_value,random_N, random_Nhom, random_q)
    df = generate_errors(df)
    df = run_bayes_analysis(df)
    df = run_method_a(df)
    df = run_method_b(df)
    return df