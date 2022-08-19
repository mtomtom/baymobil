## Functions to create simulated datasets
import random
import pandas as pd
import numpy as np
import itertools
import configparser
import json
from distutils.util import strtobool
import os
import baymobil as baymob
import shutil

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

    ## Load in the number of replicates, and SNPs per transcript
    no_reps = int(config.get("Simulation parameters","no_reps"))
    no_SNPs = int(config.get("Simulation parameters","no_SNPs"))
    min_read_thresh = int(config.get("Simulation parameters","min_read_thresh"))

    params = [N_values, q_values, N2_values, no_transcripts, constant_Nhom, constant_Nhom_value, random_N, random_Nhom, random_q, no_reps, no_SNPs]
    param_names = ["N_values", "q_values", "N2_values", "no_transcripts", "constant_Nhom", "constant_Nhom_value", "random_N", "random_Nhom", "random_q", "no_reps", "no_SNPs"]
    return params, param_names

def create_errors(row, q):
        rand_numbs = np.random.rand(int(row))
        nhom1 = (rand_numbs < q).sum()
        return int(nhom1)

def create_homograft_data(N_values, q_values, N2_values, no_transcripts, constant_Nhom, constant_Nhom_value,random_N, random_Nhom, random_q, no_reps, no_SNPs):
    """The homograft dataframe is only created once and used to evaluate all heterograft replicates
    """
    ## Create iterations of all parameter values
    settings_list = [N_values, q_values, N2_values]
    data = (list(itertools.product(*settings_list)))
    df = pd.DataFrame(data, columns = ['N','q',"N2_func"])
    df = pd.concat([df]*no_transcripts * no_SNPs)

    if constant_Nhom:
        df["Nhom1"] = constant_Nhom_value
        df["Nhom2"] = constant_Nhom_value
    else:
        df["Nhom1"] = df["N"]
        df["Nhom2"] = df["N"]

    ## Apply the random flags (if applicable)
    if random_Nhom:
        df["Nhom1"] = df["Nhom1"].apply(lambda x: np.random.randint(low = 10, high = x))
        df["Nhom2"] = df["Nhom2"].apply(lambda x: np.random.randint(low = 10, high = x))
    if random_q:
        df["q"] = df["q"].apply(lambda x: np.random.uniform(low = 0, high = x))
    
    ## Add in the sequencing errors
    df["nhom1"] = df.apply(lambda x: create_errors(x.Nhom1, x.q),axis=1)
    df["nhom2"] = df.apply(lambda x: create_errors(x.Nhom2, x.q),axis = 1)

    ## Create and add in the SNP IDs
    snp_ids = []
    ## First value is the transcript number
    for i in range(no_transcripts):
        ## Second value is the SNP no
        for j in range(no_SNPs):
            ## Third value is the condition number
            for k in range(len(data)):
                snp_ids.append(str(i) + "_" + str(j) + "_" + str(k))

    df["SNP"] = snp_ids

    return df[["Nhom1","Nhom2","nhom1","nhom2"]]

def create_heterograft_data(N_values, q_values, N2_values, no_transcripts, constant_Nhom, constant_Nhom_value,random_N, random_Nhom, random_q, no_reps, no_SNPs):
    settings_list = [N_values, q_values, N2_values]
    data = (list(itertools.product(*settings_list)))
    df = pd.DataFrame(data, columns = ['N','q',"N2_func"])
    df = pd.concat([df]*no_transcripts * no_SNPs)

    ## Define our mobile transcripts - each transcript should have one unique definition which is kept the same for all parameter values
    mobile_def = random.choices([True, False], weights=[0.5, 0.5], k=no_transcripts)
    mobile_def = np.repeat(mobile_def, len(data * no_SNPs))
    df["mobile"] = mobile_def

    ## Apply random flags (if appropriate)
    if random_N:
        df["N"] = df["N"].apply(lambda x: np.random.randint(low = 10, high = x))
    if random_q:
        df["q"] = df["q"].apply(lambda x: np.random.uniform(low = 0, high = x))

    ## Add in errors and mobile reads
    df["n"] = df.apply(lambda x: create_errors(x.N, x.q), axis = 1)
    df["N2"] = df["N2_func"] * df["mobile"] * df["q"] * df["N"]

    df["n"] = df["n"] + df["N2"]
    ## Create and add in the SNP IDs
    snp_ids = []
     ## First value is the transcript number
    for i in range(no_transcripts):
        ## Second value is the SNP no
        for j in range(no_SNPs):
            ## Third value is the condition number
            for k in range(len(data)):
                snp_ids.append(str(i) + "_" + str(j)+ "_" + str(k))

    df["SNP"] = snp_ids
    return df

def run_analysis(dfhom, dfhet):
    ## Merge the dataframes
    df = pd.concat([dfhet, dfhom], axis=1)
    ## Run Bayes analysis
    df["log10BF"] = df.apply(lambda x: baymob.fasterpostN2(x.Nhom1,x.nhom1,x.Nhom2,x.nhom2,x.N,x.n,10)[2], axis=1)

    ## Sum up the Bayes Factors across SNPs: IDs are transcript_SNP_condition. So we sum all SNPs for each transcript for each condition

    df["condition"] = df["SNP"].apply(lambda x: x.split("_")[2])
    df["transcript"] = df["SNP"].apply(lambda x: x.split("_")[0])

    ## Method A
    ## Evaluate each outcome
    df["Method_A"] = df["n"]
    df.loc[df.Method_A>0, "Method_A"] = 1

    ## Method B
    ## Universal pipeline
    df["Method_B"] = df["n"]
    df.loc[df.Method_B>0, "Method_B"] = 1
    df.loc[(df.nhom1>0)|(df.nhom2>0),"Method_B"] = 0

    df_grouped = df.groupby(["condition","transcript"]).sum().reset_index()
    
    ## As we have summed across the mobile values, these have now become counts, so we need to evaluate 0 and above 0 as False and True respectively

    ## Define how many SNPs need to be positive (default = 1). This can be changed for more conservative analysis
    snp_thresh = 1
    df_grouped.loc[df_grouped.Method_A>=snp_thresh,"Method_A"] = 1
    df_grouped.loc[df_grouped.Method_B>= snp_thresh,"Method_B"] = 1

    return df_grouped

def create_simulated_data():
    params, param_names = load_parameters()
    dfhom = create_homograft_data(*params)
    ## Create replicates for heterograft data
    no_reps = params[param_names.index('no_reps')]
    ## Create output folder
    if os.path.exists('output'):
        shutil.rmtree('output')

    if not os.path.exists('output'):
        os.makedirs('output')

    for i in range(no_reps):
        dfhet = create_heterograft_data(*params)
        df = run_analysis(dfhom, dfhet)
        df.to_csv("output/output" + str(i) + ".csv")