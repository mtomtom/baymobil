import numpy as np
import pandas as pd
import os
import warnings
from tqdm import tqdm
import multiprocessing as mp
import sys
import scipy.special
import numpy.ma as ma

## Define the main functions
def safebeta(N,n, alpha,beta):
    facterm = 1.0
    if alpha>1:
        a = np.arange(1, alpha)
        facterm = np.prod((n+a)/(N+a))
    if beta > 1.0:
        a = np.arange(1, beta)
        facterm = facterm * np.prod((N-n+a)/(N+alpha+a))
    return facterm / (N+alpha)

safebeta = np.vectorize(safebeta)

def logbeta(N,n, alpha,beta):
    facterm = 1.0
    if alpha>1:
        a = np.arange(1, alpha)
        facterm = (np.log(n+a) - np.log(N+a)).sum()
    if beta > 1.0:
        a = np.arange(1, beta)
        facterm = facterm + (np.log(N-n+a) - np.log(N+alpha+a)).sum()
    return facterm - np.log(N+alpha)

logbeta = np.vectorize(logbeta)

## Function to calculate the posterior ratio
def fasterpostN2(Nhomo1,nhomo1,Nhomo2,nhomo2,N,n,nmax):
    ## From the Team conversation 23.08.22
    """
    Idea for NANs: We can compute the error rate for the homograft data, e_hom=(n_hom+1)/(N_hom+2), and for the heterograft data, e_het=(n_het+1)/(N_het+2). We know if e_het < e_hom that we'll get negative Bayes factors and we know that if the errors were consistent this wouldn't get below -2 so we can set the NAN to -2 in this case. If e_het > e_hom, we'll get positive Bayes factors. NANs will occur here for very high numbers of N or when n is about half of N for moderately high N, but there is still information in there and would be good to capture. In this case we could set the NAN to, say, +10 (anything above +5 is so high that the actual number doesn't really matter too much)."""
    N = int(N)
    alpha1 = nhomo1+1
    beta1 = Nhomo1-nhomo1+1
    alpha2 = nhomo2+1
    beta2 = Nhomo2-nhomo2+1
    postN2 = np.zeros(N+2)
    postN2xN2 = np.zeros(N+2)
    PN2max = -10.0
    N2max = 0
    for N2 in np.arange(0,min(N+1,n+nmax+1)):
        N2 = int(N2)
        i = N2 + 1
        n2_min = max(0,N2-n)
        n2_max = min(N-n, N2)
        n2_array = np.arange(n2_min, n2_max+1)
        postN2[i] = np.sum(safebeta(N-N2,n-N2+n2_array,alpha1,beta1) * safebeta(N2,n2_array,alpha2,beta2))
        postN2[i]=postN2[i]/i
        if (N2>0) & (postN2[i]>PN2max):
            PN2max = postN2[i]
            N2max = N2
        
        postN2xN2[i]=postN2[i]*N2
    
    sumpostN2=sum(postN2)
    ## Need to handle the output if sumpostN2 has gone to zero
    if sumpostN2 == 0:
        print("Warning! Error rate in homograft higher than heterograft. Setting BF to -2.")
        results = ["nd","nd",-2]
    else:
        postN2=postN2/sumpostN2
        postN2xN2=postN2xN2/sumpostN2
        if postN2[1] == 0:
            print("Warning! Error rate in heterograft is significantly higher than in homograft. Setting BF to 10")
            logBF21N2 = 10
        else:
            logBF21N2 = np.log10(postN2[N2max+1]/postN2[1]) # +1 because of the index
        meanN2=sum(postN2xN2)
        results = [meanN2,N2max,logBF21N2]
    return results

def fasterpostN2new(Nhomo1, nhomo1, Nhomo2, nhomo2, N, n, nmax):
    """
    New version of the fastpostN2 function that removes the loops, and also uses log values to allow for higher read depths
    """

    ## Define functions to be used
    def create_n2_array(val):
        n2_min = max(0,val-n)
        n2_max = min(N-n, val)
        n2_array = np.arange(n2_min, n2_max+1)
        return(n2_array)

    def calculate_postN2(N2_val, n2_arrays_val):
        return scipy.special.logsumexp(logbeta(N-N2_val,n-N2_val+n2_arrays_val,alpha1,beta1) + logbeta(N2_val,n2_arrays_val,alpha2,beta2))
    
    ## Initialise the paramters
    alpha1 = nhomo1+1
    beta1 = Nhomo1-nhomo1+1
    alpha2 = nhomo2+1
    beta2 = Nhomo2-nhomo2+1
    N2_values = np.arange(0,min(N+1,n+nmax+1))
    N2_values = N2_values.astype(int)

    ## For each element in N2_values, create an array between n2_min and m2_max
    #n2_arrays = [create_n2_array(x) for x in N2_values] 
    n2_arrays = map(create_n2_array, N2_values)
    #postN2 = [calculate_postN2(x,y) for x,y in zip(N2_values,n2_arrays)]
    postN2 = map(calculate_postN2, N2_values, n2_arrays)

    ## postN2 is log(postN2)
    i = np.log(np.arange(0,min(N+1,n+nmax+1)) + 1)
    ## Subtract the log of the array
    postN2 = list(postN2) - i

    ## Multiply by log N2_values
    x = ma.log(N2_values)
    x = x.filled(0)
    postN2xN2 = postN2 + x

    ## Sum the values
    sumpostN2=scipy.special.logsumexp(postN2)

    ## Need to get the maximum of the posterior, ignoring N2 = 0 (as that's hypothesis 1 - no mobile reads)
    postN2_2 = postN2[N2_values > 0]
    N2_values_2 = N2_values[N2_values>0]
    N2max = N2_values_2[postN2_2==max(postN2_2)][0]

    postN2 = postN2 - sumpostN2
    postN2xN2 = postN2xN2 - sumpostN2
    logBF21N2 = (postN2[N2max] - postN2[0])/np.log(10)

    
    postN2xN2 = np.exp(postN2xN2)

    ## There may be a better way of doing this, but this will do for now. As the postN2xN2 array was supposed to have been multiplied by the N2_values, the first value should be 0. Subsequent multiplications and divisions wouldn't have changed this.  We therefore set these values to 0 here.
    postN2xN2[N2_values==0]=0
    meanN2 = sum(postN2xN2)
    ## Set N2max to N2_values where postN2 is maximised
    
    results = [meanN2,N2max,logBF21N2]
    return results

## Data can be passed in three formats: single value, dataframe, file
## Default case = single values, nmax=10
def run_bayes_analysis(data_list, nmax=10):
    """
    Data passed as a list that can contain ints, dataframes, or filenames.
    If dataframes or filenames, then headers need to be: SNP, N, n
    Files need to be in csv format.
    """
    if isinstance(data_list, list): 
        ## If 6 single values are passed, then return a single BF 
        if (len(data_list) == 6) & (all([isinstance(item, int) for item in data_list])):
            Nhomo1,nhomo1,Nhomo2,nhomo2,N,n = data_list
            [meanN2,N2max,log10BF] = fasterpostN2new(Nhomo1,nhomo1,Nhomo2,nhomo2,N,n,nmax)
            return meanN2, N2max, log10BF
        ## If 3 dataframes are passed, then process and return a single dataframe
        if (len(data_list) == 3) & (all([isinstance(item, pd.DataFrame) for item in data_list])):  
            ## Run analysis
            [df_hom_eco1, df_hom_eco2, df_het] = data_list
            results_df = run_bayes_analysis_df(df_hom_eco1, df_hom_eco2, df_het, nmax)
            return results_df
        ## If 3 strings are passed, then load in the files and process the data
        if (len(data_list) == 3) & (all([isinstance(item, str) for item in data_list])):
            ## Run analysis
            [df_hom_eco1, df_hom_eco2, het_file] = data_list
            run_bayes_analysis_files(df_hom_eco1, df_hom_eco2, het_file, nmax)
        ## If wrong length of list is passed, then print an error
        if (len(data_list)!= 6) & (len(data_list) != 3):
            print("Error: list should contain either 6 values, 3 dataframes or 3 filepaths.")
            return None
    else:
        print("Data need to be in list format")
        return False


## Function to check that the files are in the appropriate format
def check_data(data_file):
    """
    Column headings need to be "SNP","N","n
    """
    df = pd.read_csv(data_file, low_memory=False, index_col = None)
    ### Check column headings
    cols = df.columns.to_list()
    check_cols = ['SNP','N','n']
    if not set(check_cols).issubset(cols):
        print("File must include columns: SNP N n")
        raise ValueError('Files not formatted correctly.')
    ### Check values can be converted to numeric (if not already)
    cols = df.columns.drop(['SNP'])
    try:
        df[cols] = df[cols].apply(pd.to_numeric, errors='coerce')
        len_df = len(df) 
        df.dropna(inplace=True)
        if len(df) < len_df:
            print(f"Warning: {len_df - len(df)} Non numeric values were removed.")
        len_df = len(df)
        df = df[df["N"]>0]
        if len(df) < len_df:
            print(f"Warning: {len_df - len(df)} Zero values were removed.")
    except:
        raise ValueError('Unable to convert values to numeric.')

    ## Need to check whether any eco2 columns are greater than the N columns, as this will cause the code to crash
    error_values = len(df.loc[df["n"] > df["N"]])
    if error_values>0:
        raise ValueError('Error: n values cannot be greater than N')
    return df

def check_data_df(df):
    cols = df.columns.to_list()
    check_cols = ['SNP','N','n']
    if not set(check_cols).issubset(cols):
        print("File must include columns: SNP N n")
        raise ValueError('Files not formatted correctly.')
    ### Check values can be converted to numeric (if not already)
    cols = df.columns.drop(['SNP'])
    try:
        df[cols] = df[cols].apply(pd.to_numeric, errors='coerce')
        len_df = len(df) 
        df.dropna(inplace=True)
        if len(df) < len_df:
            print(f"Warning: {len_df - len(df)} Non numeric values were removed.")
        len_df = len(df)
        df = df.loc[df["N"]>0]
        if len(df) < len_df:
            print(f"Warning: {len_df - len(df)} Zero values were removed.")
    except:
        raise ValueError('Unable to convert values to numeric.')

    ## Need to check whether any eco2 columns are greater than the N columns, as this will cause the code to crash
    error_values = len(df.loc[df["n"] > df["N"]])
    if error_values>0:
        raise ValueError('Error: n values cannot be greater than N')
    return df

## Apply all the functions and calculate the BF
def calculate_evidence(df):
    df["pos"] = df.apply(lambda x:fasterpostN2new(x.Nhomo1,x.nhomo1,x.Nhomo2,x.nhomo2,x.N,x.n, x.nmax), axis=1 )
    df[['meanN2','N2max','log10BF']] = pd.DataFrame(df.pos.tolist(), index= df.index)
    df.drop(columns=["pos"], inplace=True)
    return df

def run_bayes_analysis_files(df_hom_eco1:str, df_hom_eco2:str, het_file:str, nmax):
    ## Check that all of the files are in the correct format
    df_het_file = check_data(het_file)
    df_hom_eco1 = check_data(df_hom_eco1)
    df_hom_eco2 = check_data(df_hom_eco2)
    ## Rename the columns in the homograft files
    df_hom_eco1 = df_hom_eco1.rename(columns={'N': 'Nhomo1', 'n': 'nhomo1'})
    df_hom_eco2 = df_hom_eco2.rename(columns={'N': 'Nhomo2', 'n': 'nhomo2'})
    ## Create a single dataframe with all of the values
    het_file_merged = pd.merge(df_het_file, df_hom_eco1[["SNP","Nhomo1","nhomo1"]], on = "SNP")
    het_file_merged = pd.merge(het_file_merged, df_hom_eco2[["SNP","Nhomo2","nhomo2"]], on = "SNP")
    ## Add in our nmax values
    if nmax == "max":
        het_file_merged["nmax"] = het_file_merged["N"]
    else:
        try:
            het_file_merged["nmax"] = nmax
            pd.to_numeric(het_file_merged["nmax"], errors='coerce')
        except:
            raise ValueError('Unable to convert nmax value to numeric.')
    ## Split the dataframe into chunks to parallelise
    dfs = np.array_split(het_file_merged, mp.cpu_count() - 1)
    pool = mp.Pool(processes = (mp.cpu_count() - 1))
    results = list(tqdm(pool.imap(calculate_evidence, dfs)))
    pool.close()
    pool.join()
    df = pd.concat(results)
    outfile = het_file.split(".")[0] + "_results.csv"
    df.to_csv(outfile, index = None)

def run_bayes_analysis_df(df_hom_eco1:pd.DataFrame, df_hom_eco2:pd.DataFrame, df_het:pd.DataFrame, nmax):
    ## Check that all of the dataframes are in the correct format
    df_het_file = check_data_df(df_het)
    df_hom_eco1 = check_data_df(df_hom_eco1)
    df_hom_eco2 = check_data_df(df_hom_eco2)
    ## Rename the columns in the homograft files
    df_hom_eco1 = df_hom_eco1.rename(columns={'N': 'Nhomo1', 'n': 'nhomo1'})
    df_hom_eco2 = df_hom_eco2.rename(columns={'N': 'Nhomo2', 'n': 'nhomo2'})
    ## Create a single dataframe with all of the values
    het_file_merged = pd.merge(df_het, df_hom_eco1[["SNP","Nhomo1","nhomo1"]], on = "SNP")
    het_file_merged = pd.merge(het_file_merged, df_hom_eco2[["SNP","Nhomo2","nhomo2"]], on = "SNP")
    ## Add in our nmax values
    if nmax == "max":
        het_file_merged["nmax"] = het_file_merged["N"]
    else:
        try:
            het_file_merged["nmax"] = nmax
            pd.to_numeric(het_file_merged["nmax"], errors='coerce')
        except:
            raise ValueError('Unable to convert nmax value to numeric.')
    ## Split the dataframe into chunks to parallelise
    dfs = np.array_split(het_file_merged, mp.cpu_count() - 1)
    pool = mp.Pool(processes = (mp.cpu_count() - 1))
    results = list(tqdm(pool.imap(calculate_evidence, dfs)))
    pool.close()
    pool.join()
    return pd.concat(results)