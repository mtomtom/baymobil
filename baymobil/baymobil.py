### Baymobil takes in 2 sets of homograft data, and 1 set of heterograft data, and compares the SNP specific error rates between the heterograft and homograft data sets. If these differ significantly, then the assumption is that the heterograft data contains both errors and reads from mobile transcripts, suggesting that they come from transcripts that have been transported across the graft junction (mobile)

import numpy as np
import scipy.special
import numpy.ma as ma
import pandas as pd
from multipledispatch import dispatch
from pandarallel import pandarallel
import os

## Code for 3 different cases: single values, files, and dataframes
def logbeta(N, n, alpha, beta):
    facterm = 0.0
    if alpha > 1:
        a = np.arange(1, alpha)
        facterm += np.sum(np.log(n + a) - np.log(N + a))
    if beta > 1:
        a = np.arange(1, beta)
        facterm += np.sum(np.log(N - n + a) - np.log(N + alpha + a))
    return facterm - np.log(N + alpha)

logbeta = np.vectorize(logbeta)

def fasterpostN2(Nh1, nh1, Nh2, nh2, N, n, nmax):
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
    alpha1 = nh1+1
    beta1 = Nh1-nh1+1
    alpha2 = nh2+1
    beta2 = Nh2-nh2+1
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
    
    results = [meanN2,N2max,logBF21N2]
    return results

## Check the dataframes
def check_data(df):
    ### Check column headings
    cols = df.columns.to_list()
    check_cols = ['SNP','N','n','Nh1','nh1','Nh2','nh2']
    if not set(check_cols).issubset(cols):
        print("File must include columns: SNP N n Nh1 nh1 Nh2 nh2")
        raise ValueError('Files not formatted correctly.')
    ### Check values can be converted to numeric (if not already)
    cols = ["N","n","Nh1","nh1","Nh2","nh2"]
    try:
        df[cols] = df[cols].apply(pd.to_numeric, errors='coerce')
        len_df = len(df) 
        df.dropna(inplace=True)
        if len(df) < len_df:
            print(f"Warning: {len_df - len(df)} Non numeric values were removed.")
        len_df = len(df)
        df = df[df["N"]>0].copy()
        if len(df) < len_df:
            print(f"Warning: {len_df - len(df)} Zero values were removed.")
    except:
        raise ValueError('Unable to convert values to numeric.')

    ## Need to check whether any eco2 columns are greater than the N columns, as this will cause the code to crash
    error_values = len(df.loc[df["n"] > df["N"]])
    if error_values>0:
        raise ValueError('Error: n values cannot be greater than N')
    return df

## Instead of an if statement, create overloaded methods using multiple inputs:
## If inputs are all ints, then return a single value
## If inputs are all dataframes, then return a dataframe
## If inputs are files, then return a file

## Passing single values
@dispatch(int, int, int, int, int, int, int)
def run_bayes(Nh1, nh1, Nh2, nh2, N, n, nmax):
    result = fasterpostN2(Nh1, nh1, Nh2, nh2, N, n, nmax)
    return result

# passing a single dataframe. Needs to include columns: SNP N n Nh1 nh1 Nh2 nh2
@dispatch(pd.DataFrame, int)
def run_bayes(df, nmax):
    df = check_data(df)
    ## Add in our nmax values
    if nmax == "max":
        df["nmax"] = df["N"]
    else:
        try:
            df["nmax"] = nmax
            pd.to_numeric(df["nmax"], errors='coerce')
        except:
            raise ValueError('Unable to convert nmax value to numeric.')
    ## Split the dataframe into chunks to parallelise
    pandarallel.initialize(progress_bar=True)
    df["pos"] = df.parallel_apply(lambda x: fasterpostN2(x.Nh1, x.nh1, x.Nh2, x.nh2, x.N, x.n, nmax), axis = 1)
    df[['meanN2','N2max','log10BF']] = pd.DataFrame(df.pos.tolist(), index= df.index)
    df.drop(columns=["pos"], inplace=True)
    return df


# Passing three filenames
@dispatch(str, str, str, int)
def run_bayes(df_hom_eco1:str, df_hom_eco2:str, het_file:str, nmax):
    ## Check that all of the files are in the correct format
    df_hom1 = pd.read_csv(df_hom_eco1, low_memory=False, index_col = None)
    df_hom2 = pd.read_csv(df_hom_eco2, low_memory=False, index_col = None)
    df_het = pd.read_csv(het_file, low_memory=False, index_col = None)
    df_het = check_data(df_het)
    df_hom1 = check_data(df_hom1)
    df_hom2 = check_data(df_hom2)
    ## Rename the columns in the homograft files
    df_hom1 = df_hom1.rename(columns={'N': 'Nhomo1', 'n': 'nhomo1'})
    df_hom2 = df_hom2.rename(columns={'N': 'Nhomo2', 'n': 'nhomo2'})
    ## Create a single dataframe with all of the values
    het_merged = pd.merge(df_het, df_hom1[["SNP","Nhomo1","nhomo1"]], on = "SNP")
    het_merged = pd.merge(het_merged, df_hom2[["SNP","Nhomo2","nhomo2"]], on = "SNP")
    ## Add in our nmax values
    if nmax == "max":
        het_merged["nmax"] = het_merged["N"]
    else:
        try:
            het_merged["nmax"] = nmax
            pd.to_numeric(het_merged["nmax"], errors='coerce')
        except:
            raise ValueError('Unable to convert nmax value to numeric.')
    ## Split the dataframe into chunks to parallelise
    pandarallel.initialize(progress_bar=False)
    het_merged["pos"] = het_merged.parallel_apply(lambda x: fasterpostN2(x.Nhomo1, x.nhomo1, x.Nhomo2, x.nhomo2, x.N, x.n, 10), axis = 1)
    het_merged[['meanN2','N2max','log10BF']] = pd.DataFrame(het_merged.pos.tolist(), index= het_merged.index)
    het_merged.drop(columns=["pos"], inplace=True)
    ## Get the file name - remove path first
    file_name = os.path.basename(het_file)
    outfile = file_name.split(".")[0] + "_results.csv"
    het_merged.to_csv(outfile, index = None)