import numpy as np
import scipy.special as ss
import pandas as pd
from decimal import Decimal
import os
import warnings
from tqdm import tqdm

## The stirling approximation is used for large read depths, which would return a binomial
## value of inf
def stirling_binom(N,n2):
  if (N == 0) | (n2 == 0): return np.log(1)
  elif N == n2: return np.log(1)
  else: return N*np.log(N) - n2*np.log(n2) - (N-n2) * np.log(N-n2) + 0.5 *(np.log(N) - np.log(n2) - np.log(N-n2) -np.log(2 * np.pi))

## Function to check that the files are in the appropriate format
def check_data(data_file):
    df = pd.read_csv(data_file, low_memory=False, index_col = None)
    ### Check column headings
    cols = df.columns.to_list()
    check_cols = ['SNP','N','eco1','eco2']
    if not set(check_cols).issubset(cols):
        print("File must include columns: SNP N eco1 eco2")
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
    return df

## Apply all the functions and calculate the BF
def calculate_evidence(df):
    tqdm.pandas()
    ## Functions output natural logs - convert to log10
    print("Calculating Bayes Factors...")
    df["pos"] = df.progress_apply(lambda x:fasterpostN2(x.Nhomo1,x.nhomo1,x.Nhomo2,x.nhomo2,x.N,x.eco2, x.nmax), axis=1 )
    df[['meanN2','N2max','log10BF']] = pd.DataFrame(df.pos.tolist(), index= df.index)
    df.drop(columns=["pos"], inplace=True)
    return df

def calculate_evidence_stirling(df):
    tqdm.pandas()
    ## Functions output natural logs - convert to log10
    print("Calculating Bayes Factors using Stirling approximation...")
    df.drop(columns=["meanN2","N2max","log10BF"],inplace=True)
    df["pos"] = df.progress_apply(lambda x:fasterpostN2_stirling(x.Nhomo1,x.nhomo1,x.Nhomo2,x.nhomo2,x.N,x.eco2, x.nmax), axis=1 )
    df[['meanN2','N2max','log10BF']] = pd.DataFrame(df.pos.tolist(), index= df.index)
    df.drop(columns=["pos"], inplace=True)
    return df

### This function will take in 3 float values and return the output 
def run_bayes_analysis(hom_eco1N:float, hom_eco1n:float, hom_eco2N:float, hom_eco2n:float, hetN: float, hetn: float, nmax):
    if nmax == "max": nmax = hetN
    result = fasterpostN2(hom_eco1N, hom_eco1n, hom_eco2N, hom_eco2n, hetN, hetn, nmax)
    ## Result is meanN2, N2max, log10BF
    return result

### This function will take in 3 float values and return the output, using the Stirling function
def run_bayes_analysis_stirling(hom_eco1N:float, hom_eco1n:float, hom_eco2N:float, hom_eco2n:float, hetN: float, hetn: float, nmax):
    if nmax == "max": nmax = hetN
    result = fasterpostN2_stirling(hom_eco1N, hom_eco1n, hom_eco2N, hom_eco2n, hetN, hetn, nmax)
    ## Result is meanN2, N2max, log10BF
    return result

def run_bayes_analysis_files(df_hom_eco1:str, df_hom_eco2:str, het_file:str, nmax):
    ## Check that all of the files are in the correct format
    df_het_file = check_data(het_file)
    df_hom_eco1 = check_data(df_hom_eco1)
    df_hom_eco2 = check_data(df_hom_eco2)
    ## Rename the columns in the homograft files
    df_hom_eco1 = df_hom_eco1.rename(columns={'N': 'Nhomo1', 'eco2': 'nhomo1'})
    df_hom_eco2 = df_hom_eco2.rename(columns={'N': 'Nhomo2', 'eco2': 'nhomo2'})
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
    ## If the binomial returns an inf value, this code will raise a warning. We ignore it here, because we will replace these values later
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        df = calculate_evidence(het_file_merged)
    ## Check for inf values
    df_inf = df[df["log10BF"]==np.inf]
    if len(df_inf)>0:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            df_inf_results = calculate_evidence_stirling(df_inf)
        df_no_inf = df[df["log10BF"]!=np.inf]
        df_new = pd.concat([df_no_inf, df_inf_results],axis=1)
        df = df_new
    outfile = het_file.split(".")[0] + "_results.csv"
    df.to_csv(outfile, index = None)
    
## The main function that handles the Bayesian anaylsis. This takes the parameters: Nhomo1, ## nhomo2, the reads that map to eco1 and eco2 in homograft 1, Nhomo2, nhomo2, the reads ## that map to eco2 and eco1 in homograft 2, N, n, the reads that map to eco1 and eco2 in the heterograft dataset and nmax, the reads from N over which we want to integrate.

def fasterpostN2 (Nhomo1,nhomo1,Nhomo2,nhomo2,N,n,nmax):
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
        for n2 in np.arange(n2_min, n2_max+1):
            postN2[i] = postN2[i] + ss.binom(N-N2,n-N2+n2)*ss.binom(N2,n2)* ss.beta(n-N2+n2+alpha1,N-n-n2+beta1)*ss.beta(n2+alpha2,N2-n2+beta2)
        postN2[i]=postN2[i]/i
        if (N2>0) & (postN2[i]>PN2max):
            PN2max = postN2[i]
            N2max = N2
        
        postN2xN2[i]=postN2[i]*N2
    
    sumpostN2=sum(postN2)
    postN2=postN2/sumpostN2
    postN2xN2=postN2xN2/sumpostN2
    logBF21N2 = np.log10(postN2[N2max+1]/postN2[1]) # +1 because of the index
    meanN2=sum(postN2xN2)
    results = [meanN2,N2max,logBF21N2]
    return results

### The main function using the Stirling approximation for those values with large read depths
def fasterpostN2_stirling(Nhomo1,nhomo1,Nhomo2,nhomo2,N,n,nmax):
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
        for n2 in np.arange(n2_min, n2_max+1):
            log_values = Decimal(stirling_binom(N-N2,n-N2+n2) + stirling_binom(N2,n2)+ ss.betaln(n-N2+n2+alpha1,N-n-n2+beta1)+ss.betaln(n2+alpha2,N2-n2+beta2))
            postN2[i] = Decimal(postN2[N2+1]) + log_values.exp()
        postN2[i]=postN2[i]/i
        if (N2>0) & (postN2[i]>PN2max):
            PN2max = postN2[i]
            N2max = N2
        
        postN2xN2[i]=postN2[i]*N2
    
    sumpostN2=sum(postN2)
    postN2=postN2/sumpostN2
    postN2xN2=postN2xN2/sumpostN2
    logBF21N2 = np.log10(postN2[N2max+1]/postN2[1]) # +1 because of the index
    meanN2=sum(postN2xN2)
    results = [meanN2,N2max,logBF21N2]
    return results