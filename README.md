# baymobil

To install: pip install -i https://test.pypi.org/simple/ baymobil

To run: 

`import baymobil as baymob`

The main function takes in a list of either single values, dataframes, or filenames, and depending on which objects are passed, handles the output differently.

Case 1: Passing list of single values:

`baymob.run_bayes_analysis(<values_list>, <nmax>)`

where values_list contains Nhom1, nhom1, Nhom2, nhom2, Nhet, nhet, which are all ints. If nmax is omitted, then the default value of 10 will be used. If it is set to "max", then the Nhet will be used: e.g. `bf.run_bayes_analysis(<values_list> "max")`

Case 2: Passing list of dfs

If the function is passed a list of dataframes (hom1, hom2, het), then it will return a new dataframe with new columns for maximum number of reads with the highest evidence for coming from the second ecotype (n2max), the average over all the reads that could from the second ecotype (N2mean) and the log 10 Bayes factor (log10BF). The dataframes need to include columns with the headers: SNP N eco1 eco2, where SNP is the SNP identifier, N is the total number of reads, eco1 are the reads that map to the local ecotype, eco2 are the reads that map to the distal ecotype. Any additional columns (for example, transcript ID) won't affect the code and will just be ignored, but retained for the results dataframes.

Case 3: Passing list of filenames (hom1, hom2, het)
  
Files need to be in csv format and include columns with the headers: SNP N eco1 eco2, which describe the data as above. The code will then output a results file comprising the original data plus the maximum number of reads with the highest evidence for coming from the second ecotype (n2max), the average over all the reads that could from the second ecotype (N2mean) and the log 10 Bayes factor (log10BF).

## Running the simulations

See the jupyter notebook for more how to run the simulations seen in the paper.

