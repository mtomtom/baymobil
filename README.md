# bf_calculator

To install: pip install -i https://test.pypi.org/simple/ bf-calculator

To run: 

1. Run individual values. Code examples:

  `import bf_calculator as bf`
  
  `bf.run_bayes_analysis(Nhom1, nhom1, Nhom2, nhom2, Nhet, nhet, nmax)`

if you want the function to test all values of N:

  `bf.run_bayes_analysis(<Nhom1>, <nhom1>, <Nhom2>, <nhom2>, <Nhet>, <nhet>, "max")`
  
2. Pass in files with multiple values. Files need to be in csv format and include columns with the headers: SNP N eco1 eco2, where SNP is the SNP identifier, N is the total number of reads, eco1 are the reads that map to the local ecotype, eco2 are the reads that map to the distal ecotype. Any additional columns (for example, transcript ID) won't affect the code.

The code will then output a results file comprising the original data plus the maximum number of reads with the highest evidence for coming from the second ecotype (n2max), the average over all the reads that could from the second ecotype (N2mean) and the log 10 Bayes factor (log10BF).

To use this function: bf.run_bayes_analysis_files(<hom1_file>,<hom2_file>,<het_file>,<nmax>)

To run the examples in the test data set: `bf.run_bayes_analysis_files(""test_data/C-C-root-FN.csv", "test_data/P-P-root-FN.csv", test_data/Col-FN-root-1.csv", 10)`
