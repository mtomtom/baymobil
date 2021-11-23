# bayes-snp-diff

To run: `python run_bf.py <het_file> <hom_file_eco1> <hom_file_eco1>`

Files need to be in csv format with the headers: SNP N eco1 eco2, where SNP is the SNP identifier, N is the total number of reads, eco1 are the reads that map to the local ecotype, eco2 are the reads that map to the distal ecotype.

The code will then output a results file containing the SNP, the log 10 Bayes factor, and the number of reads with the highest evidence for coming from the second ecotype.

To run the examples in the test data set: `python run_bf.py test_data/Col-FN-shoot-1.csv test_data/C-C-shoot-FN.csv test_data/P-P-shoot-FN.csv`
