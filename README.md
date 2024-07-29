# Baymobil SNP Analysis Tool

Baymobil is a tool designed to analyze SNP (Single Nucleotide Polymorphism) specific error rates in grafted plants. It compares the error rates between homograft and heterograft data sets to determine if the heterograft data contains mobile transcripts.

## Overview

Baymobil takes in two sets of homograft data and one set of heterograft data. It compares the SNP-specific error rates between the heterograft and homograft data sets. If the error rates differ significantly, it suggests that the heterograft data contains both errors and reads from mobile transcripts, indicating that they come from transcripts that have been transported across the graft junction.

## Dependencies

The following Python packages are required to run Baymobil:
- numpy
- scipy
- pandas
- multipledispatch
- pandarallel

## Usage

### Passing Single Values

```python
result = run_bayes(Nh1, nh1, Nh2, nh2, N, n, nmax)
print(result)
```

### Passing a Dataframe
```python
import pandas as pd

# Create or load a dataframe
df = pd.read_csv('snp_data.csv')

# Process the dataframe
result_df = run_bayes(df, nmax)
print(result_df)
```

### Passing Filenames
```python
# Set nmax (a good default is 10)
nmax = 10
# Process the files
run_bayes('homograft1.csv', 'homograft2.csv', 'heterograft.csv', nmax)
```

### Example DataFrame Format
Input data frames and CSV files should include the following columns:

- SNP: Identifier for the SNP
- N: Total number of reads
- n: Number of variant reads
- Nh1: Total reads in homograft 1
- nh1: Variant reads in homograft 1
- Nh2: Total reads in homograft 2
- nh2: Variant reads in homograft 2

| SNP |  N  | n  | Nh1 | nh1 | Nh2 | nh2 |
|-----|-----|----|-----|-----|-----|-----|
| rs1 | 100 | 10 | 120 | 12  | 130 | 15  |
| rs2 | 200 | 20 | 220 | 25  | 230 | 22  |

## Notes
- Ensure all input values are numeric and non-negative.
- The n values must not exceed the corresponding N values.
- The run_bayes function uses pandarallel for parallel processing. Ensure it is properly initialized in your environment.



