import bf
import argparse
import os

parser = argparse.ArgumentParser()

## User specifies heterograft file and two homograft files
parser.add_argument("het_file", metavar="het_file", type = str, help="Heterograft data file")
parser.add_argument("hom_file_eco1", metavar="hom_file_eco1", type = str, help="Heterograft data file")
parser.add_argument("hom_file_eco2", metavar="hom_file_eco2", type = str, help="Heterograft data file")

args = parser.parse_args()

het_file = args.het_file
hom_file_eco1 = args.hom_file_eco1
hom_file_eco2 = args.hom_file_eco2

## Confirm the files back to the user
print(f"Running Bayesian analysis for: het file = {het_file}, hom file eco1 = {hom_file_eco1}, and hom file eco2 = {hom_file_eco2}")

df_hom_eco1, df_hom_eco2 = bf.get_homograft_data(hom_file_eco1, hom_file_eco2)
df_het = bf.run_bayes_analysis(het_file, df_hom_eco1, df_hom_eco2)
## Write output to file for further analysis
       
if not os.path.exists('output'):
    os.makedirs('output')
path = 'output/'
    
## Output the data
filename = het_file.split("/")[-1]
filename = filename.split(".")[0:-1]
filename = ' '.join(filename)
#df_het[["SNP","Bayes factor","n2"]].to_csv(path + filename + "_results.csv", index = None)
df_het.to_csv(path + filename + "_results.csv", index = None)




    