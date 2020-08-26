# septic_shock_metaanalysis
This repository contains the code for the project: 
## Transcriptomic meta-analysis reveals up-regulation of gene expression functional in osteoclast differentiation in human septic shock
 Samanwoy Mukhopadhyay,Pravat K. Thatoi,Abhay D. Pandey,Bidyut K. Das,Balachandran Ravindran,Samsiddhi Bhattacharjee,Saroj K. Mohapatra
## Availability of data
The dataset(s) and the R code supporting the conclusions of this article are available under the project ssnibmg in figshare https://doi.org/10.6084/m9.figshare.4592800.v1.

### Step 1: Installation of the Data Package ssnibmg
1. Dowload the file ssnibmg_1.0.tar.gz from
https://figshare.com/ (search for ssnibmg)
2. Change the directory to where you saved the file. Start R.
3. At the R prompt, issue the following command:
$ install.packages(pkgs=”ssnibmg_1.0.tar.gz”, repos=NULL)
4. Now the data package ssnibmg is installed on your computer.
5. Check with the following command:
$ library(“ssnibmg”)
### Step 2: Running the analysis code
6. Dowload the file ssnibmgdoc.zip from
https://figshare.com/ (search for ssnibmg)
7. Save the file at a suitable location and extract the contents. Change into the newly created directory ssnibmgdoc.
8. The code for data analysis are listed in the vignette ssnibmg_vignette.pdf.
9. The Rcode and Metadata are subdirectories under ssnibmgdoc (your current working directory in R).
10. The vignette is a self-computable document. Run the R code in steps to reproduce the results described in the vignette.
