README file for R package supporting the paper "A joint graphical model for inferring gene networks
across multiple subpopulations and data types".


Contents of this archive
------------------------
This archive contains 
(1) pkg: subdirectory that contains the R package.
(2) JEGN-manual.pdf: reference manual.
(3) simulation: subdirectory that contains codes for carrying out simulation studies. Run the "demo_simulation_JEGN.R" file to perform simulation studies.
(4) TCGA_BRCA: subdirectory that contains codes for applying JEGN to TCGA breast cancer datasets. Run the "TCGA_BRAC_Application.R" file to perform real data analysis.

The JEGN package has the following R-package dependencies: Matrix, stats.
The dependents are automatically installed along with JEGN. You can use the following commands to install JEGN from GitHub.


# Step 1. Install the devtools package. Invoke R and then type
install.packages("devtools")

# Step 2. Load the devtools package.
library("devtools")

# Step 3. Install the JEGN package from GitHub.
install_github("Zhangxf-ccnu/JEGN", subdir="pkg")


Useage
Load the library JEGN in R console, by running
library(JEGN)

Simply run the one of diffential analysis methods on your favorite datasets. For example,
data("TCGA.BRCA")
TCGA.BRCA.JEGN = JEGN(TCGA.BRCA$X, lambda =  0.95, alpha = 0.4, model = "nonparanormal")

For more examples about simulation studies, please refer to the "simulation" subdirectory. 
For more examples about real data application, please refer to the "TCGA_BRCA" subdirectory. 

Please do not hesitate to contact Dr. Xiao-Fei Zhang at zhangxf@mail.ccnu.edu.cn to 
seek any clarifications regarding any contents or operation of the archive.