# Splice Junction Filter Workflow



# Setup

## R Requirments

* `R` - https://cran.r-project.org/src/base/R-4/
* `R Studio` - https://www.rstudio.com/products/rstudio/download/

## R Dependencies

|  |  |  |  |
| --- | --- | --- | --- |
| recount3_1.10.2 | data.table_1.14.8 | slam_0.1-50 | stringr_1.5.0 |
| readr_2.1.4 | dplyr_1.1.2 | tibble_3.2.1 | tidyr_1.3.0 |

* Please note some R packages may have their own dependencies. Please Install prior to running the workflow.

## Select your study(s) of interests

Go to [Recount3's Study Explorer](https://jhubiostatistics.shinyapps.io/recount3-study-explorer/) to browse the data and choose a study or list of studies of interest. The workflow will require the proper project name(s), as shown in the Study Explorer.

![alt text](https://github.com/shawlab-moffitt/Splice-Junction-Filter-Workflow/blob/main/Example_Data/Recount3_Demo_Pic.png?raw=true)

# R/SpliceJxn_Workflow.R

Here is a full workflow script where users can provide a project names from Recount3 and obtain normalized and filtered Splice Junction counts and meta data.

## Input

* **Project_Directory:** The path to where the results will be written out
* **Project_Name_List:** A comma seperated list of Recount3 project names. This can be one or more.
* **Expr_CutPoint & SampleN_CutPoint:** A numeric filtering criteria for the junction data (default Expr_CutPoint=1 and SampleN_CutPoint=10).
  * When applied this will filter the junction matrix and only keep junctions with an expression level above [Expr_CutPoint (1)] in [SampleN_CutPoint (10)] or more samples.
 
## Script Workflow

1. Set up environment to download Recount3 data, load in libraries and functions, and check user inputs
2. Download project data from Recount3
3. Normalize junction data
4. Filter out junctions that are found in GTEx and Bone Marrow samples
5. Filter out junctions that have an expresion of zero accross all samples
6. Filter out junctions according to the criteria given as user input
7. Generate BED files for each junction matrix output

## Script Output

* **Junction Filter Script Log:** A log file tracking the number of junctions at each filtering step and the names of the junction output files
* **Meta Data:** Meta data extracted from the Recount3 project data
* **Junction Expression Matrices:** The junction counts extracted after each filtering step
* **Junction Name Annotation:** Junction name annotation data for each junction matrices output
* **BED Files:** Bed files denoting the average junction expression and total samples expressed for each junction matrix output

# R/SpliceJxn_Subset.R

Some project data, for example TCGA, has a mix of sample types (e.g. Tumor, Normal, Metastatic). This script allows the user to subset the junction data output from [SpliceJxn_Workflow](https://github.com/shawlab-moffitt/Splice-Junction-Filter-Workflow/tree/main#rsplicejxn_workflowr) by a column of choice from the meta data. 

## Input

* **Junction_File:** A junction matrix output from [SpliceJxn_Workflow](https://github.com/shawlab-moffitt/Splice-Junction-Filter-Workflow/tree/main#rsplicejxn_workflowr)
* **Meta_File:** Meta data from the same project as the Junction_File
* **Column_Name:** A desired column name from the meta data to subset the junction data based on (e.g. "tcga.cgc_sample_sample_type")

## Script Workflow

1. Load library and read in data
2. Subset a junction and meta data table based on each unique category from the desired column input
3. Write out subset data

## Output

* **Subset Junction & Meta data file:** For each unique category and junction and meta data file will be written to the same folder the input data was from.

# R/SpliceJxn2Bed.R

This script converts a splice junction expression file to a BED file. This is already performed in the SpliceJxn_Workflow.R script.

## Input

* **JxnFile:** A junction matrix output from [SpliceJxn_Workflow](https://github.com/shawlab-moffitt/Splice-Junction-Filter-Workflow/tree/main#rsplicejxn_workflowr)
* **AnnoFile:** Junction name annotation data from the same filtering step and project as the Junction_File
* **OutFile:** A desired output file name

## Script Workflow

1. Load libraries and read in data
2. Split the junction name into separate columns (Chromosome, Start, Stop, Strand)
3. Merge annotation information from the junction name annotation file (this denotes if the junction is known or unknown)
4. Calculate the average junction expression for each junction and the sum of the total samples expressing each junction
5. Write out BED files

## Output

* **Average Expression BED File:** Junction position information, annotation information, and average expression for each junction
* **Sample Count Expressing BED File:** Junction position information, annotation information, and total samples expressing each junction



