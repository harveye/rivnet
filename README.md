# Rivnet project 

## Code and R data repository for the Harvey and Altermatt XX article

### Repository structure:
- data: Original data
- scripts: R scripts (**see below**)
- R: Custom R function called by the main code from Scripts folder
- doc: Documents that are not generated by the code in this repository 
- figs: Figures generated by the code
- output: Documents generated by the code other than figures (tables, data files)

### Scripts order: 
The scripts have to be run in a specific order to run properly
- 1- Rivnet_DatMan.R (**this script cannot be run - 'Important notice' below**)
   - input: Raw data
   - output: Data for ordination
- 2- Rivnet_ORD.R
   - input: Data for ordination
   - output: Data for SEM
- 3- Rivnet_SEM.R
   - input: Data for SEM
   - output: SEM results

Rivnet_Piemap.R and Rivnet_TimeBxplot can be used in any order. The former will generate a raw version of Figure 1 in the article, and the second will generate a simple boxplot figure to investigate the effect of 'sampling year' on patterns. 

**Please note**
Unfortunately because of copyright issues we were not allowed to include the original data from the Swiss BioDiversity Monitoring program. This means that all data in the data folder are pre-processed information from the original data and that the Rivnet_DatMan.R script is here for information purpose only but will not run properly. Inquiries to obtain the original data can be made directily to the BDM organization: http://www.biodiversitymonitoring.ch/en/home.html
