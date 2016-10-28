Data processing workflow
=======
### Introduction
1. Alison send xls file (Release Manifest)
2. Copy all the GTCs from illumina-lims to ERISONE (server mounted as /external/ppm-illumina-lims)
3. Run gtc2ped.py
   a. Concordance table, tab-delimited (BiobankID, SentrixID)
   b. Gender table, tab-delimited SentrixID, BiobankID, Gender(1/2)
1. Follow diagram below:
![Diagram of data processing workflow](img/data_processing_workflow.png)

### PART 1: SNP remapping workflow
![PED file workflow](img/ped_file_generation_workflow.png)


### PART 2: SNP mapping and snpinfo file creation
![Snpinfo file workflow](img/snpinfo_processing_workflow.png)

