Biobank
=======
Bioinformatics processing pipelines for the Partners Biobank Initiative. These scripts are designed for the Illumina Multi-Ethnic SNP array, but they could be easily adapated to accomodate all Illumina-based SNP arrays.

etsai [at] bwh [dot] harvard [dot] edu

### Files Reference
```
These python classes are used in the data processing & QC scripts below.
Python classes:
   BPM.py
   GTC.py
   SnpInfo.py
   Strand.py

Data processing and QC: 
Mapping SNP Markers
   strand_check.py
   makeMAPfile.py
   MakeQCPlots.py

Data conversion (GTC -> PED, VCF):
   gtc2ped.py
   ped2vcf.py
   collapseVCF.py

Data delivery:
   GetImputedGenotypes.py			# need to update this for combining across data sets!
   GetIndividualGenotypesPED.py
   GetIndividualGenotypesVCF.py

Readme:
   readme/???
   readme/genotype_ped.r2.readme
   readme/genotype_vcf.r3.readme

Config:
   etc/???
   etc/genotype_ped.r2.yml
   etc/genotype_vcf.r2.yml

These scripts were written for some quick one-time analysis. They might be helpful.
Supplemental:
   suppl/CompareBpmCsvFiles.py
   suppl/CompareSnpInfoFiles.py
   suppl/FindSnpDuplicates.py
```
## [Data processing workflow](docs/processdata.md)

### Data delivery workflow
#### Given:
```
- Biobank_USERID_DATETIME_x10_FILETYPE.csv file from Biobank Portal OPS email
- specification of delivery format type (PED, VCF, imputedPED)
- delivery options (ERISONE, HPCWIN3, Mac/PC)
```
#### Step 1:
PED data release:
> python GetIndividualGenotypesPED.py -y etc/genotype_ped.r2.yml -i Biobank_USERID_DATETIME_x10_FILETYPE.csv

VCF data release:
> python GetIndividualGenotypesVCF.py -y etc/genotype_vcf.r2.yml -i Biobank_USERID_DATETIME_x10_FILETYPE.csv


### Imputation data delivery workflow
#### For a full data release:

#### Step 1:
For a release on a subset of data, please use the attached sample sheet to run the following command:
> python GetImputedGenotypes.py -y etc/impute_v1 -i Biobank_USERID_DATETIME_x10_FILETYPE.csv

*Please note:* This processing style can only handle one batch at a time. Please separate out the samples prior to processing and carefully reference the appropriate impute file paths.
