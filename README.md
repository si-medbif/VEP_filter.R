# VEP_filter.R
Variant prioritization algorithm developed by Nutchavadee Vorasan for SUDS project to filter variants annotated with VEP by R

VEP annotation was performed using GRCh37 (hg19) annotation 

# VEP_filter.Chunky.template.R

Editor: K2

Synopsis: Read file in chunks (for Large file size)

Support: VEP95.2 , dbSNP 3.5a


# Oncotator_filter.R

Editor: K2

Synopsis:  for filtering results from Oncotator
Filtering default: 

1. Select variants with 'PASS' 

2. Filter with Impact (High, Moderate) 

3. Filter with LR and SVM (include "D" or "-") 

4.  Filter with CADD score (exculde values which is "< 15") 

5. Filter with Allele Frequency (Using General AF to exclude variant that greater than  10% (>0.1) )

* **able to adjust filtering steps**

Support: Oncotator1.9.8.0+
