# Publication Extension Inputs

The frozen six-cohort atlas and native-state discovery are never retrained with
these inputs. Every extension is validation-only.

## Automatically downloaded

`GSE308624` is a tumour-only CosMx cohort with 71 gastric tumour sections. The
runner downloads the public H5AD and clinical workbook from GEO and evaluates
frozen program neighbourhoods. It cannot test tumour-versus-normal depletion.

## Optional paired spatial validation

Place a raw-count H5AD at:

`data/external/validation_extensions/PAIRED_GC_SPATIAL/paired_tumor_normal_spatial.h5ad`

It must have unique cell/spot IDs, integer counts in `X` or `layers['counts']`,
spatial coordinates in `obsm['spatial']`, and explicit `patient_id`,
`sample_id`, and tumour/normal condition metadata. At least three paired
patients are required before it can be considered an inferential tumour-normal
validation.

## PRJEB25780 pembrolizumab validation

The public study is `PRJEB25780` (pembrolizumab in metastatic gastric cancer):
https://www.ebi.ac.uk/ena/browser/view/PRJEB25780

The public raw sequencing release is too large and does not itself provide an
auditable sample-response table. Do not infer response labels from FASTQ names.
Use a matched processed, gene-by-sample expression table and clinical table:

`data/external/validation_extensions/PRJEB25780/PRJEB25780_expression.csv`
`data/external/validation_extensions/PRJEB25780/PRJEB25780_clinical.csv`

The expression CSV must have a `gene`/`gene_symbol` first column and one column
per clinical sample. The clinical CSV needs `sample_id` plus a response field
encoded as CR/PR/SD/PD or 1/0. The runner tests frozen signatures only; it does
not fit or report a response-prediction model.

## TCGA molecular subtype context

Place matched tables at:

`data/external/validation_extensions/TCGA_STAD_SUBTYPES/TCGA_STAD_expression.csv`
`data/external/validation_extensions/TCGA_STAD_SUBTYPES/TCGA_STAD_clinical.csv`

The clinical table needs `sample_id` and `molecular_subtype` (or
`tcga_subtype`/`subtype`). The analysis is a frozen-signature association with
supplied EBV/MSI/CIN/genomically stable labels; it is not a subtype classifier.
