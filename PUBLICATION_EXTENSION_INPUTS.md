# Publication Extension Inputs

The frozen six-cohort atlas and native-state discovery are never retrained with
these inputs. Every extension is validation-only.

## Automatically downloaded

`GSE308624` is a tumour-only CosMx cohort with 71 gastric tumour sections. The
runner downloads the public H5AD and clinical workbook from GEO and evaluates
frozen program neighbourhoods. It cannot test tumour-versus-normal depletion.

`GSE163558` is also downloaded and reconstructed from its public 10x count
matrices. It contains three primary tumours, one adjacent non-tumoral sample,
and six metastases from six patients. It is used only for descriptive
primary/adjacent/metastatic context because it has fewer than three auditable
tumour-normal pairs. It is never added to the six-cohort atlas or the
tumour-normal random-effects meta-analysis.

`GSE189926` is an automatically downloadable public CD45-selected single-cell
gastric tumour cohort from 13 patients treated with disulfiram plus nivolumab,
with 22 pre/post samples and response labels in the GEO SOFT record. The
extension runner builds sample-level immune pseudobulk directly from the public
count matrices and projects only frozen myeloid/T/NK state signatures. It is a
small, single-arm treatment-context analysis: it does not test CAF/endothelial
states, train a response model, or produce an AUROC.

`TCGA-STAD` molecular-subtype context can also be prepared automatically from
the original open GDC 2014 publication supplements: the 291-sample RPKM matrix
and the 295-patient master table containing EBV, MSI, genomically stable, and
CIN labels. The runner restricts this to matched primary tumours and uses it
only for frozen-signature versus subtype associations.

## Optional paired spatial validation

Place a raw-count H5AD at:

`data/external/validation_extensions/PAIRED_GC_SPATIAL/paired_tumor_normal_spatial.h5ad`

It must have unique cell/spot IDs, integer counts in `X` or `layers['counts']`,
spatial coordinates in `obsm['spatial']`, and explicit `patient_id`,
`sample_id`, and tumour/normal condition metadata. At least three paired
patients are required before it can be considered an inferential tumour-normal
validation.

As of the public-data audit, neither `GSE251950` (nine primary tumours plus one
primary-metastatic pair) nor `GSE308624` (71 tumour sections) is a paired
tumour-normal spatial cohort. They support spatial architecture only and must
not be described as tumour-normal spatial replication.

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

## What the extension runner does when inputs are absent

The runner still completes Korean longitudinal/pathology context, GSE163558
descriptive context, and GSE308624 tumour-only spatial ecology. With
`--download-gse189926` and `--download-tcga`, it also completes the two public,
auditable extension branches above. It writes one preflight row and one manifest
entry for every unavailable paired-spatial or optional PRJEB component. Missing
inputs are reported as `not_run`/`pending`; no substitute cohort, inferred
response label, or synthetic clinical result is created.
