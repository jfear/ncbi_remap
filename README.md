# DrosSRA Workflow

## Objectives

- Generate a new fly annotation for release 6 of the genome
- Create a gene expression atlas across stages and tissue

## Outline

1. Download all *Drosophila melanogaster* RNA-seq data from SRA and remap to a 
   common reference using multiple pipelines.
2. Perform extensive QC to ensure metadata quality.
3. Using reference guided methods (Flybase r6.07) generate a new annotation.

## Notes

``$ENTREZ_API_KEY``

## Setup
```bash
$ conda env create --file environment.yaml
```

## Development
```bash
$ conda env update -n drosSRA_workflow --file environment_dev.yaml
```
