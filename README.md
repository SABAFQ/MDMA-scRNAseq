# MDMA-scRNAseq

## Related datasets
[__Single-cell profiling and functional characterization bridge molecular and functional features of NAc neuron subtype__](https://www.nature.com/articles/s41593-021-00938-x#data-availability)

[![GEO:GSE118020](https://img.shields.io/badge/GEO-GSE118020-green.svg)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118020)
```bash 
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE118nnn/GSE118020/suppl/GSE118020_expression_matrix.csv.gz
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE118nnn/GSE118020/suppl/GSE118020_meta_data.csv.gz
```

## Analysis Workflow

**Aim:** Preprocessing scRNA-seq datasets

**Tools**: 
- [`scanpy`](https://scanpy.readthedocs.io/en/stable/)
    - https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html#Preprocessing
        - Identify highly-variable genes

**Tasks:**
- [X] Load `adata` objects and merge them into one big dataset
- [X] Identify "dead cells"!
- [X] Data cleaning
    - In the first round of the analysis we have seen that we need to find more dead cells!
    - Run Scrublet with `batch_key='batch'` option. 
