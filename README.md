# MDMA/Meth Prosocial Effects - scRNA-seq

In this study, we applied single-cell RNA sequencing (scRNA-seq) to examine how 3,4-methylenedioxymethamphetamine (MDMA) alters transcriptional activity in the nucleus accumbens (NAc), a brain region central to social behavior and reward. MDMA produces robust prosocial behaviors in mice, whereas methamphetamine (Meth), despite structural similarity, fails to elicit such effects and instead drives addictive phenotypes. Our goal was to begin delineating RNA-level mechanisms that distinguish these divergent outcomes. By comparing the transcriptional activity elicited by each drug, we can better identify the cellular mechanisms in the NAc that promote prosocial behavior versus those that drive addiction.

Using TRAP2;Ai14 mice, in which neuronal activity–dependent Fos expression drives tdTomato labeling, we profiled transcriptional responses across five experimental groups: MDMA (male, acute male, acute female), Meth (male), and saline (male). Two mice were pooled per replicate (n = 3 per condition). NAc tissue was dissected, FACS-enriched for tdTomato+ cells (~90–95%), and processed using the 10x Genomics platform.

Preliminarily we performed systematic preprocessing including ambient RNA correction and contextual filtering to recover high quality cells and their mRNA expression. Then, we performed cell clustering to identify major subtypes of cells based on expression similarities and cross validated our cell type annotations with multiple published NAc single cell atlases.

Then, we focused on neuronal cell populations and used gene level and pathway level differential expression analysis across treatments to identify transcriptional patterns (i.e. upregulation or downregulation) that appear uniquely associated with MDMA’s prosocial effects.

We hypothesized that MDMA preferentially induces immediate early gene programs and regulates receptor and signaling transcripts in neuronal subsets implicated in social behavior, while Meth is more strongly aligned with addiction pathways. These exploratory findings provide an initial RNA-level perspective on MDMA’s action in the NAc and motivate future work to clarify molecular targets that may disentangle prosocial benefits from addictive liabilities.

## Poster Presentation

- Faghihzadeh, S., Arab, A., Wu, X., & Heifets, B. D. (2025). **Single-cell transcriptomic insights into the neural mechanisms of MDMA’s prosocial and addictive effects in the nucleus accumbens.** _3rd NYC RNA Symposium, New York, NY_

<p align="center">
  <img src="2025 NYC RNA Symposium Poster.png" alt="Poster" width="600" />
</p>

<!-- ## Related datasets
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
    - Run Scrublet with `batch_key='batch'` option.  -->
