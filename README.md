# Single-Cell RNA-seq Analysis
This workflow aims to identify neurons involved in stress-dependent activation through the exploration of a publicly available single-cell atlas. The cells were subsetted for those expressing pdyn, crhb, and/or avp. Differential expression analysis was performed and visualized through volcano plots. 

## Table of Contents 
* Installation Requirements 
* Set Up + Dataset 
* Author
* References 

## Installation Requirements
### Software 
* R 
* R studio
* Excel (refer to DESEQ.R) 

### Packages 
* Seurat 
* dplyr 
* pasilla
* DESeq2
* Enhanced Volcano

## Dataset + Set up
The data used in this workflow came from [this paper](https://doi.org/10.1016/j.neuron.2020.09.023) and can be downloaded [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE158142). However, this workflow may have flexibility for other scRNA datasets. 

The data was downloaded and set up in the following way (with a singular barcodes, genes, and matrix file per sample folder) : 
![image](https://user-images.githubusercontent.com/87867639/162547844-7c341e5d-eae6-45ab-b84a-50a3dd97763c.png)

## Author 
Amara Alexander, Undergraduate Research Assistant, VA-MD College of Veterinary Medicine.

Mentor: Dr. Xie at the Epigenetics and Computational Laboratory located within the VA-MD College of Veterinary Medicine

## References 
As I am writing this almost a year after developing this workflow, I have no recollection of the many references I used, so I will include some helpful resources for scRNA analysis. 
* [Github page with large compilation of scRNA analysis resources](https://github.com/crazyhottommy/scRNAseq-analysis-notes)






