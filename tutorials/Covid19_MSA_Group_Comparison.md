# Comparison of Covid-19 Spike Protein Sequences Using Multiple Sequence Alignment

## Introduction to the Workflow

![welcome screenshot](./screenshots/0.0_welcome.png) 

The workflow compares two sets of peptide, DNA, or RNA sequences using following steps:
1. **Multiple sequence alignment** (MSA) and calculate a distance matrix using the Bioconductor’s [msa](http://www.bioconductor.org/packages/release/bioc/html/msa.html) package.
2. Generate and plot a **phylogenetic tree** by neighbor-joining using the packages [ape](https://guangchuangyu.github.io/software/ggtree/documentation/) and [ggtree](https://guangchuangyu.github.io/software/ggtree/documentation/).

There are two required input parameters (fasta_1 and fasta_2) that represent **two sets of sequences** in fasta format. Optional parameters include the sequence type, msa method, and distance type. The workflow’s main output is a plot of a **phylogenetic tree** showing the relation of the two sets.

The workflow can be found at: [https://github.com/CompEpigen/msa_group_compare](https://github.com/CompEpigen/msa_group_compare)

For further information, please also see:

https://w3id.org/cwl/view/git/93d3f03cdd9c44bdc609a11f097a4bad9451be84/CWL/workflows/msa_group_compare.cwl 