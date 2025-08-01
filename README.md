# Single-Cell RNA Sequencing of Brain Metastases in Mouse Models
Here we re-analyse a data set [published](https://www.cell.com/cancer-cell/abstract/S1535-6108(24)00314-3) in the group of J. Massague.

The experimental setup contains sorted cancer cells from (brain-tropic variants of triple-negative breast cancer model MDA231 and HER2+ HCC1954), labeled tumor microenvironment (mouse cells in proximity of the cancer cells) and unlabelled mouse cells from the brain. Samples from both cell lines were pooled for sequencing and using antibody-based HTOs for demultiplexing. Experiment was performed in 2 batches.

The analysis contains the following steps:
- QC and removal of damaged cells
- sample demultiplexing using HTOs
- doublet detection
- batch correction
- cell annotation using multiple references (Allen brain atlas, ImmGen)

![Mouse cells](https://github.com/MikeKlocCZ/2025_scRNA_brain_mets_MassaugeJ/blob/main/figures_examples/04_TSNE_annotation-0.png "Annotated mouse cells")

We also included cell-cycle analysis (using cyclone) and in the case of Microglia (the biggest cell population) we performed pseudotime analysis using scVelo (wrapped in a velociraptor R package)

<img src="https://github.com/MikeKlocCZ/2025_scRNA_brain_mets_MassaugeJ/blob/main/figures_examples/05_TSNE_velocity.png " width="600">