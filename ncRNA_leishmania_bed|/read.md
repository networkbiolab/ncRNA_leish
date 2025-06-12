# Predicted ncRNA BED Files

This folder contains **BED files** with the coordinates of **predicted non-coding RNAs (ncRNAs)** for each species analyzed in this project.

## Description

Each `.bed` file includes the genomic coordinates of regions predicted as ncRNAs, along with essential information for visualization and further analysis using genome browsers such as IGV or UCSC Genome Browser.

These predictions were obtained using a combined approach of **ab initio** prediction and **comparative searches** against ncRNA databases, followed by filtering based on RNA-seq expression data (FPKM ≥ 1).

## Contents

- One `.bed` file per species, named using the corresponding species identifier (e.g., `Lmajor_ncRNAs.bed`, `Linfantum_ncRNAs.bed`, etc.).
- Each line in a BED file represents a predicted ncRNA with the following fields:
  - **chrom**: chromosome or contig
  - **start**: start position (0-based)
  - **end**: end position (1-based)
  - **name**: ncRNA identifier
  - **score**: optional score (if applicable)
  - **strand**: strand orientation (`+` or `-`)

## Usage

These files can be used for:
- Visualization in genome browsers
- Overlap analysis with known annotations
- Integration with other genomic or transcriptomic datasets

## Notes

- Make sure to use the correct reference genome assemblies corresponding to each species when working with these files.
- File names include the species name or abbreviation for easy identification.

---

**Maintainer:** Eduardo Martínez  
**Project:** ncRNA Prediction in *Leishmania* Species  
**mail:** jmartinezh@udla.cl or lalomartinez92@gmail.com
**paper:** # Predicted ncRNA BED Files

This folder contains **BED files** with the coordinates of **predicted non-coding RNAs (ncRNAs)** for each species analyzed in this project.

## Description

Each `.bed` file includes the genomic coordinates of regions predicted as ncRNAs, along with essential information for visualization and further analysis using genome browsers such as IGV or UCSC Genome Browser.

These predictions were obtained using a combined approach of **ab initio** prediction and **comparative searches** against ncRNA databases, followed by filtering based on RNA-seq expression data (FPKM ≥ 1).

## Contents

- One `.bed` file per species, named using the corresponding species identifier (e.g., `Lmajor_ncRNAs.bed`, `Linfantum_ncRNAs.bed`, etc.).
- Each line in a BED file represents a predicted ncRNA with the following fields:
  - **chrom**: chromosome or contig
  - **start**: start position (0-based)
  - **end**: end position (1-based)
  - **name**: ncRNA identifier
  - **score**: optional score (if applicable)
  - **strand**: strand orientation (`+` or `-`)

## Usage

These files can be used for:
- Visualization in genome browsers
- Overlap analysis with known annotations
- Integration with other genomic or transcriptomic datasets

## Notes

- Make sure to use the correct reference genome assemblies corresponding to each species when working with these files.
- File names include the species name or abbreviation for easy identification.

---

**Maintainer:** Eduardo Martínez  
**Project:** Comparative and systems analyses of _Leishmania_ spp. non-coding RNAs through developmental stages  
**mail:** jmartinezh@udla.cl or lalomartinez92@gmail.com  
**paper:** https://doi.org/10.1371/journal.pntd.0013108 
