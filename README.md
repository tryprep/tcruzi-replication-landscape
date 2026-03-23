# Trypanosoma cruzi Replication Landscape 🧬

This repository contains the analytical frameworks, pipelines, and scripts developed for the study:  
**"Genome compartmentalization is linked to the replication program and mutational outcomes in Trypanosoma cruzi"**  
*Franco, T. A., Freitas, R. P., Pires, D. S., Souza, V. C., Damasceno, J. D., McCulloch, R., & Elias, M. C. (2026).*

---

## 🔬 Project Overview

Our research investigates how the spatial organization of the *Trypanosoma cruzi* genome influences its evolution. By integrating high-resolution **D-NAscent** mapping (Nanopore sequencing) with variant calling, we delineate the contrasting replication dynamics between the **Core** (conserved) and **Disruptive** (rapidly evolving) compartments and their impact on the mutational landscape (SNPs and Loss of Heterozygosity - LOH).

### 👥 Contributors
*   **Thiago Andrade Franco** (Lead Developer/First Author)
*   **David da Silva Pires**
*   [**Vinícius Carius de Souza**](https://github.com/vcarius)
*   **Maria Carolina Elias** (Principal Investigator/Corresponding Author)

---

## 🛠️ Repository Content
This repository is organized to provide the specific workflows used for:
*   **Replication Mapping:** Implementation of Replication Fork Directionality (RFD) and Origin Efficiency Metrics (OEM).
*   **Genomic Landmarks:** Identification of Initiation (IZD) and Termination (TZD) Zone Domains.
*   **Compartment Analysis:** Classification of Polycistronic Units (PTUs) and genomic regions (Core vs. Disruptive/GpDR).
*   **Mutational Load:** Pipelines for variant calling (SNPs) and quantification of allelic imbalance (LOH index).

---

## 🛠️ Technical Stack

The project utilizes a multi-language approach and industry-standard bioinformatics tools:

*   **Languages:** Python, R, C++, Shell Script
*   **Bioinformatics Tools:**
    *   `Freebayes`: Variant calling (SNPs/Indels)
    *   `BCFtools`: SNP analysis and filtering
    *   `BEDtools`: Genomic coordinate processing
    *   `D-Nascent`: Replication fork identification
    *   `Winnowmap`: Long-read mapping
    *   `Picard` & `Samtools`: Data processing and quality control
    *   `SnpEff`: Variant annotation and dN/dS estimation

---

## ⚙️ Pipelines Description

The analysis is divided into three specialized workflows:

### 1. DNA Replication Dynamics Identification
This pipeline processes raw outputs to identify and classify replication events at base-pair resolution.
*   **Functionality:** Employs C++ programs and ad-hoc scripts to construct and classify datasets for replication origins (IZDs) and termination sites (TZDs).
*   **Key Output:** A curated, high-precision dataset of genomic replication landmarks.

### 2. Genomic Mapping & Enrichment Analysis
Focuses on the spatial distribution of replication events across genomic compartments and polycistronic units.
*   **Mapping:** Origins and terminations are intersected with genomic features using `BEDtools` and custom Shell/R scripts.
*   **Polycistron Definition:** Dedicated pipeline for defining polycistronic regions based on gene family and GC content (as detailed in the paper’s Methods).
*   **Statistics:** Calculates fold-enrichment of replication events relative to genome-wide expected frequencies.

### 3. Variant Analysis & Mutational Asymmetry
A comprehensive workflow for analyzing genetic diversity and mutational bias.
*   **Variant Calling:** High-quality read selection and mapping followed by variant discovery.
*   **Asymmetry Analysis:** Calculates strand-specific SNP density at replication landmarks and performs **Loss of Heterozygosity (LOH)** analysis.
*   **Functional Impact:** Variant annotation via `SnpEff` to provide the basis for dN/dS ratio calculations and evolutionary interpretation.

---

## 📊 Data Availability

The raw sequencing data analyzed in this study are publicly available in the **NCBI Sequence Read Archive (SRA)** under:
*   **BioProject:** [PRJNA1002335](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1002335)
*   **Accession Numbers:** SAMN36840645 to SAMN36840649.

---

## 🚀 How to Use

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/Tryp-Rep-Laboratory/tcruzi-replication-landscape.git
    ```
2.  **Environment:** Ensure all tools listed in the "Technical Stack" are installed (we recommend using **Conda** for environment management).
3.  **Scripts:** Detailed execution instructions for each script can be found within their respective directories.

---

## 📝 Citation

If you use these pipelines or data in your research, please cite:
> **Franco, T. A., et al. (2026). Genome compartmentalization is linked to the replication program and mutational outcomes in Trypanosoma cruzi. [Journal Name/DOI Link].**

---
© 2024 Tryp-rep Lab - Butantan Institute
