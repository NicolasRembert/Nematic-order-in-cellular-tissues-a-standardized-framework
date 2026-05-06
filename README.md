# Nematic-order-in-cellular-tissues-a-standardized-framework-and-anomalous-defect-dynamics

## Overview

Cellular monolayers often exhibit orientational order, with nematic alignment of cell shape and cytoskeletal structures governing tissue-scale collective dynamics. Despite extensive studies, a unified analysis framework for characterizing active nematics in living systems remains partial, and key discrepancies with theory persist.

Here, we present a systematic and comparative analysis of nematic order and tissue flow dynamics across twelve distinct cell types. We quantify the impact of analysis parameters and provide data-driven guidelines to improve reproducibility and cross-study comparability. Across all nematic systems, we uncover remarkably consistent static properties, supporting the universality of nematic behavior in living tissues.

---

## Summary figure

<p align="center">
  <img src="Figure/Summary figure.png" width="900">
</p>

The pipeline illustrates the full workflow used in this repository:
1. Extraction of orientation fields from microscopy **square** images  
2. Computation of nematic order and defect detection, Quantification of defect dynamics  
3. Velocity-field analysis using PIV and quantification of cell flow around defects   

---

## Analysis workflow

The typical workflow is:

### 1. Orientation field generation
Orientation vector fields are computed from microscopy images using **ImageJ / OrientationJ**,and then saved as csv. Multiple csvs can be generated with varying **feature size** and **Down-sampling** size.

### 2. Coherency analysis
Parameter exploration to determine optimal settings for OrientationJ analysis.

### 3. Nematic analysis
Using MATLAB:

- computation of nematic correlation length
- computation of nematic order
- detection of ±1/2 topological defects
- defect tracking in time series
- defect densities
- nearest neighbor distance
- TAMSD (time-averaged mean square displacement)

### 4. Flow analysis around defects
Velocity fields obtained with **PIVlab** are used to compute:

- velocity magnitude around defects
- velocity divergence




### Main scripts

**Python**

`Generating vector fields on stitched.py`  
Generates nematic orientation vector fields from microscopy images using outputs from **ImageJ / OrientationJ**, user can change the two main parameters "feature size" and "Down-sampling factor" by varying respectively the **tensor** and **grid** variables in the script. This code saves the csv results tables with the formatting required for further analysis with MATLAB.



**MATLAB**

`coherency_analysis.py`  
Allows visualization and quantification of the **coherency ** as a function of feature size set on OrientationJ .

`Main_tracking_defects.m`  

- detection of topological defects
- tracking of defects in time series

`Velocity_around_defect_analysis.m`  
Computes the **dynamics of cells around topological defects** using velocity fields obtained from **PIVlab**.

---

---

## Example datasets

A small dataset is included in the repository **Confluence 48h cell data** to reproduce the **defect dynamics analysis**.

### Additional datasets

Due to file size limitations, additional datasets are hosted externally.

**Coherency test dataset**

[Download dataset](https://drive.switch.ch/index.php/s/s3NtlIEZ17n7KVl)

**PIV dataset for defect dynamics**

[Download dataset](https://drive.switch.ch/index.php/s/ajCHMFG8D5xC6yB)

