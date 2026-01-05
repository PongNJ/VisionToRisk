# Vision to Risk

**From Vision to Risk: Explaining Motorcycle Hotspot Risk Using Vision–Language Attribute Extraction**

This repository contains the **data processing, analysis, and modelling code** used in the *Vision to Risk* study. The project integrates **network-based hotspot detection** with **vision-based infrastructure feature extraction** from street-level imagery to analyse motorcycle crash risk.

The workflow combines:

* Exploratory spatial analysis
* Network Kernel Density Estimation (NKDE)
* Vision–Language Model (VLM)–based road attribute extraction
* Hotspot–control comparison and factor analysis

The code is organised to reflect the **methodological pipeline of the paper**.

---

## Repository Structure

```
VisionToRisk/
│
├── exploratory/
│   ├── DOH.ipynb
│   ├── OSM_explore.ipynb
│
├── nkde/
│   ├── NKDE_TRAMS.R
│   ├── NKDE_TRAMS_LAMPANG_TAK.R
│   ├── nkde_result.ipynb
│
├── vlm/
│   ├── SVI_feature_retrival.ipynb
│
├── analysis/
│   ├── hotspot_and_control.ipynb
│   ├── hotspot_analysis.ipynb
│   ├── factor_analysis.ipynb
│
│
└── README.md
```

---

## Exploratory Analysis

**Directory:** `exploratory/`

These notebooks are used for **initial data inspection and spatial understanding** prior to formal modelling.

* `DOH.ipynb`
  Exploratory analysis of Department of Highways (DOH) road network and crash data.

* `OSM_explore.ipynb`
  Exploration and validation of OpenStreetMap (OSM) road classifications and network attributes used in NKDE.

---

## Network Kernel Density Estimation (NKDE)

**Directory:** `nkde/`

This folder contains the scripts used to identify **motorcycle crash hotspots on road networks** using Network Kernel Density Estimation.

* `NKDE_TRAMS.R`
  Core NKDE implementation using national crash data.

* `NKDE_TRAMS_LAMPANG_TAK.R`
  NKDE analysis focused on Lampang–Tak case study region.


* `nkde_result.ipynb`
  Post-processing and visualisation of NKDE outputs.

---

## Vision–Language Model (VLM) Feature Extraction

**Directory:** `vlm/`

This component extracts **road infrastructure and roadside attributes** directly from street-level imagery.

* `SVI_feature_retrival.ipynb`
  Vision–Language Model pipeline for extracting road safety attributes from Street View Imagery (SVI), aligned with the study’s infrastructure taxonomy.

---

## Hotspot and Factor Analysis

**Directory:** `analysis/`

These notebooks link **visual attributes** with **hotspot likelihood**.

* `hotspot_and_control.ipynb`
  Construction of hotspot and control samples from NKDE outputs.

* `hotspot_analysis.ipynb`
  Descriptive and spatial analysis of identified hotspot segments.

* `factor_analysis.ipynb`
  Statistical analysis examining associations between extracted visual attributes and hotspot occurrence.


---

## Data Availability

Due to data licensing and privacy restrictions, **raw crash data and imagery are not included** in this repository.
Scripts assume access to:

* Road network data (OSM / DOH)
* Crash records (processed form)
* Street View Imagery (e.g., Mapillary)

---

## Reproducibility Notes

* NKDE scripts are written in **R**
* Feature extraction and analysis notebooks are written in **Python**
* Coordinate reference systems and network segmentation choices are documented in the notebooks
* Results may vary depending on kernel bandwidths and network preprocessing

