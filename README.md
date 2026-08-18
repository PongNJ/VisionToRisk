# VisionToRisk

Code accompanying the paper:

**Street view-based infrastructure assessment of motorcycle hotspots: A national road network analysis in Thailand**  
Accident Analysis & Prevention, 2026. DOI: `10.1016/j.aap.2026.108706`

This repository contains the reproducible analysis workflow used to connect motorcycle casualty hotspots on Thailand's road network with infrastructure attributes observed from street-level imagery. The workflow combines crash data preparation, OpenStreetMap road-network processing, Network Kernel Density Estimation (NKDE), Mapillary street-view sampling, vision-language-model feature extraction, and hotspot/control factor analysis.

Raw crash records, road-network extracts, street-view images, API credentials, and generated model outputs are not included because of data licensing, privacy, and file-size constraints.

## Repository Structure

```text
VisionToRisk/
|-- exploratory/
|   |-- TRAMS.ipynb
|   `-- OSM_explore.ipynb
|-- nkde/
|   |-- NKDE_TRAMS.R
|   |-- NKDE_TRAMS_LAMPANG_TAK.R
|   `-- nkde_result.ipynb
|-- vlm/
|   |-- SVI_feature_retrival.ipynb
|   `-- index.js
|-- hotspot_analysis/
|   |-- hotspot_and_control.ipynb
|   |-- hotspot_analysis.ipynb
|   `-- factor_analysis.ipynb
|-- fig/
|   |-- SHAP.png
|   |-- SHAP_bee.png
|   `-- SHAP_noAADT.png
```

## Workflow

1. **Crash and network exploration**  
   `exploratory/TRAMS.ipynb` prepares and inspects Thailand motorcycle crash records. `exploratory/OSM_explore.ipynb` evaluates OpenStreetMap road classes and network attributes used for network-based analysis.

2. **Network hotspot detection**  
   `nkde/NKDE_TRAMS.R` and `nkde/NKDE_TRAMS_LAMPANG_TAK.R` implement NKDE for motorcycle casualty hotspots. `nkde/nkde_result.ipynb` post-processes hotspot outputs and evaluates hotspot capture using severity-weighted crash measures.

3. **Street-view infrastructure extraction**  
   `vlm/SVI_feature_retrival.ipynb` samples lixel locations, retrieves Mapillary street-level imagery, and extracts infrastructure attributes with a vision-language-model workflow aligned with the paper's road-safety feature taxonomy.

4. **Hotspot/control analysis**  
   `hotspot_analysis/hotspot_and_control.ipynb` constructs hotspot and control samples. `hotspot_analysis/hotspot_analysis.ipynb` summarizes hotspot characteristics, and `hotspot_analysis/factor_analysis.ipynb` estimates associations between extracted infrastructure features and hotspot status.

## Data Availability

The repository intentionally excludes:

- raw TRAMS/crash records;
- downloaded OpenStreetMap or other road-network data;
- Mapillary images and API credentials;
- intermediate geospatial files;
- generated SHAP outputs, model comparison tables, and scratch figures.

Users should place local inputs in ignored data/output directories and update notebook paths for their own environment before running the workflow.

## Software Notes

The analysis uses both R and Python notebooks/scripts. Core geospatial and modelling dependencies include typical scientific Python/R packages for spatial data handling, network analysis, machine learning, and visualization. Exact package versions should be recorded from the execution environment used to reproduce the paper results.

## Citation

If you use this repository, please cite the paper:

```text
Street view-based infrastructure assessment of motorcycle hotspots:
A national road network analysis in Thailand.
Accident Analysis & Prevention, 2026.
https://doi.org/10.1016/j.aap.2026.108706
```
