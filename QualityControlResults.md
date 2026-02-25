# Quality Control Results: Resting State EEG Pipeline
**Authors:** Natalia Lopez, Julia Reina, Guiomar Niso  
**Date:** January 2026  
**Institution:** Cajal Institute (CSIC), Madrid, Spain  

---

## 1. Project Description
This pipeline performs automated Quality Control (QC) and Power Spectral Density (PSD) analysis on resting-state EEG signals using **Brainstorm**. 

## 2. File Organization
The results are stored in JSON format for easy programmatic access.
- **Directory:** `/out_reports`
- **Pattern:** `QC-<subject>-<protocol>.json`

## 3. Data Dictionary
The JSON output contains the following fields:

| Field | Type | Description |
| :--- | :--- | :--- |
| `participant` | String | Subject ID from BIDS dataset. |
| `protocol`| String | Protocol name, as specified in the parameters. |
| `date`| String | Date and time of execution, in format YYYY-MM-DD HH:MM:SS |
| `avgPerChannelPre` | Vector | Mean voltage per sensor ($\mu V$) prior to processing. |
| `stdPerChannelPre` | Vector | Signal standard deviation per sensor prior to processing. |
| `peakFrequencies` | Vector | Frequencies ($ Hz $) where PSD peaks were detected. |
| `peaks` | Vector | Power amplitude at the detected frequencies ($\mu V^2/Hz$). |

---


