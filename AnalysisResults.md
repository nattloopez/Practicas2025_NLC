# Analysis Results: Resting State EEG Pipeline
**Authors:** Natalia Lopez, Julia Reina, Guiomar Niso  
**Date:** January 2026  
**Institution:** Cajal Institute (CSIC), Madrid, Spain  

---

## 1. Project Description
This pipeline performs automated frequency-based analysis on both sensors and sources of resting-state EEG signals using **Brainstorm**. 

## 2. File Organization
The results are stored in JSON format for easy programmatic access.
- **Directory:** `/out_reports`
- **Pattern:** `Analysis-<subject>-<protocol>.json`

## 3. Data Dictionary
The JSON output contains the following fields:

| Field | Type | Description |
| :--- | :--- | :--- |
| `participant` | String | Subject ID from BIDS dataset. |
| `protocol`| String | Protocol name, as specified in the parameters. |
| `date`| String | Date and time of execution, in format YYYY-MM-DD HH:MM:SS |
---



