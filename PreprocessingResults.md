# Preprocessing Results: Resting State EEG Pipeline
**Authors:** Natalia Lopez, Julia Reina, Guiomar Niso  
**Date:** January 2026  
**Institution:** Cajal Institute (CSIC), Madrid, Spain  

---

## 1. Project Description
This pipeline performs automated Preprocessing of resting-state EEG signals using **Brainstorm**. It includes filtering, re-referencing, artifact cleaning, detection of bad segments or channels, PSD analysis and source computation. 

## 2. File Organization
The results are stored in JSON format for easy programmatic access.
- **Directory:** `/out_reports`
- **Pattern:** `Preproc-<subject>-<protocol>.json`

## 3. Data Dictionary
The JSON output contains the following fields:

| Field | Type | Description |
| :--- | :--- | :--- |
| `participant` | String | Subject ID from BIDS dataset. |
| `protocol`| String | Protocol name, as specified in the parameters. |
| `date`| String | Date and time of execution, in format YYYY-MM-DD HH:MM:SS |
| `numBadChannels` | Integer | Number of detected bad channels through peak-to-peak analysis. |
| `badChannelIndexes` | Vector | Indexes of the channels that were marked as bad. |
| `numBadSegments` | Integer | Number of bad segments detected throughout the recording via peak-to-peak analysis. |
| `numBlinks` | Integer | Number of blinks detected throughout the recording. |
| `numBadSegments` | Integer | Number of heartbeats detected throughout the recording. |
| `badSegments` | Structure | Initial and ending times of each bad segment detected, arranged in two different vectors. |
| `blinkTimes`| Structure | Array of vectors containing the times where blink events were detected. | 
| `cardiacTimes`| Structure | Array of vectors containing the times where cardiac events were detected. |

---



