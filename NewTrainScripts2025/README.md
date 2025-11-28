# README

## Overview

The files used for training are organized into **three main categories**:

---

## 1. CSV Files

A single text file contains **all CSV files** that include:

* the keyword **`NpNpi`**
* the **highest energy range**

These CSV paths are grouped together for convenience.

---

## 2. JSON Files

There are **five text files**, each containing **both vector and scalar JSON configuration files**.
The number of JSON files in each is shown in parentheses:

* **AtmDUNE_FlatJson.txt** (612)
* **AtmDUNE_NatJson.txt** (596)
* **BeamDUNE_FlatAr40Json.txt** (256)
* **BeamDuneOnAxisND_NatJson.txt** (256)
* **BeamNOvAND_NatJson.txt** (256)

---

## 3. Shell Scripts for GPU Jobs

There are **five shell scripts**, one corresponding to each JSON list above:

* **AtmDUNE_Flat.sh**
* **AtmDUNE_Nat.sh**
* **BeamDUNE_FlatAr40.sh**
* **BeamDuneOnAxisND_Nat.sh**
* **BeamNOvAND_Nat.sh**

Each shell script is intended to be submitted on a **separate GPU node**.

---

## 4. Updating the JSON File Ranges in Each Shell Script

Inside every shell script, you will find a line like:

```bash
Json_files=( "${Json_files[@]:0:20}" )
```

This example selects the **first 20 JSON files**.
Each person must **change the range** to match the subset of files they are responsible for training.

### **File Assignments**

#### **AtmDUNE_FlatJson.txt (612 JSON files)**

* **Richi:** 0–203
* **Dr. Barrow:** 204–407
* **Casey:** 408–611

#### **AtmDUNE_NatJson.txt (596 JSON files)**

* **Richi:** 0–198
* **Dr. Barrow:** 199–397
* **Casey:** 398–595

#### **Beam JSON lists (each contains 256 files)**

Applies to:

* BeamDUNE_FlatAr40Json.txt
* BeamDuneOnAxisND_NatJson.txt
* BeamNOvAND_NatJson.txt

Ranges:

* **Richi:** 0–84
* **Dr. Barrow:** 85–169
* **Casey:** 170–255

---

## 5. Running the Shell Scripts

Each person must run their assigned **five shell scripts** on **five separate GPU nodes**, one script per node.

Make sure the following are done **before running**:

1. Clone or pull the **`wide` branch** of the `transformer_EE` repository from GitHub.
2. Download **all JSON, CSV, and shell script files** into the **same working directory**.
3. Ensure that each script has the correct slice of `Json_files` based on the ranges above.



