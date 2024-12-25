# **README: From producing CSV files & Running training to Plotting**  

## **Step 1: CSV and Signal Selection Files**  

### **CSV File Organization**  
The work is divided based on whether leptons are treated as **scalars** or **vectors**, as described below:  

- **Vector Leptons without NC** – Only charged leptons enter as vectors; neutral leptons are excluded.  
- **Vector Leptons with NC** – Both charged and neutral leptons enter as vectors.  
- **Vector Leptons with 0 NC** – Charged leptons enter as vectors, while neutral leptons have all kinematic values set to zero.  
- **Scalar Leptons with 0 NC** – True kinematic information for charged leptons entering as scalars; neutral leptons have kinematic values set to zero.  
- **Scalar Leptons with NC** – Both charged and neutral leptons enter as scalars.  

---

### **Signal Selection Format**  

**Sample Signal Selection File:**  

```
EventNum: All            # Use 'All' to process all events or specify a number for a subset  
Initial_Nu_Energy_min: 0.1  
Initial_Nu_Energy_max: 1  
Initial_PDG: Any         # Use 'Any' for all neutrinos or specify a PDG code for a specific neutrino  
Final_state_lepton_PDG_code_or_process: Inclusive  # Use 'CC' for charged leptons, 'NC' for neutral leptons, or 'Inclusive' for both.  
Proton_KE: 0.025         # Kinetic energy in GeV  
Pi+-_KE: 0.07  
K+-_KE: 0.03  
Muon_KE: 0.005  
Electron_KE: 0.005  
num_proton: N            # Use 'N' for any number of protons/pions or specify a value  
num_pion: N  
INPUT_ROOT_FILE: NNBarAtm_hA_BR.100000000.gtrac.root  
OUTPUT_ROOT_FILE: AnyNu_Inclusive_Thresh_p1to1_  
OUTPUT_DIR: AnyNu_Inclusive_Thresh_p1to1_  
OUTPUT_NAME: AnyNu_Inclusive_Thresh_p1to1_  
```  

---

### **Running C++ Files**  

To run a single signal selection file with the C++ code, use the following command:  
```bash
root -b -q 'ScalarLept_wNC.C("Signal_Selection1.txt")'
```  
**Note:**  
- Use the provided bash script to process multiple text files simultaneously.  

---

## **Step 2: JSON File Production**  

### **General Notes:**  
JSON files are created based on whether leptons are treated as **scalars** or **vectors**. The process involves generating configurations that define vector, scalar, and target types.  

---

### **Directory Setup**  
```python
# Define the directory containing CSV files
csv_directory = "/exp/dune/app/users/rrichi/FinalCodes"

# Path to save generated JSON configuration files
output_directory = "/home/rrichi/TransformerPolars/transformer_EE/transformer_ee/config/"
```  

---

### **Defining Vector, Scalar, and Target Types**  

**Vector Types** (Final state particles):  
```python
vector_type = [
    "Final_State_Particles_PDG",
    "Final_State_Particles_Mass",
    "Final_State_Particles_Energy",
    "Final_State_Particles_Momentum_X",
    "Final_State_Particles_Momentum_Y",
    "Final_State_Particles_Momentum_Z",
    "Final_State_Particles_CosTheta",
    "Final_State_Particles_Theta"
]
```  

**Scalar Types**:  
```python
scalar_types = [
    "Lept_PDG",
    "Lept_Mass",
    "Lept_Energy",
    "Lept_MomX",
    "Lept_MomY",
    "Lept_MomZ",
    "Lept_CosTheta",
    "Lept_Theta",
    "tot_fKE",
    "p_tot",
    "P_miss"
]
```  

**Target Types** (Variables for predictions):  
```python
target_types = [
    "Nu_Theta",
    "Nu_CosTheta"
]
```  

---

### **Naming Scheme for JSON Files**  
```python
# Convert target types to string format for filenames
target_type_str = '_'.join(target_type)

# Generate the JSON filename based on CSV name and target type
json_filetojoin = f"input_GENIEv3-0-6-Honda-Truth-hA-LFG_{os.path.splitext(csv_file)[0]}_{target_type_str}_{acronym}.json"
json_filename = os.path.join(output_directory, json_filetojoin)
```  

---

### **Saving Model and Paths**  
```python
"save_path": f"/home/rrichi/TransformerPolars/save/model/LossVars_{target_type_str}/GENIEv3-0-6-Honda-Truth-hA-LFG_{os.path.splitext(csv_file)[0]}_{acronym}",
"model_phys_name": f"GENIEv3-0-6-Honda-Truth-hA-LFG_{os.path.splitext(csv_file)[0]}_{acronym}"
```  

---

**Note:**  
- JSON files for **vector leptons** are generated in a similar way to scalar leptons.  
---

## Step 3: Training File

### Specify the Directories for CSV and JSON Files
```python
csv_directory = "/exp/dune/app/users/rrichi/FinalCodes"
json_directory = "/home/rrichi/TransformerPolars/transformer_EE/transformer_ee/config"
```

### Shortening Names for WandB
The following script reduces the name size for WandB logging:
```python
target_string = base_jsonfile.split('_' + acronym)[0].split('pi_', 1)[-1]
id_wandb = base_csvfile.replace('input_GENIEv3-0-6-Honda-Truth-hA-LFG_', '') \
    .replace('_Thresh', '') \
    .replace('_eventnum', '')
```

### Logging Information in WandB
Use the following format to view training details in WandB:
```python
my_logger = WandBLogger(
    project="GENIE-Atmo",
    entity="neutrinoenenergyestimators",
    config=input_d,
    dir="/home/rrichi/TransformerPolars/save",
    id=id_wandb
)
```

**Important Note:**  
Ensure the save path for JSON files and training files is consistent:
```python
input_dict["save_path"] = f"/home/rrichi/TransformerPolars/save/model/LossVars_{target_type_str}/GENIEv3-0-6-Honda-Truth-hA-LFG_{os.path.splitext(csv_file)[0]}_{acronym}"
```

### Running the Training File
Run the training file using the following command (or automate it using a bash script):
```bash
python3 new_train.py "$json_file" "$csv_file"
```

---

## Step 4: Producing Pictures from `.npz` Files

### Define a Function to Extract Values from `.npz` Files
Use a dictionary or any other preferred method to capture true and predicted values:
```python
def automate(filepath):
    print("Contents of the npz file:")
    with np.load(filepath) as file:
        trueval = file['trueval']
        prediction = file['prediction']
        for key in file.keys():
            print(key)
        if (prediction < 0).all():
            print("Prediction is less than 0")

    n_val, dim = trueval.shape
    return trueval, prediction, n_val, dim
```

### Specify the Path for `.npz` Files
File paths are separated based on different target types:
```python
file1path = f"/home/rrichi/TransformerPolars/save/model/LossVars_Nu_Energy/GENIEv3-0-6-Honda-Truth-hA-LFG_{os.path.splitext(csv_file)[0]}_{acronym}/model/result.npz"
```

### Key Formatting for Variables
Separate the CSV file base name and loss function type for dictionary keys:
```python
add_name = f"{os.path.splitext(csv_file)[0]}_{acronym}"

# Example
variables[f"TrueNuEn_4Mom_{add_name}"] = trueval_4Mom[:, 0]
```

**Note:**  
For plotting loss functions across different target types, use the full key names when plotting:
```python
TrueNuEn_VectorLeptWithoutNC = variables[f"TrueNuEn_4Mom_AnyNu_Inclusive_Thresh_p1to1_VectorLeptWithoutNC_eventnum_All_NpNpi_MAE"]
```

### Running the File
To run the file, use:
```bash
python3 compare_main.py
```
```
