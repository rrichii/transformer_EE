#!/bin/bash
###For MAE, MSE, MAPE, MCE for AnyNuCC_VectorLept_withNC
# Read the files into arrays

readarray -t Json_files < AnyNuCCVectorLeptwNC.txt
readarray -t CSV_files < CSV_AnyNuCC.txt

# Iterate over json files
for json_file in "${Json_files[@]}"; do
    # Extract the base name of the json file (without extension)
    json_base=$(basename "$json_file" .json)
    
    # Iterate over csv files
    for csv_file in "${CSV_files[@]}"; do
        # Extract the base name of the csv file (without extension)
        csv_base=$(basename "$csv_file" .csv)
        
        # Check if the base name of the csv file is found in the json file name
        if [[ "$json_base" == *"$csv_base"* ]]; then
            # Run the training script with the matching files
            # ./training.sh "$json_file" "$csv_file"
            echo "Matching JSON and CSV files found:"
            echo "JSON: $json_file"
            echo "CSV: $csv_file"
            python3 new_trainCopy.py "$json_file" "$csv_file" 
        fi
    done
done
