#!/bin/bash

readarray -t Json_files < AtmDUNE_NatJson.txt
readarray -t CSV_files < CSV.txt

Json_files=( "${Json_files[@]:0:20}" )

# Iterate over json files
for json_file in "${Json_files[@]}"; do
    scalar_list=$(jq -r '.scalar[]' "$json_file")
    scalar_list_cli=$(echo $scalar_list)

    vector_list=$(jq -r '.vector[]' "$json_file")
    vector_list_cli=$(echo $vector_list)
    
    json_base=$(basename "$json_file" .json)

    for csv_file in "${CSV_files[@]}"; do

        csv_base=$(basename "$csv_file" .csv)

        if [[ "$json_base" == *"$csv_base"* ]]; then

            echo "Matching JSON and CSV files found:"
            echo "JSON: $json_file"
            echo "CSV: $csv_file"
            short_base="${json_base#Numu_CC_Train_}"

            # Make output directories
            out_dir="/exp/dune/data/users/${USER}/MLProject/Training_Samples/Atmospherics_DUNE_Like/Natural_Spectra/${json_base}"
            mkdir -p "${out_dir}/wandb"
            

            # Run training
            python3 train_wide.py \
                --epochs 40 \
                --d-model 256 \
                --nhead 8 \
                --noise-scalar $scalar_list_cli \
                --noise-vector $vector_list_cli \
                --num-layers 6 \
                --optimizer adamw \
                --lr 3e-4 \
                --weight-decay 0.01 \
                --dropout 0.3 \
                --base-config "${json_file}" \
                --data-path "${csv_file}" \
                --save-path "${out_dir}" \
                --wandb-dir "${out_dir}/wandb" \
                --wandb-project GENIE-Train-Atmo-DUNE-Natural-2025 \
                --wandb-id "$short_base" \
                --dataframe-type polars
        fi
    done
done
