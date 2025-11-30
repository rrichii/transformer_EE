#!/bin/bash

readarray -t Json_files < BeamDuneOnAxisND_NatJson.txt
readarray -t CSV_files < CSV.txt

Json_files=( "${Json_files[@]:0:4}" )

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
            
            if [[ "$json_base" == *"_Lept_"* ]]; then
                short_base="$json_base"
                short_base="${short_base#Numu_CC_Train_DUNEBeam_Natural_OnAxisND_p1to10_ScalarLeptwNC_eventnum_All_NpNpi_}"

            else
              short_base="$json_base"
              short_base="${short_base#Numu_CC_Train_DUNEBeam_Natural_OnAxisND_}"

            fi

            # Make output directories
            out_dir="/exp/dune/data/users/${USER}/MLProject/Training_Samples/Beam_Like/Natural_Spectra/DUNEOnAxisND/${json_base}"
            mkdir -p "${out_dir}/wandb"

            # Run training
            python3 train_wide.py \
                --epochs 25 \
                --d-model 256 \
                --nhead 16 \
                --enable-noise \
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
                --wandb-project GENIE-Train-Beam-DUNE-Natural-2025 \
                --wandb-id "$short_base" \
                --dataframe-type polars
        fi
    done
done
