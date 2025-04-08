## This is Richi's
#This has WANDB in it, YOU CAN RUN THIS BY CALLING input_trainingfile.py "json_file" "csv_file" on the terminal
import json
import os
import argparse
from transformer_ee.train import MVtrainer
from transformer_ee.logger.wandb_train_logger import WandBLogger
import wandb
import shutil
import os

# Initialize WandB Logger


# Your code logic here...



# 
##Directory of the CSV and json files
csv_directory = "/exp/dune/app/users/rrichi/FinalCodes"
json_directory="/home/rrichi/TransformerPolars/transformer_EE/transformer_ee/config"

##Acronyms for the Loss functions to be used later on in the save path
acronyms = ["MSE", "MAE", "MAPE","MCE"]

# wandb_dir_path = "/exp/dune/app/users/rrichi/wandb"  # Change this if the path is different

def training(csv_file, input_dict, acronym, target_type_str):
    input_dict["num_workers"] = 10
    input_dict["data_path"] = os.path.join(csv_directory, csv_file)
    input_dict["model"]["kwargs"]["nhead"] = 2
    input_dict["model"]["epochs"] = 100
    input_dict["model"]["kwargs"]["num_layers"] = 5
    input_dict["optimizer"]["name"] = "sgd" #what?
    input_dict["optimizer"]["kwargs"]["lr"] = 0.01
    input_dict["optimizer"]["kwargs"]["momentum"] = 0.9
    input_dict["dataframe_type"] = "polars"
    #demo save path: home/rrichi/MLProject/save/model/LossVars_"Nu_Energy"/GENIEv3-0-6-Honda-Truth-hA-LFG_"AnyNu_CC_Thresh_p1to1_eventnum_All_NpNpi_MAE"(The name inside the quotes is changeable)
    input_dict["save_path"] =f"/home/rrichi/TransformerPolars/save/model/LossVars_{target_type_str}/GENIEv3-0-6-Honda-Truth-hA-LFG_{os.path.splitext(csv_file)[0]}_{acronym}"

    # =f"/exp/dune/app/users/rrichi/TransformerPolars/save/model/LossVars_{target_type_str}/GENIEv3-0-6-Honda-Truth-hA-LFG_{os.path.splitext(csv_file)[0]}_{acronym}"

    # 
    return input_dict
    

def main(json_file, csv_file):
    #splits the json file for using it in the target_string
    base_jsonfile = os.path.splitext(json_file)[0]
    base_csvfile=os.path.splitext(csv_file)[0]
    json_filename=os.path.join(json_directory,json_file)
    for acronym in acronyms:
        if acronym in base_jsonfile:
            with open(json_filename, encoding="UTF-8", mode="r") as f:
                input_d = json.load(f)
                #target_string separates the target variables from the json name to enter as an input in the training function for using in the save path
                target_string = base_jsonfile.split('_' + acronym)[0].split('pi_', 1)[-1]
                train_dict = training(csv_file, input_d, acronym, target_string)
                id_wandb = base_jsonfile.replace('input_GENIEv3-0-6-Honda-Truth-hA-LFG_', '') \
                                             .replace('_Thresh', '') \
                                             .replace('_eventnum', '')
                my_logger = WandBLogger(
    project="GENIE-Atmo", entity="neutrinoenenergyestimators", config=input_d, dir="/exp/dune/app/users/rrichi", id=id_wandb )
         
                my_trainer = MVtrainer(input_d, logger=my_logger)
                # my_trainer = MVtrainer(train_dict)
                my_trainer.train()
                my_trainer.eval()


# wandb.finish()

# # Remove the WandB directory safely
# wandb_dir = "/exp/dune/app/users/rrichi/wandb"
# if os.path.exists(wandb_dir):
#     shutil.rmtree(wandb_dir)

##Takes two arguments, 1st input json file and then csv file
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Train a model using specified JSON and CSV files.")
    parser.add_argument("json_file", type=str, help="The path to the JSON file.")
    parser.add_argument("csv_file", type=str, help="The path to the CSV file.")
    args = parser.parse_args()

    main(args.json_file, args.csv_file)
