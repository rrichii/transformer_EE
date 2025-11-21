import os
import json
import itertools

# Directories
csv_directory = "/exp/dune/data/users/rrichi/MLProject/Training_Samples/Atmospherics_DUNE_Like/Flat_Spectra"
output_directory = "/exp/dune/data/users/rrichi/MLProject/Training_Samples/Atmospherics_DUNE_Like/Flat_Spectra/Jsonfiles/VectorFiles"

# Vectors & Scalars
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

scalar_type = ["tot_fKE", "p_tot", "P_miss", "tot_hadronic_energy"]

# All combinations (exact ordering preserved)
combinations = [
    (["Nu_Energy"],["Topology"]),                           
    (["Nu_Mom_X"],["Topology"]),
    (["Nu_Mom_Y"],["Topology"]),
    (["Nu_Mom_Z"],["Topology"]),
    (["Nu_Theta"],["Topology"]),
    (["Nu_Baseline"],["Topology"]),
    (["Nu_Phi"],["Topology"]),
    (["Nu_Energy"],["Nu_Mom_X", "Nu_Mom_Y", "Nu_Mom_Z"],["Topology"]),
    (["Nu_Energy"],["Nu_Theta"],["Topology"]),
    (["Nu_Theta"],["Nu_Mom_X", "Nu_Mom_Y", "Nu_Mom_Z"],["Topology"]),
    (["Nu_Energy"],["Nu_Theta"],["Nu_Phi"],["Topology"]),
    (["Nu_Energy"],["Nu_Baseline"],["Topology"])
]

# Loss functions
loss_functions = ["mean squared error", "mean absolute error", "mean absolute percentage error", "mean cubic error"]
acronyms = ["MSE", "MAE", "MAPE", "MCE"]
loss_map = dict(zip(loss_functions, acronyms))

# Safe short name
def short_name(v):
    parts = v.split("_", 1)
    return parts[1] if len(parts) > 1 else v

# Build config
def generate_configuration(csv_file, target_vars, losses):
    base_loss_names = losses.copy()

    # Topology always MAE
    if "Topology" in target_vars:
        topo_index = target_vars.index("Topology")
        base_loss_names[topo_index] = "mean absolute error"

    # Build tag
    tag_parts = [short_name(v) + "_" + loss_map[base_loss_names[i]] 
                 for i, v in enumerate(target_vars)]
    combo_tag = "_".join(tag_parts)

    # Coefficients (Topology = 0)
    coefficients = [0.5] * len(target_vars)
    if "Topology" in target_vars:
        coefficients[target_vars.index("Topology")] = 0

    config = {
        "data_path": os.path.join(csv_directory, csv_file),
        "num_workers": 10,
        "dataframe_type": "polars",
        "vector": vector_type,
        "scalar": scalar_type,
        "target": target_vars,
        "max_num_prongs": 35,
        "batch_size_train": 1024,
        "batch_size_valid": 256,
        "batch_size_test": 3000,
        "test_size": 0.2,
        "valid_size": 0.05,
        "seed": 0,
        "loss": {
            "kwargs": {
                "coefficients": coefficients,
                "base_loss_names": base_loss_names
            }
        },
        "optimizer": {"name": "AdamW", "kwargs": {"lr": 0.001}},
        "model": {"name": "Transformer_EE_MV", "kwargs": {}},
        "save_path":
            f"/exp/dune/data/users/rrichi/MLProject/Training_Samples/"
            f"Atmospherics_DUNE_Like/Flat_Spectra/LossVars_{combo_tag}/"
            f"{os.path.splitext(csv_file)[0]}_{combo_tag}",
        "model_phys_name":
            f"{os.path.splitext(csv_file)[0]}_{combo_tag}"
    }

    return config, combo_tag

# Flatten combination lists
def flatten_combo(combo):
    target_vars = []
    for item in combo:
        if isinstance(item, list):
            target_vars.extend(item)
        else:
            target_vars.append(item)
    return target_vars

# Helpers
def is_momentum(v):
    return v in ["Nu_Mom_X","Nu_Mom_Y","Nu_Mom_Z"]

def is_topology(v):
    return v == "Topology"

# Get CSV files
csv_files = [f for f in os.listdir(csv_directory) if f.endswith(".csv") and "VectorLeptwNC" in f]

# -------- JSON Generation -------- #
for csv_file in csv_files:
    for combo in combinations:

        target_vars = flatten_combo(combo)

        # Identify variable groups
        momentum_vars = [v for v in target_vars if is_momentum(v)]
        other_vars = [v for v in target_vars if (not is_momentum(v) and not is_topology(v))]
        topo_present = ("Topology" in target_vars)

        # Build loss assignment loops
        mom_loss_choices = loss_functions if momentum_vars else [None]

        other_var_loss_choices = [loss_functions] * len(other_vars)

        # Cartesian product
        for mom_loss in mom_loss_choices:
            for other_losses in itertools.product(*other_var_loss_choices):

                # Build full loss list in order
                losses = []
                idx_other = 0

                for v in target_vars:

                    if is_momentum(v):
                        losses.append(mom_loss)      # same for all momentum components

                    elif is_topology(v):
                        losses.append("mean absolute error")   # ALWAYS MAE

                    else:
                        losses.append(other_losses[idx_other])
                        idx_other += 1

                # Now build config
                config, combo_tag = generate_configuration(csv_file, target_vars, losses)

                # Save JSON
                json_filename = f"{os.path.splitext(csv_file)[0]}_{combo_tag}.json"
                output_path = os.path.join(output_directory, json_filename)

                with open(output_path, "w") as jf:
                    json.dump(config, jf, indent=4)

                print("Saved:", output_path)
