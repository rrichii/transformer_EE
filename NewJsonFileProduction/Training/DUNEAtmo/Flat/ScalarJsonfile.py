import os
import json
import itertools

# Directories
csv_directory = "/exp/dune/data/users/rrichi/MLProject/Training_Samples/Atmospherics_DUNE_Like/Flat_Spectra"
output_directory = "/exp/dune/data/users/rrichi/MLProject/Training_Samples/Atmospherics_DUNE_Like/Flat_Spectra/Jsonfiles/ScalarFiles"

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

scalar_types = [
    [
        "Lept_PDG",
        "Lept_Mass",
        "Lept_Energy",
        "Lept_MomX",
        "Lept_MomY",
        "Lept_MomZ",
        "Lept_CosTheta",
        "Lept_Theta",
        "Lept_Phi_z",
        "tot_fKE",
        "p_tot",
        "P_miss","tot_hadronic_energy"
    ],
    [  
        "Lept_PDG",
        "Lept_Mass",
        "Lept_Energy",
        "Lept_MomX",
        "Lept_MomY",
        "Lept_MomZ",
        "Lept_CosTheta",
        "Lept_Theta",
        "Lept_Phi_z",
        "tot_fKE",
        "p_tot",
        "P_miss","tot_hadronic_energy"
    ],
    [
        "Lept_PDG",
        "Lept_Mass",
        "Lept_Energy",
        "Lept_MomX",
        "Lept_MomY",
        "Lept_MomZ",
        "Lept_CosTheta",
        "Lept_Theta",
        "Lept_Phi_z",
        "tot_fKE",
        "p_tot",
        "P_miss","tot_hadronic_energy"
    ],
    [
        "tot_fKE",
        "p_tot",
        "P_miss",
        "tot_hadronic_energy"
    ],
    
    [#can calculate the magnitude of the momentum from lepton energy 
        "Lept_CosTheta",
        "Lept_Theta",
        "Lept_Phi_z",
        "tot_fKE",
        "p_tot",
        "P_miss","tot_hadronic_energy"
    ],
    [
        "Lept_PDG",
        "Lept_Mass",
        "Lept_Energy",
        "Lept_MomX",
        "Lept_MomY",
        "Lept_MomZ",
        "Lept_CosTheta",
        "Lept_Theta",
        "Lept_Phi_z",
        "tot_fKE",
        "p_tot",
        "P_miss","tot_hadronic_energy"

    ],
    [
        "Lept_PDG",
        "Lept_Mass",
        "Lept_Energy",
        "Lept_MomX",
        "Lept_MomY",
        "Lept_MomZ",
        "Lept_CosTheta",
        "Lept_Theta",
        "Lept_Phi_z",
        "tot_fKE",
        "p_tot","tot_hadronic_energy"
    ],
    [
        "Lept_PDG",
        "Lept_Mass",
        "Lept_Energy",
        "Lept_MomX",
        "Lept_MomY",
        "Lept_MomZ",
        "Lept_CosTheta",
        "Lept_Theta",
        "Lept_Phi_z",
        "tot_fKE",
        "p_tot",
        "P_miss","tot_hadronic_energy"

    ],
    [
        "Lept_PDG",
        "Lept_Mass",
        "Lept_Energy",
        "Lept_MomX",
        "Lept_MomY",
        "Lept_MomZ",
        "Lept_CosTheta",
        "Lept_Theta",
        "Lept_Phi_z",
        "tot_fKE",
        "p_tot",
        "P_miss","tot_hadronic_energy"
    ],
    [
        "Lept_PDG",
        "Lept_Mass",
        "Lept_Energy",
        "Lept_MomX",
        "Lept_MomY",
        "Lept_MomZ",
        "Lept_CosTheta",
        "Lept_Theta",
        "Lept_Phi_z",
        "tot_fKE",
        "p_tot",
        "P_miss","tot_hadronic_energy"

    ],
    [
        "Lept_PDG",
        "Lept_Mass",
        "Lept_Energy",
        "Lept_MomX",
        "Lept_MomY",
        "Lept_MomZ",
        "Lept_CosTheta",
        "Lept_Theta",
        "Lept_Phi_z",
        "tot_fKE",
        "p_tot",
        "P_miss","tot_hadronic_energy"

    ],
    [
        "Lept_PDG",
        "Lept_Mass",
        "Lept_Energy",
        "Lept_MomX",
        "Lept_MomY",
        "Lept_MomZ",
        "Lept_CosTheta",
        "Lept_Theta",
        "Lept_Phi_z",
        "tot_fKE",
        "p_tot",
        "P_miss","tot_hadronic_energy"

    ],
    [
        "Lept_PDG",
        "Lept_Mass",
        "Lept_Energy",
        "Lept_MomX",
        "Lept_MomY",
        "Lept_MomZ",
        "Lept_CosTheta",
        "Lept_Theta",
        "Lept_Phi_z",
        "tot_fKE",
        "p_tot",
        "P_miss","tot_hadronic_energy"

    ],
    [   
        "Lept_PDG",
        "Lept_Mass",
        "Lept_Energy",
        "Lept_MomX",
        "Lept_MomY",
        "Lept_MomZ",
        "Lept_CosTheta",
        "Lept_Theta",
        "Lept_Phi_z",
        "tot_fKE",
        "p_tot",
        "P_miss","tot_hadronic_energy"
        
    ]

]

# All combinations (exact ordering preserved)
target_types = [["Nu_Theta","Nu_CosTheta","Topology"],
                ["Nu_Energy","Topology"],
                ["Nu_Energy","Nu_Mom_X","Nu_Mom_Y","Nu_Mom_Z","Topology"],
                ["Nu_Energy","Nu_Mom_X","Nu_Mom_Y","Nu_Mom_Z","Lept_Energy","Lept_MomX","Lept_MomY","Lept_MomZ","Topology"],
                ["Nu_Energy","Lept_Energy","Nu_Baseline","Topology"],
                ["Nu_Baseline"],["Nu_Energy","Nu_Mom_X","Nu_Mom_Y","Nu_Mom_Z","P_miss","Topology"],
                ["Nu_Energy","Nu_Mom_X","Nu_Mom_Y","Nu_Mom_Z","Nu_Baseline","Topology"],
                ["Nu_Energy","Nu_Baseline"],
                [
                    "Nu_Energy",
                    "Nu_Theta","Topology"
            
                ],
                [
                    "Nu_Theta","Topology"
                ],
                [
                    "Nu_CosTheta","Topology"
                ],
                [
                    "Nu_Mom_X",
                    "Nu_Mom_Y",
                    "Nu_Mom_Z","Topology"
                ],
                [
                    "Nu_Baseline","Topology"
                ]]

loss_functions = [
    "mean squared error",
    "mean absolute error",
    "mean absolute percentage error",
    "mean cubic error"
]

loss_map = {
    "mean squared error":       "MSE",
    "mean absolute error":      "MAE",
    "mean absolute percentage error": "MAPE",
    "mean cubic error":         "MCE"
}

# ---------------------------------- #
# HELPERS
# ---------------------------------- #
# ----------------------------------------
def group_type(v):
    if "Energy" in v: return "energy"
    if "Mom"    in v: return "momentum"
    if any(k in v for k in ["Theta","CosTheta","Phi"]): return "angle"
    if "Baseline" in v: return "baseline"
    if v == "Topology": return "topology"
    if "P_miss" in v: return "pmiss"
    return "other"

def species_prefix(v):
    if v.startswith("Nu_"): return "Nu"
    if v.startswith("Lept_"): return "Lept"
    return ""

def select_scalar_list_for_targets(index):
    return scalar_types[index]

def build_group_tag_preserve_order(target_vars, losses):
    parts = []
    used_groups = set()   # ensure each group is added ONCE

    for i, v in enumerate(target_vars):
        g = group_type(v)
        sp = species_prefix(v)
        loss_code = loss_map[losses[i]]

        # ---- ENERGY GROUPS (Nu_Energy, Lept_Energy) ----
        if g == "energy":
            key = f"{sp}_Energy"
            if key not in used_groups:
                parts.append(f"{key}_{loss_code}")
                used_groups.add(key)

        # ---- MOMENTUM GROUPS (Nu_Mom, Lept_Mom) ----
        elif g == "momentum":
            key = f"{sp}_Mom"
            if key not in used_groups:
                parts.append(f"{key}_{loss_code}")
                used_groups.add(key)

        # ---- ANGLES KEEP INDIVIDUAL TAGS ----
        elif g == "angle":
            parts.append(f"{v}_{loss_code}")

        # ---- BASELINE SINGLE TAG ----
        elif g == "baseline":
            if "Baseline" not in used_groups:
                parts.append(f"{v}_{loss_code}")
                used_groups.add("Baseline")

        # ---- P_miss SINGLE TAG ----
        elif g == "pmiss":
            if "P_miss" not in used_groups:
                parts.append(f"P_miss_{loss_code}")
                used_groups.add("P_miss")

        # ---- TOPOLOGY ALWAYS ONCE ----
        elif g == "topology":
            if "Topology" not in used_groups:
                parts.append("Topology_MAE")
                used_groups.add("Topology")

    return "_".join(parts)


def generate_configuration(csv_file, target_vars, losses, scalar_list):

    losses = list(losses)
    if "Topology" in target_vars:
        losses[target_vars.index("Topology")] = "mean absolute error"

    combo_tag = build_group_tag_preserve_order(target_vars, losses)

    coefficients = [0.5] * len(target_vars)
    if "Topology" in target_vars:
        coefficients[target_vars.index("Topology")] = 0

    config = {
        "data_path": os.path.join(csv_directory, csv_file),
        "vector": vector_type,
        "scalar": scalar_list,
        "target": target_vars,
        "loss": {
            "kwargs": {
                "coefficients": coefficients,
                "base_loss_names": losses
            }
        },
        "model_phys_name": f"{csvbase}_{combo_tag}",
        "save_path": f"{output_directory}/LossVars_{combo_tag}/{combo_tag}",
    }
    return config, combo_tag

# ----------------------------------------
# JSON GENERATION
# ----------------------------------------

csv_files = [
    f for f in os.listdir(csv_directory)
    if f.endswith(".csv") and "ScalarLeptwNC" in f
]

    
for csv_file in csv_files:

    csvbase = os.path.splitext(csv_file)[0]

    for idx, target_vars in enumerate(target_types):

        scalar_list = select_scalar_list_for_targets(idx)

        groups_present = {g: any(group_type(v)==g for v in target_vars)
                          for g in ["energy","momentum","angle","baseline","pmiss"]}

        energy_choices   = loss_functions if groups_present["energy"]   else ["mean absolute error"]
        momentum_choices = loss_functions if groups_present["momentum"] else ["mean absolute error"]
        angle_choices    = loss_functions if groups_present["angle"]    else ["mean absolute error"]
        baseline_choices = loss_functions if groups_present["baseline"] else ["mean absolute error"]
        pmiss_choices    = loss_functions if groups_present["pmiss"]    else ["mean absolute error"]

        for L_e, L_m, L_a, L_b, L_p in itertools.product(
            energy_choices, momentum_choices,
            angle_choices, baseline_choices, pmiss_choices
        ):

            losses = []
            for v in target_vars:
                g = group_type(v)
                if g == "energy":   losses.append(L_e)
                elif g == "momentum": losses.append(L_m)
                elif g == "angle":    losses.append(L_a)
                elif g == "baseline": losses.append(L_b)
                elif g == "pmiss":    losses.append(L_p)
                elif g == "topology": losses.append("mean absolute error")
                else: losses.append("mean squared error")

            config, combo_tag = generate_configuration(csv_file, target_vars, losses, scalar_list)

            json_filename = f"{csvbase}_{combo_tag}.json"
            output_path = os.path.join(output_directory, json_filename)

            with open(output_path, "w") as jf:
                json.dump(config, jf, indent=4)

            print("Saved:", json_filename)

