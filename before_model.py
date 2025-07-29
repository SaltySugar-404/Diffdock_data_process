import math
import os.path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import SDWriter

from init import *


def split_sdf(task_id: str, split_keyword: str):
    unsplit_sdf_file = os.path.join(task_id, "input_one_protein_here")
    split_sdf_dir = os.path.join(task_id, "input_split_ligands_here")

    supplier = Chem.SDMolSupplier(unsplit_sdf_file, removeHs=False)
    for idx, mol in enumerate(supplier):
        mol_id = mol.GetProp(split_keyword) if mol.HasProp(split_keyword) else f"mol_{idx}"
        writer = SDWriter(os.path.join(split_sdf_dir, f"{mol_id}.sdf"))
        writer.write(mol)
        writer.close()


def generate_model_input(task_id: str, num_chunks: int):
    protein_id = os.listdir(os.path.join(task_id, "input_one_protein_here"))[0].split(".")[0]
    all_ligand_ids = [file.split(".")[0] for file in
                      sorted(os.listdir(os.path.join(task_id, "input_split_ligands_here")))]

    all_complex_name, all_protein_path, all_ligand_description, all_protein_sequence = [], [], [], []
    for ligand_id in all_ligand_ids:
        all_complex_name.append(f"{protein_id}_{ligand_id}")
        all_protein_path.append(
            os.path.join(CODE_ROOT_DIR_NAME, task_id, "input_one_protein_here", protein_id + ".pdb"))
        all_ligand_description.append(
            os.path.join(CODE_ROOT_DIR_NAME, task_id, "input_split_ligands_here", ligand_id + ".sdf"))
    all_model_inputs = pd.DataFrame({"complex_name": all_complex_name, "protein_path": all_protein_path,
                                     "ligand_description": all_ligand_description,
                                     "protein_sequence": all_protein_sequence})

    total = len(all_model_inputs.index)
    chunk_size = math.ceil(total / num_chunks)
    split_model_inputs = [all_model_inputs[i:i + chunk_size] for i in range(0, total, chunk_size)]
    model_inputs_dir = os.path.join(task_id, "all_model_inputs")
    for i, one_model_input in enumerate(split_model_inputs):
        one_model_input.to_csv(os.path.join(model_inputs_dir, f"chunk_{i}.csv"), index=False)


def get_shell(task_id: str):
    task_dir = os.path.join(CODE_ROOT_DIR_NAME, task_id)
    all_commands = []
    for index in range(len(os.listdir(os.path.join(task_id, "all_model_inputs")))):
        current_command = ""
        current_command += f"CUDA_VISIBLE_DEVICES={index} nohup python -m inference "
        current_command += f"--config default_inference_args.yaml "
        current_command += f"--protein_ligand_csv {os.path.join(task_dir, 'all_model_inputs', f'chunk_{index}.csv')} "
        current_command += f"--out_dir {os.path.join(task_dir, 'all_model_inputs')} "
        current_command += f"{os.path.join(task_dir, 'logs')}/chunk_{index}_log.out 2>&1 &"
        all_commands.append(current_command)
    return all_commands
