import os.path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import SDWriter

from init import *


def split_sdf(input_dir: str, output_dir: str):
    unsplit_sdf_file = os.path.join(input_dir, os.listdir(input_dir)[0])
    supplier = Chem.SDMolSupplier(unsplit_sdf_file, removeHs=False)

    for idx, mol in enumerate(supplier):
        if mol is None:
            print(f"Warning: molecule {idx} is None. Skipping.")
            continue

        try:
            mol_id = mol.GetProp(SPLIT_KEYWORD) if mol.HasProp(SPLIT_KEYWORD) else f"mol_{idx}"
            output_path = os.path.join(output_dir, f"{mol_id}.sdf")

            if os.path.exists(output_path):
                print(f"File {output_path} already exists. Skipping.")
                continue

            writer = SDWriter(output_path)
            writer.write(mol)
            writer.close()
        except Exception as e:
            print(f"Error processing molecule {idx}: {e}")
            continue


def check_split(split_ligands_dir: str, remove: bool = False) -> bool:
    all_split_ligand_files = [os.path.join(split_ligands_dir, f) for f in os.listdir(split_ligands_dir)]
    all_split_ligand_files = [f for f in all_split_ligand_files if f.endswith('.sdf')]
    if len(all_split_ligand_files) == 1:
        ligand_file = all_split_ligand_files[0]
        with open(ligand_file, 'r', encoding='utf-8') as f:
            count_delimiters = 0
            line_num = 0
            for line in f:
                line_num += 1
                line = line.strip()
                if line == '$$$$':
                    count_delimiters += 1
                    if count_delimiters > 1:
                        if remove:
                            os.remove(ligand_file)
                        return False
    return True


def generate_model_input(protein_file: str, split_ligands_dir: str, model_inputs_dir: str):
    protein_id = os.path.basename(protein_file).split(".")[0]
    all_ligand_files = [os.path.join(split_ligands_dir, file) for file in sorted(os.listdir(split_ligands_dir))]
    all_ligand_ids = [file.split(".")[0] for file in sorted(os.listdir(split_ligands_dir))]

    all_complex_name, all_protein_path, all_ligand_description = [], [], []
    for index in range(len(all_ligand_ids)):
        all_complex_name.append(f"{protein_id}_{all_ligand_ids[index]}")
        all_protein_path.append(os.path.join(BASE_DIR_NAME, protein_file))
        all_ligand_description.append(os.path.join(BASE_DIR_NAME, all_ligand_files[index]))

    all_model_inputs = pd.DataFrame({
        "complex_name": all_complex_name,
        "protein_path": all_protein_path,
        "ligand_description": all_ligand_description
    })
    all_model_inputs["protein_sequence"] = pd.NA

    total = len(all_model_inputs)
    split_model_inputs = [all_model_inputs[i:i + CHUNK_SIZE] for i in range(0, total, CHUNK_SIZE)]

    for i, one_model_input in enumerate(split_model_inputs):
        one_model_input.to_csv(os.path.join(model_inputs_dir, f"chunk_{i:03d}.csv"), index=False)


def get_save_commands(model_inputs_dir: str, model_outputs_dir: str, logs_dir: str, sh_dir):
    all_files = sorted([f for f in os.listdir(model_inputs_dir) if f.startswith("chunk_") and f.endswith(".csv")])
    num_gpus = len(GPU_INDEXES)
    all_commands = []

    for index, file_name in enumerate(all_files):
        gpu_index = GPU_INDEXES[index % num_gpus]
        run_index = index // num_gpus
        current_command = ""
        current_command += f"CUDA_VISIBLE_DEVICES={gpu_index} "
        current_command += f"nohup python -m inference --config default_inference_args.yaml "
        current_command += f"--protein_ligand_csv {os.path.join(BASE_DIR_NAME, model_inputs_dir, file_name)} "
        current_command += f"--out_dir {os.path.join(BASE_DIR_NAME, model_outputs_dir)} "
        current_command += f"> {os.path.join(BASE_DIR_NAME, logs_dir, f'run_{run_index:03d}_gpu_{gpu_index}.out')} 2>&1 & "
        current_command += f"echo $! > {os.path.join(BASE_DIR_NAME, sh_dir, f'run_{run_index:03d}_gpu_{gpu_index}.pid')}"
        all_commands.append(current_command)

        for i in range(0, len(all_commands), num_gpus):
            sh_path = os.path.join(sh_dir, f"run_{(i // num_gpus):03d}.sh")
        with open(sh_path, "w") as f:
            f.write("#!/bin/bash\n\n")
            batch_cmds = all_commands[i:i + num_gpus]
            for cmd in batch_cmds:
                f.write(cmd + "\n")


def delete_pid_file(pid_dir: str):
    for filename in os.listdir(pid_dir):
        if filename.endswith('.pid'):
            file_path = os.path.join(pid_dir, filename)
            if os.path.isfile(file_path):
                os.remove(file_path)
