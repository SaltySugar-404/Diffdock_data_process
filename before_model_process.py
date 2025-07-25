import argparse
import datetime
import json
import math
import os.path
import shutil
import subprocess
from typing import List

import pandas as pd
import torch
from rdkit import Chem
from rdkit.Chem import SDWriter

from init import *

parser = argparse.ArgumentParser(description="Run this before model process.")
parser.add_argument('--sdf_split_keyword', type=str, default='DATABASE_ID', help='Keyword to split SDF files.')
parser.add_argument('--num_gpus', type=int, default=torch.cuda.device_count(), help='Number of GPUs to use.')
args = parser.parse_args()

SDF_SPLIT_KEYWORD = args.sdf_split_keyword
NUM_GPUS = args.num_gpus

# 当前任务文件夹，默认是蛋白质名
CURRENT_PROTEIN_NAME = os.listdir(OUTSIDE_PROTEIN_INPUT_DIR)[0].split(".")[0]
current_task_name = ""
all_task_name_and_index = os.listdir(ALL_TASKS)
max_task_index = -1
for task_name_and_index in all_task_name_and_index:
    if len(task_name_and_index.split(".")) == 2:
        continue
    task_name = task_name_and_index.split("_index")[0]
    task_index = int(task_name_and_index.split("_index_")[1])
    if task_name == CURRENT_PROTEIN_NAME and task_index > max_task_index:
        max_task_index = task_index
current_task_name = CURRENT_PROTEIN_NAME + f"_index_{max_task_index + 1}"
CURRENT_TASK_DIR = os.path.join(ALL_TASKS, current_task_name)
os.makedirs(CURRENT_TASK_DIR, exist_ok=True)

# 蛋白质输入文件夹
shutil.move(OUTSIDE_PROTEIN_INPUT_DIR, CURRENT_TASK_DIR)
os.makedirs(OUTSIDE_PROTEIN_INPUT_DIR, exist_ok=True)
PROTEIN_INPUT_DIR = os.path.join(CURRENT_TASK_DIR, OUTSIDE_PROTEIN_INPUT_DIR)

# 未分离配体文件夹
shutil.move(OUTSIDE_UNSPLIT_LIGAND_INPUT_DIR, CURRENT_TASK_DIR)
os.makedirs(OUTSIDE_UNSPLIT_LIGAND_INPUT_DIR, exist_ok=True)
UNSPLIT_LIGAND_DIR = os.path.join(CURRENT_TASK_DIR, OUTSIDE_UNSPLIT_LIGAND_INPUT_DIR)

# 已分离配体文件夹
shutil.move(OUTSIDE_SPLIT_LIGAND_INPUT_DIR, CURRENT_TASK_DIR)
os.makedirs(OUTSIDE_SPLIT_LIGAND_INPUT_DIR, exist_ok=True)
SPLIT_LIGAND_DIR = os.path.join(CURRENT_TASK_DIR, OUTSIDE_SPLIT_LIGAND_INPUT_DIR)

# 蛋白质索引csv
PROTEIN_INDEX_FILE = os.path.join(CURRENT_TASK_DIR, r"protein_index.csv")

# 配体索引csv
LIGAND_INDEX_FILE = os.path.join(CURRENT_TASK_DIR, r"ligand_index.csv")

# 模型输入文件夹
MODEL_INPUTS_DIR = os.path.join(CURRENT_TASK_DIR, r"all_model_inputs")
os.makedirs(MODEL_INPUTS_DIR, exist_ok=True)

# 模型日志文件夹
MODEL_LOGS_DIR = os.path.join(CURRENT_TASK_DIR, r"logs")
os.makedirs(MODEL_LOGS_DIR, exist_ok=True)

# 所有的模型输出
MODEL_OUTPUTS_DIR = os.path.join(CURRENT_TASK_DIR, "all_model_outputs")
os.makedirs(MODEL_OUTPUTS_DIR, exist_ok=True)

# 模型输出分析
MODEL_OUTPUTS_ANALYZE_FILE = os.path.join(CURRENT_TASK_DIR, "output_analyze.csv")


# 生成蛋白质索引
def generate_protein_index():
    if os.path.exists(PROTEIN_INDEX_FILE):
        print(f"Task name: {current_task_name} protein_index.csv already exist.")
        return

    all_files = [os.path.join(PROTEIN_INPUT_DIR, file_name) for file_name in sorted(os.listdir(PROTEIN_INPUT_DIR))]
    all_id, all_sequence, all_file_path = [], [], []
    for file in all_files:
        all_id.append(os.path.splitext(os.path.basename(file))[0])
        all_file_path.append(os.path.join(CODE_ROOT_DIR_NAME, file))

    data = {"id": all_id, "file_path": all_file_path}
    df = pd.DataFrame(data)
    df.to_csv(PROTEIN_INDEX_FILE, index=False)
    print(f"Task name: {current_task_name} ligand_index.csv generated")


# 分离整块sdf
def split_sdf(sdf_file: str, id_class: str):
    supplier = Chem.SDMolSupplier(sdf_file, removeHs=False)
    for idx, mol in enumerate(supplier):
        mol_id = mol.GetProp(id_class) if mol.HasProp(id_class) else f"mol_{idx}"
        writer = SDWriter(os.path.join(SPLIT_LIGAND_DIR, f"{mol_id}.sdf"))
        writer.write(mol)
        writer.close()


def generate_ligand_index():
    if os.path.exists(LIGAND_INDEX_FILE):
        print(f"Task name: {current_task_name} ligand_index.csv already exist.")
        return

    all_files = [os.path.join(SPLIT_LIGAND_DIR, file_name) for file_name in
                 sorted(os.listdir(SPLIT_LIGAND_DIR))]
    all_id, all_sequence, all_file_path = [], [], []
    for file in all_files:
        all_id.append(os.path.splitext(os.path.basename(file))[0])
        all_sequence.append(Chem.MolToSmiles(Chem.SDMolSupplier(file)[0]))
        all_file_path.append(os.path.join(CODE_ROOT_DIR_NAME, file))

    data = {"id": all_id, "sequence": all_sequence, "file_path": all_file_path}
    df = pd.DataFrame(data)
    df.to_csv(LIGAND_INDEX_FILE, index=False)
    print(f"Task name: {current_task_name} protein_index.csv generated")


def split_dataframe(df: pd.DataFrame, batch_size: int) -> List[pd.DataFrame]:
    return [df[i:i + batch_size] for i in range(0, len(df), batch_size)]


def generate_model_input(num_chunks: int, selected_protein_id: str,
                         selected_ligand_ids: List[str] = None):
    ligand_index = pd.read_csv(LIGAND_INDEX_FILE, index_col="id")
    if not selected_ligand_ids:
        selected_ligand_ids = sorted(ligand_index.index.tolist())
    protein_index = pd.read_csv(PROTEIN_INDEX_FILE, index_col="id")

    model_input = pd.DataFrame(columns=["complex_name", "protein_path", "ligand_description", "protein_sequence"])
    all_complex_name, all_protein_path, all_ligand_description = [], [], []

    for ligand_id in selected_ligand_ids:
        all_complex_name.append(f"{selected_protein_id}_{ligand_id}")
        all_protein_path.append(protein_index.loc[selected_protein_id]["file_path"])
        all_ligand_description.append(ligand_index.loc[ligand_id]["file_path"])

    model_input["complex_name"] = all_complex_name
    model_input["protein_path"] = all_protein_path
    model_input["ligand_description"] = all_ligand_description

    total = len(model_input.index)
    chunk_size = math.ceil(total / num_chunks)
    split_model_inputs = [model_input[i:i + chunk_size] for i in range(0, total, chunk_size)]

    for i, one_model_input in enumerate(split_model_inputs):
        one_model_input.to_csv(os.path.join(MODEL_INPUTS_DIR, f"chunk_{i}.csv"), index=False)
        print(f"Task name: {current_task_name} {MODEL_INPUTS_DIR}/chunk{i}.csv generated")


def generate_all_model_input(num_chunks: int, selected_protein_ids: List[str] = None,
                             selected_ligand_ids: List[str] = None):
    curren_protein_ids = [os.path.basename(file).split(".")[0] for file in sorted(os.listdir(PROTEIN_INPUT_DIR))]
    if selected_protein_ids:
        for protein_id in selected_protein_ids:
            if protein_id in curren_protein_ids:
                generate_model_input(num_chunks, protein_id, selected_ligand_ids)
    else:
        for protein_id in curren_protein_ids:
            generate_model_input(num_chunks, protein_id)


if __name__ == "__main__":
    if len(os.listdir(UNSPLIT_LIGAND_DIR)) == 1 and len(os.listdir(SPLIT_LIGAND_DIR)) == 0:
        unsplit_sdf_file = os.path.join(UNSPLIT_LIGAND_DIR, os.listdir(UNSPLIT_LIGAND_DIR)[0])
        split_sdf(sdf_file=unsplit_sdf_file, id_class=SDF_SPLIT_KEYWORD)
    generate_protein_index()
    generate_ligand_index()
    generate_model_input(num_chunks=NUM_GPUS, selected_protein_id=CURRENT_PROTEIN_NAME)

    task_config = {"create_time": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                   "num_gpus": NUM_GPUS, "sdf_split_keyword": SDF_SPLIT_KEYWORD}
    with open(os.path.join(CURRENT_TASK_DIR, "task_config.json"), "w") as json_file:
        json.dump(task_config, json_file)

    commands = []
    for index in range(NUM_GPUS):
        current_command = ""
        current_command += f"CUDA_VISIBLE_DEVICES={index} nohup python -m inference "
        current_command += f"--config default_inference_args.yaml "
        current_command += f"--protein_ligand_csv {os.path.join(MODEL_INPUTS_DIR, 'chunk_{index}.csv')} "
        current_command += f"--out_dir {MODEL_OUTPUTS_DIR} "
        current_command += f"{MODEL_LOGS_DIR}/chunk_{index}_log.out 2>&1 &"
        commands.append(current_command)

    for command in commands:
        print(command)
        subprocess.Popen(command, shell=True)
