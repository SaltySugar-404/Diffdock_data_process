import math
import os.path
from typing import List

import pandas as pd
from rdkit import Chem
from rdkit.Chem import SDWriter

from init import *


# 生成蛋白质索引
def generate_protein_index():
    if os.path.exists(PROTEIN_INDEX_FILE):
        print(f"Task name: {CURRENT_TASK_DIR} protein_index.csv already exist.")
        return

    all_files = [os.path.join(PROTEIN_INPUT_DIR, file_name) for file_name in sorted(os.listdir(PROTEIN_INPUT_DIR))]
    all_id, all_sequence, all_file_path = [], [], []
    for file in all_files:
        all_id.append(os.path.splitext(os.path.basename(file))[0])
        all_file_path.append(os.path.join(CODE_ROOT_DIR_NAME, file))

    data = {"id": all_id, "file_path": all_file_path}
    df = pd.DataFrame(data)
    df.to_csv(PROTEIN_INDEX_FILE, index=False)
    print(f"Task name: {CURRENT_TASK_DIR} ligand_index.csv generated")


# 分离整块sdf
def split_sdf(sdf_file: str, id_class: str):
    supplier = Chem.SDMolSupplier(sdf_file, removeHs=False)
    for idx, mol in enumerate(supplier):
        mol_id = mol.GetProp(id_class) if mol.HasProp(id_class) else f"mol_{idx}"
        writer = SDWriter(os.path.join(SPLIT_LIGAND_DIR, f"{mol_id}.sdf"))
        writer.write(mol)
        writer.close()


#
def generate_ligand_index():
    if os.path.exists(LIGAND_INDEX_FILE):
        print(f"Task name: {CURRENT_TASK_DIR} ligand_index.csv already exist.")
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
    print(f"Task name: {CURRENT_TASK_DIR} protein_index.csv generated")


def split_dataframe(df: pd.DataFrame, batch_size: int) -> List[pd.DataFrame]:
    return [df[i:i + batch_size] for i in range(0, len(df), batch_size)]


def generate_model_input(num_chunks: int, selected_protein_id: str,
                         selected_ligand_ids: List[str] = None):
    ligand_index = pd.read_csv(LIGAND_INDEX_FILE, index_col="id")
    if not selected_ligand_ids:
        selected_ligand_ids = sorted(ligand_index.index.tolist())
    protein_index = pd.read_csv(PROTEIN_INDEX_FILE, index_col="id")

    model_input = pd.DataFrame(columns=["complex_name", "protein_path", "ligand_description","protein_sequence"])
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
        print(f"Task name: {CURRENT_TASK_DIR} {MODEL_INPUTS_DIR}/chunk{i}.csv generated")


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
    # generate_all_model_input(num_chunks=NUM_GPUS)
