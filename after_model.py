import os.path
import re
from typing import List

import pandas as pd

from init import *


def get_rank(file_name):
    match = re.search(r'rank(\d+)', file_name)
    return int(match.group(1)) if match else float('inf')


def get_confidence(file_name: str):
    return float(file_name.split('confidence')[-1].split('.sd')[0])


def get_confidences(one_predict_result_dir: str) -> List[float]:
    predicted_names = os.listdir(one_predict_result_dir)
    confidence_files = [f for f in predicted_names if 'confidence' in f]
    sorted_confidence_files = sorted(confidence_files, key=get_rank)
    confidences = [get_confidence(file_name) for file_name in sorted_confidence_files]
    return confidences


def get_save_results(protein_id: str, all_results_dir: str, num_predicts: int, output_file: str):
    proteins, ligands = [], []
    rank_confidences = [[] for _ in range(num_predicts)]
    all_dirs = [os.path.join(all_results_dir, dir_name) for dir_name in os.listdir(all_results_dir)]
    for dir in all_dirs:
        proteins.append(protein_id)
        ligand_id = dir.split('_')[-1]
        ligands.append(ligand_id)
        if not len(os.listdir(dir)):
            print(f"Invalid dir: {dir}")
            for rank_confidence in rank_confidences:
                rank_confidence.append(-1000)
        else:
            confidences = get_confidences(dir)
            for index in range(len(rank_confidences)):
                rank_confidences[index].append(confidences[index])

    data_dict = {"proteins": proteins, "ligands": ligands}
    for index in range(len(rank_confidences)):
        data_dict[f"rank_{index + 1}_confidence"] = rank_confidences[index]
    all_results_file = pd.DataFrame(data_dict).sort_values(by='rank_1_confidence', ascending=False)
    all_results_file.to_csv(os.path.join(output_file), index=False)
