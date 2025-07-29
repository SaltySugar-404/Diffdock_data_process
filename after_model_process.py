import argparse
import os.path
import re
from typing import List

import pandas as pd

from init import *

parser = argparse.ArgumentParser(description="Run this after model process.")
parser.add_argument('--task_name', type=str, required=True, help='Select task name in dir tasks.')
args = parser.parse_args()

current_task_name = args.task_name
if current_task_name not in os.listdir(ALL_TASKS_DIR):
    print(r"Task name {current_task_name} not found.")
    exit(1)


def get_rank(file_name):
    match = re.search(r'rank(\d+)', file_name)
    return int(match.group(1)) if match else float('inf')


def get_confidence(file_name: str):
    return float(file_name.split('confidence')[-1].split('.sd')[0])


def get_confidences(one_predict_result_dir: str) -> (str, List[float]):
    predicted_names = os.listdir(one_predict_result_dir)
    sorted_predicted_names = sorted(predicted_names, key=get_rank)
    confidences = [get_confidence(file_name) for file_name in sorted_predicted_names[1:]]
    return confidences


def get_results(all_results_dir: str, num_predicts):
    proteins, ligands = [], []
    rank_confidences = [[] for _ in range(num_predicts)]
    all_dirs = [os.path.join(all_results_dir, dir_name) for dir_name in os.listdir(all_results_dir)]
    protein_id = current_task_name.split('_index_')[0]
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
    all_results_file.to_csv(os.path.join(os.path.dirname(all_results_dir), 'all_results.csv'), index=False)


if __name__ == "__main__":
    get_results(os.path.join(ALL_TASKS_DIR, current_task_name, "all_model_outputs"), NUM_PREDICTS)
