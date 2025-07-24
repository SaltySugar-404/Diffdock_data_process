import json
import os
import shutil

# current task config

with open("task_config.json") as json_file:
    current_task_config = json.load(json_file)
SDF_SPLIT_KEYWORD = current_task_config["sdf_split_keyword"]
NUM_GPUS = current_task_config["num_gpus"]

# before model process

# 项目文件夹
CODE_ROOT_DIR_NAME = r"run_diffdock"

# 所有任务
ALL_TASKS = "tasks"
os.makedirs(ALL_TASKS, exist_ok=True)

# 当前任务文件夹
CURRENT_PROTEIN_NAME = os.listdir(r"put_one_protein_here")[0].split(".")[0]
all_task_name_and_index = os.listdir(ALL_TASKS)

max_task_index = -1
for task_name_and_index in all_task_name_and_index:
    if len(task_name_and_index.split(".")) == 2:
        continue
    task_name = task_name_and_index.split("_")[0]
    task_index = int(task_name_and_index.split("_")[1])
    if task_name == CURRENT_PROTEIN_NAME and task_index > max_task_index:
        max_task_index = task_index

CURRENT_TASK_DIR = os.path.join(ALL_TASKS, CURRENT_PROTEIN_NAME + f"{max_task_index + 1}")
os.makedirs(CURRENT_TASK_DIR, exist_ok=True)
shutil.copy("task_config.json", os.path.join(CURRENT_TASK_DIR, "task_config.json"))

# 蛋白质输入文件夹
PROTEIN_INPUT_DIR = os.path.join(CURRENT_TASK_DIR, "put_one_protein_here")
if not os.path.exists(PROTEIN_INPUT_DIR):
    shutil.move(r"put_one_protein_here", CURRENT_TASK_DIR)
os.makedirs(r"put_one_protein_here", exist_ok=True)

# 配体输入文件夹
UNSPLIT_LIGAND_DIR = os.path.join(CURRENT_TASK_DIR, "put_unsplit_ligands_here")
if not os.path.exists(UNSPLIT_LIGAND_DIR):
    shutil.move(r"put_unsplit_ligands_here", CURRENT_TASK_DIR)
os.makedirs(r"put_unsplit_ligands_here", exist_ok=True)

# 配体分离文件夹
SPLIT_LIGAND_DIR = os.path.join(CURRENT_TASK_DIR, "put_split_ligands_here")
if not os.path.exists(SPLIT_LIGAND_DIR):
    shutil.move(r"put_split_ligands_here", CURRENT_TASK_DIR)
os.makedirs(r"put_split_ligands_here", exist_ok=True)

# 蛋白质索引csv
PROTEIN_INDEX_FILE = os.path.join(CURRENT_TASK_DIR, "protein_index.csv")

# 配体索引csv
LIGAND_INDEX_FILE = os.path.join(CURRENT_TASK_DIR, "ligand_index.csv")

# 模型输入文件夹
MODEL_INPUTS_DIR = os.path.join(CURRENT_TASK_DIR, "all_model_inputs")
os.makedirs(MODEL_INPUTS_DIR, exist_ok=True)

# during model process

# 单个输入生产过的预测结构数量
NUM_PREDICTS = 10

# after model process

# 所有的模型输出
MODEL_OUTPUTS_DIR = os.path.join(CURRENT_TASK_DIR, "all_model_outputs")
os.makedirs(MODEL_OUTPUTS_DIR, exist_ok=True)

# 模型输出分析
MODEL_OUTPUTS_ANALYZE_FILE = os.path.join(CURRENT_TASK_DIR, "output_analyze.csv")
