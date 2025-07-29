import os

# 项目文件夹
CODE_ROOT_DIR_NAME = r"run_diffdock"

# 所有任务
ALL_TASKS_DIR = "all_tasks"
os.makedirs(ALL_TASKS_DIR, exist_ok=True)

# 任务记录
TASK_MANAGER_FILE = "all_tasks.json"

# # 外部蛋白质输入
# OUTSIDE_PROTEIN_INPUT_DIR = r"put_one_protein_here"
# os.makedirs(OUTSIDE_PROTEIN_INPUT_DIR, exist_ok=True)
#
# # 外部已分离配体输入
# OUTSIDE_SPLIT_LIGAND_INPUT_DIR = r"put_split_ligands_here"
# os.makedirs(OUTSIDE_SPLIT_LIGAND_INPUT_DIR, exist_ok=True)
#
# # 外部未分离配体输入
# OUTSIDE_UNSPLIT_LIGAND_INPUT_DIR = r"put_unsplit_ligands_here"
# os.makedirs(OUTSIDE_UNSPLIT_LIGAND_INPUT_DIR, exist_ok=True)
#
# # 单个输入生产过的预测结构数量
# NUM_PREDICTS = 10
