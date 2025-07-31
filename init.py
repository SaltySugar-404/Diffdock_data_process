import os

# 项目文件夹
CODE_ROOT_DIR_NAME = r"run_diffdock"

OUTSIDE_PROTEIN_INPUT_DIR = "input_one_protein_here"
os.makedirs(OUTSIDE_PROTEIN_INPUT_DIR, exist_ok=True)

OUTSIDE_UNSPLIT_LIGANDS_DIR = "input_unsplit_ligands_here"
os.makedirs(OUTSIDE_UNSPLIT_LIGANDS_DIR, exist_ok=True)

OUTSIDE_SPLIT_LIGANDS_DIR = "input_split_ligands_here"
os.makedirs(OUTSIDE_SPLIT_LIGANDS_DIR, exist_ok=True)

# 所有任务
ALL_TASKS_DIR = "all_tasks"
os.makedirs(ALL_TASKS_DIR, exist_ok=True)

# 任务记录
TASK_MANAGER_FILE = "all_tasks.json"

# configs
# 可用gpu数量
GPU_INDEXES = [0, 1, 2, 3]

# 用于划分原始sdf的关键词
SPLIT_KEYWORD = "zinc_id"

# 每个模型输入区块的结合物数量，用以减小显存占用，提高处理速度
CHUNK_SIZE = 100

if __name__ == "__main__":
    print("Initializing diffdock model...")
