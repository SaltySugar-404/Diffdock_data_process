import os

# dirs

# 当前目录
BASE_DIR_NAME = os.path.basename(os.getcwd())

# 上一级路径
PARENT_DIR_PATH = os.path.dirname(os.getcwd())

BASE_PROTEIN_INPUT_DIR = "input_one_protein_here"
os.makedirs(BASE_PROTEIN_INPUT_DIR, exist_ok=True)

BASE_UNSPLIT_LIGANDS_DIR = "input_unsplit_ligands_here"
os.makedirs(BASE_UNSPLIT_LIGANDS_DIR, exist_ok=True)

BASE_SPLIT_LIGANDS_DIR = "input_split_ligands_here"
os.makedirs(BASE_SPLIT_LIGANDS_DIR, exist_ok=True)

# 所有任务
BASE_TASKS_DIR = "all_tasks"
os.makedirs(BASE_TASKS_DIR, exist_ok=True)

# sh脚本存放目录
BASE_SH_DIR = "sh"
os.makedirs(BASE_SH_DIR, exist_ok=True)

# configs

# 可用gpu数量
GPU_INDEXES = [0, 1, 2, 3]

# 用于划分原始sdf的关键词
SPLIT_KEYWORD = "zinc_id"

# 每个模型输入区块的结合物数量，用以减小显存占用，提高处理速度
CHUNK_SIZE = 1000

if __name__ == "__main__":
    print("Initializing diffdock model...")
