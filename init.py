import os

# 项目文件夹
CODE_ROOT_DIR_NAME = r"run_diffdock"

# 所有任务
ALL_TASKS_DIR = "all_tasks"
os.makedirs(ALL_TASKS_DIR, exist_ok=True)

# 任务记录
TASK_MANAGER_FILE = "all_tasks.json"

# 可用gpu数量
GPU_INDEXES = [0, 1, 2, 3]
