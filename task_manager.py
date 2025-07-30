import json
from typing import List

from before_model import *


def get_tasks():
    with open(TASK_MANAGER_FILE, "r") as f:
        all_tasks = json.load(f)
    return all_tasks["finished_tasks"], all_tasks["running_tasks"]


def update_tasks(finished_tasks: List[str], running_tasks: List[str]):
    with open(TASK_MANAGER_FILE, "w") as f:
        json.dump({"finished_tasks": finished_tasks, "running_tasks": running_tasks}, f)


def new_task(task_id: str):
    finished_tasks, running_tasks = get_tasks()
    all_tasks = finished_tasks + running_tasks
    if task_id in all_tasks:
        print(f'Task {task_id} already exists, adding as {task_id}_new')
        task_id += "_new"

    task_dir = os.path.join(ALL_TASKS_DIR, task_id)
    input_one_protein_dir = os.path.join(task_dir, "input_one_protein_here")
    input_unsplit_ligands_dir = os.path.join(task_dir, "input_unsplit_ligands_here")
    input_split_ligands_dir = os.path.join(task_dir, "input_split_ligands_here")
    all_model_inputs_dir = os.path.join(task_dir, "all_model_inputs")
    all_model_outputs_dir = os.path.join(task_dir, "all_model_outputs")
    logs_dir = os.path.join(task_dir, "logs")
    os.makedirs(task_dir, exist_ok=True)
    os.makedirs(input_one_protein_dir, exist_ok=True)
    os.makedirs(input_unsplit_ligands_dir, exist_ok=True)
    os.makedirs(input_split_ligands_dir, exist_ok=True)
    os.makedirs(all_model_inputs_dir, exist_ok=True)
    os.makedirs(all_model_outputs_dir, exist_ok=True)
    os.makedirs(logs_dir, exist_ok=True)


def before_model_process(task_id: str, split_key_word: str):
    task_dir = os.path.join(ALL_TASKS_DIR, task_id)
    input_one_protein_dir = os.path.join(task_dir, "input_one_protein_here")
    input_unsplit_ligands_dir = os.path.join(task_dir, "input_unsplit_ligands_here")
    input_split_ligands_dir = os.path.join(task_dir, "input_split_ligands_here")
    all_model_inputs_dir = os.path.join(task_dir, "all_model_inputs")
    all_model_outputs_dir = os.path.join(task_dir, "all_model_outputs")
    logs_dir = os.path.join(task_dir, "logs")
    if len(os.listdir(all_model_inputs_dir)):
        return
    if len(os.listdir(input_one_protein_dir)) == 0:
        return
    if not len(os.listdir(input_one_protein_dir)) == 1:
        return
    if not len(os.listdir(input_unsplit_ligands_dir)) + len(os.listdir(input_split_ligands_dir)):
        return

    if len(os.listdir(input_split_ligands_dir)) == 0 and len(os.listdir(input_unsplit_ligands_dir)):
        split_sdf(input_unsplit_ligands_dir, input_unsplit_ligands_dir, split_key_word)
    generate_model_input(task_id, len(GPU_INDEXES))
    commands = get_commands(task_id)
    for command in commands:
        print(command)


def show_tasks():
    finished_tasks, running_tasks = get_tasks()
    print(f"finished tasks:{finished_tasks}")
    print(f"running tasks:{running_tasks}")


def insert_task(task_id: str, index: int):
    finished_tasks, running_tasks = get_tasks()
    if not (1 <= index <= len(running_tasks)):
        print(f"Fail to insert task: {task_id}, index:{index} is out of valid range [1, {len(running_tasks)}]")
        return
    running_tasks.insert(index, task_id)
    update_tasks(finished_tasks, running_tasks)


if __name__ == "__main__":
    show_tasks()
    before_model_process(task_id="test_3", split_key_word="zinc_id")
