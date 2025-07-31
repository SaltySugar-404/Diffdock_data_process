import argparse
import os.path
import shutil

from before_model import *


def get_tasks():
    return [file for file in os.listdir(ALL_TASKS_DIR) if os.path.isdir(os.path.join(ALL_TASKS_DIR, file))]


def list_tasks():
    print(f"All tasks: {get_tasks()}")


def new_task(task_id: str):
    # check split keyword
    if not len(SPLIT_KEYWORD):
        print(f"Split keyword empty")
        exit(1)
    # check task ids
    all_tasks = get_tasks()
    while task_id in all_tasks:
        print(f"Task id: {task_id} already in tasks list, trying {str(task_id + '_new')}")
        task_id = str(task_id + "_new")
    # check protein inputs
    all_root_protein_files = [os.path.join(OUTSIDE_PROTEIN_INPUT_DIR, f) for f in os.listdir(OUTSIDE_PROTEIN_INPUT_DIR)]
    all_root_protein_files = [f for f in all_root_protein_files if f.endswith('.pdb')]
    if len(all_root_protein_files) == 0:
        print(f"No protein file found in {OUTSIDE_PROTEIN_INPUT_DIR}")
        exit(1)
    if len(all_root_protein_files) > 1:
        print(f"Too mach protein files in {OUTSIDE_PROTEIN_INPUT_DIR}")
        exit(1)

    # check ligand inputs
    all_root_unsplit_ligand_files = [os.path.join(OUTSIDE_UNSPLIT_LIGANDS_DIR, f) for f in os.listdir(OUTSIDE_UNSPLIT_LIGANDS_DIR)]
    all_root_unsplit_ligand_files = [f for f in all_root_unsplit_ligand_files if f.endswith('.sdf')]
    all_root_split_ligand_files = [os.path.join(OUTSIDE_SPLIT_LIGANDS_DIR, f) for f in os.listdir(OUTSIDE_SPLIT_LIGANDS_DIR)]
    all_root_split_ligand_files = [f for f in all_root_split_ligand_files if f.endswith('.sdf')]
    if len(all_root_split_ligand_files) == 0:
        if len(all_root_unsplit_ligand_files) == 0:
            print(f"No ligand file found in {OUTSIDE_UNSPLIT_LIGANDS_DIR}")
            exit(1)
        if len(all_root_unsplit_ligand_files) > 1:
            print(f"Too much ligand files found in {OUTSIDE_UNSPLIT_LIGANDS_DIR}")
            exit(1)
        split_sdf(OUTSIDE_UNSPLIT_LIGANDS_DIR, OUTSIDE_SPLIT_LIGANDS_DIR)
        # check split
        if not check_split(OUTSIDE_SPLIT_LIGANDS_DIR):
            print(f"Error split keyword: {SPLIT_KEYWORD}")
            exit(1)

    # make dirs
    task_dir = os.path.join(ALL_TASKS_DIR, task_id)
    all_model_inputs_dir = os.path.join(task_dir, "all_model_inputs")
    all_model_outputs_dir = os.path.join(task_dir, "all_model_outputs")
    logs_dir = os.path.join(task_dir, "logs")
    sh_dir = os.path.join(task_dir, "sh")
    os.makedirs(task_dir, exist_ok=True)
    os.makedirs(all_model_inputs_dir, exist_ok=True)
    os.makedirs(all_model_outputs_dir, exist_ok=True)
    os.makedirs(logs_dir, exist_ok=True)
    os.makedirs(sh_dir, exist_ok=True)

    # process
    root_protein_file = all_root_protein_files[0]
    task_protein_file = os.path.join(task_dir, os.path.basename(root_protein_file))
    if not os.path.isfile(task_protein_file):
        shutil.copy(root_protein_file, task_dir)
    generate_model_input(task_protein_file, OUTSIDE_SPLIT_LIGANDS_DIR, all_model_inputs_dir)
    get_save_commands(all_model_inputs_dir, all_model_outputs_dir, logs_dir, sh_dir)
    print(f"Task id: {task_id} created")


def delete_task(task_id: str):
    all_tasks = get_tasks()
    if task_id not in all_tasks:
        print(f"Task id: {task_id} not found")
    else:
        shutil.rmtree(os.path.join(ALL_TASKS_DIR, task_id))
        print(f"Task id: {task_id} deleted")


def run_task(task_id: str, command_index: int):
    sh_dir = os.path.join(ALL_TASKS_DIR, "commands", task_id)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Task Manager")

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--list", action="store_true", help="List all tasks")
    group.add_argument("--new", metavar="TASK_ID", type=str, help="Create new task with TASK_ID")
    group.add_argument("--delete", metavar="TASK_ID", type=str, help="Delete task with TASK_ID")
    group.add_argument("--run", metavar="TASK_ID", type=str, help="Run task with TASK_ID")

    args = parser.parse_args()

    if args.list:
        list_tasks()
    elif args.new:
        new_task(args.new)
    elif args.delete:
        delete_task(args.delete)
