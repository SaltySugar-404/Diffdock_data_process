import argparse
import os.path
import shutil

from after_model import *
from before_model import *
from during_model import *


def get_tasks():
    return [file for file in os.listdir(BASE_TASKS_DIR) if os.path.isdir(os.path.join(BASE_TASKS_DIR, file))]


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
    all_root_protein_files = [os.path.join(BASE_PROTEIN_INPUT_DIR, f) for f in os.listdir(BASE_PROTEIN_INPUT_DIR)]
    all_root_protein_files = [f for f in all_root_protein_files if f.endswith('.pdb')]
    if len(all_root_protein_files) == 0:
        print(f"No protein file found in {BASE_PROTEIN_INPUT_DIR}")
        exit(1)
    if len(all_root_protein_files) > 1:
        print(f"Too mach protein files in {BASE_PROTEIN_INPUT_DIR}")
        exit(1)

    # check ligand inputs
    all_root_unsplit_ligand_files = [os.path.join(BASE_UNSPLIT_LIGANDS_DIR, f) for f in os.listdir(BASE_UNSPLIT_LIGANDS_DIR)]
    all_root_unsplit_ligand_files = [f for f in all_root_unsplit_ligand_files if f.endswith('.sdf')]
    all_root_split_ligand_files = [os.path.join(BASE_SPLIT_LIGANDS_DIR, f) for f in os.listdir(BASE_SPLIT_LIGANDS_DIR)]
    all_root_split_ligand_files = [f for f in all_root_split_ligand_files if f.endswith('.sdf')]
    if len(all_root_split_ligand_files) == 0:
        if len(all_root_unsplit_ligand_files) == 0:
            print(f"No ligand file found in {BASE_UNSPLIT_LIGANDS_DIR}")
            exit(1)
        if len(all_root_unsplit_ligand_files) > 1:
            print(f"Too much ligand files found in {BASE_UNSPLIT_LIGANDS_DIR}")
            exit(1)
        split_sdf(BASE_UNSPLIT_LIGANDS_DIR, BASE_SPLIT_LIGANDS_DIR)
        # check split
        if not check_split(BASE_SPLIT_LIGANDS_DIR):
            print(f"Error split keyword: {SPLIT_KEYWORD}")
            exit(1)

    # make dirs
    task_dir = os.path.join(BASE_TASKS_DIR, task_id)
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
    generate_model_input(task_protein_file, BASE_SPLIT_LIGANDS_DIR, all_model_inputs_dir)
    get_save_commands(all_model_inputs_dir, all_model_outputs_dir, logs_dir, sh_dir)
    print(f"Task id: {task_id} created")
    print(f"To run task id: {task_id}, do as follows:")
    print(f"cd ..")
    print(f"nohup bash {os.path.join(BASE_SH_DIR, 'run_all.sh')} {os.path.join(BASE_DIR_NAME, sh_dir)} > run_all.log 2>&1 &")


def delete_task(task_id: str):
    if task_id not in get_tasks():
        print(f"Task id: {task_id} not found")
        return

    shutil.rmtree(os.path.join(BASE_TASKS_DIR, task_id))
    print(f"Task id: {task_id} deleted")


def clear_pid(task_id: str):
    if task_id not in get_tasks():
        print(f"Task id: {task_id} not found")
        return

    task_dir = os.path.join(BASE_TASKS_DIR, task_id)
    delete_pid_file(os.path.join(task_dir, "sh"))
    print(f"Task id: {task_id} cached pid files cleared")


def check_task_status(task_id: str):
    if task_id not in get_tasks():
        print(f"Task id: {task_id} not found")
        return

    task_dir = os.path.join(BASE_TASKS_DIR, task_id)
    sh_dir = os.path.join(task_dir, "sh")
    latest_pid_status = get_latest_pid_status(sh_dir, len(GPU_INDEXES))
    progress_now = get_progress(sh_dir)
    time_used = get_time_used(sh_dir)
    if True in latest_pid_status:
        print(f"Task id: {task_id} running, progress: {progress_now}, time used: {time_used} minutes")
    else:
        print(f"Task id: {task_id} not running, progress: {progress_now}, time used: {time_used} minutes")


def save_results(task_id: str):
    if task_id not in get_tasks():
        print(f"Task id: {task_id} not found")
        return

    task_dir = os.path.join(BASE_TASKS_DIR, task_id)
    protein_id = [f for f in os.listdir(task_dir) if f.endswith('.pdb')][0].split('.')[0]
    all_model_outputs_dir = os.path.join(task_dir, "all_model_outputs")
    output_file = os.path.join(task_dir, "results.csv")
    get_save_results(protein_id, all_model_outputs_dir, 10, output_file)
    print(f"Task id: {task_id} {output_file} saved")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Task Manager")

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--list", action="store_true", help="List all tasks")
    group.add_argument("--new", metavar="TASK_ID", type=str, help="Create new task with TASK_ID")
    group.add_argument("--clear", metavar="TASK_ID", type=str, help="Clear cached pid with TASK_ID")
    group.add_argument("--status", metavar="TASK_ID", type=str, help="Show status with TASK_ID")
    group.add_argument("--result", metavar="TASK_ID", type=str, help="Generate result.csv with TASK_ID")
    group.add_argument("--delete", metavar="TASK_ID", type=str, help="Delete task with TASK_ID")

    args = parser.parse_args()

    if args.list:
        list_tasks()
    elif args.new:
        new_task(args.new)
    elif args.clear:
        clear_pid(args.clear)
    elif args.status:
        clear_pid(args.status)
    elif args.result:
        save_results(args.result)
    elif args.delete:
        delete_task(args.delete)
