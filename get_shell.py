import subprocess

from init import *

commands = []

for index in range(NUM_GPUS):
    current_command = ""
    current_command += f"CUDA_VISIBLE_DEVICES={index} nohup python -m inference "
    current_command += f"--config default_inference_args.yaml "
    current_command += f"--protein_ligand_csv /home/user0/DiffDock/{CODE_ROOT_DIR_NAME}/tasks/{CURRENT_PROTEIN_NAME}/all_model_inputs/chunk_{index}.csv "

for index in range(NUM_GPUS):
    command = (
        f"CUDA_VISIBLE_DEVICES={index} nohup python -m inference "
        f"--config /home/user0/DiffDock/default_inference_args.yaml "
        f"--protein_ligand_csv /home/user0/DiffDock/{CODE_ROOT_DIR_NAME}/tasks/{CURRENT_PROTEIN_NAME}/all_model_inputs/chunk_{index}.csv "
        f"--out_dir /home/user0/DiffDock/{CODE_ROOT_DIR_NAME}/tasks/{CURRENT_PROTEIN_NAME}/all_model_outputs "
        f"> /home/user0/DiffDock/{CODE_ROOT_DIR_NAME}/tasks/{CURRENT_PROTEIN_NAME}/chunk_{index}_log.out 2>&1 &"
    )
    print(command)
    subprocess.Popen(command, shell=True)
