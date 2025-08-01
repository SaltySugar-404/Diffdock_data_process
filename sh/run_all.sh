#!/bin/bash

sh_dir="$1"

if [ -z "$sh_dir" ]; then
    echo "$(date +"%Y-%m-%d %H:%M:%S") Usage: $0 <sh_dir>"
    exit 1
fi

for run_script in $(ls "$sh_dir"/run_*.sh | sort); do
    echo "$(date +"%Y-%m-%d %H:%M:%S") Running $run_script"
    bash "$run_script"

    script_base=$(basename "$run_script" .sh)
    pids_to_wait=()

    for i in {1..30}; do
        pid_files=($(ls "$sh_dir"/"${script_base}"_gpu_*.pid 2>/dev/null))
        if [ ${#pid_files[@]} -eq 4 ]; then
            break
        fi
        sleep 1
    done

    for pid_file in "$sh_dir"/"${script_base}"_gpu_*.pid; do
        if [ -f "$pid_file" ]; then
            pid=$(cat "$pid_file")
            if [ -n "$pid" ]; then
                echo "$(date +"%Y-%m-%d %H:%M:%S") Waiting for PID $pid from $pid_file"
                pids_to_wait+=("$pid")
            fi
        fi
    done

    for pid in "${pids_to_wait[@]}"; do
        while kill -0 "$pid" 2>/dev/null; do
            sleep 5
        done
        echo "$(date +"%Y-%m-%d %H:%M:%S") PID $pid finished."
    done

    echo "$(date +"%Y-%m-%d %H:%M:%S") $run_script done."
done

echo "$(date +"%Y-%m-%d %H:%M:%S") All scripts finished."
