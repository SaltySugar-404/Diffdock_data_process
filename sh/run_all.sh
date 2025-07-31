#!/bin/bash

sh_dir="sh"

for run_script in $(ls "$sh_dir"/run_*.sh | sort); do
    echo "Running $run_script"
    bash "$run_script"

    for pid_file in $(grep -oP '>\s*\K\S+\.pid' "$run_script"); do
        pid_path="$sh_dir/$pid_file"
        if [ -f "$pid_path" ]; then
            pid=$(cat "$pid_path")
            if [ -n "$pid" ]; then
                while kill -0 "$pid" 2>/dev/null; do
                    sleep 10
                done
            fi
        fi
    done
    echo "$run_script done."
done

echo "All scripts finished."