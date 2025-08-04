from init import *
import time

def get_pid_status(pid: int):
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    else:
        return True


def get_latest_pid_status(pid_dir: str, num_pids: int):
    all_pid_files = [f for f in os.listdir(pid_dir) if f.endswith(".pid")]
    if not all_pid_files:
        return []

    all_pid_files = [os.path.join(pid_dir, f) for f in all_pid_files]
    all_pid_files.sort(key=lambda x: os.path.getmtime(x), reverse=True)

    statuses = []
    for f in all_pid_files[:num_pids]:
        try:
            with open(f, "r") as fp:
                pid_str = fp.read().strip()
                pid = int(pid_str)
                statuses.append(get_pid_status(pid))
        except Exception as e:
            print(f"Failed to read or parse PID file {f}: {e}")
            statuses.append(False)
    return statuses


def get_progress(sh_dir: str):
    all_pid_files = [f for f in os.listdir(sh_dir) if f.endswith(".pid")]
    all_sh_files = [f for f in os.listdir(sh_dir) if f.endswith(".sh")]
    return f"{int(len(all_pid_files) / len(GPU_INDEXES))}/{len(all_sh_files)}"


def get_time_used(sh_dir: str):
    all_pid_files = [f for f in os.listdir(sh_dir) if f.endswith(".pid")]
    if not len(all_pid_files):
        return 0

    all_pid_files = [os.path.join(sh_dir, f) for f in all_pid_files]
    all_pid_files.sort(key=lambda x: os.path.getmtime(x), reverse=True)
    earliest_time = os.path.getmtime(all_pid_files[-1])
    all_sh_files = [f for f in os.listdir(sh_dir) if f.endswith(".sh")]
    if len(all_pid_files) < len(all_sh_files):
        latest_time = time.time()
    else:
        if True in get_latest_pid_status(sh_dir, len(GPU_INDEXES)):
            latest_time = time.time()
        else:
            latest_time = os.path.getmtime(all_pid_files[0])

    time_diff_sec = latest_time - earliest_time
    time_diff_minutes = round(time_diff_sec / 60.0, 1)
    return time_diff_minutes
