import os


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
    if not len(all_pid_files):
        return []

    return [get_pid_status(f) for f in all_pid_files.sort(key=lambda x: os.path.getmtime(x), reverse=True)[:num_pids]]


def get_progress(sh_dir: str):
    all_pid_files = [f for f in os.listdir(sh_dir) if f.endswith(".pid")]
    all_sh_files = [f for f in os.listdir(sh_dir) if f.endswith(".sh")]
    return f"{len(all_pid_files)}/{len(all_sh_files)}"


def get_time_used(pid_dir: str):
    all_pid_files = [f for f in os.listdir(pid_dir) if f.endswith(".pid")]
    if not len(all_pid_files):
        return 0

    all_pid_files.sort(key=lambda x: os.path.getmtime(x), reverse=True)
    latest_time = os.path.getmtime(all_pid_files[0])
    earliest_time = os.path.getmtime(all_pid_files[-1])
    time_diff_sec = latest_time - earliest_time
    time_diff_minutes = time_diff_sec / 60.0
    return time_diff_minutes
