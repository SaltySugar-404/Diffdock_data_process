This project implements a task-splitting mechanism for DiffDock, which helps reduce GPU memory usage and prevents performance degradation caused by memory leaks.

This project supports multi-GPU parallel execution and real-time progress monitoring, making large-scale docking tasks more efficient and manageable. 

## Quick Start:

Place this project under path/to/Diffdock/. You can modify the contents of init.py according to your specific use case, such as setting global configurations, paths, or default parameters.

The entry point of the project is task_manager.py.

example:

```python task_manager.py --new Task_id```