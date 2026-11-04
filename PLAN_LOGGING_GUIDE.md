# Plan Logging Guide

This guide explains how to use the plan logging system for tracking `run_*` plan executions.

## Overview

The plan logging system automatically tracks when plans start and end, saving records to JSONL files organized by plan name and month/year. This enables recovery if plans are interrupted.

- **Log Location**: `/nsls2/data/hxn/legacy/users/json_plans_log/`
- **Log File Naming**: `{plan_name}_{MMYYYY}.jsonl` (e.g., `run_diff_json_112026.jsonl` for November 2026)
- **Format**: JSONL (JSON Lines) - one JSON object per line for easy parsing and recovery

## Setup

The logging system is initialized automatically when `startup/99-plan-logger.py` is loaded. Ensure this file is in your startup folder.

## Usage in Plans

### Basic Usage

At the beginning of your `run_*` plan:

```python
from plan_logger import PlanLogger

def run_diff_json(path_to_json, tracking_file=None, do_confirm=True):
    # Log the start
    PlanLogger.log_plan_start(
        plan_name="run_diff_json",
        json_path=path_to_json,
        extra_params={"do_confirm": do_confirm}
    )
    
    try:
        # Your plan code here
        # ... perform scans ...
        
        # Log successful completion
        PlanLogger.log_plan_end(
            plan_name="run_diff_json",
            status="completed"
        )
    except Exception as e:
        # Log failure
        PlanLogger.log_plan_end(
            plan_name="run_diff_json",
            status="failed",
            error_msg=str(e)
        )
        raise
```

### API Reference

#### `PlanLogger.log_plan_start(plan_name, json_path, extra_params=None)`

Log the start of a plan execution.

**Arguments:**
- `plan_name` (str): Name of the plan, e.g., "run_diff_json"
- `json_path` (str): Path to the JSON configuration file
- `extra_params` (dict, optional): Additional parameters to log (e.g., number of repeats, angles, etc.)

**Example:**
```python
PlanLogger.log_plan_start(
    plan_name="run_repeat_2d_json",
    json_path="/path/to/diff_params.json",
    extra_params={"n_repeat": 10, "scan_angle": 45.5}
)
```

#### `PlanLogger.log_plan_end(plan_name, status="completed", error_msg=None, extra_data=None)`

Log the completion or failure of a plan.

**Arguments:**
- `plan_name` (str): Name of the plan
- `status` (str): "completed" or "failed"
- `error_msg` (str, optional): Error message if status is "failed"
- `extra_data` (dict, optional): Additional data like scan IDs, number of images collected, etc.

**Example:**
```python
PlanLogger.log_plan_end(
    plan_name="run_diff_json",
    status="completed",
    extra_data={
        "scan_ids": [100, 101, 102],
        "angles_completed": 45,
        "total_angles": 50
    }
)
```

## Log File Format

Each line in a JSONL file is a complete JSON object. Example:

```json
{"timestamp": "2026-11-04T15:30:45.123456", "event": "started", "json_path": "/data/users/current_user/diff_params.json", "params": {"do_confirm": true}}
{"timestamp": "2026-11-04T15:45:30.654321", "event": "completed", "data": {"scan_ids": [100, 101, 102]}}
```

## Querying Logs

### View recent log entries (command line)

```bash
tail -20 /nsls2/data/hxn/legacy/users/json_plans_log/run_diff_json_112026.jsonl
```

### Parse and analyze logs (Python)

```python
import json
from pathlib import Path

log_path = Path("/nsls2/data/hxn/legacy/users/json_plans_log/run_diff_json_112026.jsonl")

with open(log_path) as f:
    for line in f:
        record = json.loads(line)
        print(f"{record['timestamp']}: {record['event']}")
        if 'error' in record:
            print(f"  ERROR: {record['error']}")
```

### Recover from interrupted plans

If a plan fails, check the logs to see:
1. When it started
2. What parameters were used
3. Where/why it failed
4. How many iterations completed (if applicable)

Then re-run the plan or manually resume from the last completed angle/iteration.

## Log Directory Setup

The log directory is created automatically on first use. If you need to set it up manually:

```bash
mkdir -p /nsls2/data/hxn/legacy/users/json_plans_log
chmod 775 /nsls2/data/hxn/legacy/users/json_plans_log
```

## Examples for Existing Plans

### `run_tomo_json`

```python
def run_tomo_json(path_to_json, tracking_file=None):
    PlanLogger.log_plan_start(
        plan_name="run_tomo_json",
        json_path=path_to_json
    )
    
    try:
        # ... existing plan code ...
        PlanLogger.log_plan_end(
            plan_name="run_tomo_json",
            status="completed",
            extra_data={"final_angle": angles[-1]}
        )
    except Exception as e:
        PlanLogger.log_plan_end(
            plan_name="run_tomo_json",
            status="failed",
            error_msg=str(e)
        )
        raise
```

### `run_diff_json`

```python
def run_diff_json(path_to_json, tracking_file=None, do_confirm=True):
    PlanLogger.log_plan_start(
        plan_name="run_diff_json",
        json_path=path_to_json,
        extra_params={"do_confirm": do_confirm}
    )
    
    try:
        # ... existing plan code ...
        PlanLogger.log_plan_end(
            plan_name="run_diff_json",
            status="completed",
            extra_data={"angles_completed": len(angles)}
        )
    except Exception as e:
        PlanLogger.log_plan_end(
            plan_name="run_diff_json",
            status="failed",
            error_msg=str(e)
        )
        raise
```

### `run_repeat_2d_json`

```python
def run_repeat_2d_json(path_to_json, n_repeat=10, tracking_file=None, do_confirm=True):
    PlanLogger.log_plan_start(
        plan_name="run_repeat_2d_json",
        json_path=path_to_json,
        extra_params={"n_repeat": n_repeat}
    )
    
    try:
        # ... existing plan code ...
        PlanLogger.log_plan_end(
            plan_name="run_repeat_2d_json",
            status="completed",
            extra_data={"repeats_completed": n_repeat}
        )
    except Exception as e:
        PlanLogger.log_plan_end(
            plan_name="run_repeat_2d_json",
            status="failed",
            error_msg=str(e)
        )
        raise
```

## Troubleshooting

**Q: Why aren't logs being created?**
- Check if `/nsls2/data/hxn/legacy/users/json_plans_log/` directory exists and has write permissions
- Check if startup file is being loaded (look for initialization message in startup log)

**Q: Can I change the log location?**
- Yes, modify `PlanLogger.LOG_BASE_PATH` in `startup/99-plan-logger.py`

**Q: How do I prevent log files from getting too large?**
- Files are automatically organized by month/year, creating a new file each month
- If a month's file gets too large, you can archive old logs manually
