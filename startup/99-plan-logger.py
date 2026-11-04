print(f"Loading {__file__!r} ...")

import json
import os
from datetime import datetime
from pathlib import Path


class PlanLogger:
    """Logger for tracking plan executions to JSONL files organized by month/year."""
    
    LOG_BASE_PATH = "/nsls2/data/hxn/legacy/users/json_plans_log"
    
    @classmethod
    def ensure_log_directory(cls):
        """Create log directory if it doesn't exist."""
        try:
            log_dir = Path(cls.LOG_BASE_PATH)
            log_dir.mkdir(parents=True, exist_ok=True)
            return True
        except Exception as e:
            print(f"Warning: Could not create log directory {cls.LOG_BASE_PATH}: {e}")
            return False
    
    @classmethod
    def get_log_filepath(cls, plan_name):
        """Get the JSONL filepath for a plan, named with month/year.
        
        Examples:
            run_diff_json -> run_diff_json_112026.jsonl (for Nov 2026)
            run_tomo_json -> run_tomo_json_112026.jsonl
        """
        now = datetime.now()
        month_year = now.strftime("%m%Y")  # Format: MMYYYY
        filename = f"{plan_name}_{month_year}.jsonl"
        return os.path.join(cls.LOG_BASE_PATH, filename)
    
    @classmethod
    def log_plan_start(cls, plan_name, json_path, extra_params=None):
        """Log plan execution start."""
        if not cls.ensure_log_directory():
            return
        
        log_entry = {
            "timestamp": datetime.now().isoformat(),
            "event": "started",
            "json_path": str(json_path),
        }
        if extra_params:
            log_entry["params"] = extra_params
        
        try:
            log_filepath = cls.get_log_filepath(plan_name)
            with open(log_filepath, 'a') as f:
                f.write(json.dumps(log_entry) + '\n')
        except Exception as e:
            print(f"Error logging plan start: {e}")
    
    @classmethod
    def log_plan_end(cls, plan_name, status="completed", error_msg=None, extra_data=None):
        """Log plan execution completion or failure.
        
        Args:
            plan_name: Name of the plan (e.g., 'run_diff_json')
            status: 'completed' or 'failed'
            error_msg: Error message if status is 'failed'
            extra_data: Optional dict with additional data to log
        """
        if not cls.ensure_log_directory():
            return
        
        log_entry = {
            "timestamp": datetime.now().isoformat(),
            "event": status,
        }
        if error_msg:
            log_entry["error"] = error_msg
        if extra_data:
            log_entry["data"] = extra_data
        
        try:
            log_filepath = cls.get_log_filepath(plan_name)
            with open(log_filepath, 'a') as f:
                f.write(json.dumps(log_entry) + '\n')
        except Exception as e:
            print(f"Error logging plan end: {e}")


# Initialize logging directory on startup
if PlanLogger.ensure_log_directory():
    print(f"Plan logger initialized. Logs will be saved to: {PlanLogger.LOG_BASE_PATH}")
else:
    print(f"Warning: Plan logger could not initialize directory at {PlanLogger.LOG_BASE_PATH}")
