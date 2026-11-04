# Plan Logging Integration Checklist

This checklist guides the integration of the plan logging system into the `run_*` plans.

## ✅ System Setup (COMPLETED)

- [x] Integrated `PlanLogger` class into `startup/72-json-plans.py`
- [x] Log directory path: `/nsls2/data/hxn/legacy/users/json_plans_log/`
- [x] JSONL format with monthly rotation (plan_name_MMYYYY.jsonl)
- [x] Complete documentation as docstrings in the code

## 📝 Integration Steps (TO DO)

### Step 1: Import the Logger

The `PlanLogger` class is already available in `startup/72-json-plans.py`. No additional import needed.

### Step 2: Add Logging to `run_tomo_json`

Location: `startup/72-json-plans.py`, line ~575

- [ ] Add at function start (after docstring):
  ```python
  PlanLogger.log_plan_start(
      plan_name="run_tomo_json",
      json_path=path_to_json
  )
  ```

- [ ] Add at end of function (before return or final operations):
  ```python
  PlanLogger.log_plan_end(
      plan_name="run_tomo_json",
      status="completed"
  )
  ```

- [ ] Wrap in try/except to capture failures:
  ```python
  except Exception as e:
      PlanLogger.log_plan_end(
          plan_name="run_tomo_json",
          status="failed",
          error_msg=str(e)
      )
      raise
  ```

### Step 3: Add Logging to `run_diff_json`

Location: `startup/72-json-plans.py`, line ~1094

- [ ] Add at function start (after docstring):
  ```python
  PlanLogger.log_plan_start(
      plan_name="run_diff_json",
      json_path=path_to_json,
      extra_params={"do_confirm": do_confirm}
  )
  ```

- [ ] Add at end of function:
  ```python
  PlanLogger.log_plan_end(
      plan_name="run_diff_json",
      status="completed",
      extra_data={"angles_completed": len(angles)}
  )
  ```

- [ ] Add error handling:
  ```python
  except Exception as e:
      PlanLogger.log_plan_end(
          plan_name="run_diff_json",
          status="failed",
          error_msg=str(e)
      )
      raise
  ```

### Step 4: Add Logging to `run_repeat_2d_json`

Location: `startup/72-json-plans.py`, line ~1445

- [ ] Add at function start (after docstring):
  ```python
  PlanLogger.log_plan_start(
      plan_name="run_repeat_2d_json",
      json_path=path_to_json,
      extra_params={"n_repeat": n_repeat}
  )
  ```

- [ ] Add at end of function:
  ```python
  PlanLogger.log_plan_end(
      plan_name="run_repeat_2d_json",
      status="completed",
      extra_data={"repeats_completed": n_repeat}
  )
  ```

- [ ] Add error handling:
  ```python
  except Exception as e:
      PlanLogger.log_plan_end(
          plan_name="run_repeat_2d_json",
          status="failed",
          error_msg=str(e)
      )
      raise
  ```

## 🧪 Testing

### Syntax Verification

- [ ] Verify modified file compiles:
  ```bash
  python3 -m py_compile startup/72-json-plans.py
  ```

### Integration Testing

- [ ] Create test JSON configuration files
- [ ] Run a test plan (with minimal angles/repeats)
- [ ] Verify log file is created at:
  ```bash
  ls -la /nsls2/data/hxn/legacy/users/json_plans_log/
  ```

- [ ] Check log file content:
  ```bash
  head -5 /nsls2/data/hxn/legacy/users/json_plans_log/run_diff_json_*.jsonl
  ```

- [ ] Verify log entries are valid JSON:
  ```bash
  python3 << 'EOF'
  import json
  with open('/nsls2/data/hxn/legacy/users/json_plans_log/run_diff_json_*.jsonl') as f:
      for line in f:
          json.loads(line)
  print("✓ All log entries are valid JSON")
  EOF
  ```

### Failure Recovery Testing

- [ ] Interrupt a running plan mid-execution
- [ ] Verify "failed" event is logged with error message
- [ ] Modify JSON config to resume from last completed state
- [ ] Re-run plan and verify new log entries are appended

## 📊 Documentation

The complete PlanLogger documentation is available as docstrings in:
- `startup/72-json-plans.py` - Search for "class PlanLogger"

View the docstrings:
```python
from startup.file_72-json-plans import PlanLogger
help(PlanLogger)
help(PlanLogger.log_plan_start)
help(PlanLogger.log_plan_end)
```

## 🚀 Deployment

- [ ] Commit changes:
  ```bash
  git add startup/72-json-plans.py
  git commit -m "Integrate plan logging for run_* functions"
  ```

- [ ] Test in beamline environment
- [ ] Monitor first week of logs
- [ ] Document any issues or improvements needed

## 📋 Usage Reference

### Log Start Event

```python
PlanLogger.log_plan_start(
    plan_name="run_diff_json",
    json_path=path_to_json,
    extra_params={"param_name": param_value}  # optional
)
```

### Log End Event - Success

```python
PlanLogger.log_plan_end(
    plan_name="run_diff_json",
    status="completed",
    extra_data={"scan_ids": [100, 101], "angles_completed": 45}  # optional
)
```

### Log End Event - Failure

```python
PlanLogger.log_plan_end(
    plan_name="run_diff_json",
    status="failed",
    error_msg="Beam dump at angle 30.5"
)
```

## 📂 Log File Organization

```
/nsls2/data/hxn/legacy/users/json_plans_log/
├── run_tomo_json_112026.jsonl         (November 2026)
├── run_diff_json_112026.jsonl
├── run_repeat_2d_json_112026.jsonl
├── run_tomo_json_122026.jsonl         (December 2026)
└── ... (one file per plan per month)
```

## 📋 Maintenance

### Regular Tasks

- [ ] Monthly: Review log file sizes and archive if needed
- [ ] Quarterly: Review recovery procedures
- [ ] As needed: Add additional tracked metrics

### Troubleshooting

- If logs not created: Check directory permissions at `/nsls2/data/hxn/legacy/users/json_plans_log/`
- If JSON parsing fails: Validate log file format with `python3 -m json.tool <logfile>`
- The logging system is non-intrusive - failures in logging won't break plans

---

**Last Updated:** 2026-11-04
**Status:** Code integrated and ready for plan function integration
