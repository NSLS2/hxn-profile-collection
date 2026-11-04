# Plan Logging Integration Checklist

This checklist guides the integration of the plan logging system into the existing `run_*` plans.

## ✅ System Setup (COMPLETED)

- [x] Created `startup/99-plan-logger.py` with `PlanLogger` class
- [x] Log directory path: `/nsls2/data/hxn/legacy/users/json_plans_log/`
- [x] JSONL format with monthly rotation (plan_name_MMYYYY.jsonl)
- [x] Created comprehensive documentation in `PLAN_LOGGING_GUIDE.md`

## 📝 Integration Steps (TO DO)

### Step 1: Update `startup/72-json-plans.py`

- [ ] Add import at the top:
  ```python
  from plan_logger import PlanLogger
  ```

### Step 2: Add Logging to `run_tomo_json`

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

### Unit Testing

- [ ] Verify `startup/99-plan-logger.py` syntax:
  ```bash
  python3 -m py_compile startup/99-plan-logger.py
  ```

- [ ] Verify modified `startup/72-json-plans.py` syntax:
  ```bash
  python3 -m py_compile startup/72-json-plans.py
  ```

### Integration Testing

- [ ] Create test JSON configuration files
- [ ] Run a test plan (e.g., with minimal angles/repeats)
- [ ] Verify log file is created at:
  ```bash
  /nsls2/data/hxn/legacy/users/json_plans_log/
  ```

- [ ] Check log file format and content:
  ```bash
  head -5 /nsls2/data/hxn/legacy/users/json_plans_log/run_diff_json_*.jsonl
  ```

- [ ] Verify each log entry is valid JSON:
  ```bash
  python3 -c "
  import json
  with open('/nsls2/data/hxn/legacy/users/json_plans_log/run_diff_json_*.jsonl') as f:
      for line in f:
          json.loads(line)  # Will raise if invalid
  print('✓ All log entries are valid JSON')
  "
  ```

### Failure Recovery Testing

- [ ] Interrupt a running plan mid-execution
- [ ] Verify "failed" event is logged with error message
- [ ] Modify JSON config to resume from last completed state
- [ ] Re-run plan and verify new log entries are appended

## 📊 Documentation

- [x] Main guide: `PLAN_LOGGING_GUIDE.md`
- [ ] Add code comments in `startup/72-json-plans.py` explaining logging
- [ ] Update beamline wiki/documentation if applicable

## 🚀 Deployment

- [ ] Commit changes:
  ```bash
  git add startup/72-json-plans.py
  git commit -m "Integrate plan logging for run_* functions"
  ```

- [ ] Test in beamline environment
- [ ] Monitor first week of logs
- [ ] Document any issues or improvements needed

## 📋 Maintenance

### Regular Tasks

- [ ] Monthly: Review log file sizes and archive if needed
- [ ] Quarterly: Review recovery procedures and document lessons learned
- [ ] As needed: Update logging parameters or add additional tracked metrics

### Troubleshooting Checklist

- [ ] If logs not created: Check directory permissions at `/nsls2/data/hxn/legacy/users/json_plans_log/`
- [ ] If import fails: Verify `startup/99-plan-logger.py` is in startup folder
- [ ] If JSON parsing fails: Validate log file format with `python3 -m json.tool`

## Notes

- The logging system is designed to be non-intrusive and won't break plans if logging fails
- Logs are organized by month to prevent files from becoming too large
- Each month creates new files (e.g., `run_diff_json_112026.jsonl`)
- Log entries can be easily parsed and analyzed programmatically

---

**Last Updated:** 2026-11-04
**Status:** Ready for integration
