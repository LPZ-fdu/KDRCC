# Server execution

1. Copy the source project to the server. Set `code_root` and `work_root` in `05_submit.R` as needed; use a distinct work directory for this study run.
2. Run `Rscript --vanilla 05_submit.R install`, then `Rscript --vanilla 05_submit.R check`.
3. Review `Rscript --vanilla 05_submit.R plan A`.
4. Submit `nohup Rscript --vanilla 05_submit.R run A > A_run.txt 2>&1 &`.
5. After completion, submit `nohup Rscript --vanilla 05_submit.R report A > A_report.txt 2>&1 &`.
6. Review A's parameter comparisons, set `freeze$config_id` in the submission file, and run `Rscript --vanilla 05_submit.R freeze`.
7. Review `Rscript --vanilla 05_submit.R plan B,C`, then submit `nohup Rscript --vanilla 05_submit.R run B,C > BC_run.txt 2>&1 &`.
8. Report completed B and C separately with `nohup Rscript --vanilla 05_submit.R report B > B_report.txt 2>&1 &` and the corresponding C command.

The words after `05_submit.R` are command-line arguments: the action and the module list. The submission script reads them and records the effective action. All configuration changes remain in that one script.

Default workers=32 launches independent R processes. Each uses one internal numerical-library thread. Active PIDs are saved in `worker_processes.csv`; console output gives timestamps, total completed fits, elapsed time, estimated remaining time and active workers' major stages. Resource records retain component timings and measured memory/cache information where supported.

Use `tail -f A_run.txt` or `tail -f BC_run.txt` to inspect progress. Ctrl+C exits tail. The run itself remains in the background. Each fit has an independent checkpoint. Reporting can be run again from the saved result directory after presentation changes.

Set `common_fixed`, `B_overrides`, `C_overrides` and `D_overrides` to override frozen parameters explicitly. Empty lists inherit the frozen values. Default C null counts are 8,16,40,80; `null_counts` permits a manually specified list. The parameter file is checked against the generating-mechanism identity.

## Module D submission

Point `work_root` to the directory containing the reviewed parameter file. Submit `nohup Rscript --vanilla 05_submit.R plan D > D_plan.txt 2>&1 &` and review the 1,600-new-fit, 200-kernel default plan and the reused B reference. Then submit `nohup Rscript --vanilla 05_submit.R run D > D_run.txt 2>&1 &`. After fitting finishes, run `nohup Rscript --vanilla 05_submit.R report D > D_report.txt 2>&1 &`.

The D command selects only that module. Existing modules retain their own result directories and checkpoints. Full submission settings, correlation definitions, output paths and parameter override examples are in MODULE_D.md.

## Module E submission

Keep the work directory containing the reviewed `config/fixed_parameters.R`. Select E explicitly on the command line. Run `nohup Rscript --vanilla 05_submit.R check E > E_check.txt 2>&1 &` and review the completed checks. Then run `nohup Rscript --vanilla 05_submit.R plan E > E_plan.txt 2>&1 &` and review its scenarios and candidates. Submit `nohup Rscript --vanilla 05_submit.R run E > E_run.txt 2>&1 &`. After completion, submit `nohup Rscript --vanilla 05_submit.R report E > E_report.txt 2>&1 &`.

E writes new result and cache directories. It performs its own paired fixed-configuration comparison on fresh data and records the frozen configuration's provenance. `tail -f E_run.txt` shows final-fit completion, validation progress within each active repetition, elapsed time and remaining-time estimates. Settings and override examples are in MODULE_E.md.


## Saved-model contribution evaluation

Contribution shares evaluate the eight prespecified active exposures at the reference background, normalized within those eight variables. Use `05_submit.R contributions C,E` followed by `05_submit.R report C,E` to calculate contributions from saved fitted functions and write the standard reports. Execution settings remain in `05_submit.R`. See [Contribution evaluation](CONTRIBUTIONS.md).
