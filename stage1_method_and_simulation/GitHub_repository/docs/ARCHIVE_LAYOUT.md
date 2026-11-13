# Archive layout and restoration

## Git repository

The repository root contains the executable simulation project. `manuscript/` contains the editable appendix and its document source. `records/results/` contains all saved publication reports and module-level configuration, plan and session records. `records/submissions/` contains actual submitted settings and effective configuration records. `archive/` contains file-level and asset-level SHA-256 inventories.

## Release assets

Upload every `stage1_runtime_*.zip` file as an asset of the GitHub Release carrying the `stage1-20261107` source tag. Each ZIP is independently readable and contains paths rooted at `KDRCC_formal_20261027/`. Extracting all assets into the same parent restores the retained work directory. Each asset is bounded by approximately 1.25 GiB of file payload and fits the GitHub Release per-file limit.

The retained work directory includes all A–E fit records, learned kernels, auxiliary estimates, evaluation inputs, calibration objects, candidate scores, choices, task records, logs, source snapshots and publication reports. It includes successful and failed records. Snapshot files preserve their recorded names and contents.

The retention inventory classifies deterministic data/work caches and small-sample check output in `archive/excluded_files.csv`. The simulation generator, semantic random seeds, task plans and configuration records define regeneration of those caches. New simulation runs create their own writable work directory. Completed report reconstruction reads the retained fits and evaluation objects.

## Commands

Run from the source repository:

```bash
python tools/verify_stage1.py
python tools/verify_stage1.py --assets /path/to/release_assets
python tools/restore_records.py --assets /path/to/release_assets --destination /path/to/parent
```

The restoration command verifies asset SHA-256 hashes and ZIP CRCs, validates archive paths, and checks existing files against their archived bytes before reuse. Set `work_root` in `05_submit.R` to the restored directory when regenerating reports. Set `report_run_directories` to individual restored module directories as needed.

## Upload procedure

1. Upload the contents of `GitHub_repository/` as the Git repository and create the source tag `stage1-20261107`.
2. Create the corresponding GitHub Release and attach all runtime ZIP assets and `SHA256SUMS.csv`.
3. Publish the source and the matching runtime assets together. Record the final repository and Release links in the project publication metadata.

GitHub blocks ordinary Git files above 100 MiB and recommends keeping repositories small. Browser file upload supports files up to 25 MiB. GitHub Releases support up to 1,000 assets per release with each file under 2 GiB. This layout places detailed runtime objects in Release assets and keeps directly browsable material in Git.

References: [GitHub large-file guidance](https://docs.github.com/en/repositories/working-with-files/managing-large-files/about-large-files-on-github), [GitHub Release limits](https://docs.github.com/en/repositories/releasing-projects-on-github/about-releases).
