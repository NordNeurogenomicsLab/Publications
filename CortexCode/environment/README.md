# R environment records

- [sessionInfo.txt](sessionInfo.txt): actual R session captured for the reconstruction plotting environment on 2026-11-06. It lists R, operating system, loaded packages and namespaces. No scientific model was fitted to capture this record.
- [package_versions.csv](package_versions.csv): installed versions and availability of packages referenced by the notebooks and helpers, including optional analysis dependencies. This is a version inventory, not an `renv` lockfile or proof that every optional analysis was executed.
- [archived_sessions/](archived_sessions/): unmodified session records from earlier reconstruction exports. Their numerical prefixes refer to the notebook numbering at capture time; in particular old `09` is the archived author-plot notebook, not the current `09_manuscript_graphics.Rmd`.
- `sessions/`: new notebook-specific session records. Main notebooks write these through `finish_notebook()`; pipeline notebooks call `write_session_info()`. These files can be committed independently of ignored figure reports.

These describe the reconstruction environment. They do not establish the original authors' historical package/database versions. Each notebook also prints `sessionInfo()` in its rendered report.
