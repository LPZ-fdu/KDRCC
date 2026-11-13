# This is the only file to edit when submitting or reporting an experiment.
# It may remain here or be copied elsewhere. If copied, set code_root explicitly.
# Linux: Rscript --vanilla /path/to/05_submit.R plan
#        Rscript --vanilla /path/to/05_submit.R run A
#        Rscript --vanilla /path/to/05_submit.R report A
#        Rscript --vanilla /path/to/05_submit.R freeze
#        Rscript --vanilla /path/to/05_submit.R run B,C
#        Rscript --vanilla /path/to/05_submit.R report B,C

# ===================== EDITABLE SUBMISSION SETTINGS =====================
# code_root <- NULL  # NULL: source directory containing this submission script.
#                   # Linux example: "/home/user19/KDRCC_source"
# work_root <- NULL  # NULL: sibling directory "KDRCC_server_work".
#                   # Linux example: "/home/user19/KDRCC_work"
#                   # Keep this separate from the source directory.

code_root <- "/d2/home/user19/LPZ/KDRCC/simulated"
work_root <- "/d2/home/user19/LPZ/KDRCC/KDRCC_formal_20261027"

submission <- list(
  action = "run",                # "install", "check", "plan", "run", "report", "freeze"
  modules = "A",                 # After freezing: c("B", "C")
  seed = 27182818L,
  workers = 32L,                 # Independent R processes, not internal BLAS threads
  repetitions = list(A=100L, B=100L, C=100L),
  null_counts = c(8L,16L,40L,80L),
  resume = TRUE,
  retry_failed = FALSE,
  clear_stale_lock = FALSE,      # TRUE only after confirming previous processes stopped

  # All five values are N, L, M, kappa, lambda. Empty means inherit frozen values.
  # These common overrides apply to BOTH B and C; module-specific values win.
  common_fixed = list(),
  B_overrides = list(),          # Example: list(N=20000L, L=20000L, M=144L, kappa=1e-4, lambda=10)
  C_overrides = list(),

  # Optional module A changes; all other defaults remain unchanged.
  # Example: list(grid=list(N=c(1000L,5000L,20000L,50000L)))
  A_overrides = list(),

  # For less common settings, use the same nested structure as 00_parameters.R.
  # Example: list(execution=list(cache_parent="/scratch/user19/kdrcc_cache"))
  # Explicit settings above take precedence over this list.
  advanced_overrides = list(),

  # Used only for action="freeze". NULL ID accepts the eligible A recommendation.
  # freeze = list(
  #   source_run = NULL,
  #   config_id = NULL,
  #   manual_overrides = list()
  # ),
  freeze = list(
    source_run = "module_A_9ed1fac7995a",
    config_id = "A035",
    manual_overrides = list()
  ),
  # Report settings do not affect fitting. modules selects latest runs to report.
  report_run_directories = NULL, # Or full paths to specific saved module directories
  report_overrides = list(
    font="Times New Roman", dpi=300L, width_inches=12, base_cex=1.12,
    scatter_points=300L, save_pointwise_tables=TRUE,
    C_null_counts=NULL,          # NULL plots every completed null-count setting
    C_panels_per_page=6L
  )
)
# =================== END OF EDITABLE SETTINGS ===========================

# Convenience entry scripts may request one action without editing this file.
if(exists(".submission_action",inherits=FALSE))submission$action<-.submission_action
if(exists(".submission_modules",inherits=FALSE))submission$modules<-.submission_modules

.submission_file <- local({
  paths<-Filter(Negate(is.null),lapply(sys.frames(),function(x)x$ofile))
  arg<-grep("^--file=",commandArgs(FALSE),value=TRUE)
  path<-if(length(paths))tail(paths,1)[[1]] else if(length(arg))sub("^--file=","",arg[1]) else "05_submit.R"
  normalizePath(path,winslash="/",mustWork=TRUE)
})
if(is.null(code_root))code_root<-dirname(.submission_file)
code_root<-normalizePath(path.expand(code_root),winslash="/",mustWork=TRUE)
if(!file.exists(file.path(code_root,"load.R")))stop("Set code_root to the directory containing load.R.")
if(is.null(work_root))work_root<-file.path(dirname(code_root),"KDRCC_server_work")

# Command-line action/modules override the two corresponding settings only.
.file_arg<-grep("^--file=",commandArgs(FALSE),value=TRUE)
if(length(.file_arg) && identical(normalizePath(sub("^--file=","",.file_arg[1]),winslash="/",mustWork=FALSE),.submission_file)) {
  .args<-commandArgs(trailingOnly=TRUE)
  if(length(.args)>2L)stop("Usage: Rscript 05_submit.R [action] [A or B,C]")
  if(length(.args)>=1L)submission$action<-.args[1]
  if(length(.args)>=2L)submission$modules<-strsplit(.args[2],",",fixed=TRUE)[[1]]
}
source(file.path(code_root,"load.R"),encoding="UTF-8")
load_kdrcc_simulation(code_root)
submission_result <- execute_submission(submission,code_root,work_root,.submission_file)
