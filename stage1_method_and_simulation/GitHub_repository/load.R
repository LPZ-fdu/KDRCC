load_kdrcc_simulation <- function(root) {
  if (getRversion() < "3.6.3") stop("R 3.6.3 or newer is required.")
  # Match R >= 4's character-column default on R 3.6, including PSOCK workers.
  if (getRversion() < "4.0.0") options(stringsAsFactors = FALSE)
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  for (f in sort(list.files(file.path(root, "R"), pattern = "[.]R$", full.names = TRUE)))
    source(f, local = .GlobalEnv, encoding = "UTF-8", keep.source = FALSE)
  source(file.path(root, "config", "mechanism.R"), local = .GlobalEnv, encoding = "UTF-8")
  options(kdr.root = root)
  invisible(root)
}
