load_data <- function(n, p, save_dir = "spatial_data") {
  exp_dir <- file.path(save_dir, sprintf("n%d_p%d", n, p))

  # Check whether the directory exists
  if (!dir.exists(exp_dir)) {
    stop(sprintf("Data directory %s does not exist", exp_dir))
  }

  # Load all data
  load(file.path(exp_dir, "train_data.RData"))
  load(file.path(exp_dir, "test_data.RData"))
  load(file.path(exp_dir, "grid_data.RData"))
  load(file.path(exp_dir, "meta_info.RData"))

  # Return all data
  return(list(
    train = train_data,
    test = test_data,
    grid = grid_data,
    meta = meta_info
  ))
}
