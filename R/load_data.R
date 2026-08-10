load_data <- function(n, p, save_dir = "spatial_data") {
  exp_dir <- file.path(save_dir, sprintf("n%d_p%d", n, p))
  
  # 检查目录是否存在
  if (!dir.exists(exp_dir)) {
    stop(sprintf("数据目录 %s 不存在", exp_dir))
  }
  
  # 加载所有数据
  load(file.path(exp_dir, "train_data.RData"))
  load(file.path(exp_dir, "test_data.RData"))
  load(file.path(exp_dir, "grid_data.RData"))
  load(file.path(exp_dir, "meta_info.RData"))
  
  # 返回所有数据
  return(list(
    train = train_data,
    test = test_data,
    grid = grid_data,
    meta = meta_info
  ))
}