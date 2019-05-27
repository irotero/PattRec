list.of.packages <- c("outliers", "xlsx", "WriteXLS", "RMySQL", "gdata", "tools", "random")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cloud.r-project.org/")