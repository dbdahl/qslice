
dte <- 260518
data_now <- "db40"
data_now <- "db"

files_all <- list.files("output")
files_use <- files_all[grep(paste0(data_now, ".*_dte", dte), x = files_all)]

temp <- lapply(paste0("output/", files_use), read.csv)
out <- do.call(rbind, temp)
str(out)

## clean up
# system(paste("rm ", paste0("output/*.csv")))
# system(paste("rm ", paste0("logs/*.txt")))

write.csv(file = paste0("output/combined_", data_now, "_", dte, ".csv"),
          out,
          row.names = FALSE)

system("ls output")
