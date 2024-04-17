target <- "normal"
target <- "gamma"
target <- "gammalog"
target <- "igamma"
target <- "igammalog"

rnd <- 2
dte <- 240308
dte <- 240309

files_all <- list.files("output")
# files_use <- files_all[grep(paste0("round", rnd, "_target", target, ".*_dte", dte), x = files_all)]
files_use <- files_all[grep(paste0("round", rnd, ".*_dte", dte), x = files_all)]

temp <- lapply(paste0("output/", files_use), read.csv)
out <- do.call(rbind, temp)
str(out)

# write.csv(file = paste0("output/combined_round", rnd, "_target", target, "_", dte, ".csv"),
#           out,
#           row.names = FALSE)
write.csv(file = paste0("output/combined_round", rnd, "_all_", dte, ".csv"),
          out,
          row.names = FALSE)

## clean up
system(paste("rm", paste0("output/", files_use, collapse = " ")))
system(paste("rm", paste0("logs/*.txt")))
