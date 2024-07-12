rnd <- 2

dte <- 240627
dte <- 240709

files_all <- list.files("output")
files_use <- files_all[grep(paste0("round", rnd, ".*_dte", dte), x = files_all)]

temp <- lapply(paste0("output/", files_use), read.csv)
out <- do.call(rbind, temp)
str(out)

## clean up
# system(paste("rm", paste0("output/*.csv")))
# system(paste("rm", paste0("logs/*.txt")))

## write output
(outfile <- paste0("output/combined_round", rnd, "_all_", dte, ".csv"))
write.csv(file = outfile,
          out,
          row.names = FALSE)

list.files("output")
