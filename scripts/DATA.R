###

message("importing series...")

bench <- TRUE
if (bench == TRUE){
    bench_1 <- Sys.time()
}

data_meta <- read.csv(
    file = file.path(paste0(directory_data, "/meta.csv")),
    header = TRUE,
    sep = ";"
)
data_meta$name <- as.character(data_meta$name)
data_meta$site <- as.character(data_meta$site)

date_start      <- as.Date("2008-01-01", format = "%Y-%m-%d")
date_end        <- as.Date("2017-12-31", format = "%Y-%m-%d")
date_seq        <- seq(date_start, date_end, by = "day")

bh_depths       <- c(0, 5, 10, 15, 20)
bh_selection    <- sort(data_meta$name)
bh_locations    <- rep(NA, length(bh_selection))

for (i in 1:length(bh_selection)){
    bh_locations[i] <- data_meta$site[data_meta$name == bh_selection[i]]
}
bh_param <- data.frame(
    ID   = bh_selection,
    name = bh_locations
)

data_all <- list()
n_files  <- length(list.files(path = paste0(directory_data, "/boreholes/")))
n_series <- 0

for (iter_depth in 1:length(bh_depths)){
    depth <- bh_depths[iter_depth]
    for (iter_bh in 1:nrow(data_meta)){
        borehole <- data_meta$name[iter_bh]
        bh <- read.csv(
            file = file.path(paste0(directory_data, "/boreholes/", borehole, ".csv")),
            sep = ",", 
            header = TRUE,
            check.names = FALSE
        )
        loggers     <- as.numeric(colnames(bh)[-1])
        depthlog    <- as.numeric(colnames(bh)[as.numeric(which.min(abs(loggers-depth)))+1])
        if (depthlog > depth-2.5 && depthlog < depth+2.5){
            message(paste0("   uploading: ", borehole, " at depth: ~", depth, "m [", depthlog, "m] from: ", borehole, ".csv"))
            df <- data.frame(
                loc = as.character(paste(borehole)),
                depth = depthlog,
                date = as.Date(
                    bh$Time,
                    format = "%Y-%m-%d"
                ),
                gst = bh[, as.character(depthlog)]
            )
            df$loc <- as.character(df$loc)
            rm(loggers, depthlog)
            df <- df[!duplicated(df$date), ]
            df <- subset(
                x = df,
                date >= date_start & date <= date_end,
                resetRownames = TRUE
            )
            t <- data.frame(
                date = seq.Date(
                    from = df$date[1],
                    to = df$date[nrow(df)],
                    by = "days"
                )
            )
            df <- merge(
                t, df, 
                by    = "date",
                all.x = TRUE
            )
            if (iter_bh == 1){
                data <- df
            } else {
                data <- rbind(data, df)
            }
            n_series <- n_series +1
        }
    }
    data_all[[iter_depth]] <- data
}
rm(df, bh, data, iter_bh, iter_depth, borehole, t)

# message(cat(green("\n   • time span:")))
huxbh_timespans <- huxtable::hux(
    names = c("Start date", "End date"),
    dates = c(date_start, date_end),
    add_colnames = FALSE
)
# huxtable::print_screen(
#     huxbh_timespans,
#     compact = TRUE,
#     colnames = FALSE
# )
rm(huxbh_timespans)
# message(cat(green("\n   • selected boreholes:")))
huxbh_param <- huxtable::hux(
    ID = bh_selection,
    name = bh_locations,
    add_colnames = FALSE
)
huxtable::bold(huxbh_param)[1, ]             <- FALSE
huxtable::bottom_border(huxbh_param)[1, ]    <- FALSE
# huxtable::print_screen(
#     huxbh_param,
#     compact = TRUE,
#     colnames = FALSE
# )
rm(huxbh_param)
# message(cat(green("\n   • at depths:")))
huxbh_depths <- huxtable::hux(
    depth = bh_depths,
    unit = rep("m", length(bh_depths)),
    add_colnames = FALSE
)
huxtable::align(huxbh_depths) <-  'left'
# huxtable::print_screen(
#     huxbh_depths,
#     compact = TRUE,
#     colnames = FALSE
# )
rm(huxbh_depths)
# message("\n")

if (bench == TRUE){
    bench_2 <- Sys.time()
    message(paste0(n_series, " series imported in ", round(difftime(bench_2, bench_1,units = "mins"), 3), " min"))
    rm(bench, bench_1, bench_2)
}
rm(n_files, n_series, bh_param)

###