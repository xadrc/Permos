###

update_GFA <- 0

if (update_GFA == TRUE){
    for (iter_depth in 1:length(bh_depths)){
        depth <- bh_depths[iter_depth]
        message(paste0("running gap-filling algorithm at depth ~", depth, "m..."))
        data <- data_all[[iter_depth]]
        data <- as_data_frame(data)
        directory_GF <- paste0(directory_outputs, "/gapfilling/depth_", bh_depths[iter_depth], "m")
        dir.create(
            path = directory_GF,
            showWarnings = FALSE,
            recursive = TRUE
        )
        source(file.path(paste0(directory_scripts, "/GFAL.R")), echo = FALSE)
    }
    setwd(directory_main)
}

data_gapsfilled <- list()

for (iter_depth in 1:length(bh_depths)){
    gapsfilled <- utils::read.delim(
        file = file.path(paste0(directory_outputs, "/gapfilling/depth_", bh_depths[iter_depth], "m/output/gst_filled.txt"))
    )
    colnames(gapsfilled) <- c("loc", "date", "sim", "errnorm", "errsyst")
    data_gapsfilled[[paste(bh_depths[iter_depth])]] <- gapsfilled
}

rm(iter_depth, update_GFA, data, df, gaps, directory_GF, gapsfilled)

###

update_NAS <- TRUE

if (update_NAS == TRUE){
    message("completing series with  10% < NAs < 60%...")
    bench <- TRUE
    if (bench == TRUE){
        bench_1 <- Sys.time()
    }
    directory_tables <- paste0(directory_outputs, "/filledseries/tables")
    dir.create(paste0(directory_tables), recursive = TRUE, showWarnings = FALSE)
    do.call(file.remove, list(list.files(directory_tables, full.names = TRUE)))
    na_methods  <- c(
        "NULL",     # NO FILL
        "MED",     # MEDIAN
        # "MEAN",   # MEAN
        "LINT",     # LINEAR INTERPOLATION
        "STIN",     # STINEMAN'S ALGORITHM
        # "PINT",   # POLYNOMIAL INTERPOLATION
        "GFA"       # GAP FILLING ALGORITHM®
    )
    n_series <- 0
    for (iter_depth in 1:length(bh_depths)){
        data  <- data_all[[iter_depth]]
        depth <- bh_depths[iter_depth]
        for (iter_bh in 1:length(unique(data$loc))){
            borehole <- unique(data$loc)[iter_bh]
            df <- subset(
                data,
                data$loc == borehole,
                resetRownames = TRUE
            )
            df       <- df[!duplicated(df$date), ]
            date_min <- na.omit(df)$date[1]
            date_max <- na.omit(df)$date[nrow(na.omit(df))]
            if (!is.na(date_min) && !is.na(date_max)){
                n_nas    <- length(subset(df, df$date > date_min & df$date < date_max)$gst) - length(na.omit(subset(df, df$date > date_min & df$date < date_max)$gst))
                n_obs    <- length(subset(df, df$date > date_min & df$date < date_max)$gst)
            } else {
                n_nas    <- -1
                n_obs    <- -1
            }
            # fillNA methods
            if (n_nas/n_obs>.1 && n_nas/n_obs<.6){
                n_series <- n_series +4
                df$isNA  <- is.na(df$gst)
                message(paste0("   completing: ", borehole, " at depth: ~", depth, "m - NA method: MED"))
                df$MED     <- stats::median(df$gst, na.rm = TRUE)                 # MEAN
                # message(paste0("   completing: ", borehole, " at depth: ~", depth, "m - NA method: MEAN"))
                # df$MEAN    <- base::mean(df$gst, na.rm = TRUE)                 # MEAN
                message(paste0("   completing: ", borehole, " at depth: ~", depth, "m - NA method: LINT"))
                df$LINT    <- zoo::na.approx(df$gst, na.rm = FALSE)            # LINEAR INTERPOLATION
                message(paste0("   completing: ", borehole, " at depth: ~", depth, "m - NA method: STIN"))
                df$STIN    <- stinepack::na.stinterp(df$gst, na.rm = FALSE)    # STINEMAN'S ALG
                # message(paste0("   completing: ", borehole, " at depth: ~", depth, "m - NA method: PINT"))
                # df$PINT    <- zoo::na.spline(df$gst, na.rm = FALSE)            # POLYNOMIAL INTERPOLATION
                # fill GFA outputs
                message(paste0("   completing: ", borehole, " at depth: ~", depth, "m - NA method: GFA"))
                gfa      <- data.frame(
                    loc  = data_gapsfilled[[as.character(depth)]]$loc,
                    date = data_gapsfilled[[as.character(depth)]]$date,
                    GFA  = data_gapsfilled[[as.character(depth)]]$sim
                )
                gfa$loc  <- as.character(gfa$loc)
                gfa$date <- as.Date(gfa$date, format = "%Y-%m-%d")
                gfa <- subset(
                    x = gfa,
                    loc == borehole,
                    resetRownames = TRUE
                )
                gfa <- gfa[!duplicated(gfa$date), ]
                gfa$loc <- NULL
                df <- merge(
                    df, gfa,
                    by.x = "date",
                    all.x = TRUE
                )
                rm(t, gfa)
                data.frame()
                for (i in 1:nrow(df)){
                    if (df$isNA[i] != TRUE){
                        df$MED[i]  <- NA
                        # df$MEAN[i] <- NA
                        df$LINT[i] <- NA
                        df$STIN[i] <- NA
                        # df$PINT[i] <- NA
                        df$GFA[i]  <- NA
                    }
                    if (df$date[i] < date_min | df$date[i] > date_max){
                        df$MED[i]  <- NA
                        # df$MEAN[i] <- NA
                        df$LINT[i] <- NA
                        df$STIN[i] <- NA
                        # df$PINT[i] <- NA
                        df$GFA[i]  <- NA
                    }
                }
                saveRDS(
                    object = df,
                    file = file.path(paste0(directory_tables, "/", borehole, "_", depth, ".rds"))
                )
                rm(df)
            }
            rm(df, borehole)
        }
        rm(data)
    }
    if (bench == TRUE){
        bench_2 <- Sys.time()
        message(paste0(n_series, " series processed in ", round(difftime(bench_2, bench_1,units = "mins"), 3), " min"))
        rm(bench, bench_1, bench_2)
    }
    rm(depth, n_series)
}

rm(directory_tables, update_NAS, iter_bh, iter_depth, i, data_gapsfilled, n_nas, n_obs, date_min, date_max)

###