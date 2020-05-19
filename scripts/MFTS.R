###

directory_tables <- paste0(directory_outputs, "/forecasts/tables")
dir.create(paste0(directory_tables), recursive = TRUE, showWarnings = FALSE)

date_limit  <- as.Date("2015-12-31", format = "%Y-%m-%d")
date_seqnew <- seq.Date(date_seq[1], date_limit, by = "days")

###

update_TSM  <- TRUE  # !! long !!

if (update_TSM == TRUE){
    do.call(file.remove, list(list.files(directory_tables, full.names = TRUE)))
    # Bench
    bench <- TRUE
    if (bench == TRUE){
        bench_1 <- Sys.time()
    }
    # Output RMSD
    data_RMSD <- data.frame(x = NA)[,-1]
    for (i in 1:length(na_methods)){
        data_RMSD[ ,as.character(paste0(na_methods[i]))] <- NA
    }
    data_RMSD <- data_RMSD[-1,]
    message("modelling time-series...")
    n_series <- 0
    for (iter_depth in 1:length(bh_depths)){ # iterate through depths
        depth <- bh_depths[iter_depth]
        for (iter_bh in 1:length(data_meta$name)){ # iterate through boreholes at defined depth
            borehole    <- data_meta$name[iter_bh]
            if (file.exists(file.path(paste0(directory_outputs, "/filledseries/tables/", borehole, "_", depth, ".rds")))){
                df       <- readRDS(file = file.path(paste0(directory_outputs, "/filledseries/tables/", borehole, "_", depth, ".rds")))
                df       <- df[!duplicated(df$date), ]
                df_train <- subset(x = df, date <= date_limit)
                if (nrow(df_train)>0){
                    while (anyNA(tail(df_train$gst, 30))){
                        df_train <- df_train[-nrow(df_train), ]
                    }
                    df_test  <- subset(x = df, date > df_train$date[nrow(df_train)])  
                }
                if (length(df_train$gst)>1825){ # run if training set > 3yrs
                    n_obs       <- length(df_train$gst)
                    n_nas       <- as.numeric(summary(is.na(df_train$gst))[3])
                    if (n_nas/n_obs>.1 && n_nas/n_obs<.6){ # run if NAs > 10% && < 60% of the series
                        for (iter_method in 1:length(na_methods)){ # iterate through na methods
                            na_method <- na_methods[iter_method]
                            message(paste0("   retrodicting: ", borehole, " at depth: ~", depth, "m from: ", df_train$date[nrow(df_train)]+1, " to: ", df_train$date[nrow(df_train)]+1+730," - Model: TBATS - NA Method: ", na_method))
                            if (na_method == "NULL"){
                                rec  <- df$gst
                            } else {
                                rec  <- dplyr::coalesce(df[na_method][,1], df$gst)
                            }
                            rec <- data.frame(
                                date = df$date,
                                absdate = seq(1, nrow(df), 1),
                                rec  = rec
                            )
                            rec         <- na.omit(rec)
                            # Separate training et testing sets
                            rec_train   <- rec
                            rec_test    <- df
                            rec_train   <- subset(
                                x = rec,
                                date <= df_train$date[nrow(df_train)],
                                resetRownames = TRUE
                            )
                            rec_test    <- subset(
                                x = df,
                                date > df_train$date[nrow(df_train)],
                                resetRownames = TRUE
                            )
                            # Convert to time-series with f/ frequency
                            series      <- ts(rec_train$rec, frequency = 365,25)
                            # model
                            model       <- forecast::tbats(
                                y                = series,
                                use.box.cox      = FALSE,
                                use.parallel     = TRUE,
                                seasonal.periods = 365.25
                            )
                            # retrodiction / forecast
                            forecast    <- forecast::forecast(
                                object = model,
                                h      = length(seq.Date(date_limit+1, date_seq[length(date_seq)], by = "days"))
                            )
                            forecast      <- as.data.frame(forecast)
                            forecast$date <- seq.Date(
                                from = rec_train$date[nrow(rec_train)]+1, 
                                to   = rec_train$date[nrow(rec_train)]+1+730, #length(seq.Date(date_limit+1, date_seq[length(date_seq)], by = "days")),  # +2yrs from last
                                by   = "days"
                            )
                            rownames(forecast) <- NULL
                            # save model and forecast
                            base::saveRDS(
                                object = forecast,
                                file   = file.path(paste0(directory_tables, "/forecast_", unique(df$loc)[1], "_", depth, "m_", na_method, ".rds"))
                            )
                            base::save(
                                object = model,
                                file   = file.path(paste0(directory_tables, "/model_", unique(df$loc)[1], "_", depth, "m_", na_method, ".rda"))
                            )
                            # Backtesting with Root Mean Squared Deviation
                            bcktest <- merge(
                                x = data.frame(
                                    date     = forecast$date,
                                    forecast = forecast$`Point Forecast`
                                ),
                                y = data.frame(
                                    date     = rec_test$date,
                                    gst      = rec_test$gst
                                ),
                                by    = "date",
                                all.x = TRUE
                            )
                            bcktest <- na.omit(bcktest)
                            x       <- sqrt(mean((bcktest$forecast-bcktest$gst)^2, na.rm = TRUE))
                            data_RMSD[paste0(borehole, "_", depth), paste0(na_method)] <- x
                            n_series <- n_series+1
                        }
                    }
                }
            }
        }
    }
    data_RMSD["MEAN", ]   <- NA
    for (i in 1:length(na_methods)){
        data_RMSD["MEAN",   as.character(paste0(na_methods[i]))] <- base::mean(as.vector(data_RMSD[ ,paste0(na_methods[i])])[1:nrow(data_RMSD)-1], na.rm = TRUE)
    }
    data_RMSD["MEDIAN", ] <- NA
    for (i in 1:length(na_methods)){
        data_RMSD["MEDIAN", as.character(paste0(na_methods[i]))] <- stats::median(as.vector(data_RMSD[ ,paste0(na_methods[i])])[1:(nrow(data_RMSD)-2)], na.rm = TRUE)
    }
    base::saveRDS(
        object = data_RMSD,
        file = file.path(paste0(directory_tables, "/RMSD.rds"))
    )
    if (exists("bench")){
        if (bench == TRUE){
            bench_2 <- Sys.time()
            message(paste0(n_series, " series forecasted in ", round(difftime(bench_2, bench_1,units = "mins"), 3), " min"))
            rm(bench, bench_1, bench_2)
        }
    }
    rm(
        iter_depth, iter_bh, iter_method, update_TSM, depth, borehole, df, 
        n_series, n_obs, n_nas, n_days, df_train, date_seqnew, na_method, 
        rec, rec_test, rec_train, series, model, forecast, bcktest, x, i
    )
}

rm(directory_tables)

###