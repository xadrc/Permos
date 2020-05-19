###

directory_tables <- paste0(directory_outputs, "/forecasts/tables")
directory_plots  <- paste0(directory_outputs, "/forecasts/plots")
dir.create(paste0(directory_plots),  recursive = TRUE, showWarnings = FALSE)

date_limit  <- as.Date("2015-12-31", format = "%Y-%m-%d")
date_seqnew <- seq.Date(date_seq[1], date_limit, by = "days")

do.call(file.remove, list(list.files(directory_plots, full.names = TRUE)))

n_plot <- 0

###

plot_MFTS <- TRUE

if (plot_MFTS == TRUE){
    message("plotting forecasted series...")
    if (!exists("data_RMSD")){
        data_RMSD <- readRDS(
            file = file.path(paste0(directory_tables, "/RMSD.rds"))
        )
    }
    for (iter_depth in 1:length(bh_depths)){
        depth <- bh_depths[iter_depth]
        for (iter_bh in 1:length(data_meta$name)){
            borehole <- data_meta$name[iter_bh]
            for (iter_method in 1:length(na_methods)){
                na_method <- na_methods[iter_method]
                RMSD      <- data_RMSD[paste0(borehole, "_", depth), paste0(na_method)]
                if (!is.na(RMSD) && file.exists(file.path(paste0(directory_tables, "/forecast_", borehole, "_", depth, "m_", na_method, ".rds")))){
                    # reload data sets & models
                    forecast <- base::readRDS(
                        file = file.path(paste0(directory_tables, "/forecast_", borehole, "_", depth, "m_", na_method, ".rds"))
                    )
                    colnames(forecast) <- c("forecast", "Lo80", "Hi80", "Lo95", "Hi95", "date")
                    load(file = file.path(paste0(directory_tables, "/model_", borehole, "_", depth, "m_", na_method, ".rda")))
                    df <- base::readRDS(file = file.path(paste0(directory_outputs, "/filledseries/tables/", borehole, "_", depth, ".rds")))
                    # sep train & test sets
                    df       <- df[!duplicated(df$date), ]
                    df_train <- subset(x = df, date <= date_limit)
                    while (is.na(df_train$gst[nrow(df_train)])){
                        df_train <- df_train[-nrow(df_train),]
                    }
                    df_test  <- subset(x = df, date > df_train$date[nrow(df_train)])
                    # plot
                    message(paste0("   plotting: ", borehole, " at depth: ~", depth, "m - NA method: ", na_method))
                    {
                        quartz(width = 5, height = 3, pointsize = 11)
                        par(omi = c(rep(0, 4)))
                        # width param
                        xlim <- c(date_seq[1], date_seq[length(date_seq)]+1)
                        ylim <- extendrange(
                            c(min(c(min(forecast$`Lo 95`, na.rm = TRUE), min(df$gst, na.rm = TRUE))), max(c(max(forecast$`Hi 95`, na.rm = TRUE), max(df$gst, na.rm = TRUE)))),
                            f = .05
                        )
                        ylim[1] <- ylim[1] - abs(4/2 - abs(diff(ylim))/2)
                        ylim[2] <- ylim[2] + abs(4/2 - abs(diff(ylim))/2)
                        ticks   <- seq.Date(from = as.Date(date_start), to = as.Date(date_end)+1, by = "year")
                        ticks   <- as.Date(ticks, format = "%Y")
                        # titles
                        title1  <- paste0(unique(df$loc)[1], " | Depth: ", unique(df$depth)[1], "m | NA_method: ", na_method)
                        title2  <- paste0("RMSD: ", round(RMSD, 3))
                        sub     <- paste0(model[["method"]], "(", ifelse(model[["parameters"]][["control"]][["use.box.cox"]] != FALSE, 0, 1), ", {", ifelse(is.null(model[["p"]]), "-", model[["p"]]), ",", ifelse(is.null(model[["q"]]), "-", model[["q"]]), "}, ", ifelse(is.null(model[["damping.parameter"]]), "-", round(model[["damping.parameter"]], 3)), ", {<", ifelse(is.null(model[["seasonal.periods"]]), "-", model[["seasonal.periods"]]), ",", ifelse(is.null(model[["k.vector"]]), "-", model[["k.vector"]]), ">})")
                        # draw window
                        par(mar = c(2, 2, 2, 0), pty = "m")
                        plot(
                            0, 0,
                            type = "n",
                            xlim = xlim, ylim = ylim,
                            axes = FALSE,
                            xlab = "", ylab = ""
                        )
                        title(
                            main = c(title1, sub),
                            adj = 0, cex.main = c(.8),
                            col = c("black", "grey80")
                        )
                        title(
                            main = c(NA, title2),
                            adj = 1, cex.main = .9
                        )
                        axis.Date(
                            side = 1,
                            date_seq,
                            format = "%Y",
                            at = ticks,
                            cex.axis = .8
                        )
                        axis(
                            side = 2,
                            cex.axis = .8
                        )
                        abline(
                            v=ticks,
                            h= pretty(extendrange(ylim)),
                            lty = 3, col = "lightgrey"
                        )
                        box(lwd = .15)
                        # plot observations
                        lines(
                            df_train$date, df_train$gst,
                            type = "l", lwd = 1, 
                            col = "black"
                        )
                        # plot replaced missing values
                        if (na_method != "NULL"){
                            NAs <- df_train[, paste0(na_method)]
                            lines(
                                df_train$date, NAs,
                                type = "l", lwd = 1,
                                col = "orange"
                            )
                            rm(NAs)
                        }
                        # plot forecasts intervals
                        segments(
                            x0 = forecast$date, y0 = forecast$Lo95,
                            x1 = forecast$date, y1 = forecast$Hi95,
                            lwd = .1,
                            col = grDevices::rgb(0/255, 120/255, 255/255, alpha = .1)
                        )
                        segments(
                            x0 = forecast$date, y0 = forecast$Lo80,
                            x1 = forecast$date, y1 = forecast$Hi80,
                            lwd = .1,
                            col = grDevices::rgb(0/255, 120/255, 255/255, alpha = .2)
                        )
                        # plot observations
                        lines(
                            df_test$date, df_test$gst,
                            type = "l", lwd = 1,
                            col = "red"
                        )
                        # plot forecasts
                        lines(
                            forecast$date, forecast$forecast,
                            type = "l", lwd = 1,
                            col = "blue"
                        )
                        # legend
                        legend(
                            x = "bottomleft",
                            legend = c(
                                "observations",
                                "interpolated values"
                            ),
                            pch = rep(NA, 2),
                            lty = rep(1,  2),
                            col = c(
                                "black",
                                "orange"
                            ),
                            bty = "n", 
                            horiz = FALSE,
                            cex = .8
                        )
                        legend(
                            x = "bottom",
                            legend = c(
                                "retrodiction", 
                                "observations"
                            ),
                            pch = rep(NA, 2),
                            lty = rep(1,  2),
                            col = c(
                                "blue",
                                "red"
                            ),
                            bty = "n", 
                            horiz = FALSE,
                            cex = .8
                        )
                        legend(
                            x = "bottomright",
                            legend = c(
                                "80% confidence interval",
                                "95% confidence interval"
                            ),
                            pch = rep(15, 2),
                            col = c(
                                grDevices::rgb(0/255, 120/255, 255/255, alpha = .3),
                                grDevices::rgb(0/255, 120/255, 255/255, alpha = .2)
                            ),
                            bty = "n", 
                            horiz = FALSE,
                            cex = .8
                        )
                    }
                    # save plots
                    quartz.save(
                        file = file.path(paste0(directory_plots, "/", borehole, "_", depth, "_", na_method, ".jpg")),
                        type = "jpg", dpi = 300
                    )
                    n_plot <- n_plot +1
                    dev.off()
                }
            }
        }
    }
    rm(forecast, df, df_test, df_train, RMSD, na_method, depth, 
       borehole, iter_depth, iter_bh, iter_method, data_RMSD, model)
    rm(xlim, ylim, ticks, title1, title2, sub, i)
}

###

plot_RMSD <- TRUE

if (plot_RMSD == TRUE){
    if (!exists('data_RMSD')){
        data_RMSD <- readRDS(file = file.path(paste0(directory_tables, "/RMSD.rds")))
    }
    data_RMSD <- na.omit(data_RMSD)
    quartz(
        width = 5, height = 6, pointsize = 11,
        title = "RMSD"
    )
    par(
        mfrow = c(length(na_methods), length(bh_depths)),
        omi = c(2/3, .2, .2, 0.05)
    )
    for(iter_method in 1:length(na_methods)){
        na_method <- na_methods[iter_method]
        for (iter_depth in 1:length(bh_depths)){
            depth <- bh_depths[iter_depth]
            par(mar = c(0, 2, 1, 0), pty = "m")
            if (iter_method !=1 && iter_depth !=1){
                par(new = FALSE)
            }
            data_RMSD$select <- stringr::str_sub(rownames(data_RMSD), -2) == stringr::str_sub(paste0("_", depth), -2)
            data <- subset(
                x = data_RMSD,
                data_RMSD$select == TRUE,
                resetRownames = TRUE
            )
            # data$select <- NULL
            ticks_values <- seq(1, nrow(data), 1)
            ticks_names  <- stringr::str_sub(rownames(data), end = -str_count(paste0("_", depth))-1)
            xlim = range(ticks_values)
            ylim = c(0, ceiling(max(data)))
            # if (ylim[2]%%2 !=0){
            #     ylim[2] <- ylim[2]+1
            # }
            {
                plot(
                    0, 0,
                    type = "n",
                    xlim = xlim, ylim = ylim,
                    axes = FALSE,
                    xlab = "", ylab = ""
                )
                axis(
                    side = 2,
                    cex.axis = .8
                )
                if (iter_method == length(na_methods)){
                    axis(
                        side = 1,
                        cex.axis = .8,
                        labels = ticks_names, 
                        at = seq(xlim[1], xlim[length(xlim)], 1),
                        las = 2
                    )
                } else {
                    axis(
                        side = 1,
                        cex.axis = 0.0001,
                        labels = NULL, 
                        at = seq(xlim[1], xlim[length(xlim)], 1)
                    )
                }
                abline(
                    h = pretty(extendrange(ylim)),
                #   v = ticks_values,
                    lty = 3, col = "lightgrey"
                )
                box(lwd = .3)
                lines(
                    x = ticks_values, y = data[,paste0(na_method)],
                    type = "h"
                )
                abline(
                    h = median(data[,paste0(na_method)], na.rm = TRUE),
                    lty = 2, lwd = 1, col = "red"
                )
                par(new = TRUE)
                plot(
                    0, 0,
                    type = "n",
                    xlim = c(0,1), ylim = c(0,1),
                    axes = FALSE,
                    xlab = "", ylab = ""
                )
                # text(
                #     x = 1,
                #     y = .95,
                #     labels = paste0("med = ", round(median(data[,paste0(na_method)], na.rm = TRUE), 3)),
                #     adj = 1,
                #     col = "red"
                # )
                if (iter_depth == 1){
                    mtext(
                        text = paste0("method: ", na_method),
                        side = 2,
                        padj = 0,
                        line = 2,
                        cex  = 3/4
                    )
                }
                if (iter_method == 1){
                    mtext(
                        text = paste0("depth: ~", depth, "m"),
                        side = 3,
                        padj = 0,
                        line = 1/2,
                        cex  = 3/4
                    )
                }
            }
        }
        data_RMSD["select"] <- NULL
    }
    quartz.save(
        file = file.path(directory_outputs, "REPORT_RMSD.jpg")),
        type = "jpg", dpi = 300
    )
    # dev.off()
    n_plot <- n_plot +1
    data_RMSD <- readRDS(file = file.path(paste0(directory_tables, "/RMSD.rds")))
}

rm(plot_MFTS, data, depth, na_method, ticks_names, ticks_values, xlim, ylim, iter_depth, iter_method)

###

message("plots saved in ", directory_plots)
rm(n_plot, directory_plots)

###