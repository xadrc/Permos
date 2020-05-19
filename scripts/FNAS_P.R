###

directory_tables <- paste0(directory_outputs, "/filledseries/tables")
directory_plots  <- paste0(directory_outputs, "/filledseries/plots")
dir.create(paste0(directory_plots),  recursive = TRUE, showWarnings = FALSE)

###

plot_fNA <- TRUE

if (plot_fNA == TRUE){
    do.call(file.remove, list(list.files(directory_plots, full.names = TRUE)))
    message("plotting filled series...")
    n_plot <- 0
    for (iter_depth in 1:length(bh_depths)){
        depth <- bh_depths[iter_depth]
        for (iter_bh in 1:length(data_meta$name)){
            borehole <- data_meta$name[iter_bh]
            if (file.exists(file.path(paste0(directory_tables, "/", borehole, "_", depth, ".rds")))){
                df       <- base::readRDS(file = file.path(paste0(directory_tables, "/", borehole, "_", depth, ".rds")))
                for (iter_method in 1:length(na_methods)){
                    na_method <- na_methods[iter_method]
                    # data prep
                    {
                        if (na_method == "NULL"){
                            rec <- df$gst
                        } else {
                            gst <- df$gst
                            nas <- df[na_method][,1]
                            rec <- dplyr::coalesce(nas, gst)
                        }
                        rec <- data.frame(
                            date = df$date,
                            rec  = rec
                        )
                        rec <- na.omit(rec)
                        ymean <- data.frame(
                            date  = lubridate::year(as.Date(row.names(as.data.frame(xts::apply.yearly(zoo::zoo(rec$rec, rec$date), base::mean))))),
                            ymean = as.data.frame(xts::apply.yearly(zoo::zoo(rec$rec, rec$date), base::mean))[,1]
                        )
                        ymean <- merge(
                            data.frame(date = seq(ymean$date[1], ymean$date[nrow(ymean)], 1)), ymean,
                            by = "date", 
                            all.x = TRUE
                        )
                        utils::write.csv(
                            x = data.frame(
                                year = ymean$date,
                                mean = round(ymean$ymean, 3)
                            ),
                            row.names = FALSE, col.names = TRUE,
                            file = file.path(paste0(directory_tables, "/", borehole, "_", depth, "_", na_method, ".csv"))
                        )
                        t <- data.frame(date = seq.Date(rec$date[1], rec$date[nrow(rec)], by = "day"))
                        ymean$date <- as.Date(lubridate::ymd(sprintf("%d-01-01", ymean$date)), format = "%Y-%m-%d")
                        ymean[is.na(ymean)] <- 0
                        ymean <- merge(
                            t, ymean,
                            by = "date",
                            all.x = TRUE
                        )
                        ymean$ymean <- zoo::na.locf(ymean$ymean, na.rm = FALSE, fromLast = FALSE)
                    }
                    # plot
                    message(paste0("   plotting: ", borehole, " at depth: ~", depth, "m - NA method: ", na_method))
                    {
                        quartz(width = 5, height = 3, pointsize = 11)
                        par(omi = c(rep(0, 4)))
                        # width param
                        xlim <- c(date_seq[1], date_seq[length(date_seq)]+1)
                        ylim <- extendrange(c(min(c(min(df$gst, na.rm = TRUE), min(rec$rec, na.rm = TRUE))), max(c(max(df$gst, na.rm = TRUE), max(rec$rec, na.rm = TRUE)))), f = .05)
                        ylim[1] <- ylim[1] - abs(4/2 - abs(diff(ylim))/2)
                        ylim[2] <- ylim[2] + abs(4/2 - abs(diff(ylim))/2)
                        ticks   <- seq.Date(from = as.Date(date_start), to = as.Date(date_end)+1, by = "year")
                        ticks   <- as.Date(ticks, format = "%Y")
                        # titles
                        title = paste0(unique(df$loc)[1], " | Depth: ", unique(df$depth)[1], "m | NA_method: ", na_method)
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
                            main = title,
                            adj = 0, cex.main = .8,
                            col = "black"
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
                        lines(
                            ymean$date, ymean$ymean,
                            type = "h", lwd = .2,
                            col = "grey90"
                        )
                        abline(
                            v=ticks,
                            h= pretty(extendrange(ylim)),
                            lty = 3, col = "lightgrey"
                        )
                        box(lwd = .15)
                        # plot observations
                        lines(
                            df$date, df$gst,
                            type = "l", lwd = 1,
                            col = "black"
                        )
                        # plot replaced missing values
                        if (na_method != "NULL"){
                            lines(
                                df$date, nas,
                                type = "l", lwd = 1, 
                                col = "orange"
                            )
                        }
                        # legend
                        legend(
                            x = "bottomleft",
                            legend = c(
                                "observations",
                                "replaced values"
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
                            x = "bottomright",
                            legend = "yearly means",
                            pch = 15,
                            col = "grey90",
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
    message("plots saved in ", directory_plots)
    rm(depth, borehole, df, na_method, rec, gst, nas, ymean, t)
    rm(xlim, ylim, ticks, title, n_plot)
}

rm(directory_tables, directory_plots, plot_fNA)

###