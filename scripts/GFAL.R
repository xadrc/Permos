
# add output directories
dir.create(paste0(directory_GF,"/test/"), showWarnings = FALSE)
dir.create(paste0(directory_GF,"/output/plots/gaps/"), showWarnings = FALSE, recursive = TRUE)

# define output folder for results (gap-filled values and summary plots)
outputfolder <- paste0(directory_GF,"/output/")

# The following source file will load all libraries and functions, please adapt the FILE accordingly
#  source file 'GST_gapfilling_functions.R' in supplementary material
# source('functions.R')

# load sample data set with selected PERMOS GST time series from supplementary material folder
# load('sample_data.RData')

setwd(outputfolder)
bench_1 <- Sys.time()

 
# ---- Define parameters -------------------------------------------------------

#   date_start   <- as.Date("2005-01-01")    # min-date for data set
#   date_end     <- as.Date("2015-12-31")    # max-date for data set
date_start <- as.Date(date_start)
date_end   <- as.Date(date_end)

maxlen_shortgap <- 3 # max. gap-duration to treat as shortgap (relevant also for removing slivers for long gaps)


# ---- Basic data processing ---------------------------------------------------

df <- fillNA(
 	df = data,
 	loc = data$loc,
 	date = data$date,
 	tolerance = 150
 ) 

# ---- Detect gaps and create indices ------------------------------------------

gaps <- gapdetection(df=df,maxgapdur = 500)
  
for (i in unique(gaps$index) ){
    g <- gaps[gaps$index==i,]
    if(exists('gaptable')) {
        gaptable <- rbind(gaptable,data.frame(loc=g$loc,index=g$index,date=seq(g$start,g$end,by="day")))
    } else {
        gaptable <- data.frame(loc=g$loc,index=g$index,date=seq(g$start,g$end,by="day"))
    }
}
rm(i)

df <- as.data.frame(left_join(df,gaptable,by=c("loc","date"))) 
df$index <- as.character(df$index)

# remove unused things
# rm(g,i,gaptable,gst,gt)
rm(g,i,gaptable)

n_gap <- 0

for(i in 1:length(gaps$index)){
    gapname <- gaps$index[i]
    g <- gaps[gaps$index == gapname, ]  # gap-metadata
    y <- df %>% filter(loc == g$loc) %>% dplyr::select(loc, date, gst)  # entire target logger
    # GAPS
    if (g$duration <= maxlen_shortgap) {
        # this is a short gap -> LI_NN
        t <- data.frame(fill_LI(x = y, info = g))
        t$error_norm <- t$error
        t$error_syst <- 0
        sem_gap <- round(mean(t$error, na.rm = T), 2)
        gapstats <-
            data.frame(
                g,
                method = "LI",
                win = NA,
                otherloc = NA,
                SEM = sem_gap
            )
    } else {
        # this is a long gap
        t <- fill_QM(
                df = df,
                x = y,
                info = g,
                win = 15,
                mindur = 30,
                minobs = 3,
                minreg = 5
            )
        t$error_norm <- ifelse(t$normdist == 1, t$error, 0)
        t$error_syst <- ifelse(t$normdist  < 1, t$error, 0)
        sem_gap <- ifelse(
            sum(t$normdist) < 1,
            round(mean(t$error_syst, na.rm = T), 3),
            round(sqrt((mean(t$error_norm, na.rm = T) ^ 2 / sum(t$normdist)) + # stochastic errors
                        mean(t$error_syst, na.rm = T) ^ 2 # systematic errors
                   ), 3))
        gapstats <- data.frame(
                g,
                method = "QM",
                win = 15,
                otherloc = first(t$otherloc),
                SEM = sem_gap
            )
    }
    if (is.finite(first(t$sim))) {
        gapdat <-
            data.frame(
                loc = g$loc,
                date = t$date,
                sim = t$sim,
                error_norm = t$error_norm,
                error_syst = t$error_syst
            )
        #write.table(gapdat,file=gsub(" ","",paste(outputfolder,"gst_filled.txt")),quote=FALSE,sep="\t",na="",row.names=FALSE,col.names=FALSE,append=TRUE)
        write.table(
            gapdat,
            file = "gst_filled.txt",
            quote = FALSE,
            sep = "\t",
            na = "",
            row.names = FALSE,
            col.names = FALSE,
            append = TRUE
        )
        #write.table(gapstats,file=gsub(" ","",paste(outputfolder,"gapstats.txt")),quote=FALSE,sep="\t",na="",row.names=FALSE,col.names=FALSE,append=TRUE)
        write.table(
            gapstats,
            file = "gapstats.txt",
            quote = FALSE,
            sep = "\t",
            na = "",
            row.names = FALSE,
            col.names = FALSE,
            append = TRUE
        )
        title <- paste("Gap:", i, "| Method:", gapstats$method, "| otherloc:", gapstats$otherloc, "| SEM:", sem_gap, "K | Bias rMAGST:", round(g$duration / 365 * sem_gap, 2), "K")
        merge <- left_join(y, gapdat[, c("date", "sim", "error_norm", "error_syst")], by = "date")
        extend_plot <- ifelse(
            g$duration < 30,
            ceiling(30 - g$duration) / 2,
            round(g$duration * 2.5, 0)
        )
        # PLOTS
        if (g$duration <= maxlen_shortgap) {
            # make plot for short gaps
            min_y <- floor(min(y$gst, na.rm = T) - 5)
            max_y <- ceiling(max(y$gst, na.rm = T) + 5)
            p <-
                suppressWarnings(
                    qplot(
                        data = merge,
                        x = date,
                        y = gst,
                        xlab = "",
                        ylab = "GST (?C)",
                        main = title,
                        geom = "line"
                    ) +
                    geom_errorbar(
                        aes(
                            ymin = merge$sim - merge$error_syst - merge$error_norm,
                            ymax = merge$sim + merge$error_syst + merge$error_norm
                        ),
                        colour = "orange"
                    ) +
                    geom_point(y = merge$sim, colour =
                                   "red4") +
                    scale_y_continuous(
                        limits = c(min_y, max_y),
                        minor_breaks = NULL
                    ) +
                    scale_x_date(
                        limits = c(
                            as.Date(g$start - extend_plot),
                            as.Date(g$end + extend_plot)
                        ),
                        minor_breaks = NULL
                    ) +
                    theme_bw()
                )
        } else {
            # for lager gaps
            min_y <-
                floor(min(gapdat$sim, na.rm = T) - max(c(
                    gapdat$error_norm, gapdat$error_syst
                ), na.rm = T) - 2)
            max_y <-
                ceiling(max(gapdat$sim, na.rm = T) + max(c(
                    gapdat$error_norm, gapdat$error_syst
                ), na.rm = T) + 2)
            p <-
                suppressWarnings(
                    qplot(
                        data = merge,
                        x = date,
                        y = gst,
                        xlab = "",
                        ylab = "GST (?C)",
                        main = title,
                        geom = "line"
                    ) +
                    geom_ribbon(
                        aes(
                            ymin = merge$sim - merge$error_syst,
                            ymax = merge$sim + merge$error_syst
                        ),
                        fill = "blue",
                        alpha = 0.5
                    ) +
                    geom_ribbon(
                        aes(
                            ymin = merge$sim - merge$error_norm,
                            ymax = merge$sim + merge$error_norm
                        ),
                        fill = "orange",
                        alpha = 0.3
                    ) +
                    geom_line(y = merge$sim, colour =
                                  "red4") +
                    scale_y_continuous(
                        limits = c(min_y, max_y),
                        minor_breaks = NULL
                    ) +
                    # scale_x_date(limits=c(as.Date(g$start-extend_plot),as.Date(g$end+extend_plot)), minor_breaks=NULL ) +
                    scale_x_date(
                        limits = c(g$start - extend_plot, g$end + extend_plot),
                        minor_breaks = NULL
                    ) +
                    theme_bw()
                )
        }
        suppressWarnings(ggsave(
            p,
            file = paste0(outputfolder, "/plots/gaps/", gapname, ".pdf"),
            width = 16,
            height = 8
        ))
        n_gap <- n_gap +1
    } # else {
    #    message(paste0("   gap: ", gapname, " not filled!"))
    # }
    #gapstats <- gapstats ; saveRDS(object = gapstats, file = "gapstats.rds")
    #gapdat <- gapdat ; saveRDS(object = gapdat, file = "gapdat.rds")
    rm(g, y, t, title, merge, gapstats, gapdat)
}

bench_2 <- Sys.time()
message(paste0(n_gap, " gaps filled in ", round(difftime(bench_2, bench_1,units = "mins"), 3), " min"))

setwd(directory_main)
rm(outputfolder, extend_plot, i, max_y, maxlen_shortgap, min_y, sem_gap, df, gaps, p, gapname, bench_1, bench_2, n_gap)
