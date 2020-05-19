message("importing functions...")

# ---- BASIC FUNCTIONS ----

# Function that returns Root Mean Squared Error
RMSE <- function(error) sqrt( mean(error^2, na.rm=TRUE) ) 

# computes the standard error of the mean for the difference vector error 
SEM <- function(error) sqrt( var(error,na.rm=TRUE)/length(na.omit(error)) )

# to make regular time-series with NA where values are missing
fillNA <- function(df=df,loc=loc,date=date,tolerance=30) {  # tolerance: difference in days from October 1st
  df2 <- df %>% group_by(loc) %>% summarize(mindate=min(date),maxdate=max(date))
  for(i in df2$loc) {
    x <- df2[df2$loc==i,]
    mindate <- as.Date( ifelse( month(x$mindate) == 10 & day(x$mindate)-tolerance > 1, 
                                as.Date(gsub(" ","",paste(year(x$mindate),"-10-01"))), x$mindate),origin='1970-01-01' )
    maxdate <- as.Date( ifelse( month(x$maxdate) == 9 & day(x$maxdate)+tolerance < 30, 
                                as.Date(gsub(" ","",paste(year(x$maxdate),"-09-30"))), x$maxdate),origin='1970-01-01' )
    if(exists('dummydates')) {
      dummydates <- rbind(dummydates,data.frame(loc=as.character(i),date=seq(mindate,maxdate,by="day")))
    } else {
      dummydates <- data.frame(loc=as.character(i),date=seq(mindate,maxdate,by="day"))
    }
  }
  return(left_join(dummydates,df,by=c("loc","date")))
  rm(dummydates,df2,x,i)
}

fill_NA <- function(sample_data){
    sample_data = as.data.frame(sample_data)
    boreholes <- unique(sample_data$loc)
    for (i in 1:length(boreholes)){
        df <- subset(
            x = sample_data,
            loc == boreholes[i],
            resetRownames = TRUE
        )
        t <- seq(df$date[1], df$date[nrow(df)], by = "day")
        t <- data.frame(date = t)
        df <- merge(
            t, df,
            by.x = "date",
            all.x = TRUE
        )
        if (i == 1){
            df2 <- df
        } else {
            df2 <- rbind(df2, df)
        }
    }
    df <- data.frame(
        loc = df2$loc,
        date = df2$date,
        gst = df2$gst
    )
    df$loc <- as.character(df$loc)
    df$date <- as.Date(df$date, format = "%Y-%m-%d")
    rm(boreholes, t, df2)
    return(df)
}

# ---- AUTOMATIC GAP DETECTION ----

gapdetection <- function(df=df, maxgapdur=500){
  gaps <- df[!is.finite(df$gst),c("loc","date")]  %>% 
    group_by(loc)%>% arrange(date) %>% mutate(
      start = as.Date( ifelse(date==min(date) | (date-lag(date,1))>1 , date, NA), origin='1970-01-01' ),
      end = as.Date( ifelse(date==max(date) | (lead(date,1)-date)>1 , date, NA), origin='1970-01-01' )
    ) %>% filter(start>1|end>1) 
  gaps <- gaps %>% dplyr:::select(loc,start,end) 
  
  gaps <- data.frame( 
    loc=gaps[is.finite(gaps$start),1], 
    start=gaps[is.finite(gaps$start),2], 
    end=gaps[is.finite(gaps$end),3] )
  
  gaps <- gaps %>% group_by(loc) %>% arrange(start) %>% mutate(
    duration = as.numeric(1+(end-start)),
    index = gsub(' ','_',paste(loc,seq(1,100,1)[1:length(duration)]))
  ) %>% filter(duration<=maxgapdur)
  
  return(gaps)  # return data frame with gap-metadata
  rm(g,i,gaps)
}

# ---- LINEAR INTERPOLATION (SHORT GAPS) ----

fill_LI <- function(x=x,info=info) {  
  x <- x %>% filter((date>=(info$start-ceiling(info$duration/2)) & 
                    date<(info$start+ceiling(info$duration/2)+info$duration)) ) # subset to relevant part of df  
  dates <- seq(info$start,(info$start+info$duration-1),1)  # dates with missing values
  left_gst <- mean(x[x$date<info$start,'gst'],na.rm=T)  # first interpolation value
  right_gst <- mean(x[x$date>=(info$start+info$duration),'gst'],na.rm=T)  # last interpolation value
  
  if(is.finite(left_gst) & is.finite(right_gst) & abs(left_gst-right_gst)>0 & info$duration>1) {
    gst <- round( seq( left_gst, right_gst, (right_gst-left_gst)/(info$duration+1) ) ,2)[1:info$duration+1]  # LI
  } else {
    gst <- round( mean(x$gst,na.rm=T) ,2)  # handle special situations 
  }
  error <- round( sd(x$gst,na.rm=T) ,2)  # uncertainty calculation based on the mean standard variation
  
  return(data.frame(date=dates,sim=gst,error=error,normdist=1))
  rm(x,info,left_gst,right_gst,left_sd,right_gst,dates)
}

# ---- QUANTILE MAPPING ----

fill_QM <- function(df=df, x=x, info=info, win=15, mindur=30, minobs=3, minreg=5){
    # parameters
    # win <- 15  # used for filtering and smoothing - uneven for computational reasons
    # mindur <- 30  # minimal duration of QM-fitting window
    # minobs <- 3  # minimal factor by which the common observations must surpass the gap duration
    # minreg <- 5  # minimal number of regressors for the pre-selection to continue the procedure

  gap_dates <- seq((info$start),(info$start+info$duration-1),1)
  
  if(info$duration<mindur) {  # extend very short gaps to reach at least 30 days
    extension <- round(ceiling(mindur-info$duration)/2)
  } else {
    extension <- 0
  }
  
  ydays <- yday(seq((info$start-extension),(info$start+info$duration-1+extension),1))
  x2 <- x %>% filter(date %in% gap_dates) %>% dplyr:::select(date,gst)  # target-gap only
  x <- x %>% filter(yday(date) %in% ydays & is.finite(gst)) %>% dplyr:::select(date,gst)  # target with all finite data
  
  # search best regressor logger: pre-selection
  otherlocs <- df %>% filter(date %in% gap_dates & loc!=g$loc & is.finite(gst)) %>% 
    group_by(loc) %>% summarise( n=n() ) %>% filter(n>=g$duration) %>% dplyr:::select(loc)  # covers the entire gap and ...
  reg <- df %>% filter(loc %in% otherlocs$loc & is.finite(gst) & yday(date) %in% ydays) %>% 
    dplyr:::select(loc,date,gst) 
  otherlocs <- as.data.frame( left_join(x[,c('date','gst')],reg,by="date") ) %>% 
    filter(is.finite(gst.x) & is.finite(gst.y)) %>% group_by(loc) %>% summarize( 
      n=n() ) %>% filter(n>=minobs*length(ydays) ) %>% dplyr:::select(loc)  # at least 3*ydays of common data
  
  if(length(unique(otherlocs$loc)) >= minreg) {  # only proceed if at least 5 potential regressors are avilable
    o <- reg %>% filter(loc %in% otherlocs$loc) %>% left_join(.,x,by="date") %>% 
      filter(is.finite(gst.x) & is.finite(gst.y)) %>%
      group_by(otherloc=loc) %>% summarize( 
        sd_finite = ifelse(sd(gst.y,na.rm=T)!=0 & sd(gst.x,na.rm=T!=0),1,0),
        pearson = ifelse(sd_finite==1, cor(x=gst.y, y=gst.x, use="complete.obs", method="pearson"), NA), 
        spearman = ifelse(sd_finite==1, cor(x=gst.y, y=gst.x, use="complete.obs", method="spearman"), NA),
        SEMD = SEM(gst.y-gst.x),
        n = sum(gst.x*gst.y*0+1,na.rm=T) ) 
    o <- o %>% filter(n>=median(o$n)) 
    o <- o %>% filter(pearson>=quantile(o$pearson,na.rm=T,probs=0.95) | spearman>=quantile(o$spearman,na.rm=T,probs=0.95) |
                        SEMD<=quantile(o$SEMD,na.rm=T,probs=0.05)) %>% dplyr:::select(otherloc) 
    # select regressor by minimizing QM residuals
    for(oloc in unique(o$otherloc)) {
      r <- reg %>% filter(loc==oloc) %>% dplyr:::select(date,gst)  # selected regressor vector
      # merge <- na.omit( merge(r, x, by.x="date", by.y="date") ) 
      merge <- na.omit( inner_join(r, x, by="date" ) ) 
      if(g$duration>30 & g$duration<330 & length(unique(year(merge$date)))>5) {  # reduce to complete years
        merge2 <- merge %>% group_by(year=year(date)) %>% summarize(n=n()) 
        merge2 <- merge2 %>% filter(n>=median(merge2$n,na.rm=T))
        merge <- merge %>% filter(year(date) %in% merge2$year)
      }
      # fit quantiles with error handling
      qm.fit <- try( fitQmapQUANT(obs=merge$gst.y, mod=merge$gst.x, qstep=0.01, wet.day=F) , silent=T) # standard res.
      # apply bias-correction for each quantile
      fit <- try( data.frame(date=r$date, sim=round( doQmapQUANT( r$gst, qm.fit, type ="linear") ,2)), silent=T)
      if(class(qm.fit)=="try-error" | class(fit)=="try-error"){  # to avoid unsolvable situations
        qm.fit <- fitQmapQUANT(obs = merge$gst.y + rnorm(n=length(merge$date), mean=0,sd=0.00001), 
                               mod = merge$gst.x + rnorm(n=length(merge$date), mean=0,sd=0.00001), qstep=0.01, wet.day=F)
        fit <- data.frame(date=r$date, sim=round( doQmapQUANT( r$gst, qm.fit, type ="linear") ,2))
      }    
      if(exists('stat') && class(stat)=="data.frame"){  
        stat <- rbind(stat, data.frame(oloc, merge(fit, x, by.x="date", by.y="date") %>% 
                            summarize( SEMD = SEM(abs(sim-gst)), nullsum=sum(sim*gst,na.rm=T), n=n() ) ) )
      } else {
        stat <- data.frame(oloc, merge(fit, x, by.x="date", by.y="date") %>% 
                summarize( SEMD = SEM(abs(sim-gst)), nullsum=sum(sim*gst,na.rm=T), n=n() ) )
      }
      stat <- stat %>% filter(nullsum>0)  # contains performance stats for each potential regressor
    }
    o <- stat %>% filter(SEMD==min(stat$SEMD,na.rm=T)) %>% head(1) %>% dplyr:::select(oloc)  # final selection of best regressor
    if(!is.na(o$oloc)) {
      r <- reg %>% filter(loc==o$oloc) %>% dplyr:::select(date,gst)  # selected regressor vector
      merge <- na.omit( inner_join(r, x, by="date" ) )   # ...merged with target vector
      
      # fit quantiles with error handling
      qm.fit <- try( fitQmapQUANT(obs=merge$gst.y, mod=merge$gst.x, qstep=0.01, wet.day=F) , silent=T) # standard res.
      
      if(class(qm.fit)=="try-error")  {
        qm.fit <- try( fitQmapQUANT(obs=merge$gst.y, mod=merge$gst.x, qstep=0.005, wet.day=F) , silent=T) # finer res.
      }
      
      if(class(qm.fit)=="try-error")  {
        qm.fit <- try( fitQmapQUANT(obs=merge$gst.y, mod=merge$gst.x, qstep=0.05, wet.day=F) , silent=T) # coarser res.
      }
      # apply bias-correction for each quantile
      fit <- try( data.frame(date=r$date, sim=round( doQmapQUANT( r$gst, qm.fit, type ="linear") ,2)) , silent=T)
      if(class(qm.fit)=="try-error" | class(fit)=="try-error") {  # to avoid unsolvable solutions
        qm.fit <- fitQmapQUANT(obs = merge$gst.y + rnorm(n=length(merge$date), mean=0,sd=0.00001), 
                               mod = merge$gst.x + rnorm(n=length(merge$date), mean=0,sd=0.00001), qstep=0.01, wet.day=F)
        fit <- data.frame(date=r$date, sim=round( doQmapQUANT( r$gst, qm.fit, type ="linear") ,2))    
      }
      # estimate uncertainty (error=estimated bias, normdist 1:=potentially normally distributed)
      if(info$duration>mindur) {
        error <- try ( na.omit( inner_join(fit, x, by="date") ) %>% 
          filter( is.finite(sim) & is.finite(gst) & yday(date) %in% yday(gap_dates) ) %>% 
          group_by( yday=yday(date)) %>% summarize(
            error = round( sqrt( (1/(n()-1))*sum(abs(sim-gst)^2,na.rm=T)) ,3),
            normdist = ifelse( median(gst,na.rm=T)>=1 & median(sim,na.rm=T) >=1 , 1, 0)  ) %>% ungroup() %>% 
          dplyr:::select(yday,error,normdist) %>% arrange(yday) %>% mutate(
            error = ifelse(is.finite(error),error,median(error,na.rm=T)),
            normdist = ifelse(is.finite(normdist), normdist, 0),
            normdist = round(runmean(x=normdist,k=win,alg="C",endrule="mean",align="center"),0) ) ,silent=T)
      }
      if(info$duration<=mindur) {
        error <- try( na.omit( inner_join(fit, x, by="date") ) %>% arrange(as.numeric(date)) %>% mutate(
          rbias = runmean(x=sim-gst,k=win,alg="C",endrule="mean",align="center"),
          rbias_abs = runmean(x=abs(sim-gst),k=win,alg="C",endrule="mean",align="center"),
          # check if is this day of the year might be snow covered or if the melt-out is less than 14 days back...
          snow = ifelse(runsd(gst,k=win,endrule="sd",align="center")<0.25 | runsd(sim,k=win,endrule="sd",align="center")<0.25 | 
                          runsd(lag(gst,14),k=win,endrule="sd",align="center")<=0.1 | runsd(lag(sim,14),k=win,endrule="sd",align="center")<=0.1, 1, 0),
          yday = yday(date) ) %>% 
            filter( is.finite(sim) & is.finite(gst) & is.finite(rbias) & yday %in% yday(gap_dates) ) %>% 
            group_by( yday) %>% summarize( 
              error = round( sqrt( (1/(n()-1))*sum(rbias_abs^2,na.rm=T)) ,3),
              n = n(),
              normdist = ifelse(n<3 | median(snow,na.rm=T)>=0.5, 0, 1) ) %>% ungroup() %>% 
            dplyr:::select(yday,error,normdist) %>% arrange(yday) %>% mutate(
              error = ifelse(is.finite(error),error,median(error,na.rm=T)),
              normdist = ifelse(is.finite(normdist), normdist, 0),
              normdist = round(runmean(x=normdist,k=win,alg="C",endrule="mean",align="center"),0) ) ,silent=T)
        } 
      if(last(error$yday)==365) {  # add yday 366 (=365) just in case...
        error <- rbind(error, data.frame(yday=366,error=last(error$error),normdist=last(error$normdist)))
      }
      # extract GST estimate for gap
      fit <- fit[fit$date %in% gap_dates,] %>% mutate( yday=yday(date) )
      fit <- na.omit( inner_join(fit, error, by="yday") ) %>% 
        dplyr:::select(date, sim, error, normdist) %>% arrange(as.numeric(date))
      return( data.frame(fit, otherloc=o$oloc) )
    } else {
      return( data.frame(date=gap_dates, sim=NA, error=NA, normdist=NA, otherloc=NA) )
    }
  } else {  # not enough potential regressors to proceed
    return( data.frame(date=gap_dates, sim=NA, error=NA, normdist=NA, otherloc=NA) )
  }
  rm(x,info,x2,gap_dates,ydays,otherlocs,reg,o,r,merge,qm.fit,fit,error,stat,win)
}
