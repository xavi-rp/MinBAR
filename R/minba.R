###########################################################################
########                                                       ############
########             Minimum Background Area for SDMs          ############
########                          MinBA                        ############
########                                                       ############
###########################################################################

#' Downloading species ocurrence data (presences) from Bioatles
#'
#' minba aims at (1) defining what is the minimum or optimal background extent 
#' necessary to fit good partial SDMs and/or (2) determining if the background 
#' area used to fit a partial SDM is reliable enough to extract ecologically 
#' relevant conclusions from it.
#'
#' @author Xavier Rotllan-Puig
#' @description It aims at (1) defining what is the minimum or optimal background extent necessary to fit good partial SDMs and/or (2) determining if the background area used to fit a partial SDM is reliable enough to extract ecologically relevant conclusions from it.
#' @param occ Data set with presences (occurrences). A csv file with 3 columns: long, lat and species name
#' @param climdat An .RData file containing a RasterStack with the independent variables
#' @param num_bands Number of bandwidths
#' @param n_times Number of models
#' @param BI_part 
#' @param BI_tot 
#' @param SD_BI_part 
#' @param SD_BI_tot 
#' @return \code{selfinfo_mod_}, \code{info_mod_} and \code{info_mod_means_} (all followed by the name of the species). The first two tables are merely informative about how the modelling process has been developed and the results of each model. Whereas \code{info_mod_means_} shows the means of the n models run for each bandwidth
#' @name minba()
#'
# Created on: Summer 2018
#
# Created by: Xavier Rotllan-Puig (xavi.rotllan.puig@gmail.com)
# 
# 

minba <- function(occ = NULL, climdat = NULL, 
                  num_bands = 10, n_times = 3,
                  BI_part = NULL,BI_tot = NULL,
                  SD_BI_part = NULL, SD_BI_tot = NULL
                  ){
  #### Settings ####
  #setwd("/Users/xavi/Google Drive/MinBAR/data")
  wd <- getwd()
  dir2save <- paste0(wd, "/minba_", format(Sys.Date(), format="%Y%m%d"))
  if(!file.exists(dir2save)) dir.create(dir2save)
  
  #### Retrieving Presence Records ####
  #occ <- "sp_records_bioatles.csv"
  presences <- read.csv(occ, header = TRUE)
  head(presences)
  presences$sp2 <- tolower(paste(substr(presences$species, 1, 3), substr(sub(".* ", "", presences$species), 1, 3), sep = "_"))
  colnames(presences)[1:2] <- c("lat", "lon") 
  presences <- presences[, c(2,1,3,4)]

  #### Climatic Data ####
  #climdat <- "wc05.RData"
  load(paste0(wd, "/wc0.5/", climdat), verbose = TRUE)
  bioclim
  
  
  #Downloading data from worldclim.org
  if(tolower(climdat) != "no"){
    if (data_rep == "gbif"){
      bioclim <- getData('worldclim', var='bio', res = resol, path = paste0(wd))
      save(bioclim, file = "wc5/wc5.RData")
    }
    if (data_rep == "bioatles"){
      bioclim_16 <- getData('worldclim', var='bio', res = resol, lon=3, lat=39,
                            path = paste0(wd))  # importing tile 16
      bioclim <- bioclim_16
      save(bioclim, file = "wc0.5/wc05.RData")
      rm(bioclim_16) ; gc()
    }
  }
  if (data_rep == "gbif")  load(paste0(wd, "/wc5/wc5.RData"), verbose = FALSE)
  if (data_rep == "bioatles")  
  
  
  #Making a mask (if necessary)
  mskng <- "no"
  if(tolower(mskng) != "no"){
    msk <- bioclim$bio12
    msk <- reclassify(msk, c(0, msk@data@max, 1)) # reclassify to 1
  }
  
  #### Modelling per each species ####
  specs <- unique(presences$sp2)
  #specs <- unique(presences$sp2)[12]
  
  best2_bnd_2exp <- as.data.frame(matrix(ncol = 0, nrow = 0)) # a table to export rankings of best and 2nd best bandwidth
  
  for(sps in specs){
    pres <- presences[presences$sp2 %in% sps, ] # selecting for species
    specs_long <- as.character(unique(pres$species)) # complete name of the species
    pres <- pres[, c(1,2,4)]
    coordinates(pres) <- c("lon", "lat")  # setting spatial coordinates
    pres@proj4string <- CRS("+init=EPSG:4326")
    #plot(pres)
    
    #### Calculating the centre of the population, its most distant point and "bands" ####
    #pop_cent <- c(mean(presences$x, na.rm =TRUE), mean(presences$y, na.rm =TRUE))
    #pop_cent <- as.data.frame(matrix(pop_cent, ncol = 2))
    #names(pop_cent) <- c("x", "y")
    
    geocntr <- as.data.frame(geomean(pres))  #mean location for spherical (longitude/latitude) coordinates that deals with the angularity
    #geocntr1 <- geocntr
    #coordinates(geocntr1) <- c("x", "y") 
    #plot(geocntr1, col = "green", add = TRUE)
    
    pres$dist2centr <- distGeo(pres, geocntr) #in meters
    pres$dist2centr <- pres$dist2centr/1000   #in km
    furthest <- max(pres$dist2centr) 
    
    # by % of presences equally distributed
    bndwidth <- as.vector(quantile(pres$dist2centr, probs = seq(0, 1, 1/num_bands), names = TRUE))
    bndwidth <- c(bndwidth[2:(length(bndwidth)-1)], furthest) # Not defined by distance, but by % of presences equally distributed
    # This is particularly useful for very discontinuous distributions (e.g. introduced or invasive species),
    # while not affecting more aggregated populations
    # by equally distant bandwidths
    #bndwidth <- as.vector(seq(0, furthest, furthest/num_bands))[-1]
    
    #### Croping variables to pres extent  + 5%
    ext1 <- pres@bbox
    incr1 <- apply(ext1, 1, function(x) x[2] - x[1]) * 0.05
    ext1[1, 1] <- pres@bbox[1, 1] - incr1[1]
    ext1[1, 2] <- pres@bbox[1, 2] + incr1[1]
    ext1[2, 1] <- pres@bbox[2, 1] - incr1[2]
    ext1[2, 2] <- pres@bbox[2, 2] + incr1[2]
    varbles1 <- stack(crop(bioclim, ext1))
    
    # number of background points (see Guevara et al, 2017)
    num_bckgr1 <- (varbles1@ncols * varbles1@nrows) * 50/100 
    
    
    #### Making models for each bandwidth ####
    #table with info to be exported
    dt2exp <- as.data.frame(matrix(ncol = 12, nrow = 0))
    selfinfo2exp <- as.data.frame(matrix(ncol = 9, nrow = 0))
    dt2exp_mean <- as.data.frame(matrix(ncol = 7, nrow = 0))
    
    for (bdw in 1:length(bndwidth)) { # for each bandwidth
      x <- 1
      # set of presences for modeling within the bandwidth
      pres4model <- pres[pres$dist2centr <= bndwidth[bdw], ] 
      
      # croping variables to pres4model extent  + 5%  <-- to fit the model
      ext <- pres4model@bbox
      incr <- apply(ext, 1, function(x) x[2] - x[1]) * 0.05
      ext[1, 1] <- pres4model@bbox[1, 1] - incr[1]
      ext[1, 2] <- pres4model@bbox[1, 2] + incr[1]
      ext[2, 1] <- pres4model@bbox[2, 1] - incr[2]
      ext[2, 2] <- pres4model@bbox[2, 2] + incr[2]
      varbles <- stack(crop(bioclim, ext))
      #dev.off()
      #plot(varbles$bio2_16)
      #plot(pres, add = TRUE)
      
      # number of background points (see Guevara et al, 2017)
      num_bckgr <- (varbles@ncols * varbles@nrows) * 50/100 
      #if(num_bckgr<100) {
      #  pres4model1 <- sample(1:nrow(pres4model), nrow(pres4model)*0.1) 
      #  pres4model <- pres4model[pres4model1,]
      #}
      
      
      # sampling presences for calibrating and testing (70-30%) within the bandwidth
      folds <- sample(1:nrow(pres4model), nrow(pres4model)*0.7)  
      samp <- as.numeric(unlist(folds))
      pres4cali <- pres4model[samp, 1]
      pres4test <- pres4model[-samp, 1]
      rm(pres4model); gc()
      
      # sampling presences for testing on the whole extent (30% of total presences except those for calibrating)
      pres1 <- pres[-samp, 1]
      folds1 <- sample(1:nrow(pres1), nrow(pres1)*0.3)  
      samp1 <- as.numeric(unlist(folds1))
      pres4test_tot <- pres1[-samp1, 1]
      
      
      repeat{   # maybe it can be done directly with maxent; if so, we would also have the "average-model"
        t1 <- Sys.time()
        cat("\r","modelling for",specs_long,"- bandwidth #",bdw,"_",x)
        
        # Running maxent from dismo 
        if(!file.exists(paste0(dir2save,"/results_", sps))) dir.create(paste0(dir2save,"/results_", sps))
        if(!file.exists(paste0(dir2save,"/results_", sps, "/model_", sps, "_", bdw, "_", x))) dir.create(paste0(dir2save,"/results_", sps, "/model_", sps, "_", bdw, "_", x))
        path <- paste0(dir2save,"/results_", sps,"/model_", sps, "_", bdw, "_", x)
        
        dir_func <- function(varbles, pres4cali, num_bckgr, path){ # to avoid stop modelling if low number of background points or other errors
          res <- tryCatch(
            {
              modl <- maxent(varbles, pres4cali, removeDuplicates = TRUE,
                             nbg=num_bckgr)
              if(exists("modl")) save(modl, file = paste0(path, "/model.RData"))
            },
            error = function(con){
              message(con)
              return(NULL)
            }
          )
          if(exists("modl")){ return(modl) }else{ return(NULL) } 
        } #end of dir_func
        
        modl <- dir_func(varbles, pres4cali, num_bckgr, path)
        if(is.null(modl)){ break }
        
        #making predictions on the same extent
        preds <- predict(modl, varbles, filename = paste0(path, "/predictions"), progress = '', overwrite = TRUE)
        
        #make evaluations (on the same extent with 30% to test)
        bg <- randomPoints(varbles, num_bckgr) # background points
        evs <- evaluate(modl, p=pres4test, a=bg, x=varbles)
        save(evs, file = paste0(path, "/evaluations.RData"))
        #load(paste0(path, "/evaluations.RData"), verbose = TRUE)
        
        #Computing Boyce Index (on the same extent with 30% to test)
        byce <- ecospat.boyce(fit = preds, obs = pres4test@coords, nclass=0, window.w="default", res=100, PEplot = TRUE)
        byce$Spearman.cor
        save(byce, file = paste0(path, "/boyce.RData"))
        #load(paste0(path, "/boyce.RData"), verbose = TRUE)
        
        # In the last bandwidth it would make no sense repeating predictions/evaluations on the same extent
        # to assess for transferability. However, for the sake of consistency in the execution time, they are 
        # calculated.
        # To return to the version where they are not calculated, check commit e6e0040 of 25/08/2018
        
        ## making predictions on the whole species extent
        preds1 <- predict(modl, varbles1, filename=paste0(path, "/predictions_tot"), progress = '', overwrite = TRUE)
        
        #make evaluations
        bg1 <- randomPoints(varbles1, num_bckgr1) # background points
        evs1 <- evaluate(modl, p=pres4test_tot, a=bg1, x=varbles1)
        save(evs1, file = paste0(path, "/evaluations_tot.RData"))
        #load(paste0(path, "/evaluations_tot.RData"), verbose = T)
        
        #Computing Boyce Index
        byce1 <- ecospat.boyce(fit = preds1, obs = pres4test_tot@coords, nclass=0, window.w="default", res=100, PEplot = TRUE)
        byce1$Spearman.cor
        save(byce1, file = paste0(path, "/boyce_tot.RData"))
        #load(paste0(path, "/boyce_tot.RData"), verbose = TRUE)
        
        # gathering info to be exported
        t2 <- Sys.time() - t1
        if(attr(t2, "units") == "hours") {t2 <- t2*60; attr(t2, "units") <- "mins"}
        if(attr(t2, "units") == "secs") {t2 <- t2/60; attr(t2, "units") <- "mins"}
        dt2exp_2 <- as.data.frame(matrix(c(specs_long, paste(bdw, x, sep="_"), bndwidth[bdw], nrow(modl@presence), evs@np, num_bckgr, evs@auc, byce$Spearman.cor, evs1@np, num_bckgr1, evs1@auc, byce1$Spearman.cor), 1, 12, byrow = TRUE))
        dt2exp <- rbind(dt2exp, dt2exp_2)  
        selfinfo2exp_2 <- as.data.frame(matrix(c(specs_long, paste(bdw, x, sep="_"), nrow(pres4cali), nrow(modl@presence), nrow(pres4test), evs@np, num_bckgr, evs@na, t2), 1, 9, byrow = TRUE))
        selfinfo2exp <- rbind(selfinfo2exp, selfinfo2exp_2)
        
        if (x == n_times){ break }else{ x <- x +1 }
      } #end of repeat n times
      
      if(is.null(modl)){ print("jumping to next bandwidth"); next }
      
      cat("\n","computing average for",specs_long,"- bandwidth #",bdw,"\n")
      dt2exp[,-c(1:6)] <- data.frame(lapply(dt2exp[-c(1:6)], function(x) as.numeric(as.character(x))))
      dt2exp_m <- mean(dt2exp[(nrow(dt2exp)-n_times+1):nrow(dt2exp), (ncol(dt2exp)-4)], na.rm = TRUE) #mean Boyce partial area
      dt2exp_m2 <- mean(dt2exp[(nrow(dt2exp)-n_times+1):nrow(dt2exp), ncol(dt2exp)], na.rm = TRUE) #mean Boyce whole area
      if (bdw > 3){
        dt2exp_sd1 <- sd(c(dt2exp_mean$V3[(bdw-3):(bdw-1)], dt2exp_m))
        dt2exp_sd2 <- sd(c(dt2exp_mean$V4[(bdw-3):(bdw-1)], dt2exp_m2))
      }else{
        dt2exp_sd1 <- NA
        dt2exp_sd2 <- NA
      }
      selfinfo2exp[,-c(1:8)] <- data.frame(lapply(selfinfo2exp[-c(1:8)], function(x) as.numeric(as.character(x))))
      dt2exp_m1 <- mean(selfinfo2exp[(nrow(selfinfo2exp)-n_times+1):nrow(selfinfo2exp), ncol(selfinfo2exp)], na.rm = TRUE)
      dt2exp_mean_2 <- as.data.frame(matrix(c(specs_long, bndwidth[bdw], dt2exp_m, dt2exp_m2, dt2exp_sd1, dt2exp_sd2, dt2exp_m1), 1, 7, byrow = TRUE))
      dt2exp_mean <- rbind(dt2exp_mean, dt2exp_mean_2)
      dt2exp_mean[, c(2:7)] <- data.frame(lapply(dt2exp_mean[c(2:7)], function(x) as.numeric(as.character(x))))
      
      rm(modl, preds, preds1, bg, evs, byce); gc()
      
      #Conditions to stop the process
      if(!is.null(BI_part) | !is.null(BI_tot) | !is.null(SD_BI_part) | !is.null(SD_BI_tot)){
        brk <- 0
        if (!is.null(BI_part) && BI_part <= dt2exp_m){ print("minimum BI_part has been reached"); brk <- 1}
        if (!is.null(BI_tot) && BI_tot <= dt2exp_m2){ print("minimum BI_tot has been reached"); brk <- 1}
        if (!is.null(SD_BI_part) && !is.na(dt2exp_sd1) && SD_BI_part >= dt2exp_sd1){ print("minimum SD_BI_part has been reached"); brk <- 1}
        if (!is.null(SD_BI_tot) && !is.na(dt2exp_sd2) && SD_BI_tot >= dt2exp_sd2){ print("minimum SD_BI_tot has been reached"); brk <- 1}
        if (brk == 1)  break
      }
      
    } # end of for each bandwidth
    
    names(dt2exp_mean) <- c("Species", "Bandwidth", "BoyceIndex_part", "BoyceIndex_tot", "SD_part", "SD_tot", "ExecutionTime")
    names(dt2exp) <- c("Species", "ModelNum", "Bandwidth", "numPresencesCalib", "numPresencesTest", "numBackground", "AUC_part", "BoyceIndex", "numPresencesTest_tot", "numBackground_tot", "AUC_tot", "BoyceIndex_tot")
    names(selfinfo2exp) <- c("Species", "ModelNum", "num_pres_calib", "num_pres_calib_used", "num_pres_test", "num_pres_test_used", "num_background", "num_bckgrnd_used", "exec_time")
    
    computing_ranks <- 1
    if (computing_ranks == 1){
      dt2exp_mean$rankBI_part <- rank(-dt2exp_mean$BoyceIndex_part, ties.method = "first")
      dt2exp_mean$rankBI_tot <- rank(-dt2exp_mean$BoyceIndex_tot, ties.method = "first")
      dt2exp_mean$rankTime <- rank(dt2exp_mean$ExecutionTime, ties.method = "first")
      dt2exp_mean$rankFinalNoTime <- rank((dt2exp_mean$rankBI_part + dt2exp_mean$rankBI_tot), ties.method = "first")
      dt2exp_mean$rankFinalWithTime <- rank((dt2exp_mean$rankBI_part + dt2exp_mean$rankBI_tot + dt2exp_mean$rankTime), ties.method = "first")
      
      best2_bnd <- c(sps, 
                     row.names(dt2exp_mean[dt2exp_mean$rankFinalNoTime == 1, ]), 
                     row.names(dt2exp_mean[dt2exp_mean$rankFinalNoTime == 2, ]),
                     row.names(dt2exp_mean[dt2exp_mean$rankFinalWithTime == 1, ]), 
                     row.names(dt2exp_mean[dt2exp_mean$rankFinalWithTime == 2, ]))
      best2_bnd <- as.data.frame(t(best2_bnd))
      names(best2_bnd) <- c("Species", "Best_Bandwidth_NoTime", "SecondBest_bandwidth_NoTime", "Best_Bandwidth_WithTime", "SecondBest_bandwidth_WithTime")
      
      best2_bnd_2exp <- rbind(best2_bnd_2exp, best2_bnd)
      write.csv(best2_bnd_2exp, paste0(dir2save, "/rankingBestBandwidth.csv"), row.names = FALSE)
    }
    
    write.csv(dt2exp_mean, paste0(dir2save, "/results_", sps, "/info_mod_means_", sps, ".csv"), row.names = FALSE)
    write.csv(dt2exp, paste0(dir2save, "/results_", sps, "/info_mod_", sps, ".csv"), row.names = FALSE)
    write.csv(selfinfo2exp, paste0(dir2save, "/results_", sps, "/selfinfo_mod_", sps, ".csv"), row.names = FALSE)
    
    #### Making a plot ####
    graphics.off()
    #dt2exp_mean[,-1] <- data.frame(lapply(dt2exp_mean[-1], function(x) as.numeric(as.character(x))))
    dt2exp_mean[,names(dt2exp_mean) %in% c("BoyceIndex_part", "BoyceIndex_tot")] <- round(dt2exp_mean[,names(dt2exp_mean) %in% c("BoyceIndex_part", "BoyceIndex_tot")], 3)
    pdf(paste0(dir2save, "/results_", sps, "/boyce_bandwidth_", sps, "_part_tot.pdf"))
    if(nrow(dt2exp_mean) < 5){ tp <- c("p") }else{ tp <- c("p", "smooth") }
    plt <- xyplot(BoyceIndex_part ~ Bandwidth, dt2exp_mean,
                  #scales = list(y = list(log = 10)),
                  #xavi180825: type = c("p", "smooth"),
                  type = tp,
                  span = 0.8,  
                  #type = c("p", "l"), 
                  #xavi180825: ylim = c(0.8, 1.05),
                  ylim = c(0.45, 1.05),
                  col = "blue",
                  main = bquote(Boyce~Index~(mean~of~.(n_times)~models)~-~italic(.(specs_long))),
                  ylab = "Boyce Index", xlab = "Bandwidth (km)",
                  #par.settings = list(par.ylab.text = list(col = "black")),
                  #par.settings = simpleTheme(col = 1),
                  key=list(#space = "right",
                    x=0.5,y=0.2,
                    lines = list(col=c("blue", "green", "magenta")),
                    text = list(c("Boyce Index Partial","Boyce Index Total", "Execution Time"))
                  ))
    plt1 <- xyplot(ExecutionTime ~ Bandwidth, dt2exp_mean,
                   #type = c("p", "smooth"), 
                   type = c("p", "r"), 
                   #ylim = c(0.8, 1.05),
                   ylab = "Execution Time (min)",
                   col = "magenta"
                   #,par.settings = simpleTheme(col = "magenta")
    )
    dbl_plt <- doubleYScale(plt, plt1, add.ylab2 = TRUE)
    #plot(dbl_plt)
    plt2 <- xyplot(BoyceIndex_tot ~ Bandwidth, dt2exp_mean,
                   #xavi081825:  type = c("p", "smooth"),
                   type = tp,
                   span = 0.8,
                   #type = c("p", "l"), 
                   #ylim = c(0.8, 1.05),
                   col = "green")
    
    plot(dbl_plt + as.layer(plt2))
    dev.off()
    
  } # end of loop for sps
  
  #
  
  #### Frequences of Best bandwidth #### 
  
  best2_bnd_2exp <- read.csv(paste0(dir2save, "/rankingBestBandwidth.csv"), header = TRUE)
  frec_best_NoTime <- as.data.frame(table(best2_bnd_2exp$Best_Bandwidth_NoTime))
  frec_best_WithTime <- as.data.frame(table(best2_bnd_2exp$Best_Bandwidth_WithTime))
  
  frec_best <- as.data.frame(matrix(seq(1:num_bands), nrow = num_bands, ncol = 1))
  frec_best <- merge(frec_best, frec_best_NoTime, by.x = "V1", by.y = "Var1", all = TRUE)
  frec_best <- merge(frec_best, frec_best_WithTime, by.x = "V1", by.y = "Var1", all = TRUE)
  frec_best[is.na(frec_best)] <- 0 
  names(frec_best) <- c("BandwidthNum", "Frec_Best_NoTime", "Frec_Best_WithTime")
  
  perc1 <- round((prop.table(frec_best$Frec_Best_NoTime)*100), 0)
  perc2 <- round((prop.table(frec_best$Frec_Best_WithTime)*100), 0)
  perc <- paste0(c(perc1, perc2), "%")
  maxY <- max(frec_best_NoTime$Freq, frec_best_WithTime$Freq) + 1
  posX <- c(1:num_bands, (num_bands + 2):((2 * (num_bands)) + 1)) - 0.3
  posY <- (maxY / 10 * 0.2) + c(frec_best$Frec_Best_NoTime, frec_best$Frec_Best_WithTime)
  if (maxY <= 5) legY <- -0.4 else legY <- -1.5
  
  palte <- colorRampPalette(colors = c("darkgreen", "green", "yellow", "orange", "red", "darkred"))(num_bands)
  
  pdf(paste0(dir2save, "/BestBandwidths.pdf"))
  par(xpd = TRUE, mar = par()$mar + c(3,1,0,0))
  bpl <- barplot(as.matrix(frec_best[,c(2:3)]),
                 main = "Best Bandwidth With and Without Execution Time", ylab = paste0("Frequencies (n = ", nrow(best2_bnd_2exp), ")"), 
                 ylim = c(0, maxY),
                 beside = TRUE, space = c(0, 1),
                 col = palte,
                 #names.arg = c(1:10, 1:10),
                 names.arg = c("Without Execution Time", "With Execution Time"),
                 cex.names = 1, las = 1)
  lg <- legend((num_bands/2), legY,
               legend = frec_best$BandwidthNum, 
               fill = palte,
               title = "Bandwidth Number", cex = 1, ncol = (num_bands/2))
  txt <- text(x = posX, y = posY, perc, cex = 0.8, pos = 4, srt = 45)
  
  dev.off()
  # 
  
  
  #### Plot with all bandwidth plots ####
  
  
  
  #### BI by groups ####
  
  dt2exp_mean_all <- as.data.frame(matrix(ncol = 0, nrow = 0))
  
  for (sps in specs){
    dt2exp_mean <- read.csv(paste0(dir2save, "/results_", sps, "/info_mod_means_", sps, ".csv"), header = TRUE)[, 1:4]
    dt2exp_mean_all <- rbind(dt2exp_mean_all, dt2exp_mean)
  }
  
  BI_SD_bySpecies <- as.data.frame(dt2exp_mean_all %>% group_by(Species) %>% summarise(Mean_BI_part = mean(BoyceIndex_part),
                                                                                       SD_BI_part = sd(BoyceIndex_part),
                                                                                       Mean_BI_tot = mean(BoyceIndex_tot),
                                                                                       SD_BI_tot = sd(BoyceIndex_tot)))
  BI_SD_bySpecies
  range(BI_SD_bySpecies$Mean_part)
  range(BI_SD_bySpecies$SD_part)
  range(BI_SD_bySpecies$Mean_tot)
  range(BI_SD_bySpecies$SD_tot)
  
  
  
  
  
}
  
