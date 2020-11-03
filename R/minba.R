###########################################################################
########                                                       ############
########             Minimum Background Area for SDMs          ############
########                          MinBA                        ############
########                                                       ############
###########################################################################
#'
#' @author Xavier Rotllan-Puig & Anna Traveset
#' @title Determining the Minimal Background Area for Species Distribution Models
#' @description A versatile tool that aims at aims at (1) defining the minimum background extent necessary to fit SDMs reliable enough to extract ecologically relevant conclusions from them and (2) optimizing the modelling process in terms of computation demands. See Rotllan-Puig, X. & Traveset, A. (2021)
#' @details Please check the article 'Determining the Minimal Background Area for Species Distribution Models: MinBAR Package' for further details on how to use this package, examples, etc.
#' @import "dismo" "maxnet"
#' @importFrom grDevices dev.off graphics.off pdf
#' @importFrom graphics plot
#' @importFrom stats quantile sd
#' @importFrom utils read.csv write.csv
#' @param occ Data frame or character. Data set with presences (occurrences). A data frame with 3 columns: long, lat and species name (in this order)
#' @param varbles Raster* object. A raster brick of the independent variables, or a directory where the rasters are. It will use all the rasters in the folder. Supported: .tif and .bil
#' @param wd Character. A directory to save the results
#' @param prj Numeric. Coordinates system (e.g. "4326" is WGS84; check \url{https://spatialreference.org/} )
#' @param num_bands Numeric. Number of buffers (default is 10)
#' @param n_rep Numeric. Number of replicates (default is 15)
#' @param occ_prop_test Numeric. Proportion of presences (occurrences) set aside for testing (default is 0.3)
#' @param maxent_tool Character. Either "dismo" or (default) "maxnet"
#' @param BI_part Numeric. Maximum Boyce Index Partial to stop the process if reached
#' @param BI_tot Numeric. Maximum Boyce Index Total to stop the process if reached
#' @param SD_BI_part Numeric. Minimum SD of the Boyce Index Partial to stop the process if reached (last 3 buffers)
#' @param SD_BI_tot Numeric. Minimum SD of the Boyce Index Total to stop the process if reached (last 3 buffers)
#' @return \code{selfinfo_mod_}, \code{info_mod_} and \code{info_mod_means_} (all followed by the name of the species). The first two tables are merely informative about how the modelling process has been developed and the results of each model. Whereas \code{info_mod_means_} shows the means of the n models run for each buffer
#' @name minba()
#' @references Rotllan-Puig, X. & Traveset, A. 2021. Determining the Minimal Background Area for Species Distribution Models: MinBAR Package. Ecological Modelling. 439:109353. https://doi.org/10.1016/j.ecolmodel.2020.109353
#' @export
#' @examples
#' \dontrun{
#' minba(occ = sprecords, varbles = bioscrop,
#'       wd = tempdir(), prj = 4326, num_bands = 3, n_rep = 3,
#'       maxent_tool = "maxnet")
#' }
#'
# Created on: Summer 2018 - Winter 2019 (Updated: Summer 2020)
#

minba <- function(occ = NULL, varbles = NULL,
                  wd = NULL,
                  prj = NULL,
                  num_bands = 10, n_rep = 15,
                  occ_prop_test = 0.3,
                  maxent_tool = "maxnet",
                  BI_part = NULL, BI_tot = NULL,
                  SD_BI_part = NULL, SD_BI_tot = NULL){
  #### Settings ####
  if(is.null(wd)) stop("Please, indicate a directory (wd) to save results")
  dir2save <- paste0(wd, "/minba_", format(Sys.Date(), format="%Y%m%d"))
  if(!file.exists(dir2save)) dir.create(dir2save)
  graphics.off()
  if(occ_prop_test <= 0 | occ_prop_test >= 1) stop("Please, specify a proportion of occurences to set aside for testing the models (default: 0.3)")


  #### Retrieving Presence Records ####
  if(is.vector(occ)){
    presences <- read.csv(occ, header = TRUE)
  }else if(!exists("occ") | ncol(occ) > 3){
    stop("Please provide a 3-columns data frame with the presences coordinates (long, lat and species name)")
  }else{
    presences <- occ
  }
  presences$sp2 <- tolower(paste(substr(presences$species, 1, 3), substr(sub(".* ", "", presences$species), 1, 3), sep = "_"))
  colnames(presences)[1:2] <- c("lon", "lat")
  presences <- presences[, c(2,1,3,4)]

  #### Climatic Data ####
  if(!exists("varbles") | is.vector(varbles)){
    rstrs <- list.files(varbles, pattern = c(".bil$"), full.names = T)
    rstrs <- c(rstrs, list.files(varbles, pattern = c(".tif$"), full.names = T))

    #vrbles <- stack()
    for(rst in 1:length(rstrs)){
      temp <- raster::raster(rstrs[rst])
      if(rst == 1){
        vrbles <- raster::stack(temp)
      }else{
        vrbles <- raster::stack(vrbles, temp)
      }
    }
  }else{
    vrbles <- varbles
  }
  if (!is.null(prj)) vrbles@crs <- sp::CRS(paste0("+init=EPSG:", prj))
  if (is.na(vrbles@crs)) stop("Please provide variables with a coordinates system or a
                              coordinate system reference through 'prj'")

  #### Modelling per each species ####
  specs <- unique(presences$sp2)

  best2_bnd_2exp <- as.data.frame(matrix(ncol = 0, nrow = 0)) # a table to export rankings of best and 2nd best buffer

  for(sps in specs){
    pres <- presences[presences$sp2 %in% sps, ] # selecting for species
    specs_long <- as.character(unique(pres$species)) # complete name of the species
    pres <- pres[, c(2, 1, 4)]
    sp::coordinates(pres) <- c("lon", "lat")  # setting spatial coordinates
    if (!is.null(prj)) pres@proj4string <- sp::CRS(paste0("+init=EPSG:", prj))
    if (is.na(pres@proj4string)) stop("Make sure the coordinates system of occurrences and variabes are the same")


    #### Calculating the centre of the population, its most distant point and "bands" ####
    if(!grepl("WGS84", vrbles@crs@projargs)){
      vrbles <- raster::projectRaster(vrbles, crs = sp::CRS(paste0("+init=EPSG:", 4326)))
      pres <- sp::spTransform(pres, CRSobj = sp::CRS(paste0("+init=EPSG:", 4326)))
    }

    geocntr <- as.data.frame(geosphere::geomean(pres))  #mean location for spherical (longitude/latitude) coordinates that deals with the angularity

    pres$dist2centr <- geosphere::distGeo(pres, geocntr) #in meters
    pres$dist2centr <- pres$dist2centr/1000   #in km
    furthest <- max(pres$dist2centr)

    # by % of presences equally distributed
    bndwidth <- as.vector(quantile(pres$dist2centr, probs = seq(0, 1, 1/num_bands), names = TRUE))
    bndwidth <- c(bndwidth[2:(length(bndwidth)-1)], furthest) # Not defined by distance, but by % of presences equally distributed
    # This is particularly useful for very discontinuous distributions (e.g. introduced or invasive species),
    # while not affecting more aggregated populations
    # by equally distant buffers
    #bndwidth <- as.vector(seq(0, furthest, furthest/num_bands))[-1]

    #### Croping variables to pres extent  + 5%
    ext1 <- pres@bbox
    incr1 <- apply(ext1, 1, function(x) x[2] - x[1]) * 0.05
    ext1[1, 1] <- pres@bbox[1, 1] - incr1[1]
    ext1[1, 2] <- pres@bbox[1, 2] + incr1[1]
    ext1[2, 1] <- pres@bbox[2, 1] - incr1[2]
    ext1[2, 2] <- pres@bbox[2, 2] + incr1[2]
    if(all(ext1[, 1] > as.vector(vrbles@extent)[c(1, 3)]) & all(ext1[, 2] < as.vector(vrbles@extent)[c(2, 4)])){
      varbles1 <<- raster::stack(raster::crop(vrbles, ext1))
    }else{
      varbles1 <<- raster::stack(raster::crop(vrbles, pres@bbox))
    }

    # number of background points (see Guevara et al, 2017)
    num_bckgr1 <- floor((varbles1@ncols * varbles1@nrows) * 50/100)
    # background points
    bckgr_pts1 <- dismo::randomPoints(varbles1[[1]], num_bckgr1, pres)


    #### Making models for each buffer ####
    #tables with info to be exported
    dt2exp <- as.data.frame(matrix(ncol = 12, nrow = 0))
    selfinfo2exp <- as.data.frame(matrix(ncol = 9, nrow = 0))
    dt2exp_mean <- as.data.frame(matrix(ncol = 7, nrow = 0))

    for (bdw in 1:length(bndwidth)) { # for each buffer
      x <- 1
      # set of presences for modeling within the buffer
      pres4model <- pres[pres$dist2centr <= bndwidth[bdw], ]

      # croping variables to pres4model extent  + 5%  <-- to fit the model
      ext <- pres4model@bbox
      incr <- apply(ext, 1, function(x) x[2] - x[1]) * 0.05
      ext[1, 1] <- pres4model@bbox[1, 1] - incr[1]
      ext[1, 2] <- pres4model@bbox[1, 2] + incr[1]
      ext[2, 1] <- pres4model@bbox[2, 1] - incr[2]
      ext[2, 2] <- pres4model@bbox[2, 2] + incr[2]
      if(all(ext[, 1] > as.vector(vrbles@extent)[c(1, 3)]) & all(ext[, 2] < as.vector(vrbles@extent)[c(2, 4)])){
        varbles2 <- raster::stack(raster::crop(vrbles, ext))
      }else{
        varbles2 <- raster::stack(raster::crop(vrbles, pres4model@bbox))
      }

      # number of background points (see Guevara et al, 2017)
      num_bckgr <- floor((varbles2@ncols * varbles2@nrows) * 50/100)
      #if(num_bckgr<100) {
      #  pres4model1 <- sample(1:nrow(pres4model), nrow(pres4model)*0.1)
      #  pres4model <- pres4model[pres4model1,]
      #}

      # sampling presences for calibrating and testing (70-30%) within the buffer
      folds <- sample(1:nrow(pres4model), nrow(pres4model)*(1-occ_prop_test))
      samp <- as.numeric(unlist(folds))
      pres4cali <- pres4model[samp, 1]
      pres4test <- pres4model[-samp, 1]

      # background points
      bckgr_pts <- dismo::randomPoints(varbles2, num_bckgr, pres4model)

      #rm(pres4model); gc()

      # sampling presences for testing on the whole extent (30% of total presences except those for calibrating)
      pres1 <- pres[-samp, 1]
      folds1 <- sample(1:nrow(pres1), nrow(pres1)*(occ_prop_test))
      samp1 <- as.numeric(unlist(folds1))
      pres4test_tot <- pres1[-samp1, 1]

      repeat{   # maybe it can be done directly with maxent; if so, we would also have the "average-model"
        t1 <- Sys.time()
        message("\r", "modelling for ", specs_long," - buffer #", bdw, "-", x)

        # Running maxent from dismo or maxnet
        if(!file.exists(paste0(dir2save,"/results_", sps))) dir.create(paste0(dir2save,"/results_", sps))
        if(!file.exists(paste0(dir2save,"/results_", sps, "/model_", sps, "_", bdw, "_", x))) dir.create(paste0(dir2save,"/results_", sps, "/model_", sps, "_", bdw, "_", x))
        path <- paste0(dir2save,"/results_", sps,"/model_", sps, "_", bdw, "_", x)

        dir_func <- function(varbles2, pres4cali, num_bckgr, bckgr_pts, path){ # to avoid stop modelling if low number of background points or other errors
          res <- tryCatch(
            {
              data_train <- raster::extract(varbles2, pres4cali)
              if(maxent_tool == "dismo"){
                #modl_dismo <- dismo::maxent(varbles2, pres4cali, nbg = num_bckgr)
                modl <- dismo::maxent(varbles2, pres4cali, a = bckgr_pts)
              }else if(maxent_tool == "maxnet"){
                bckgr_train <- raster::extract(varbles2, bckgr_pts)
                pres_abs <- c(rep(1, nrow(data_train)), rep(0, nrow(bckgr_train)))
                data_model <- data.frame(cbind(pres_abs, rbind(data_train, bckgr_train)))
                data_model <- data_model[stats::complete.cases(data_model), ]
                modl <- maxnet::maxnet(data_model[, 1], data_model[, - 1],
                                       f = maxnet::maxnet.formula(p = data_model[, 1],
                                                                  data = data_model[, - 1],
                                                                  classes = "default"))
              }

              if(exists("modl")) save(modl, file = paste0(path, "/model.RData"))
            },
            error = function(con){
              message(con)
              return(NULL)
            }
          )
          if(exists("modl")){ return(list(modl, data_train, varbles2)) }else{ return(NULL) }
        } #end of dir_func

        modl <- dir_func(varbles2, pres4cali, num_bckgr, bckgr_pts, path)
        data_train <- modl[[2]]
        varbles2 <- modl[[3]]
        modl <- modl[[1]]
        if(is.null(modl)){ break }

        #making predictions on the same extent
        if(maxent_tool == "dismo"){
          preds <- dismo::predict(modl, varbles2, args = 'outputformat=logistic',
                                  filename = paste0(path, "/predictions"), progress = '',
                                  overwrite = TRUE)
        }else if(maxent_tool == "maxnet"){
          pres4test$tovalidate <- 1
          varbles2test <- raster::rasterize(pres4test, varbles2[[1]], pres4test$tovalidate)
          #varbles2predict <- as.data.frame(matrix(nrow = length(varbles2[[1]]@data@values), ncol = dim(varbles2)[3]))
          varbles2predict <- as.data.frame(matrix(nrow = length(raster::getValues(varbles2[[1]])), ncol = dim(varbles2)[3]))
          names(varbles2predict) <- colnames(data_train)
          for(i in 1:dim(varbles2)[3]){
            #varbles2predict[, i] <- varbles2[[i]]@data@values
            varbles2predict[, i] <- raster::getValues(varbles2[[i]])
          }
          #varbles2predict$tovalidate <- varbles2test@data@values
          varbles2predict$tovalidate <- raster::getValues(varbles2test)
          varbles2predict$tovalidate[is.na(varbles2predict$tovalidate)] <- 0
          varbles2predict <- varbles2predict[stats::complete.cases(varbles2predict), ]
          varbles2predict$preds_maxnet <- predict(modl,
                                                  varbles2predict[, - length(varbles2predict)],
                                                  clamp = TRUE,
                                                  type = c("logistic"))
        }

        #make evaluations (on the same extent with 30% to test)
        evs <- dismo::evaluate(modl, p = pres4test, a = bckgr_pts, x = varbles2)
        save(evs, file = paste0(path, "/evaluations.RData"))
        #rm(varbles2)
        #gc()
        graphics.off()

        #Computing Boyce Index (on the same extent with 30% to test)
        if(maxent_tool == "dismo"){
          byce <- ecospat::ecospat.boyce(fit = preds, obs = pres4test@coords, nclass=0, window.w="default", res=100, PEplot = TRUE)
          byce$Spearman.cor
        }else if(maxent_tool == "maxnet"){
          byce <- ecospat::ecospat.boyce(fit = varbles2predict$preds_maxnet, obs = varbles2predict[varbles2predict$tovalidate == 1, ]$preds_maxnet, nclass=0, window.w="default", res=100, PEplot = TRUE)
          byce$Spearman.cor
        }

        save(byce, file = paste0(path, "/boyce.RData"))

        # In the last buffer it would make no sense repeating predictions/evaluations on the same extent
        # to assess for transferability. However, for the sake of consistency in the execution time, they are
        # calculated.
        # To return to the version where they are not calculated, check commit e6e0040 of 25/08/2018

        ## making predictions on the whole species extent
        if(maxent_tool == "dismo"){
          preds1 <- dismo::predict(modl, varbles1, args = 'outputformat=logistic',
                                   filename = paste0(path, "/predictions_tot"), progress = '',
                                   overwrite = TRUE)
        }else if(maxent_tool == "maxnet"){
          pres4test_tot$tovalidate <- 1
          varbles1 <- varbles1
          varbles2test <- raster::rasterize(pres4test_tot, varbles1[[1]], pres4test_tot$tovalidate)
          #varbles2predict <- as.data.frame(matrix(nrow = length(varbles1[[1]]@data@values), ncol = dim(varbles1)[3]))
          varbles2predict <- as.data.frame(matrix(nrow = length(raster::getValues(varbles1[[1]])), ncol = dim(varbles1)[3]))
          names(varbles2predict) <- colnames(data_train)
          for(i in 1:dim(varbles1)[3]){
            #varbles2predict[, i] <- varbles1[[i]]@data@values
            varbles2predict[, i] <- raster::getValues(varbles1[[i]])
          }
          #varbles2predict$tovalidate <- varbles2test@data@values
          varbles2predict$tovalidate <- raster::getValues(varbles2test)
          varbles2predict$tovalidate[is.na(varbles2predict$tovalidate)] <- 0
          varbles2predict <- varbles2predict[stats::complete.cases(varbles2predict), ]
          varbles2predict$preds_maxnet <- predict(modl,
                                                  varbles2predict[, - length(varbles2predict)],
                                                  clamp = TRUE,
                                                  type = c("logistic"))
        }

        #Make evaluations
        evs1 <- dismo::evaluate(modl, p = pres4test_tot, a = bckgr_pts1, x = varbles1)
        save(evs1, file = paste0(path, "/evaluations_tot.RData"))

        #MESS map
        reference_points <- raster::extract(varbles2, pres4model)
        mss <- dismo::mess(x = varbles1, v = reference_points, full = FALSE, filename = paste0(path, "/MESS_map.tif"))

        #Computing Boyce Index
        if(maxent_tool == "dismo"){
          byce1 <- ecospat::ecospat.boyce(fit = preds1, obs = pres4test_tot@coords, nclass=0, window.w="default", res=100, PEplot = TRUE)
          byce1$Spearman.cor
        }else if(maxent_tool == "maxnet"){
          byce1 <- ecospat::ecospat.boyce(fit = varbles2predict$preds_maxnet, obs = varbles2predict[varbles2predict$tovalidate == 1, ]$preds_maxnet, nclass=0, window.w="default", res=100, PEplot = TRUE)
          byce1$Spearman.cor
        }
        save(byce1, file = paste0(path, "/boyce_tot.RData"))

        # gathering info to be exported
        t2 <- Sys.time() - t1
        if(attr(t2, "units") == "hours") {t2 <- t2*60; attr(t2, "units") <- "mins"}
        if(attr(t2, "units") == "secs") {t2 <- t2/60; attr(t2, "units") <- "mins"}
        if(maxent_tool == "dismo"){
          num_pres_calib_used <- nrow(modl@presence)
        }else if(maxent_tool == "maxnet"){
          num_pres_calib_used <- nrow(data_train)
        }
        dt2exp_2 <- as.data.frame(matrix(c(specs_long, paste(bdw, x, sep="_"), bndwidth[bdw], num_pres_calib_used, evs@np, num_bckgr, evs@auc, byce$Spearman.cor, evs1@np, num_bckgr1, evs1@auc, byce1$Spearman.cor), 1, 12, byrow = TRUE))
        dt2exp <- rbind(dt2exp, dt2exp_2)
        selfinfo2exp_2 <- as.data.frame(matrix(c(specs_long, paste(bdw, x, sep="_"), nrow(pres4cali), num_pres_calib_used, nrow(pres4test), evs@np, num_bckgr, evs@na, t2), 1, 9, byrow = TRUE))
        selfinfo2exp <- rbind(selfinfo2exp, selfinfo2exp_2)

        if (x == n_rep){ break }else{ x <- x +1 }
      } #end of repeat n times

      if(is.null(modl)){ message("\r", "jumping to next buffer"); next }

      message("computing average for ", specs_long, " - buffer #", bdw, "\n")
      dt2exp[,-c(1:6)] <- data.frame(lapply(dt2exp[-c(1:6)], function(x) as.numeric(as.character(x))))
      dt2exp_m <- mean(dt2exp[(nrow(dt2exp)-n_rep+1):nrow(dt2exp), (ncol(dt2exp)-4)], na.rm = TRUE) #mean Boyce partial area
      dt2exp_m2 <- mean(dt2exp[(nrow(dt2exp)-n_rep+1):nrow(dt2exp), ncol(dt2exp)], na.rm = TRUE) #mean Boyce whole area
      if (bdw > 3){
        dt2exp_sd1 <- sd(c(dt2exp_mean$V3[(bdw-3):(bdw-1)], dt2exp_m), na.rm = TRUE)
        dt2exp_sd2 <- sd(c(dt2exp_mean$V4[(bdw-3):(bdw-1)], dt2exp_m2), na.rm = TRUE)
      }else{
        dt2exp_sd1 <- NA
        dt2exp_sd2 <- NA
      }
      selfinfo2exp[,-c(1:8)] <- data.frame(lapply(selfinfo2exp[-c(1:8)], function(x) as.numeric(as.character(x))))
      dt2exp_m1 <- mean(selfinfo2exp[(nrow(selfinfo2exp)-n_rep+1):nrow(selfinfo2exp), ncol(selfinfo2exp)], na.rm = TRUE)
      dt2exp_mean_2 <- as.data.frame(matrix(c(specs_long, bndwidth[bdw], dt2exp_m, dt2exp_m2, dt2exp_sd1, dt2exp_sd2, dt2exp_m1), 1, 7, byrow = TRUE))
      dt2exp_mean <- rbind(dt2exp_mean, dt2exp_mean_2)
      dt2exp_mean[, c(2:7)] <- data.frame(lapply(dt2exp_mean[c(2:7)], function(x) as.numeric(as.character(x))))

      rm(modl, evs, byce); gc()

      #Conditions to stop the process
      if(!is.null(BI_part) | !is.null(BI_tot) | !is.null(SD_BI_part) | !is.null(SD_BI_tot)){
        brk <- 0
        if (!is.null(BI_part) && BI_part <= dt2exp_m){ message("\r", "minimum BI_part has been reached"); brk <- 1}
        if (!is.null(BI_tot) && BI_tot <= dt2exp_m2){ message("\r", "minimum BI_tot has been reached"); brk <- 1}
        if (!is.null(SD_BI_part) && !is.na(dt2exp_sd1) && SD_BI_part >= dt2exp_sd1){ message("\r", "minimum SD_BI_part has been reached"); brk <- 1}
        if (!is.null(SD_BI_tot) && !is.na(dt2exp_sd2) && SD_BI_tot >= dt2exp_sd2){ message("\r", "minimum SD_BI_tot has been reached"); brk <- 1}
        if (brk == 1)  break
      }

    } # end of for each buffer

    names(dt2exp_mean) <- c("Species", "Buffer", "BoyceIndex_part", "BoyceIndex_tot", "SD_part", "SD_tot", "ExecutionTime")
    names(dt2exp) <- c("Species", "ModelNum", "Buffer", "numPresencesCalib", "numPresencesTest", "numBackground", "AUC_part", "BoyceIndex", "numPresencesTest_tot", "numBackground_tot", "AUC_tot", "BoyceIndex_tot")
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
      names(best2_bnd) <- c("Species", "Best_Buffer_NoTime", "SecondBest_Buffer_NoTime", "Best_Buffer_WithTime", "SecondBest_Buffer_WithTime")

      best2_bnd_2exp <- rbind(best2_bnd_2exp, best2_bnd)
      write.csv(best2_bnd_2exp, paste0(dir2save, "/rankingBestBuffer.csv"), row.names = FALSE)
    }

    write.csv(dt2exp_mean, paste0(dir2save, "/results_", sps, "/info_mod_means_", sps, ".csv"), row.names = FALSE)
    write.csv(dt2exp, paste0(dir2save, "/results_", sps, "/info_mod_", sps, ".csv"), row.names = FALSE)
    write.csv(selfinfo2exp, paste0(dir2save, "/results_", sps, "/selfinfo_mod_", sps, ".csv"), row.names = FALSE)

    #### Making a plot ####
    #graphics.off()
    dt2exp_mean[,names(dt2exp_mean) %in% c("BoyceIndex_part", "BoyceIndex_tot")] <- round(dt2exp_mean[,names(dt2exp_mean) %in% c("BoyceIndex_part", "BoyceIndex_tot")], 3)
    pdf(paste0(dir2save, "/results_", sps, "/boyce_buffer_", sps, "_part_tot.pdf"))
    #if(nrow(dt2exp_mean) < 5){ tp <- c("p") }else{ tp <- c("p", "smooth") }
    if(nrow(dt2exp_mean) < 5){ tp <- c("p") }else{ tp <- c("p", "l") }
    plt <- lattice::xyplot(BoyceIndex_part ~ Buffer, dt2exp_mean,
                           type = tp,
                           #span = 0.8,
                           ylim = c(0.45, 1.05),
                           col = "blue",
                           main = bquote(Boyce~Index~(mean~of~.(n_rep)~models)~-~italic(.(specs_long))),
                           ylab = "Boyce Index", xlab = "Buffer (km)",
                           key=list(#space = "right",
                             x=0.5,y=0.2,
                             lines = list(col=c("blue", "green", "magenta")),
                             text = list(c("Boyce Index Partial","Boyce Index Total", "Execution Time"))))
    plt1 <- lattice::xyplot(ExecutionTime ~ Buffer, dt2exp_mean,
                            type = c("p", "l"),
                            ylab = "Execution Time (min)",
                            col = "magenta")
    dbl_plt <- latticeExtra::doubleYScale(plt, plt1, add.ylab2 = TRUE)
    plt2 <- lattice::xyplot(BoyceIndex_tot ~ Buffer, dt2exp_mean,
                            type = tp,
                            #span = 0.8,
                            col = "green")
    plot(dbl_plt + latticeExtra::as.layer(plt2))
    dev.off()
  } # end of loop for sps
}
