getCNVInRegion<-function(covs_db, 
                                    region=NULL,
                                    line="",
                                    max_gap = 2,
                                    window = "", 
                                    sd_factor = 3, 
                                    het=TRUE
                                  ){
    
    windows <- getRegionFromDB(covs_db, region=region)
    covs    <- getWindowsInRange(covs_db, line=line, region=region)
    windows$limit <- windows$sd * sd_factor  
    windows$min_sd <- 1 - windows$limit 
    windows$max_sd <- 1 + windows$limit

    covs$max_sd <- windows$max_sd
    covs$out_of_range <- covs$norm_cov > windows$max_sd | covs$norm_cov < windows$min_sd
    covs$stiched <- FALSE
    current_stretch <- c()
    current_gap <- 0
    for(i in seq(from=1, to=length(covs), by=1)){
      if( covs[i]$out_of_range == FALSE ){
        if(length((current_stretch)) > 0){
          current_gap <- current_gap + 1
          if(current_gap > max_gap){            
            covs[
              seq(from=current_stretch[1], to = current_stretch[length(current_stretch)], by=1)
            ]$stiched <- TRUE
            current_stretch <- c()
            current_gap <- 0
          }
        }
        next
      }
      current_stretch <- append(current_stretch, i)

    }
    if(current_gap > max_gap){            
      covs[
        seq(from=current_stretch[1], to = current_stretch[length(current_stretch)], by=1)
      ]$stiched <- TRUE
    }
    covs <- covs[covs$stiched] 
    reduced <- reduce(covs)
    end(reduced) <- end(reduced) - 1
    end(covs) <- end(covs) - 1
    reduced<-binnedAverage(reduced, mcolAsRleList(covs, "norm_cov"), "norm_cov")
    if(het){
      reduced$cnv_level <- round(reduced$norm_cov * 2) / 2
    }else {
      reduced$cnv_level <- round(reduced$norm_cov ) 
    }
    reduced
}