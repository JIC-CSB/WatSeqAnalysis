getCNVInRegion<-function(covs_db, 
                                    region=NULL,
                                    line="",
                                    max_gap = 2,
                                    window = "", 
                                    sd_factor = 2.5, 
                                    het=TRUE,
                                    max_cov_del = 0.1
                                  ){
    
    windows <- getRegionFromDB(covs_db, region=region)
    covs    <- getWindowsInRange(covs_db, line=line, region=region, as.gr=TRUE)
    windows$limit <- windows$sd * sd_factor  
    windows$min_sd <- 1 - windows$limit 
    windows$max_sd <- 1 + windows$limit

    covs$max_sd <- windows$max_sd
    covs$out_of_range <- covs$norm_cov > windows$max_sd | covs$norm_cov < windows$min_sd
    covs$stiched <- FALSE
    covs$noisy <- windows$min_sd <= max_cov_del
    current_stretch <- c()
    current_gap <- 0
    for(i in seq(from=1, to=length(covs), by=1)){
      if(covs[i]$noisy == TRUE){
        next
      }
      if( covs[i]$out_of_range == FALSE ){
        if(length((current_stretch)) > 0){
          current_gap <- current_gap + 1
          if(current_gap > max_gap){  
             
            to_stich <-  seq( 
            from=current_stretch[1], 
            to = current_stretch[length(current_stretch)], 
            by=1)
            covs[to_stich]$stiched <- TRUE
            # print("stiching...")
            # print(current_stretch)
            # print(to_stich)
            current_stretch <- c()
            current_gap <- 0
          }
        }
        next
      }
      current_stretch <- append(current_stretch, i)
      current_gap <- 0

    }
    if(length(current_stretch) > 0){  
      # print("Last stiching...")
      # print(current_stretch)
      to_stich <-  seq( 
            from=current_stretch[1], 
            to = current_stretch[length(current_stretch)], 
            by=1)          
      # print(to_stich)
      covs[to_stich]$stiched <- TRUE
    }
    #print(as.data.frame(covs))
    #covs_all <- covs
    covs <- covs[covs$stiched] 
    reduced <- reduce(covs)
    end(reduced) <- end(reduced) - 1
    end(covs) <- end(covs) - 1
    # print(reduced)
    reduced<-binnedAverage(reduced, mcolAsRleList(covs, "norm_cov"), "norm_cov")
    if(het){
      reduced$cnv_level <- round(reduced$norm_cov * 2) / 2
    }else {
      reduced$cnv_level <- round(reduced$norm_cov ) 
    }
    # print("Now reduced...")
    # print(reduced)
    reduced$line <- line
    reduced$max_gap <- max_gap
    # c(reduced=reduced, covs=covs_all)
    reduced

}