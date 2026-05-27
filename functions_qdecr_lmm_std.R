#######################################################################################
############## FUNCTIONS for running vertex-wise linear mixed models ##################
#######################################################################################

# NOTE: This version provides the standardized betas and SE

################ load the data   
loadvw <- function(df_long, subDir, hemi, measure){
  # get dimensions
  nSubs <- nrow(df_long) 
  nVert <- QDECR::load.mgh(file.path(subDir, df_long$folders[1], 'surf', paste(hemi, measure, 'fwhm10.fsaverage.mgh', sep='.')))$ndim1 #08-11 ; still need to make this flexible
  
  cat("Number of vertices: ", nVert, '\n')
  cat("Number of subjects: ", nSubs, '\n')
  
  # set up empty matrix for storing vertices
  inputData <- matrix(NA, nrow = nVert, ncol=nSubs)
  # progress bar for loading vertex data
  pb <- txtProgressBar(min = 1,      # Minimum value of the progress bar
                       max = nSubs,  # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar
  
  cat("Loading vertex data:", as.character(Sys.time()), '\n')
  
  for (j in 1:nSubs){
    iFile <- file.path(subDir, df_long$folders[j], 'surf', paste(hemi, measure, 'fwhm10.fsaverage.mgh', sep='.'))

    if (!file.exists(iFile)){
      cat('cannot find:  ', iFile, '\n')
      inputData[,j] <- NA
    } else {
      inputData[,j] <- QDECR::load.mgh(iFile)$x
      setTxtProgressBar(pb, j)
    }
    
  }
  
  cat("Finished Loading vertex data:", as.character(Sys.time()), '\n')
  cat("input array is:")
  cat('ARRAY IS:   ', object.size(inputData), units = 'GB', '\n')
  
  pb <- txtProgressBar(min = 1,      # Minimum value of the progress bar
                       max = nVert, # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar
  return(inputData)
}


################lme functions for mice objects

run_lme_mice_quick <- function(x) {
    if (!all(x == 0)){
      #empty list for qhat and se
      qhat <- list()
      se <- list()
      residual <- list()
      for (i in 1:length(implist)) {  
        df <- implist[[i]]
        df$y <- x
        fm <- lmer(as.formula(model), data = df)
        fmsum <- summary(fm)
        fm_residual <- as.vector(residuals(fm))
        
        if (is.null(stackofinterest)) { 
          qhat[[i]] <- c(coef(fmsum)[, 'Estimate'])
          se[[i]] <- c(coef(fmsum)[,'Std. Error'])
          residual[[i]] <- fm_residual
        } else { 
          qhat[[i]] <- c(coef(fmsum)[stackofinterest, 'Estimate'])
          se[[i]] <- c(coef(fmsum)[stackofinterest, 'Std. Error'])
          residual[[i]] <- fm_residual
        }
      }
      # average the residual across imputed datasets
      residual_average <- colMeans(do.call(rbind, residual))
      #grab out the stats you care about
      Stats <- quick_pool2(qhat, se, residual_average)
    }
    else{ 
      Stats <- array(0, dim=c(nBricks, 1))
    }
    return(Stats)
  }
 
 
  
######## pooling the results of mice objects

quick_pool2 <- function (qhat, se, residual_average) {

    ### REWRITE OF mice::pool AND mice::mice_df to also handle multiple outcomes
    #set empty array to store output
    Stats <- array(NA, dim=c(nBricks, 1))
    qhat <- lapply(qhat, as.matrix)
    se <- lapply(se, as.matrix)
    
    eps <- 1e-100
    m <- length(qhat)
    dd <- dim(qhat[[1]])
    im <- (1 + 1/m)
    um <- Reduce("+", lapply(se, `^`, 2)) / m
    qbar <- Reduce("+", qhat) / m 
    qhat2 <- array(unlist(qhat), dim = c(dd, m))
    e <- sweep(qhat2, c(1,2),  qbar, `-`) 
    bm <- apply(e^2, c(1,2), sum)/(m - 1 + eps)
    t <- um + im * bm
    se2 <- sqrt(t)
    tval <- qbar/se2
    lambda <- im * (bm/t)
    eps2 <- 1e-04
    dfcom <- 1e+07
    lambda[lambda < eps2] <- eps2
    dfold <- (m - 1 + eps2)/lambda^2
    dfobs <- (dfcom + 1)/(dfcom + 3) * dfcom * (1 - lambda)
    df <- dfold * dfobs/(dfold + dfobs)
    
    pval <- 2 * stats::pt(-abs(tval), df = df)
    #Stats[1:(nStacks), 1] <- qbar[,1]
    #Stats[(nStacks+1):(nStacks*2), 1]  <- se2[,1]  
    #Stats[((nStacks*2)+1):(nStacks*3), 1] <- pval[,1]
    Stats[, 1] <- c(qbar[,1], se2[, 1], pval[, 1], residual_average)
    return(Stats)
  
  }


########## lme functions for a single data frame

run_lme_single_quick <- function(x) {

       if (!all(x == 0)){
        df_long$y <- x
        dim(df_long)
        
        fm <- lmer(as.formula(model), data = df_long)
        
        fm_residual <- as.vector(residuals(fm))
        
        # standardized coef and se
        stdCoef.merMod <- function(object) {
               sdy <- sd(lme4::getME(object,"y"))
               sdx <- apply(lme4::getME(object,"X"), 2, sd)
               sc <- fixef(object)*sdx/sdy
               se.fixef <- coef(summary(object))[,"Std. Error"]
               se <- se.fixef*sdx/sdy
               return(data.frame(stdcoef=sc, stdse=se))
        }
        
        anovaInfo <- Anova(fm, type=3)
      
        if (is.null(stackofinterest)) { 
          res <-  c(stdCoef.merMod(fm)[, "stdcoef"], 
                    stdCoef.merMod(fm)[, "stdse"], 
                    -log10(anovaInfo[,'Pr(>Chisq)']),
                    fm_residual)
        } else { 
          res <-  c(stdCoef.merMod(fm)[stackofinterest,"stdcoef"], 
                    stdCoef.merMod(fm)[stackofinterest, "stdse"], 
                    -log10(anovaInfo[,'Pr(>Chisq)']),
                    fm_residual)
        }
     
        Stats <- array(res, dim=c(nBricks, 1))
     } else {
        Stats <- array(0, dim=c(nBricks, 1))
     }
    
    return(Stats)
}



############  save out residuals as mgh file

residual2mgh <-function(residual, fname, filter = NULL) {
  MRI.UCHAR <-  0
  MRI.INT <-    1
  MRI.LONG <-   2
  MRI.FLOAT <-  3
  MRI.SHORT <-  4
  MRI.BITMAP <- 5
  MRI.TENSOR <- 6
  slices <- c(1:256)

  fid <- file(fname, open = "wb", blocking = TRUE)

  if (!is.null(filter)){
    if (max(filter) > ncol(residual)) stop("In bsfbm2mgh, the defined filter exceeds the bounds.")
    if (min(filter) < 1) stop("In bsfbm2mgh, the defined filter contains non-positive numbers.")
  } else {
    filter <- seq_len(ncol(residual))
  }
  width <- nrow(residual)
  height <- 1
  depth <- 1
  nframes <- length(filter)
  writeBin(as.integer(1), fid, size = 4, endian = "big")
  writeBin(as.integer(width), fid, size = 4, endian = "big")
  writeBin(as.integer(height), fid, size = 4, endian = "big")
  writeBin(as.integer(depth), fid, size = 4, endian = "big")
  writeBin(as.integer(nframes), fid, size = 4, endian = "big")
  writeBin(as.integer(MRI.FLOAT), fid, size = 4, endian = "big")
  writeBin(as.integer(1), fid, size = 4, endian = "big")
  UNUSED.SPACE.SIZE <- 256
  USED.SPACE.SIZE <- (3 * 4 + 4 * 3 * 4)
  unused.space.size <- UNUSED.SPACE.SIZE - 2
  writeBin(as.integer(0), fid, size = 2, endian = "big")
  writeBin(as.integer(rep.int(0, unused.space.size)), fid, size = 1)
  bpv <- 4
  nelts <- width * height
  for (i in filter) writeBin(residual[, i], fid, size = 4, endian = "big")
  close(fid)
  NULL
}

############### calculate the smoothness

calc_fwhm <- function(final_mask_path, est_fwhm_path, hemi, eres, mask = NULL, target = "fsaverage", verbose = TRUE) {
  cmdStr <- paste("mris_fwhm", "--i", eres, "--hemi", hemi, "--subject", target, "--prune", "--cortex", "--dat", est_fwhm_path, "--out-mask", final_mask_path)
  if (!is.null(mask)) paste(cmdStr, "--mask", mask)
  system(cmdStr, ignore.stdout = !verbose)
  fwhm <- utils::read.table(est_fwhm_path)
  fwhm <- round(fwhm)
  return(fwhm)
}


############### cluster wise correction
make_mri_surf_cluster <- function(hemi, pval, fwhm, mask_path = NULL, cwp_thr = 0.025, mcz_thr = 30, csd_sign = "abs", verbose = FALSE, var_name, measure) {

  mcz_thr2 <- paste0("th", mcz_thr)
  if (fwhm < 10) fwhm <- paste0("0", fwhm)
  
  csd <- file.path(Sys.getenv("FREESURFER_HOME"), "average/mult-comp-cor/fsaverage", hemi, paste0("cortex/fwhm", fwhm), "abs", mcz_thr2, "mc-z.csd")
  
  cmd_str <- paste("mri_surfcluster",
                   "--in", pval,
                   "--csd", csd, 
                   "--cwsig", paste0(var_name, "/", hemi,".", measure,".", var_name,".cluster.mgh"),
                   #"--vwsig", paste0(var_name, ".voxel.mgh"),
                   "--sum", paste0(var_name, "/", hemi,".", measure,".",var_name,".cluster.summary"),
                   "--ocn", paste0(var_name, "/", hemi,".", measure,".",var_name,".ocn.mgh"),
                   "--oannot", paste0(var_name, "/", hemi,".", measure,".",var_name,".ocn.annot"),
                   "--annot aparc",
                   "--cwpvalthresh", cwp_thr,
                   "--o", paste0(var_name, "/", hemi,".",measure,".",var_name,".masked.mgh"),
                   "--no-fixmni",
                   "--surf", "white")
                   
  cmd_str <- if (!is.null(mask_path)) paste(cmd_str, "--mask", mask_path) else paste(cmd_str, "--cortex")
  
  system(cmd_str)
  NULL
}



######################################################################################################################################
######################################  merged all the functions  into one single function ###########################################
######################################################################################################################################

qdecr_quicklmm <- function(input_data, subDir, hemi, measure, model, n_cores, stackofinterest=NULL) {

      # if mids object get single set, needed for obtaining info
      if (class(input_data) == "mids") {
        df_long <- complete(input_data, 1)
        implist <- lapply(seq_len(input_data$m), function(y) mice::complete(input_data, y))
      } else {
        df_long <- input_data
      }

      # read the data with dimension: nVert (number of vertices) * nSubj (number of subjects)
      inputData <- loadvw(df_long, subDir, hemi, measure)
      nVert <- nrow(inputData)
      nSubs <- ncol(inputData)
      mask <- QDECR::load.mgh(paste0('/gpfs/work2/0/einf1049/scratch/bxu/PGS_brain/',hemi, '.fsaverage.cortex.mask.mgh')) # change the path to mask
      mask <- mask$x
      inputData <- array(apply(inputData, 2, function(x) x*mask), dim=c(nVert, nSubs))

      
      ### STEP 1: Get all dimensions for storing output correctly
      # run one model on one set to get dimensions
      i=1  
      df_long$y <- inputData[i, ]
      fmTest <- lmer(model, data = df_long)
      fmFullSummary <- summary(fmTest)
      
      anovaTest <- Anova(fmTest, type=3)
      n_anova <- length(anovaTest[,'Pr(>Chisq)'])
      
      coefinfo <- data.frame(number = 1:length(row.names(coef(fmFullSummary))), coefname = row.names(coef(fmFullSummary)))
                               
      # number of coefficients
      if (is.null(stackofinterest)) { 
        nStacks <- length(row.names(coef(fmFullSummary)))
      } else { 
        nStacks <- length(stackofinterest)
      }
      # n of vertices
      nVert <- nrow(inputData)
      # n of subjects
      nSubs <- nrow(df_long)
      
      # split into 20 segments for efficiency
      nSeg <- 20
      # break into 20 segments
      dim_n <- nVert%/%nSeg + 1
      # number of datasets need to be filled
      to_fill <- nSeg - nVert%%nSeg #reminder of the nVert     
      # pad with extra 0s
      inData <- rbind(inputData, array(0, dim=c(to_fill, nSubs)))
      # break into 20 chunks
      dim(inData) <- c(dim_n, nSeg, nSubs)
      # number of coeff * number of results (b, se and p in this case)
      
      # if mids object get single set, needed for obtaining info
      if (class(input_data) == "mids") {
        nBricks <- nStacks * 3 + nSubs
      } else {
        nBricks <- nStacks * 2 + n_anova + nSubs # n of betas + n of p + n of residuals
      }
      
      #create an empty array to hold the stats
      Stat <- array(NA, dim=c(dim_n, nSeg, nBricks))
      
      
      ### STEP 2: Parallelization and running models
                               
      cl <- makeCluster(n_cores, type = "SOCK")
      # parellization 
      if (class(input_data) == "mids") {
        var_to_export <- list("nStacks","nBricks","dim_n","implist", "run_lme_mice_quick", "quick_pool2", 
                              "stackofinterest", "model")
      } else {
        var_to_export <- list("nStacks","nBricks","dim_n","df_long","run_lme_single_quick","stackofinterest", 
                              "model")
      }
      environment_var_to_export <- new.env()
      environment_var_to_export$var_to_export <- var_to_export
      clusterExport(cl, var_to_export, envir = environment_var_to_export)
      clusterEvalQ(cl, library("lme4"))
      clusterEvalQ(cl, library("car"))
      
       for(s in 1:nSeg) {
        # run the models in each chunk
        # if mids object get single set, needed for obtaining info
            if (class(input_data) == "mids") {   
              Stat[,s,] <- aperm(parApply(cl, inData[,s,], 1, run_lme_mice_quick), c(2,1)) 
            } else {
              Stat[,s,] <- aperm(parApply(cl, inData[,s,], 1, run_lme_single_quick), c(2,1))  
            }
   
        cat("Computation done ", 100*s/nSeg, "%: ", format(Sys.time(), "%D %H:%M:%OS3"), "\n", sep='')   
      } 
      stopCluster(cl)
      
      # combine all the chunks
      dim(Stat) <- c(dim_n*nSeg, nBricks) 
      Stat <- Stat[-c((dim_n*nSeg - to_fill + 1):(dim_n*nSeg)), drop=F,]
      
      # project date and related measures
      project <- paste0(hemi,"_", measure, "_", Sys.Date())
      
      if (!dir.exists(project)) {
        dir.create(project)
      }
           
      # estimate smoothness
      setwd(project)
      
      if (!dir.exists("fwhm")) {
          dir.create("fwhm")
        }
      
 
      if (class(input_data) == "mids") {
      
      Stats <- Stat[,1:(nBricks-nSubs)]
      Residual <- Stat[,(nBricks-nSubs+1):nBricks]
      
      # convert residual to mgh file
      residual2mgh(Residual, paste(hemi, measure, 'residual', 'mgh', sep = '.'), filter = NULL) 
      
      est_fwhm_path <- paste0("fwhm/","fwhm.dat")
      final_mask_path <- paste0("fwhm/","Finalmask.mgh")
      
      fwhm_est <- calc_fwhm(final_mask_path, est_fwhm_path, hemi, eres=paste(hemi, measure, 'residual', 'mgh', sep = '.'), 
                            mask=NULL, target="fsaverage", verbose=FALSE)
                            
      file.remove(paste(hemi, measure, 'residual', 'mgh', sep = '.'))
           
      # collect mgh files and create folder to store
      testMgh <- QDECR::load.mgh(file.path(subDir, df_long$folders[1], 'surf', paste(hemi, measure, 'fwhm10.fsaverage.mgh', sep ='.')))
      dim(Stats) <- c(nVert, nStacks, 3)
      
      # Iterate over the number of stacks
      for (p in 1:nStacks) {
        estMgh <- testMgh  
        seMgh <- testMgh
        pMgh <- testMgh
        estMgh$x <- Stats[, p, 1]
        seMgh$x <- Stats[, p, 2]
        pMgh$x <- Stats[, p, 3]
        
        # create new folders for each stack
        if (is.null(stackofinterest)) { 
        var_name <- row.names(coef(fmFullSummary))[p]
        } else { 
        var_name <- stackofinterest[p]
        }
        
        if (var_name == "(Intercept)") {
        var_name <- "intercept"
        }
        
        if (!dir.exists(var_name)) {
          dir.create(var_name)
        }
       
        # Define the file names
        estMghFile <- paste(hemi, measure, 'est', var_name, 'mgh', sep = '.')
        seMghFile <- paste(hemi, measure, 'se', var_name, 'mgh', sep = '.')
        pMghFile <- paste(hemi, measure, 'p', var_name, 'mgh', sep = '.')
        
        # Create file paths including the output folder
        estMghFilePath <- file.path(var_name, estMghFile)
        seMghFilePath <- file.path(var_name, seMghFile)
        pMghFilePath <- file.path(var_name, pMghFile)
        
        # Save the files to the specified folder
        QDECR::save.mgh(estMgh, estMghFilePath)
        QDECR::save.mgh(seMgh, seMghFilePath)
        QDECR::save.mgh(pMgh, pMghFilePath)
        
        saveRDS(Stats, paste0("stats.", hemi,".", measure, ".rds"))
        saveRDS(coefinfo, paste0("coefinfo.", hemi,".", measure, ".rds"))
        
        # cluster-wise correction of p values
        make_mri_surf_cluster(hemi,pval=paste0(var_name,"/",pMghFile), fwhm=fwhm_est, 
                              mask_path = NULL, cwp_thr = 0.025, mcz_thr = 30, csd_sign = "abs", 
                              verbose = FALSE, var_name, measure)    
      }
      cat(sprintf("Files have been saved to the '%s' folder.\n", project))
      
      
      } else { # for a single data frame
      
      Stats <- Stat[,1:(nStacks*2)]
      Pvals <- Stat[,(nStacks*2+1):(nStacks*2 + n_anova)]
      Residual <- Stat[,(nBricks-nSubs+1):nBricks]
      
      
      # convert residual to mgh file
      residual2mgh(Residual, paste(hemi, measure, 'residual', 'mgh', sep = '.'), filter = NULL) 
      
      est_fwhm_path <- paste0("fwhm/","fwhm.dat")
      final_mask_path <- paste0("fwhm/","Finalmask.mgh")
      
      fwhm_est <- calc_fwhm(final_mask_path, est_fwhm_path, hemi, eres=paste(hemi, measure, 'residual', 'mgh', sep = '.'), 
                            mask=NULL, target="fsaverage", verbose=TRUE)
     
      
      file.remove(paste(hemi, measure, 'residual', 'mgh', sep = '.'))
           
      # collect mgh files and create folder to store
      testMgh <- QDECR::load.mgh(file.path(subDir, df_long$folders[1], 'surf', paste(hemi, measure, 'fwhm10.fsaverage.mgh', sep ='.')))
      dim(Stats) <- c(nVert, nStacks, 2)
      
    
      # Iterate over the number of stacks
      for (p in 1:nStacks) {
        estMgh <- testMgh  
        seMgh <- testMgh

        estMgh$x <- Stats[, p, 1]
        seMgh$x <- Stats[, p, 2]
        
        # create new folders for each stack
        if (is.null(stackofinterest)) { 
        var_name <- row.names(coef(fmFullSummary))[p]
        } else { 
        var_name <- stackofinterest[p]
        }
        
        if (var_name == "(Intercept)") {
        var_name <- "intercept"
        }
        
        if (!dir.exists(var_name)) {
          dir.create(var_name)
        }
       
        # Define the file names
        estMghFile <- paste(hemi, measure, 'est', var_name, 'mgh', sep = '.')
        seMghFile <- paste(hemi, measure, 'se', var_name, 'mgh', sep = '.')
        
        # Create file paths including the output folder
        estMghFilePath <- file.path(var_name, estMghFile)
        seMghFilePath <- file.path(var_name, seMghFile)
        
        # Save the files to the specified folder
        QDECR::save.mgh(estMgh, estMghFilePath)
        QDECR::save.mgh(seMgh, seMghFilePath)
        
        saveRDS(Stats, paste0("stats.", hemi,".", measure, ".rds"))
        saveRDS(coefinfo, paste0("coefinfo.", hemi,".", measure, ".rds"))
        
       
      }
      
      # create a new folder for all the p values
      # Note: we seperate this because the length of p values from avona is not the same with the lme output
      # which makes the p values cannot correspond to the variables of the interest
      
      if (!dir.exists("P_values")) {
            dir.create("P_values")
      }
      
      setwd("P_values")
      # save out residuals and do cluster wise correction
      for (v in 1:n_anova) {
      
      p_name <- rownames(anovaTest)[v]
      
      if (p_name == "(Intercept)") {
        p_name <- "intercept"
      }  
      
      if (!dir.exists(p_name)) {
          dir.create(p_name)
      }
      
      pMgh <- testMgh
      pMgh$x <- Pvals[, v]
       
      pMghFile <- paste(hemi, measure, 'p', p_name, 'mgh', sep = '.')
        
      # Create file paths including the output folder
      pMghFilePath <- file.path(p_name, pMghFile)
        
      # Save the files to the specified folder
      QDECR::save.mgh(pMgh, pMghFilePath)
      # cluster-wise correction of p values
      make_mri_surf_cluster(hemi,pval=paste0(p_name,"/",pMghFile), fwhm=fwhm_est, 
                            mask_path = NULL, cwp_thr = 0.025, mcz_thr = 30, csd_sign = "abs", 
                            verbose = FALSE, p_name, measure)      
      }
      cat(sprintf("Files have been saved to the '%s' folder.\n", project))  
      }
                          
}










