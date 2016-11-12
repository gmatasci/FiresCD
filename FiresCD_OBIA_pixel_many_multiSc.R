#### CODE INFOS -------------------------------------------------------------

## Project Name: Fire severity classification 
## Authors: Giona Matasci (giona.matasci@gmail.com), Ignacio San Miguel (ignaciosanmiguel86@gmail.com)   
## File Name:                            
## Objective: 

#### TO DO -------------------------------------------------------------------

# PRIOR TO ACTUAL RUN:
# - 

## STILL TO DO:

## SOLVED:
# -V issue with images/projection -- OTB does not accept Lambert Conical when producing a shp, we have to save a sqlite file and then convert it to shp with ogr2ogr
# -V fill ypred df based in indices and not on rbind -- done
# -V read about OTB meanshift minsize -- the 2 other parameters control the size and shape of the segments: minsize is just used as a post-processing to merge segments smaller than a threshold
# -V split aspect in 2 components sin/cos -- new variables: "aspect.sin_mean", "aspect.sin_sd", "aspect.cos_mean", "aspect.cos_sd"
# -V check error: Error in `colnames<-`(`*tmp*`, value = c(NA, NA, NA, "elev_sd")) : length of 'dimnames' [2] not equal to array extent -- not appearing anymore
# -V new pixel level GT with majority class by surface within pixel -- done by Nacho
# -V integrate multiscale approach using seg.ids.df as a starting point -- done, params$multisc allows to specify scale parameters
# -V base parameter selection and assessment on a mixed metric: 0.5*Kappa + 0.5*Fmeasure_PartialMortality -- we now use Kappa
# -V segment all the fires beforehand and store segID in a table to query in the script -- done, we save a lot processing time (when it was inside big loop lots of redundant segmentations)
# -V include stratified RF -- done, controlled through params$use.strata 
# -V uniformize fire names (caps, etc.) -- not needed, with toppuer() everything works fine
# -V add lines to keep order of rows for single.scale.obj.dt: equivalent of "multi.scale.pix.out.dt[, pixIDout:=1:nrow(allpixels.out.dt)]" and "setkey(multi.scale.pix.out.dt, pixIDout)" -- done
# -V profiler to see which part has to go parallel -- parallel RF looses OOB confusion matrix, so not to be used to optimize parameters
# -V introduce LOOCV within meanshift loop instead of OOB? (think if need to redo segmentation or just RF in 2nd LOOCV) -- 


#### INIT --------------------------------------------------------------------

rm(list=ls())

start.message <- sprintf("OBIA for fires CD, started running on %s", Sys.time())
print(start.message)

#### PARAMETERS ---------------------------------------------

params <- list()

## General
                                                      
# params$fires <- c("Abraham", "Cone", "Flett", "Keane", "Leggo", "Levellers", "Liege",
#                   "Mcarther", "Overflow", "Perry", "Rail", "Rainbow", "Steephill", "Tanghe")  ## 14 fires to consider
# params$fires <- c("Abraham", "Tanghe", "Keane")

params$fires <- c("Leggo", "Abraham", "Tanghe")

params$mort.class.names <- as.factor(c("0-5%", "6-25%", "95-100%"))
params$mort.class.labels <- as.factor(c(1, 2, 3))
params$critical.class.label <- 2  ## NOT USED: corresponding to class "6-25%", used to compute F-measure on this class (as an alternative, more detailed measure to Kappa)

## Meanshift
params$meanshift.run <- F    ## if set to T it runs the MS segmentation and saves segIDs, otherwise it loads them from previous run
params$meanshift.save.shp <- T    ## if set to T saves in segments.shp.dir each segmentation result (shp) with parameter values to visually inspect the segments
params$write.maps <- T   ## wheter or not to write raster containing the pixel level prediction for each approach
params$nr.bands.seg <- 3   ## number of bands to use for segmentation (points to a different folder)
params$multisc$ms1 <- data.frame(ranger=c(50), spatialr=c(10), minsize=c(5))  ## parameters for multiscale: each row is a combination of the 3
params$multisc$ms2 <- data.frame(ranger=c(100), spatialr=c(10), minsize=c(5))  ## parameters for multiscale: each row is a combination of the 3
params$multisc$ms3 <- data.frame(ranger=c(200), spatialr=c(10), minsize=c(5))  ## parameters for multiscale: each row is a combination of the 3
# params$multisc$ms4 <- data.frame(ranger=c(300), spatialr=c(10), minsize=c(5))  ## parameters for multiscale: each row is a combination of the 3

# params$multisc$ms5 <- data.frame(ranger=c(50), spatialr=c(5), minsize=c(5))  ## parameters for multiscale: each row is a combination of the 3
# params$multisc$ms6 <- data.frame(ranger=c(100), spatialr=c(5), minsize=c(5))  ## parameters for multiscale: each row is a combination of the 3
# params$multisc$ms7 <- data.frame(ranger=c(200), spatialr=c(5), minsize=c(5))  ## parameters for multiscale: each row is a combination of the 3
# params$multisc$ms8 <- data.frame(ranger=c(300), spatialr=c(5), minsize=c(5))  ## parameters for multiscale: each row is a combination of the 3
# 
# params$multisc$ms9 <- data.frame(ranger=c(50), spatialr=c(20), minsize=c(5))  ## parameters for multiscale: each row is a combination of the 3
# params$multisc$ms10 <- data.frame(ranger=c(100), spatialr=c(20), minsize=c(5))  ## parameters for multiscale: each row is a combination of the 3
# params$multisc$ms11 <- data.frame(ranger=c(200), spatialr=c(20), minsize=c(5))  ## parameters for multiscale: each row is a combination of the 3
# params$multisc$ms12 <- data.frame(ranger=c(300), spatialr=c(20), minsize=c(5))  ## parameters for multiscale: each row is a combination of the 3

# params$multisc$ms2 <- data.frame(ranger=c(100, 200), spatialr=c(10, 10), minsize=c(5, 5)) 
# params$multisc$ms3 <- data.frame(ranger=c(100, 150, 200), spatialr=c(10, 10, 10), minsize=c(5, 5, 5)) 
# params$multisc$ms4 <- data.frame(ranger=c(50, 100, 150, 200), spatialr=c(10, 10, 10, 10), minsize=c(5, 5, 5, 5)) 

## Fill vector with approach names
params$approaches <- c("Pixel")
for (i in 1:length(params$multisc)) {
  params$approaches[1+i] <- sprintf("PixelMultiSc%d", i) 
}

## RF
params$pixel.predictors <- c("db1", "db2", "db3", "db4", "db5", "db7", 
                               "dNBR", "dTCG", "dTCB", "dTCW", "dNDVI", "dNDWI", 
                               "elev", "slope", "aspect.sin", "aspect.cos", 
                               "EOSD")   ## list of starting predictors (for the classification part, as the segmentation is run on 3-band or 6-band difference image only)
params$obj.predictors <- c("db1_mean", "db1_sd", "db2_mean", "db2_sd", "db3_mean", "db3_sd", "db4_mean", "db4_sd", "db5_mean", "db5_sd", "db7_mean", "db7_sd", 
                              "dNBR_mean", "dNBR_sd", "dTCG_mean", "dTCG_sd", "dTCB_mean", "dTCB_sd", "dTCW_mean", "dTCW_sd", "dNDVI_mean", "dNDVI_sd", "dNDWI_mean", "dNDWI_sd",
                              "elev_mean", "elev_sd", "slope_mean", "slope_sd", "aspect.sin_mean", "aspect.sin_sd", "aspect.cos_mean", "aspect.cos_sd",
                              "EOSD_maj", "EOSD_majpct",
                              "nrpixseg")    ## list of final predictors computed at the object level
params$targ <- "CLASS"     ## target variable of the classification (column name of dataframe)
params$seed <- 2016        ## seed to have same RF result
params$use.strata <- T    ## wheter to use strata or not in training the RF, if set to True, RF training will not be run in parallel
params$ntree <- 100    ## RF nr of trees
params$mtry <- 'sqrt_nr_var'  ## how to set RF mtry: 'sqrt_nr_var' or 'nr_var_div_3'
params$nodesize <- 1   ## RF nodesize: default for classification is 1
params$plot.importance <- T  ## whether to plot RF variable importance


## Directories
OTB.dir <- "C:/Users/gmatasci/Downloads/OTB-5.4.0-win64/bin"   ## directory in which OTB is located (input to function meanShiftOTB() )
OGR.dir <- "C:/OSGeo4W64/bin"    ## directory in which OGR is located (input to function meanShiftOTB() ), uses ogr2ogr command for shp conversion from sqlite
base.dir <- 'D:/Research/ANALYSES/FiresCD'    ## base working directory

#### LOAD PACKAGES ----------------------------------------------------------

list.of.packages <- c("lazyeval",
                      "psych",
                      "caret",
                      "e1071",
                      "raster",
                      "rgeos",
                      "rgdal",
                      "sp",
                      "spdep",
                      "spatstat",
                      "gplots",
                      "ggplot2",
                      "plyr",
                      "dplyr", ## to be loaded before foreach to avoid "assertion failed" errors
                      "magrittr",
                      "rlist",
                      "randomForest",
                      "rgl",
                      "vegan",
                      "snow",
                      "lubridate",
                      "doParallel", 
                      "foreach",
                      "data.table"
)
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]   ## named vector members whose name is "Package"
if(length(new.packages)) install.packages(new.packages)   ## install only unavailable packages
for (pack in list.of.packages){
  library(pack, character.only=TRUE)  ## call all the packages in list.of.packages
}

#### FUNCTIONS ----------------------------------------------------------

## runs MeanShift using OTB 
meanShiftOTB <- function(ranger, spatialr, minsize, OTB.dir, OGR.dir, data.dir, segments.shp.dir, img.file.name, seg.file.name, write.shp) {
  
  img.file <- file.path(data.dir, img.file.name)   ## image to segment
  clean.img.file <- file.path(segments.shp.dir, sprintf("ZeroPadded_%s", img.file.name))
  segment.file.sqlt <- file.path(segments.shp.dir, sprintf("%s.sqlite", seg.file.name))  ## initial segmentation result as a sqlite file (only way to work with Lambert conformal conic projection)
  segment.file <- file.path(segments.shp.dir, sprintf("%s.shp", seg.file.name))  ## segmentation result converted to shp via ogr2ogr
  
  if (file.exists(segment.file)) { ## if the file is already in the folder there are 2 options:
    if (write.shp) {  ## if we want full names to be used ("Segment_%s_ranger%d_spatial%d_minsize%d") it means the segmentation is already done
      segments <- readOGR(dsn=segments.shp.dir, layer=seg.file.name)  ## so just read back the shp with segmentation result...
      return(segments)   ## ...and return it so that it can be appended to seg.ids.df and exit function
    } else {   ## or, if we use the generic segment file name ("Segment")...
      unlink(segment.file, recursive = T, force = T)    ## delete the shp files each time bc they will be recreated each time 
      unlink(segment.file.sqlt, recursive = T, force = T)
    }
  }
  
  ## change no-data values to 0 to avoid the 1-segment-per-pixel situation outside the image (and get only a big polygon outside)
  if (!file.exists(clean.img.file)) {  ## run only if padded image does not already exist in temp folder
    command.text <- sprintf("%s -in %s -out %s -mode changevalue -mode.changevalue.newv 0",  ## create text for command to run OTB (equivalent to saving this text in a batch file .bat and running it with system2() )
                            file.path(OTB.dir, "otbcli_ManageNoData"),   ## complete path to OTB bat file doing the value change
                            img.file,
                            clean.img.file
    )
    system(command.text)  ## run command
  }
  
  ## actual segmentation
  command.text <- sprintf("%s -in %s -filter meanshift -filter.meanshift.ranger %d -filter.meanshift.spatialr %d -filter.meanshift.minsize %d -mode vector -mode.vector.out %s -mode.vector.layername %s & %s %s %s", 
                      file.path(OTB.dir, "otbcli_Segmentation"),   ## complete path to OTB bat file doing the segmentation
                      clean.img.file,  
                      ranger, 
                      spatialr,
                      minsize,
                      segment.file.sqlt,
                      seg.file.name,
                      file.path(OGR.dir, "ogr2ogr"),   ## complete path to OGR bat file doing the conversion from sqlite to shp
                      segment.file,
                      segment.file.sqlt
  )
  system(command.text)  ## run command
  
  segments <- readOGR(dsn=segments.shp.dir, layer=seg.file.name)  ## read back shp with segmentation result
  return(segments)
  
} 

## summarizes numerical variables
summarize.custom.num <- function(df, gvar1, gvar2, uniq.var, oper) {
  dt <- df %>%
      group_by_(gvar1, gvar2) %>%
      summarise_(npixel = interp(~get(oper)(number), number=as.name(uniq.var)))
  return(as.data.table(dt))
  
}

## summarizes factor variables
summarize.custom.fac <- function(df, gvar1, gvar2, gvar3, uniq.var, oper) {
  dt <- df %>%
      group_by_(gvar1, gvar2, gvar3) %>%
      summarise_(npixel = interp(~get(oper)(number), number=as.name(uniq.var))) %>%  
      slice(which.max(npixel))
  return(as.data.table(dt))
}

## summarize pixel values at the object level
summarize.all <- function(pixels.dt, cols.ID, cols.coord) {
  
  pixels.dt[, npixel:= 1]   ## add column of ones to be used to count the nr of pixels per segment
  
  ## fill object level df with per segment mean and std dev
  obj.dt <- data.table()  ## initialize df to store object-level features
  firstIter <- T  ## set flag for 1st iteration as true 

  ## loop over the columns to summarize only
  for (col in colnames(pixels.dt)[!(colnames(pixels.dt) %in% c(cols.ID, cols.coord))]){
    if ( is.factor(pixels.dt[[col]]) ){ ## if the column is a factor apply the summarizing function for factors with majority as summarizer
      target <- summarize.custom.fac(pixels.dt, cols.ID[1], cols.ID[2], col, "npixel", "sum")
      colnames(target)[4] <- paste(col, "_npixel", sep="")
    } else if (col %in% c("npixel")){   ## else if it is the columns of ones apply the summarzing function for continuous variables with sum as summarizer
      target <- summarize.custom.num(pixels.dt, cols.ID[1], cols.ID[2], col, "sum")
    } else {  ## else if it is a continuous column apply the summarzing function for continuous variables with mean (final column ending with "_mean") and standard deviation ("_sd") as summarizers
      target1 <- summarize.custom.num(pixels.dt, cols.ID[1], cols.ID[2], col, "mean")
      colnames(target1)[3] <- paste(col, "_mean", sep="") 
      target2 <- summarize.custom.num(pixels.dt, cols.ID[1], cols.ID[2], col, "sd")
      colnames(target2)[3] <- paste(col, "_sd", sep="") 
      target <- data.table(target1, target2[,ncol(target2), with=FALSE])
    }
    if (firstIter){  ## if first round of the loop obj.dt takes the values of the df target
      obj.dt <- rbind(obj.dt, target)
      firstIter <- F  ## set flag as false to go to 2nd part of if-else block from now on
    } else {   ## then, just append the newly computed columns columnwise (cbind)
      newcols <- !(colnames(target) %in% colnames(obj.dt) )
      obj.dt <- cbind(obj.dt, target[, newcols, with=FALSE])
    }
  }
  colnames(obj.dt)[ncol(obj.dt)] <- "nrpixseg"   ## last column is the number of pixels per segment
  obj.dt[, CLASS_npixel:= CLASS_npixel / nrpixseg ] ## compute percentages of the majority class for both... 
  colnames(obj.dt)[colnames(obj.dt) == "CLASS_npixel"] <- "CLASS_majpct"  ## ...the ground truth...
  obj.dt[, EOSD_npixel := EOSD_npixel / nrpixseg]
  colnames(obj.dt)[colnames(obj.dt) == "EOSD_npixel"] <- "EOSD_majpct"  ## ...and the EOSD layer
  colnames(obj.dt)[colnames(obj.dt) == "EOSD"] <- "EOSD_maj"  ## ...and the EOSD layer
  
  ## replace with 0s NA values resulting from single pixel polygons (not possible to compute std dev)
  for(j in seq_along(obj.dt)){
    set(obj.dt, i=which(is.na(obj.dt[[j]])), j=j, value=0)
  }
  
  setkeyv(obj.dt, cols.ID)  ## setkey that accepts vector with strings of column names
  
  pixels.dt[, npixel:= NULL]  ## remove npixel column bc otherwise it remains in the allpixels.dt 
  
  return(obj.dt)

}

## returns classification accuracy metrics
classif.metrics <- function(arg1, arg2) {
  
  if (nargs() == 2) {   ## if the number of arguments is 2, take them as the predicted and observed vectors
    predicted <- arg1
    observed <- arg2
  } else if (nargs() == 1) {  ## if the number of arguments is just 1, it is a df with predicted and observed as columns in it (to apply function to groups with plyr package)
    predicted <- arg1$predicted
    observed <- arg1$observed 
  }
  
  if ( length(levels(observed)) > 2 & length(unique(observed)) > 2) {  ## handle the multi-class case
    RES <- confusionMatrix(predicted, observed)
    PA <- RES$byClass[, "Sensitivity"]  ## Producer's accuracy, 1 - omission error, TP/(observed total)
    UA <- RES$byClass[, "Pos Pred Value"]  ## User's accuracy, 1 - commission error, TP/(predicted total)
  } else if ( length(levels(observed)) == 2 | length(unique(observed)) == 2 ) {  ## handle the binary class case
    RES <- confusionMatrix(predicted, observed, positive=levels(observed)[1])  ## first run assessment setting first class as positive
    PA <- RES$byClass[["Sensitivity"]]  ## save producers accuracy...
    UA <- RES$byClass[["Pos Pred Value"]]  ## ...and users accuracy
    RES <- confusionMatrix(predicted, observed, positive=levels(observed)[2])  ## then repeat with second class
    PA <- c(PA, RES$byClass[["Sensitivity"]])  ## and complete PA vector
    UA <- c(UA, RES$byClass[["Pos Pred Value"]])  ## and UA vector
  } else {
    stop("Observed values: less than 2 levels or levels to be updated")
  }
  
  return(list( Kappa=as.vector(RES$overall[2]), OA=as.vector(RES$overall[1]), Fmeas=(2*PA*UA)/(PA+UA), ConfMat=RES$table))
  
}

#### START --------------------------------------------------------------

tic <- proc.time() ## start clocking global time

data.dir <- file.path(base.dir, "Data", fsep = .Platform$file.sep)   ## data directory (nothing should be written here)
results.dir <- file.path(base.dir, "Results", fsep = .Platform$file.sep)   ## results directory (outputs go here)
figures.dir <- file.path(base.dir, "Figures", fsep = .Platform$file.sep)   ## figures directory (figures go here)
segments.shp.dir <- file.path(data.dir, "shp_segments")  ## directory for temporary files like the segmentation shps (overwritten each time)
if (!file.exists(segments.shp.dir)) {dir.create(segments.shp.dir, showWarnings=F, recursive=T)}  ## create it

PARAMS.file <- file.path(results.dir, 'PARAMS.Rdata', fsep = .Platform$file.sep)
save(params, file=PARAMS.file)

load(file.path(data.dir, "allpixels.Rdata"))
allpixels.dt <- dfp
rm(dfp)

pix.fires.OK <- allpixels.dt[["NAME"]] %in% toupper(params$fires) ## get indices of pixels beloning to subset of fires

allpixels.dt <- allpixels.dt %>% 
             filter(NAME %in% toupper(params$fires)) %>%   ## to keep only fires actually specified in params
             select_(.dots = c("NAME", "CLASS", "X", "Y", params$pixel.predictors))  ## keep only predictors of interest

allpixels.dt <- as.data.table(allpixels.dt)

RESULTS <- list()   ## object containing the results

RESULTS$ranked.scales <- list()

#### OTB SEGMENTATION -------------------------------------------------------

## Get rid of redundant scales
MS.params.combs.redund <- ldply(params$multisc, data.frame)
MS.params.combs.redund$.id <- NULL
MS.params.combs <- MS.params.combs.redund[!duplicated(MS.params.combs.redund), ]

## Run OTB MeanShift segmentation for all combinations of parameters beforehand and save segIds for each fire (dt seg.ids.dt)
if (params$meanshift.run) {

  ## df storing segment IDs for all the pixels kept in (as many columns as combinations of parameters)
  seg.ids.df <- data.frame( matrix(nrow=nrow(allpixels.dt), ncol=nrow(MS.params.combs)))  

  ## loop over the Meanshift parameters
  for (idx.params.comb in 1:nrow(MS.params.combs)) {
    
    ranger <- MS.params.combs[idx.params.comb, "ranger"]
    spatialr <- MS.params.combs[idx.params.comb, "spatialr"]
    minsize <- MS.params.combs[idx.params.comb, "minsize"]
        
    ## loop to segment all the fires
    for (fire in params$fires) {
      
      img.file.name <- sprintf("%s_%sbands.tif", fire, params$nr.bands.seg)
      
      if (params$meanshift.save.shp) {  ## if TRUE save segmentation result with a telling name 
        seg.file.name <- sprintf("Segment_%s_ranger%d_spatial%d_minsize%d", fire, ranger, spatialr, minsize)
      } else {   ## otherwise just overwrite it in a file called Segment
        seg.file.name <- "Segment"
      }
      segments <- meanShiftOTB(ranger, spatialr, minsize, OTB.dir, OGR.dir, file.path(data.dir, sprintf("%s_bands", params$nr.bands.seg)), segments.shp.dir, img.file.name, seg.file.name, params$meanshift.save.shp)  ## main Meanshift segmentation
      
      fire.image.proj <- CRS(proj4string(segments))   ## get projection string from segments
      
      ## create spatial points object with coordinates of all the pixels in current fire.in to extract segments
      coords <- SpatialPoints(allpixels.dt[NAME==toupper(fire), .(X, Y)], proj4string=fire.image.proj)  
      id.segment <- over(coords, segments)$dn  ## get segment IDs contained in dn column of the output of function over()
      
      seg.ids.df[allpixels.dt[, NAME]==toupper(fire), idx.params.comb] <- id.segment ## dynamically fill the part of the column number idx.params.comb of the df of segment IDs corresponding to fire.in
    
    }
    
    colnames(seg.ids.df)[idx.params.comb] <- sprintf("rr%ssr%sms%s", ranger, spatialr, minsize)  ## name that same column with a string containing Meanshift parameter values
    
  }  ## end for on nrow(MS.params.combs)
  
  seg.ids.dt <- as.data.table(seg.ids.df)
  rm(seg.ids.df)
  
  save(seg.ids.dt, file=file.path(segments.shp.dir, "segIDs.Rdata"))

} else {
  
  load(file.path(segments.shp.dir, "segIDs.Rdata"))
  
  ## Subset dt based on fires actually included in params$fires only if seg.ids.dt contains segment IDs for all the fires
  if (nrow(seg.ids.dt) == length(pix.fires.OK)) {
    seg.ids.dt <- seg.ids.dt[pix.fires.OK, ] 
  }
  
}

print("Finished segmentation")

#### START LOOCV OVER FIRES ---------------------------------------------------

## initialize empty factor vector with appropriate levels to store final class predictions at each round of the Leave-one-out cross-validation loop
Y.predicted.dt <- do.call("cbind", replicate(length(params$approaches), as.data.table(factor(rep(NA, nrow(allpixels.dt)), levels=levels(allpixels.dt[,CLASS]))), simplify = FALSE))
setnames(Y.predicted.dt, params$approaches)

for (fire.out in params$fires) {  ## LOO-CV loop over the fires to leave out
  
  print(sprintf("Leaving out %s", fire.out))

  fires.in <- params$fires[!params$fires %in% fire.out]   ## fires to keep in at this round of the loop 
  idx.pix.in <- allpixels.dt[,NAME] %in% toupper(fires.in)    ## logical indices of the pixels to keep in
  idx.pix.out <- !allpixels.dt[,NAME] %in% toupper(fires.in)   ## logical indices of the pixels to leave out
  
  allpixels.in.dt <- allpixels.dt[idx.pix.in,]   ## dt with pixels in
  allpixels.out.dt <- allpixels.dt[idx.pix.out,]  ## dt with pixels out

  idx.ms <- 1  ## index for multiscale approaches, different from idx.approach
  
  ## loop over the Meanshift parameters
  
  ## TO UNCOMMENT TO RUN WITH FOREACH ON params$approaches
  # nr.clusters <- length(params$approaches)  ## to uncomment when running in parallel, after successful debugging
  # cl <- makeCluster(nr.clusters)
  # registerDoParallel(cl)
  # Y.predicted.appr <- foreach (idx.approach = 1:length(params$approaches), .combine=cbind, .packages=list.of.packages) %dopar% {
  ## TO UNCOMMENT TO RUN WITH FOREACH ON params$approaches
    
  ## TO UNCOMMENT TO RUN WITH NORMAL FOR ON params$approaches
  for (idx.approach in 1:length(params$approaches)) {
  ## TO UNCOMMENT TO RUN WITH NORMAL FOR ON params$approaches
    
    appr <- params$approaches[idx.approach]
    
    print(sprintf("Approach: %s", appr))
    
    ## Set the default dataset as the allpixels data.tables (in & out), will be augmented with new OBIA features for approaches other than "Pixel" 
    multi.scale.pix.in.dt <- allpixels.in.dt[, c("NAME", "CLASS", params$pixel.predictors), with=FALSE]  ## initialize dt to store object-level features for all the kept-in pixels
    multi.scale.pix.in.dt[, pixIDin:=1:nrow(allpixels.in.dt)]   ## add unique ID to keep order of pixels
    multi.scale.pix.out.dt <- allpixels.out.dt[, c("NAME", "CLASS", params$pixel.predictors), with=FALSE]  ## initialize dt to store object-level features for all the kept-in pixels
    multi.scale.pix.out.dt[, pixIDout:=1:nrow(allpixels.out.dt)]
    predictors <- params$pixel.predictors    ## base pixel-based predictors
    
    if (appr != "Pixel") {
      
      multisc.predictors <- params$pixel.predictors  ## initialize multi scale predictor names as the base pixel-level predictors
      
      cmd <- sprintf("multisc.params.df <- params$multisc$ms%d", idx.ms)
      eval(parse(text=cmd))
      
      ## loop over the Meanshift parameters
      for (idx.params.comb in 1:nrow(multisc.params.df)) {
        
        ranger <- multisc.params.df[idx.params.comb, "ranger"]
        spatialr <- multisc.params.df[idx.params.comb, "spatialr"]
        minsize <- multisc.params.df[idx.params.comb, "minsize"]
            
        scale.name <- sprintf("rr%ssr%sms%s", ranger, spatialr, minsize)  ## name that same column with a string containing Meanshift parameter values
        
        multisc.names <- sprintf("Sc%s_%s", idx.params.comb, params$obj.predictors)
        
        multisc.predictors <- c(multisc.predictors, multisc.names)  ## growing list of multiscale predictor names 
        
#### MULTI-SCALE PRED RETRIEVAL FIRE IN -----------------------------------------------------------
        
        allpixels.in.dt[, segID:= seg.ids.dt[idx.pix.in, scale.name, with=FALSE]]  ## temporarily add segment IDs to allpixel.in.dt
        single.scale.obj.in.dt <- summarize.all(allpixels.in.dt, c("NAME", "segID"), c("X", "Y"))
        
        multi.scale.pix.in.dt[, segID:=allpixels.in.dt[,segID]]
        setkey(multi.scale.pix.in.dt, NAME, segID)
        cols <- c("NAME", "segID", params$obj.predictors)  ## colnames of columns to be joined
        
        multi.scale.pix.in.dt <- multi.scale.pix.in.dt[single.scale.obj.in.dt[, cols, with=F]]   ## merge data.tables by the defined keys (NAME and segID)
        setnames(multi.scale.pix.in.dt, old=params$obj.predictors, new=multisc.names)

        setkey(multi.scale.pix.in.dt, pixIDin) ## re-arrange multi.scale.pix.out.dt to have same order as allpixels.out (setkey() has messed up the order of the rows)
        
        
#### MULTI-SCALE PRED RETRIEVAL FIRE OUT -----------------------------------------------------------
        
        allpixels.out.dt[, segID:= seg.ids.dt[idx.pix.out, scale.name, with=FALSE]]  # add segment IDs to the df with the left-out pixels
        single.scale.obj.out.dt <- summarize.all(allpixels.out.dt, c("NAME", "segID"), c("X", "Y"))  ## overwrite single.scale.obj.out.dt bc it is not needed
        
        multi.scale.pix.out.dt[, segID:=allpixels.out.dt[,segID]]
        setkey(multi.scale.pix.out.dt, NAME, segID)
        cols <- c("NAME", "segID", params$obj.predictors)  ## colnames of columns to be joined
        
        multi.scale.pix.out.dt <- multi.scale.pix.out.dt[single.scale.obj.out.dt[, cols, with=F]]   ## merge data.tables by the defined keys (NAME and segID)
        setnames(multi.scale.pix.out.dt, old=params$obj.predictors, new=multisc.names)
        
        setkey(multi.scale.pix.out.dt, pixIDout) ## re-arrange multi.scale.pix.out.dt to have same order as allpixels.out (setkey() has messed up the order of the rows)
        
        
      } ## end for on nrow(multisc.params.df)
      
      predictors <- multisc.predictors ## multiscale predictors (base predictors + different scales set in params$multisc)
      
      idx.ms <- idx.ms+1 
      
    } ## end if appr != "Pixel"
    
#### MULTI-SCALE & PIX-BASED RF -----------------------------------------------------------

    nr.vars <- length(predictors)
    if (params$mtry == 'sqrt_nr_var') {
      mtries <- floor(sqrt(nr.vars))
    } else if (params$mtry == 'nr_var_div_3') {
      mtries <- floor(nr.vars/3)
    }
    
    ## set sample size for stratified RF training
    target.vec <- multi.scale.pix.in.dt[[params$targ]]
    if (params$use.strata) {
      nmin <- min(table(target.vec))
      ncl <- length(params$mort.class.labels)
      sampsize.vect <- rep(nmin, ncl)
    } else {
      sampsize.vect <- table(target.vec)
    }
    
    ## train RF on in pixels and predict on out pixels, always with a given set of "predictors"
    if (!params$use.strata) {  ## if we do not use strata (thus have a huge amount of class 1 samples), run RF in parallel with foreach()
      set.seed(params$seed)   ## set the same seed every time so that we obtain the same results
      nr.clusters <- min(params$ntree, detectCores())  ## set as minimum between the nr of trees and the nr of cores
      tot.nrtrees <- params$ntree      ## total nr of trees to be shared across cores (clusters)
      nrtrees.clust <- tot.nrtrees%/%(nr.clusters-1)   ## nr of trees per cluster (computed over the total nr of cluster minus 1)
      remaind <- tot.nrtrees%%(nr.clusters-1)   ## remainder of trees for last cluster
      cl <- makeCluster(nr.clusters)    ## initialize cores
      registerDoParallel(cl)   ## regirster them
      ## actual loop that will aggregate the results (output of all trees in the RF object rf.RF, the output of foreach()), use .multicombine=T for increased speed
      RF <- foreach (ntrees=c(rep(nrtrees.clust, nr.clusters-1), remaind), .combine=combine, .multicombine=T, .packages=c('randomForest', 'data.table')) %dopar% {
        randomForest(x=multi.scale.pix.in.dt[, predictors, with=FALSE], y=target.vec,
                     strata=target.vec, sampsize=sampsize.vect,
                     ntree=ntrees, mtry=mtries, nodesize=params$nodesize)
      }
      stopCluster(cl)
    } else {   ## or classical (sequential)
      set.seed(params$seed)
      RF <- randomForest(x=multi.scale.pix.in.dt[, predictors, with=FALSE], y=target.vec,
                         strata=target.vec, sampsize=sampsize.vect,
                         ntree=params$ntree, mtry=mtries, nodesize=params$nodesize, importance=params$plot.importance)
    }
    
    ## Predict on dt with all the pixel based variables (at each iteration with different predictors)
    predict(RF, multi.scale.pix.out.dt[, predictors, with=FALSE], type="response", predict.all=F, nodes=F)
    
    ## TO UNCOMMENT IF RUNNING WITH NORMAL FOR ON params$approaches
    ## Predict on dt with all the pixel based variables (at each iteration with different predictors)
    Y.predicted <- predict(RF, multi.scale.pix.out.dt[, predictors, with=FALSE], type="response", predict.all=F, nodes=F)
    ## Fill correspnding elements of dt containing final predicted values
    Y.predicted.dt[idx.pix.out, appr] <- Y.predicted
    ## TO UNCOMMENT IF RUNNING WITH NORMAL FOR ON params$approaches
    
  } ## end for on params$approaches
  
  ## TO UNCOMMENT IF RUNNING WITH FOREACH ON params$approaches
  # Y.predicted.dt[idx.pix.out, ] <- Y.predicted.appr
  ## TO UNCOMMENT IF RUNNING WITH FOREACH ON params$approaches
  

} ## end for on fire.out

#### ASSESSMENT & MAPS -----------------------------------------------------------

for (idx.approach in 1:length(params$approaches)) {
  
  appr <- params$approaches[idx.approach]
  
  Y.predicted <- Y.predicted.dt[[idx.approach]]  ## get column of data.table by index (subsetting a list)
  
  metrics <- list()
  
  ## overall assessment
  metrics$OVERALL <- classif.metrics(Y.predicted, allpixels.dt$CLASS)
  
  ## by fire assessment
  temp.df <- data.frame(predicted=Y.predicted, observed=allpixels.dt$CLASS, NAME=allpixels.dt$NAME)
  metrics <- c(metrics, dlply(temp.df, .(NAME), classif.metrics) )
  
  cmd <- sprintf("RESULTS$%s <- metrics", appr)  ## final RESULTS object as a list of lists Approach -> Fires -> Metrics
  eval(parse(text=cmd))
  
  if (params$write.maps) {
    
    maps.dir <- file.path(results.dir, 'Maps', fsep = .Platform$file.sep) 
    if (!file.exists(maps.dir)) {
      dir.create(maps.dir, showWarnings=F, recursive=T)  
    }
    
    fire.idx <- 1
    for (fire in params$fires) {
      
      idx.pix <- allpixels.dt[,NAME] == toupper(fire)   
      
      img.file.name <- sprintf("%s_polygon_rasterized.tif", fire)  ## file with which to initialize the empty raster to be filled with predicted values
      img.file.path <- file.path(data.dir, "rasterized_API_polygons")
      newraster <- raster(file.path(img.file.path, img.file.name))
      newraster[1:ncell(newraster)] <- NA     ## reset the values of each cell of newrester to NA (alternatively this could be used: newraster <- raster(newraster))

      xco <- allpixels.dt[NAME==toupper(fire), "X", with=FALSE]  ## get x coordinates for the pixels for which we produced a prediction
      yco <- allpixels.dt[NAME==toupper(fire), "Y", with=FALSE]  ## get y coordinates for the pixels for which we produced a prediction
      
      newraster[cellFromXY(newraster,cbind(xco,yco))] <- Y.predicted[idx.pix]   ## fill in raster at the cell indices specified as the output of cellFromXY()
      
      cmd <- sprintf("map.kappa <- RESULTS$%s$%s$Kappa", appr, toupper(fire))  ## retrieve Kappa for that map
      eval(parse(text=cmd))
      map.kappa <- gsub(".", "p", sprintf("%.3f", map.kappa), fixed=T)
      
      map.name <- sprintf("%s_ClassMap_%s_Kappa%s", fire, appr, map.kappa)
      
      writeRaster(newraster, filename=file.path(maps.dir, map.name, fsep = .Platform$file.sep) , format="GTiff", overwrite=TRUE,  datatype='INT1U')
      
      fire.idx <- fire.idx+1
    }
  }  ## end if write.maps
  
}  ## end for on params$approaches

RES.file = file.path(results.dir, 'RESULTS.Rdata', fsep = .Platform$file.sep) 
save(RESULTS, file = RES.file)

#### PRINT LOGS ---------------------------------------------------------

## clock global time
toc <- proc.time()-tic[3]
end.message <- sprintf("Total elapsed time: %s, finished running on %s", seconds_to_period(toc[3]), Sys.time())

## write log file
log.file <- file.path(results.dir, "log.txt", sep = '')
unlink(log.file, recursive = T, force = T)
write(end.message, file=log.file, append=T)

print(end.message)




