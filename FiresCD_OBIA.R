#### CODE INFOS -------------------------------------------------------------

## Project Name: Fire severity classification 
## Authors: Giona Matasci (giona.matasci@gmail.com), Ignacio San Miguel (ignaciosanmiguel86@gmail.com)   
## File Name:                            
## Objective: 

#### TO DO -------------------------------------------------------------------

# PRIOR TO ACTUAL RUN:
# - 

## STILL TO DO:
# - introduce LOOCV within meanshift loop instead of OOB? (think if need to redo segmentation or just RF in 2nd LOOCV)
# - run as it is on complete set of new fire images (no masks): big fires will impact more the process than the current test on Tanghe and Liege

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


#### INIT --------------------------------------------------------------------

rm(list=ls())

start.message <- sprintf("OBIA for fires CD, started running on %s", Sys.time())
print(start.message)

#### PARAMETERS ---------------------------------------------

params <- list()

## General
# params$approach <- c("Pixel", "ObjectSingleSc", "PixelMultiSc")  ## Type of approach:
                                                      ## Pixel --> classic pixel-based approach 
                                                      ## ObjectSingleSc --> optimizes the parameters of a Meanshift segmentation with an object-based approach  
                                                      ## PixelMultiSc --> features computes for different segmentations are stacked in the same df with a pixel-based approach

params$approach <- c("Pixel", "ObjectSingleSc", "PixelMultiSc")  ## Type of approach:

params$fires <- c("Abraham", "Cone", "Flett", "Keane", "Leggo", "Levellers", "Liege",
                  "Mcarther", "Overflow", "Perry", "Rail", "Rainbow", "Steephill", "Tanghe")  ## 14 fires to consider
# params$fires <- c("Abraham", "Flett", "Mcarther")

params$mort.class.names <- as.factor(c("0-5%", "6-25%", "95-100%"))
params$mort.class.labels <- as.factor(c(1, 2, 3))
params$critical.class.label <- 2  ## NOT USED: corresponding to class "6-25%", used to compute F-measure on this class (as an alternative, more detailed measure to Kappa)

## Meanshift
params$meanshift.run <- T    ## if set to T it runs the MS segmentation and saves segIDs, otherwise it loads them from previous run
params$meanshift.save.shp <- T    ## if set to T saves in temp.dir each segmentation result (shp) with parameter values to visually inspect the segments
params$write.maps <- T   ## wheter or not to write raster containing the pixel level prediction for each approach
params$nr.bands.seg <- 3   ## number of bands to use for segmentation (points to a different folder)
# params$ranger.vect <- c(10, 50, 200) ## On SF imagebest with 50
# params$spatialr.vect <- c(5, 20 , 50)  ## best with 5  
# params$minsize.vect <- c(10, 100, 1000) ## best with 100
# params$ranger.vect <- c(10, 20, 200, 500) ## On SF imagebest with 50
# params$spatialr.vect <- c(3, 10, 20, 50)  ## best with 5
# params$minsize.vect <- c(10) ## best with 100
# params$ranger.vect <- c(50, 100, 200, 300) ## Range radius: on SF imagebest with 50, Liege with 25
params$ranger.vect <- c(25, 50, 100, 150, 200, 250, 300) ## Range radius: on SF imagebest with 50, Liege with 25
params$spatialr.vect <- c(5, 10, 20, 30)  ## Spatial radius: on SF best with 5, Liege with 5
params$minsize.vect <- c(5) ## Minimum object size: on SF best with 100, for fireswith 5, setting it to 0 just produces 1-pixel segments
# params$multisc <- data.frame(ranger=c(25, 100, 200), spatialr=c(5, 10, 10), minsize=c(5, 5, 5))  ## parameters for multiscale: each row is a combination of the 3
params$multisc <- data.frame(ranger=c(100, 200, 250), spatialr=c(10, 10, 5), minsize=c(5, 5, 5))  ## parameters for multiscale: each row is a combination of the 3
params$OOB.samples <- "pixel"   ## wheter to run the MS parameter optimization on the OOB samples at the "object" or "pixel" level

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
params$use.strata <- T    ## wheter to use strata or not in training the RF
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
meanShiftOTB <- function(ranger, spatialr, minsize, OTB.dir, OGR.dir, data.dir, temp.dir, img.file.name, seg.file.name, write.shp) {
  
  img.file <- file.path(data.dir, img.file.name)   ## image to segment
  clean.img.file <- file.path(temp.dir, sprintf("ZeroPadded_%s", img.file.name))
  segment.file.sqlt <- file.path(temp.dir, sprintf("%s.sqlite", seg.file.name))  ## initial segmentation result as a sqlite file (only way to work with Lambert conformal conic projection)
  segment.file <- file.path(temp.dir, sprintf("%s.shp", seg.file.name))  ## segmentation result converted to shp via ogr2ogr
  
  if (file.exists(segment.file)) { ## if the file is already in the folder there are 2 options:
    if (write.shp) {  ## if we want full names to be used ("Segment_%s_ranger%d_spatial%d_minsize%d") it means the segmentation is already done
      segments <- readOGR(dsn=temp.dir, layer=seg.file.name)  ## so just read back the shp with segmentation result...
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
  
  segments <- readOGR(dsn=temp.dir, layer=seg.file.name)  ## read back shp with segmentation result
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
temp.dir <- file.path(data.dir, "temp")  ## directory for temporary files like the segmentation shps (overwritten each time)
if (!file.exists(temp.dir)) {dir.create(temp.dir, showWarnings=F, recursive=T)}  ## create it

load(file.path(data.dir, "allpixels.Rdata"))
allpixels.dt <- dfp
rm(dfp)

pix.fires.OK <- allpixels.dt[["NAME"]] %in% toupper(params$fires) ## get indices of pixels beloning to subset of fires

allpixels.dt <- allpixels.dt %>% 
             filter(NAME %in% toupper(params$fires)) %>%   ## to keep only fires actually specified in params
             select_(.dots = c("NAME", "CLASS", "X", "Y", params$pixel.predictors))  ## keep only predictors of interest

allpixels.dt <- as.data.table(allpixels.dt)

RESULTS <- list()   ## object containing the results

RESULTS$best.scale.names <- c()
RESULTS$ranked.scales <- list()

#### OTB SEGMENTATION -------------------------------------------------------

## Run OTB MeanShift segmentation for all combinations of parameters beforehand and save segIds for each fire (dt seg.ids.dt)
if (params$meanshift.run) {

  ## df storing segment IDs for all the pixels kept in (as many columns as combinations of parameters)
  seg.ids.df <- data.frame( matrix(nrow=nrow(allpixels.dt), ncol=length(params$ranger.vect)*length(params$spatialr.vect)*length(params$minsize.vect)))  
  ms.par.col <- 1  ## start at 1 the index over columns of dataframe of segment IDs
  
  ## 3 nested loops over the Meanshift parameters
  for (ranger in params$ranger.vect) {
    for (spatialr in params$spatialr.vect) {
      for (minsize in params$minsize.vect) {
        
        ## loop to segment all the fires
        for (fire in params$fires) {
          
          img.file.name <- sprintf("%s_%sbands.tif", fire, params$nr.bands.seg)
          
          if (params$meanshift.save.shp) {  ## if TRUE save segmentation result with a telling name 
            seg.file.name <- sprintf("Segment_%s_ranger%d_spatial%d_minsize%d", fire, ranger, spatialr, minsize)
          } else {   ## otherwise just overwrite it in a file called Segment
            seg.file.name <- "Segment"
          }
          segments <- meanShiftOTB(ranger, spatialr, minsize, OTB.dir, OGR.dir, file.path(data.dir, sprintf("%s_bands", params$nr.bands.seg)), temp.dir, img.file.name, seg.file.name, params$meanshift.save.shp)  ## main Meanshift segmentation
          
          fire.image.proj <- CRS(proj4string(segments))   ## get projection string from segments
          
          ## create spatial points object with coordinates of all the pixels in current fire.in to extract segments
          coords <- SpatialPoints(allpixels.dt[NAME==toupper(fire), .(X, Y)], proj4string=fire.image.proj)  
          id.segment <- over(coords, segments)$dn  ## get segment IDs contained in dn column of the output of function over()
          
          seg.ids.df[allpixels.dt[, NAME]==toupper(fire), ms.par.col] <- id.segment ## dynamically fill the part of the column number ms.par.col of the df of segment IDs corresponding to fire.in
        }
        
        colnames(seg.ids.df)[ms.par.col] <- sprintf("rr%ssr%sms%s", ranger, spatialr, minsize)  ## name that same column with a string containing Meanshift parameter values
        
        ms.par.col <- ms.par.col+1  ## increment index over columns of dataframe of IDs
        
      }  ## end for on ranger
    }  ## end for on spatialr
  }  ## end for on minsize
  
  seg.ids.dt <- as.data.table(seg.ids.df)
  rm(seg.ids.df)
  
  save(seg.ids.dt, file=file.path(temp.dir, "segIDs.Rdata"))

} else {
  
  load(file.path(temp.dir, "segIDs.Rdata"))
  seg.ids.dt <- seg.ids.dt[pix.fires.OK, ]
  
}

print("Finished segmentation")

#### START LOOCV OVER FIRES ---------------------------------------------------

## initialize empty factor vector with appropriate levels to store final class predictions at each round of the Leave-one-out cross-validation loop
Y.predicted.pixel.raw <- Y.predicted.object.multi <- Y.predicted.object.single <- factor(rep(NA, nrow(allpixels.dt)), levels=levels(allpixels.dt[,CLASS]))

for (fire.out in params$fires) {  ## LOO-CV loop over the fires to leave out
  
  print(sprintf("Leaving out %s", fire.out))

  fires.in <- params$fires[!params$fires %in% fire.out]   ## fires to keep in at this round of the loop 
  idx.pix.in <- allpixels.dt[,NAME] %in% toupper(fires.in)    ## logical indices of the pixels to keep in
  idx.pix.out <- !allpixels.dt[,NAME] %in% toupper(fires.in)   ## logical indices of the pixels to leave out
  
  allpixels.in.dt <- allpixels.dt[idx.pix.in,]   ## dt with pixels in
  allpixels.out.dt <- allpixels.dt[idx.pix.out,]  ## dt with pixels out

#### START OPTIMIZATION OF MS PARAMS ON KEPT-IN FIRES --------------------------------------
  
  Kappas.OOB.meanshift.dt <- data.table(matrix(nrow=0, ncol=4))   ## initialize empty matrix to store Out-of-Bag Overall Accuracies and associated parameters for each segmentation
  colnames(Kappas.OOB.meanshift.dt) <- c("ranger", "spatialr", "minsize", "Kappa")
  
  ms.par.col <- 1  ## start at 1 the index over columns of dataframe of segment IDs
  
  ## fill object level dt with per segment mean and std dev
  multi.scale.pix.in.dt <- allpixels.in.dt[, c("NAME", "CLASS", params$pixel.predictors), with=FALSE]  ## initialize dt to store object-level features for all the kept-in pixels
  multi.scale.pix.in.dt[, pixIDin:=1:nrow(allpixels.in.dt)]   ## add unique ID to keep order of pixels
  multisc.predictors <- params$pixel.predictors  ## initialize multi scale predictor names as the base pixel-level predictors
  multisc.idx <- 1  ## index iterating over the sets of multiscale parameters

  ## 3 nested loops over the Meanshift parameters
  for (ranger in params$ranger.vect) {
    for (spatialr in params$spatialr.vect) {
      for (minsize in params$minsize.vect) {
        
        scale.name <- sprintf("rr%ssr%sms%s", ranger, spatialr, minsize)  ## name that same column with a string containing Meanshift parameter values
        
        allpixels.in.dt[, segID:= seg.ids.dt[idx.pix.in, scale.name, with=FALSE]]  ## temporarily add segment IDs to allpixel.in.dt
        
        single.scale.obj.in.dt <- summarize.all(allpixels.in.dt, c("NAME", "segID"), c("X", "Y"))
        
        
#### MULTI-SCALE PRED RETRIEVAL FIRE IN -----------------------------------------------------------
        
        ## if segmentation scale is among the set specified by params$multisc, store for every pixel the summarized predictor values at the object level
        if (ranger %in% params$multisc$ranger[multisc.idx] & 
            spatialr %in% params$multisc$spatialr[multisc.idx] & 
            minsize %in% params$multisc$minsize[multisc.idx]) {
          
          multi.scale.pix.in.dt[, segID:=allpixels.in.dt[,segID]]
          setkey(multi.scale.pix.in.dt, NAME, segID)
          cols <- c("NAME", "segID", params$obj.predictors)  ## colnames of columns to be joined
          multi.scale.pix.in.dt <- multi.scale.pix.in.dt[single.scale.obj.in.dt[, cols, with=F]]   ## merge data.tables by the defined keys (NAME and segID)
          multisc.names <- sprintf("Sc%s_%s", multisc.idx, params$obj.predictors)
          setnames(multi.scale.pix.in.dt, old=params$obj.predictors, new=multisc.names)
          multisc.predictors <- c(multisc.predictors, multisc.names)  ## growing list of multiscale predictor names 
          multisc.idx <- multisc.idx + 1  
          
        }

        
#### SINGLE-SCALE RF TO GET OOB KAPPA -----------------------------------------------------------
        
        if (params$OOB.samples == "object") {  ## if we want to run at the object level keep segments
          
          dataset.MS.optimiz <- single.scale.obj.in.dt
          target.vec <- single.scale.obj.in.dt[[params$targ]]  ## select dt column containing the class labels on which the RF will be trained
          
        } else if (params$OOB.samples == "pixel") {    ## otherwise assign the value at the object level to each corresponding pixel
         
          dataset.MS.optimiz <- allpixels.in.dt[, c("NAME", "CLASS"), with=FALSE]  ## initialize dt to store object-level features for all the kept-in pixels
          dataset.MS.optimiz[, pixIDin:=1:nrow(allpixels.in.dt)]   ## add unique ID to keep order of pixels
          dataset.MS.optimiz[, segID:=allpixels.in.dt[,segID]]  ## add segment ID as a column
          setkey(dataset.MS.optimiz, NAME, segID)
          cols <- c("NAME", "segID", params$obj.predictors)  ## colnames of columns to be joined
          dataset.MS.optimiz <- dataset.MS.optimiz[single.scale.obj.in.dt[, cols, with=F]]   ## merge data.tables by the defined keys (NAME and segID)
          setkey(dataset.MS.optimiz, pixIDin)
          
          target.vec <- dataset.MS.optimiz[[params$targ]]
          
        }
        
        ## set mtry parameter according to params$mtry
        nr.vars <- length(params$obj.predictors) 
        if (params$mtry == 'sqrt_nr_var') {
          mtries <- floor(sqrt(nr.vars))
        } else if (params$mtry == 'nr_var_div_3') {
          mtries <- floor(nr.vars/3)
        }
        
        ## set sample size for stratified RF training
        if (params$use.strata) {
          nmin <- min(table(target.vec))
          ncl <- length(params$mort.class.labels)
          sampsize.vect <- rep(nmin, ncl)
        } else {
          sampsize.vect <- table(target.vec)
        }
        set.seed(params$seed)
        
        ## apply RF on df with object-level values using as predictors the columns listed in params$obj.predictors and with response variable the column specified in params$targ
        RF <- randomForest(x=dataset.MS.optimiz[, params$obj.predictors, with=FALSE], y=target.vec,
                                 strata=target.vec, sampsize=sampsize.vect,
                                 ntree=params$ntree, mtry=mtries, nodesize=params$nodesize, importance=params$plot.importance)

        OOB.Kappa <- cohen.kappa(RF$confusion[1:length(params$mort.class.labels), 1:length(params$mort.class.labels)])[[1]]

        Kappas.OOB.meanshift.dt <- rbind( Kappas.OOB.meanshift.dt, data.table(ranger, spatialr, minsize, Kappa=OOB.Kappa) )   ## grow dt with results
        
        ms.par.col <- ms.par.col+1  ## increment index over columns of dataframe of IDs
        
        
      }  ## end for on ranger
    }  ## end for on spatialr
  }  ## end for on minsize
  
  Kappas.OOB.meanshift.sorted.dt <- arrange(Kappas.OOB.meanshift.dt, desc(Kappa))   ## sort df with results by decreasing Kappa
  
  RESULTS$ranked.scales[[length(RESULTS$ranked.scales)+1]] <- Kappas.OOB.meanshift.sorted.dt
  
#### START LOOP TO COMPUTE PRED FOR LEFT-OUT FIRE --------------------------------------
  
  ## segment the left-out fire with the multiscale parameters and, if not already included, with the best parameters
  ranger.vect.LO <- unique(c(params$multisc$ranger, Kappas.OOB.meanshift.sorted.dt$ranger[1]))
  spatialr.vect.LO <- unique(c(params$multisc$spatialr, Kappas.OOB.meanshift.sorted.dt$spatialr[1]))
  minsize.vect.LO <- unique(c(params$multisc$minsize, Kappas.OOB.meanshift.sorted.dt$minsize[1]))
  
  multi.scale.pix.out.dt <- allpixels.out.dt[, c("NAME", "CLASS", params$pixel.predictors), with=FALSE]  ## initialize dt to store object-level features for all the kept-in pixels
  multi.scale.pix.out.dt[, pixIDout:=1:nrow(allpixels.out.dt)]
  multisc.idx <- 1  ## index iterating over the sets of multiscale parameters
  for (ranger in ranger.vect.LO) {
    for (spatialr in spatialr.vect.LO) {
      for (minsize in minsize.vect.LO) {
        
        scale.name <- sprintf("rr%ssr%sms%s", ranger, spatialr, minsize)  ## name that same column with a string containing Meanshift parameter values
        allpixels.out.dt[, segID:= seg.ids.dt[idx.pix.out, scale.name, with=FALSE]]  # add segment IDs to the df with the left-out pixels
      
#### SINGLE-SCALE RF -----------------------------------------------------------
        
        ## if best parameters build singlescale RF model on segments.in and apply on segments.out
        if (ranger==Kappas.OOB.meanshift.sorted.dt$ranger[1] &
            spatialr==Kappas.OOB.meanshift.sorted.dt$spatialr[1] &
            minsize==Kappas.OOB.meanshift.sorted.dt$minsize[1]) {
        
          RESULTS$best.scale.names <- c(RESULTS$best.scale.names, scale.name)
          
          ## retrieve segment IDs for kept-in fires associated with the best combination of parameters
          allpixels.in.dt[, segID:= seg.ids.dt[idx.pix.in, scale.name, with=FALSE]]  ## add IDs to the kept-in pixels dt
          
          single.scale.obj.in.dt <- summarize.all(allpixels.in.dt, c("NAME", "segID"), c("X", "Y"))  ## overwrite single.scale.obj.in.dt bc it is not needed
          single.scale.obj.out.dt <- summarize.all(allpixels.out.dt, c("NAME", "segID"), c("X", "Y"))  ## overwrite single.scale.obj.out.dt bc it is not needed

          nr.vars <- length(params$obj.predictors)
          if (params$mtry == 'sqrt_nr_var') {
            mtries <- floor(sqrt(nr.vars))
          } else if (params$mtry == 'nr_var_div_3') {
            mtries <- floor(nr.vars/3)
          }
         
          ## set sample size for stratified RF training
          target.vec <- single.scale.obj.in.dt[[params$targ]]
          if (params$use.strata) {
            nmin <- min(table(target.vec))
            ncl <- length(params$mort.class.labels)
            sampsize.vect <- rep(nmin, ncl)
          } else {
            sampsize.vect <- table(target.vec)
          }
          ## train RF on in segments and predict on out segments, always with params$obj.predictors
          set.seed(params$seed)
          RF <- randomForest(x=single.scale.obj.in.dt[, params$obj.predictors, with=FALSE], y=target.vec,
                             strata=target.vec, sampsize=sampsize.vect,
                             ntree=params$ntree, mtry=mtries, nodesize=params$nodesize, importance=params$plot.importance)
          
          Y.predicted.segments.out <- predict(RF, single.scale.obj.out.dt[, params$obj.predictors, with=FALSE], type="response", predict.all=F, nodes=F)
          
          ## build a dt with the predicted class for each segment (Y.predicted.segments.out) and the associated segment ID (single.scale.obj.out.df[, "segID"])
          y.pred.segID.dt <- data.table(segID=single.scale.obj.out.dt[, segID], ypred=Y.predicted.segments.out)
          
          ## Join the predicted labels for the segments to the corresponding segID in the complete vector of all images
          setkey(y.pred.segID.dt, segID) ## set key as the segment IDs column
          segID.allpixels.out.dt <- data.table(pixIDout=1:nrow(allpixels.out.dt), segID=allpixels.out.dt[,segID])  ## add pixel IDs so that we can retrieve their original order
          setkey(segID.allpixels.out.dt, "segID")
          merged.pix.predictions.dt <- segID.allpixels.out.dt[y.pred.segID.dt, ]  ## join to data.table based on a common key with this command allpixels.out.segID.dt[y.pred.segID.dt], then select only ypred as a column
          setkey(merged.pix.predictions.dt, "pixIDout")  ## resort back to the original order 
          
          Y.predicted.object.single[idx.pix.out] <- merged.pix.predictions.dt$ypred
          
        } ## end if best parameters

#### MULTI-SCALE PRED RETRIEVAL FIRE OUT -----------------------------------------------------------
        
        ## if segmentation scale is among the set specified by params$multisc, store for every pixel the summarized predictor values at the object level
        if (ranger %in% params$multisc$ranger[multisc.idx] & 
            spatialr %in% params$multisc$spatialr[multisc.idx] & 
            minsize %in% params$multisc$minsize[multisc.idx]) {
          
          single.scale.obj.out.dt <- summarize.all(allpixels.out.dt, c("NAME", "segID"), c("X", "Y"))  ## overwrite single.scale.obj.out.dt bc it is not needed
          
          multi.scale.pix.out.dt[, segID:=allpixels.out.dt[,segID]]
          setkey(multi.scale.pix.out.dt, NAME, segID)
          cols <- c("NAME", "segID", params$obj.predictors)  ## colnames of columns to be joined
          multi.scale.pix.out.dt <- multi.scale.pix.out.dt[single.scale.obj.out.dt[, cols, with=F]]   ## merge data.tables by the defined keys (NAME and segID)
          setnames(multi.scale.pix.out.dt, old=params$obj.predictors, new=sprintf("Sc%s_%s", multisc.idx, params$obj.predictors))
          multisc.idx <- multisc.idx + 1  
          
        } ## end if multisc parameters
        
      }
    }
  } ## end for on MeanShift parameters for Left-out fire
  
  ## re-arrange multi.scale.pix.out.dt to have same order as allpixels.out (setkey() has messed up the order of the rows)
  setkey(multi.scale.pix.out.dt, pixIDout)
  setkey(multi.scale.pix.in.dt, pixIDin)
  
  
#### MULTI-SCALE & PIX-BASED RF -----------------------------------------------------------
  
  for (appr in params$approach[params$approach %in% c("Pixel", "PixelMultiSc")]) {
    
    ## depending on the type of approach set on which predictors the RF will run
    if (appr == "Pixel") {
      predictors <- params$pixel.predictors    ## base pixel-based predictors
    } else if (appr == "PixelMultiSc") {
      predictors <- multisc.predictors    ## multiscale predictors (base predictors + different scales set in params$multisc)
    }

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
    
    ## predict on dt with all the pixel based variables (at each iteration with different predictors)
    Y.predicted <- predict(RF, multi.scale.pix.out.dt[, predictors, with=FALSE], type="response", predict.all=F, nodes=F)
    
    if (appr == "Pixel") {
      Y.predicted.pixel.raw[idx.pix.out] <- Y.predicted  
    } else if (appr == "PixelMultiSc") {
      Y.predicted.object.multi[idx.pix.out] <- Y.predicted  
    }

  }  ## end for on params$approach

} ## end for on fire.out

#### ASSESSMENT & MAPS -----------------------------------------------------------

for (appr in params$approach) {
  
  if (appr == "Pixel") {
    Y.predicted <- Y.predicted.pixel.raw  
  } else if (appr == "ObjectSingleSc") {
    Y.predicted <- Y.predicted.object.single  
  } else if (appr == "PixelMultiSc") {
    Y.predicted <- Y.predicted.object.multi  
  }

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
      
      if (appr == "ObjectSingleSc") {
        map.name <- sprintf("%s_ClassMap_%s_%s_Kappa%s", fire, appr, RESULTS$best.scale.names[fire.idx], map.kappa)
      } else {
        map.name <- sprintf("%s_ClassMap_%s_Kappa%s", fire, appr, map.kappa)
      }
      
      writeRaster(newraster, filename=file.path(maps.dir, map.name, fsep = .Platform$file.sep) , format="GTiff", overwrite=TRUE,  datatype='INT1U')
      
      fire.idx <- fire.idx+1
    }
  }  ## end if write.maps
  
}  ## end for on params$approach

RES.file = file.path(results.dir, 'RESULTS.Rdata', fsep = .Platform$file.sep) 
save(RESULTS, file = RES.file)

#### PRINT LOGS ---------------------------------------------------------

## clock global time
toc <- proc.time()-tic[3]
end.message <- sprintf("Total elapsed time: %s, finished running on %s", seconds_to_period(toc[3]), Sys.time())
print(end.message)




