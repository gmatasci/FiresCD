#### CODE INFOS --------------------------------------------------------------

## Project Name: Fire severity classification 
## Authors: Giona Matasci (giona.matasci@gmail.com), Ignacio San Miguel (ignaciosanmiguel86@gmail.com)   
## File Name:                            
## Objective: 

#### TO DO -------------------------------------------------------------------

# PRIOR TO ACTUAL RUN:
# - 

## STILL TO DO:
# - segment all the fires beforehand and store segID in a table to query in the script
# - include stratified RF
# - introduce LOOCV within meanshift loop instead of OOB? (think if need to redo segmentation or just RF in 2nd LOOCV)
# - uniformize fire names (caps, etc.)
# - run as it is on complete set of new fire images (no masks): big fires will impact more the process than the current test on Tanghe and Liege

## SOLVED:
# -V issue with images/projection -- OTB does not accept Lambert Conical when producing a shp, we have to save a sqlite file and then convert it to shp with ogr2ogr
# -V fill ypred df based in indices and not on rbind -- done
# -V read about OTB meanshift minsize -- the 2 other parameters control the size and shape of the segments: minsize is just used as a post-processing to merge segments smaller than a threshold
# -V split aspect in 2 components sin/cos -- new variables: "aspect.sin_mean", "aspect.sin_sd", "aspect.cos_mean", "aspect.cos_sd"
# -V check error: Error in `colnames<-`(`*tmp*`, value = c(NA, NA, NA, "elev_sd")) : length of 'dimnames' [2] not equal to array extent -- not appearing anymore
# -V new pixel level GT with majority class by surface within pixel -- done by Nacho
# -V integrate multiscale approach using seg.ids.df as a starting point -- now we 
# -V base parameter selection and assessment on a mixed metric: 0.5*Kappa + 0.5*Fmeasure_PartialMortality


#### INIT --------------------------------------------------------------------

rm(list=ls())

start.message <- sprintf("OBIA for fires CD, started running on %s", Sys.time())
print(start.message)

#### PARAMETERS ---------------------------------------------

params <- list()

## General
params$approach <- c("PixelBased", "ObjectSingleScale", "ObjectMultiScale")  ## Type of approach (so far only ObjectSingleScale is implemented):
                                                      ## PixelBased --> classic pixel-based approach 
                                                      ## ObjectSingleScale --> optimizes the parameters of a Meanshift segmentation with an object-based approach  
                                                      ## ObjectMultiScale --> features computes for different segmentations are stacked in the same df with a pixel-based approach

params$fires <- c("Liege", "Tanghe", "Overflow")
# params$fires <- c("Abraham", "Cone", "Flett", "Keane", "Leggo", "Levellers", "Liege",
#                      "Mcarther", "Overflow", "Perry", "Rail", "Rainbow", "Steephill", "Tanghe")  ## fires to consider
# params$fires <- c("Liege", "Rail", "Tanghe") 

# params$subsetting <- T    ## to subset data for development purposes (still to implement)
params$mort.class.names <- as.factor(c("0-5%", "6-25%", "95-100%"))
params$mort.class.labels <- as.factor(c(1, 2, 3))
params$critical.class.label <- 2  ## corresponding to class "6-25%", used to compute F-measure on this class (as an alternative, more detailed measure to Kappa)

## Meanshift
params$nr.bands.seg <- 3
params$meanshift.testmode <- F    ## if set to T saves each segmentation result (shp) with parameter values to visually inspect the segments
# params$ranger.vect <- c(10, 50, 200) ## On SF imagebest with 50
# params$spatialr.vect <- c(5, 20 , 50)  ## best with 5  
# params$minsize.vect <- c(10, 100, 1000) ## best with 100
# params$ranger.vect <- c(10, 20, 200, 500) ## On SF imagebest with 50
# params$spatialr.vect <- c(3, 10, 20, 50)  ## best with 5
# params$minsize.vect <- c(10) ## best with 100
# params$ranger.vect <- c(50, 100, 200, 300) ## Range radius: on SF imagebest with 50, Liege with 25
params$ranger.vect <- c(100, 200) ## Range radius: on SF imagebest with 50, Liege with 25
params$spatialr.vect <- c(10)  ## Spatial radius: on SF best with 5, Liege with 5
params$minsize.vect <- c(5) ## Minimum object size: on SF best with 100, for fireswith 5, setting it to 0 just produces 1-pixel segments
params$multisc <- data.frame(ranger=c(100, 200), spatialr=c(10, 10), minsize=c(5, 5))

## RF
params$base.predictors <- c("db1", "db2", "db3", "db4", "db5", "db7", 
                               "dNBR", "dTCG", "dTCB", "dTCW", "dNDVI", "dNDWI", 
                               "elev", "slope", "aspect.sin", "aspect.cos", 
                               "EOSD")   ## list of starting predictors (for the classification part, as the segmentation is run on 6-band difference image only)
params$obj.predictors <- c("db1_mean", "db1_sd", "db2_mean", "db2_sd", "db3_mean", "db3_sd", "db4_mean", "db4_sd", "db5_mean", "db5_sd", "db7_mean", "db7_sd", 
                              "dNBR_mean", "dNBR_sd", "dTCG_mean", "dTCG_sd", "dTCB_mean", "dTCB_sd", "dTCW_mean", "dTCW_sd", "dNDVI_mean", "dNDVI_sd", "dNDWI_mean", "dNDWI_sd",
                              "elev_mean", "elev_sd", "slope_mean", "slope_sd", "aspect.sin_mean", "aspect.sin_sd", "aspect.cos_mean", "aspect.cos_sd",
                              "EOSD_maj", "EOSD_majpct",
                              "nrpixseg")    ## list of final predictors computed at the object level

params$targ <- "CLASS"     ## target variable of the classification (column name of dataframe)
params$seed <- 2016        ## seed to have same RF result
params$parallel.RF <- T    ## whether to run RF in parallel or not
params$ntree <- 100     ## RF nr of trees
params$mtry <- 'sqrt_nr_var'  ## how to set RF mtry: 'sqrt_nr_var' or 'nr_var_div_3'
params$nodesize <- 1   ## RF nodesize: default for classification is 1
params$plot.importance <- F  ## whether to plot RF variable importance

OTB.dir <- "C:/Users/gmatasci/Downloads/OTB-5.4.0-win64/bin"   ## directory in which OTB is located (input to function meanShiftOTB() )
OGR.dir <- "C:/OSGeo4W64/bin"    ## directory in which OGR is located (input to function meanShiftOTB() ), uses ogr2ogr command for shp conversion from sqlite
base.dir <- 'D:/Research/ANALYSES/FiresCD'    ## base working directory

#### LOAD PACKAGES ----------------------------------------------------------

list.of.packages <- c("caret",
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
                      "lazyeval",
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
meanShiftOTB <- function(ranger, spatialr, minsize, OTB.dir, OGR.dir, data.dir, temp.dir, img.file.name, seg.file.name, test.mode) {
  
  img.file <- file.path(data.dir, img.file.name)   ## image to segment
  clean.img.file <- file.path(temp.dir, sprintf("ZeroPadded_%s", img.file.name))
  segment.file.sqlt <- file.path(temp.dir, sprintf("%s.sqlite", seg.file.name))  ## initial segmentation result as a sqlite file (only way to work with Lambert conformal conic projection)
  segment.file <- file.path(temp.dir, sprintf("%s.shp", seg.file.name))  ## segmentation result converted to shp via ogr2ogr
  
  ## delete existing segmentation files
  if (file.exists(segment.file)){ 
    unlink(segment.file, recursive = T, force = T) 
    unlink(segment.file.sqlt, recursive = T, force = T)
  }
  
  ## change no-data values to 0 to avoid the 1-segment-per-pixel situation outside the image (and get only a big polygon outside)
  command.text <- sprintf("%s -in %s -out %s -mode changevalue -mode.changevalue.newv 0",  ## create text for command to run OTB (equivalent to saving this text in a batch file .bat and running it with system2() )
                          file.path(OTB.dir, "otbcli_ManageNoData"),   ## complete path to OTB bat file doing the value change
                          img.file,
                          clean.img.file
  )
  system(command.text)  ## run command
  
  
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

## returns F-measure of the critical class and Kappa statistic
fires.classif.metrics <- function(predicted, observed, critical.class) {
  RES <- confusionMatrix(predicted, observed)
  sens <- RES$byClass[critical.class, "Sensitivity"]
  spec <- RES$byClass[critical.class, "Specificity"]
  return(data.frame(FmeasCritClass=(2*sens*spec)/(sens+spec), Kappa=as.vector(RES$overall[2])))
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

allpixels.dt <- allpixels.dt %>% 
             filter(NAME %in% toupper(params$fires)) %>%   ## to keep only fires actually specified in params 
             select_(.dots = c("NAME", "CLASS", "X", "Y", params$base.predictors))  ## keep only predictors of interest

allpixels.dt <- as.data.table(allpixels.dt)

#### LOOCV OVER FIRES ---------------------------------------------------

## initialize empty factor vector with appropriate levels to store final class predictions at each round of the Leave-one-out cross-validation loop
Y.predicted.pixel.raw <- Y.predicted.object.multi <- Y.predicted.object.single <- factor(rep(NA, nrow(allpixels.dt)), levels=levels(allpixels.dt[,CLASS]))

for (fire.out in params$fires) {  ## LOO-CV loop over the fires to leave out

  fires.in <- params$fires[!params$fires %in% fire.out]   ## fires to keep in at this round of the loop 
  idx.pix.in <- allpixels.dt[,NAME] %in% toupper(fires.in)    ## logical indices of the pixels to keep in
  idx.pix.out <- !allpixels.dt[,NAME] %in% toupper(fires.in)   ## logical indices of the pixels to leave out
  
  allpixels.in.dt <- allpixels.dt[idx.pix.in,]   ## dt with pixels in
  allpixels.out.dt <- allpixels.dt[idx.pix.out,]  ## dt with pixels out

#### SEGMENTATION OF KEPT-IN FIRES --------------------------------------
  
  OAs.OOB.meanshift.dt <- data.table(matrix(nrow=0, ncol=4))   ## initialize empty matrix to store Out-of-Bag Overall Accuracies and associated parameters for each segmentation
  colnames(OAs.OOB.meanshift.dt) <- c("ranger", "spatialr", "minsize", "OA")
  
  ## df storing segment IDs for all the pixels kept in (as many columns as combinations of parameters)
  seg.ids.df <- data.frame(matrix(nrow=nrow(allpixels.in.dt), ncol=length(params$ranger.vect)*length(params$spatialr.vect)*length(params$minsize.vect)))  
  ms.par.col <- 1  ## start at 1 the index over columns of dataframe of segment IDs
  
  ## fill object level dt with per segment mean and std dev
  multi.scale.pix.in.dt <- allpixels.in.dt[, c("NAME", "CLASS", params$base.predictors), with=FALSE]  ## initialize dt to store object-level features for all the kept-in pixels
  multisc.predictors <- params$base.predictors
  multisc.idx <- 1  ## index iterating over the sets of multiscale parameters

  ## 3 nested loops over the Meanshift parameters
  for (ranger in params$ranger.vect) {
    for (spatialr in params$spatialr.vect) {
      for (minsize in params$minsize.vect) {
        
        ## loop to segment all the fires kept in
        for (fire.in in fires.in) {

          img.file.name <- sprintf("%s_%sbands.tif", fire.in, params$nr.bands.seg)

          if (params$meanshift.testmode) {  ## if TRUE save segmentation result with a telling name 
            seg.file.name <- sprintf("Segment_%s_ranger%d_spatial%d_minsize%d", fire.in, ranger, spatialr, minsize)
          } else {   ## otherwise just overwrite it in a file called Segment
            seg.file.name <- "Segment"
          }
          segments <- meanShiftOTB(ranger, spatialr, minsize, OTB.dir, OGR.dir, file.path(data.dir, sprintf("%s_bands", params$nr.bands.seg)), temp.dir, img.file.name, seg.file.name, params$meanshift.testmode)  ## main Meanshift segmentation

          fire.image.proj <- CRS(proj4string(segments))   ## get projection string from segments
          
          ## create spatial points object with coordinates of all the pixels in current fire.in to extract segments
          coords <- SpatialPoints(allpixels.in.dt[NAME==toupper(fire.in), .(X, Y)], proj4string=fire.image.proj)  
          id.segment <- over(coords, segments)$dn  ## get segment IDs contained in dn column of the output of function over()
          
          seg.ids.df[allpixels.in.dt[,NAME]==toupper(fire.in), ms.par.col] <- id.segment  ## dynamically fill the part of the column number ms.par.col of the df of segment IDs corresponding to fire.in
          
        }
        
        colnames(seg.ids.df)[ms.par.col] <- sprintf("rr%ssr%sms%s", ranger, spatialr, minsize)  ## name that same column with a string containing Meanshift parameter values
        
        allpixels.in.dt[, segID:= seg.ids.df[, ms.par.col]]  ## temporarily add segment IDs to allpixel.in df
        
        single.scale.obj.dt <- summarize.all(allpixels.in.dt, c("NAME", "segID"), c("X", "Y"))
        
        ## if segmentation scale is among the set specified by params$multisc, store for every pixel the summarized predictor values at the object level
        if (ranger %in% params$multisc$ranger[multisc.idx] & 
            spatialr %in% params$multisc$spatialr[multisc.idx] & 
            minsize %in% params$multisc$minsize[multisc.idx]) {
          
          multi.scale.pix.in.dt[, segID:=allpixels.in.dt[,segID]]
          setkey(multi.scale.pix.in.dt, NAME, segID)
          cols <- c("NAME", "segID", params$obj.predictors)  ## colnames of columns to be joined
          multi.scale.pix.in.dt <- multi.scale.pix.in.dt[single.scale.obj.dt[, cols, with=F]]   ## merge data.tables by the defined keys (NAME and segID)
          multisc.names <- sprintf("Sc%s_%s", multisc.idx, params$obj.predictors)
          setnames(multi.scale.pix.in.dt, old=params$obj.predictors, new=multisc.names)
          multisc.predictors <- c(multisc.predictors, multisc.names)  ## save multiscale predictors names in a growing list
          multisc.idx <- multisc.idx + 1  
          
        }
        
        ## set mtry parameter according to params$mtry
        nr.vars <- length(params$obj.predictors) 
        if (params$mtry == 'sqrt_nr_var') {
          mtries <- floor(sqrt(nr.vars))
        } else if (params$mtry == 'nr_var_div_3') {
          mtries <- floor(nr.vars/3)
        }
        set.seed(params$seed)   
        ## apply RF on df with object-level values using as predictors the columns listed in params$obj.predictors and with response variable the column specified in params$targ
        RF <- randomForest(x=single.scale.obj.dt[,params$obj.predictors, with=FALSE], y=single.scale.obj.dt[[params$targ]], 
                           ntree=params$ntree, mtry=mtries, nodesize=params$nodesize, importance=params$plot.importance)
        OAs.OOB.meanshift.dt <- rbind( OAs.OOB.meanshift.dt, data.table(ranger, spatialr, minsize, OA=1-RF$err.rate[params$ntree, "OOB"]) )   ## grow dt with results
      
        ms.par.col <- ms.par.col+1  ## increment index over columns of dataframe of IDs
        
      }  ## end for on ranger
    }  ## end for on spatialr
  }  ## end for on minsize
  
  OAs.OOB.meanshift.sorted.dt <- arrange(OAs.OOB.meanshift.dt, desc(OA))   ## sort df with results by decreasing OA
  
#### SEGMENTATION OF LEFT-OUT FIRE ---------------------------------------------
  
  ## segment the left-out fire with the multiscale parameters and, if not already included, with the best parameters
  ranger.vect.LO <- unique(c(params$multisc$ranger, OAs.OOB.meanshift.sorted.dt$ranger[1]))
  spatialr.vect.LO <- unique(c(params$multisc$spatialr, OAs.OOB.meanshift.sorted.dt$spatialr[1]))
  minsize.vect.LO <- unique(c(params$multisc$minsize, OAs.OOB.meanshift.sorted.dt$minsize[1]))
  
  multi.scale.pix.out.dt <- allpixels.out.dt[, c("NAME", "CLASS", params$base.predictors), with=FALSE]  ## initialize dt to store object-level features for all the kept-in pixels
  multisc.idx <- 1  ## index iterating over the sets of multiscale parameters
  for (ranger in ranger.vect.LO) {
    for (spatialr in spatialr.vect.LO) {
      for (minsize in minsize.vect.LO) {
        
        ## segment left out fire
        img.file.name <- sprintf("%s_%sbands.tif", fire.out, params$nr.bands.seg)
        if (params$meanshift.testmode) {
          seg.file.name <- sprintf("Segment_%s_ranger%d_spatial%d_minsize%d", fire.out, ranger, spatialr, minsize)
        } else {
          seg.file.name <- "Segment"
        }
        segments <- meanShiftOTB(ranger, spatialr, minsize, OTB.dir, OGR.dir, file.path(data.dir, sprintf("%s_bands", params$nr.bands.seg)), temp.dir, img.file.name, seg.file.name, params$meanshift.testmode)
        fire.image.proj <- CRS(proj4string(segments)) 
        coords <- SpatialPoints(allpixels.out.dt[, .(X, Y)], proj4string=fire.image.proj)  ## create spatial points for extrating segments from OTB
        id.segment.out <- over(coords, segments)$dn

        allpixels.out.dt[, segID:= id.segment.out] # add segment IDs to the df with the left-out pixels
        
#### CLASSIFICATION OF LEFT-OUT FIRE ----------------------------------------
        
        ## if best parameters
        if (ranger==OAs.OOB.meanshift.sorted.dt$ranger[1] &
            spatialr==OAs.OOB.meanshift.sorted.dt$spatialr[1] &
            minsize==OAs.OOB.meanshift.sorted.dt$minsize[1]) {
        
          ## retrieve segment IDs for kept-in fires associated with the best combination of parameters
          best.colname <- sprintf("rr%ssr%sms%s", OAs.OOB.meanshift.sorted.dt$ranger[1], OAs.OOB.meanshift.sorted.dt$spatialr[1], OAs.OOB.meanshift.sorted.dt$minsize[1])
          id.segment.in <- seg.ids.df[, best.colname]   
          allpixels.in.dt[, segID:= id.segment.in]  ## add IDs to the kept-in pixels dt
          
          allpixels.final.dt <- rbind(allpixels.in.dt, allpixels.out.dt)  ## stack together the in and out pixels (unique segments are identified by fire NAME and segID) 
          
          single.scale.obj.dt <- summarize.all(allpixels.final.dt, c("NAME", "segID"), c("X", "Y"))  ## overwrite single.scale.obj.dt bc it is not needed
          
          nr.vars <- length(params$obj.predictors)
          if (params$mtry == 'sqrt_nr_var') {
            mtries <- floor(sqrt(nr.vars))
          } else if (params$mtry == 'nr_var_div_3') {
            mtries <- floor(nr.vars/3)
          }
          segments.in <- single.scale.obj.dt$NAME %in% toupper(fires.in)   ## retrieve indices of kept-in segments in the object-level df
          segments.out <- single.scale.obj.dt$NAME == toupper(fire.out)    ## retrieve indices of left-out segments in the object-level df
          
          ## train RF on in segments and predict on out segments, always with params$obj.predictors
          set.seed(params$seed)
          RF <- randomForest(x=single.scale.obj.dt[segments.in, params$obj.predictors, with=FALSE], y=single.scale.obj.dt[[params$targ]][segments.in],   ## y has to be a vector and the syntax for data.table is first getting the vector with [[]] then subsetting it from outside by adding [segments.in] 
                             ntree=params$ntree, mtry=mtries, nodesize=params$nodesize, importance=params$plot.importance)
          
          Y.predicted.segments.out <- predict(RF, single.scale.obj.dt[segments.out, params$obj.predictors, with=FALSE], type="response", predict.all=F, nodes=F)
          
          ## Save shp with predicted class
          ## build a dt with the predicted class for each segment (Y.predicted.segments.out) and the associated segment ID (single.scale.obj.df[segments.out, "segID"])
          y.pred.segID.dt <- data.table(segID=single.scale.obj.dt[segments.out, segID], ypred=Y.predicted.segments.out)
          
          ## assign the predicted class to the "data" df of the segments shp by merging by segment ID ("dn" or "segID"), all.x=T is used to keep the segments for which there is no prediction (outside of fire)
          segments@data <- merge(segments@data, y.pred.segID.dt, by.x="dn", by.y="segID", all.x=T)  
          segments@data$ypred[is.na(segments@data$ypred)] <- params$mort.class.labels[1]    ## if NA are assigned to polygons outside of fire, assign the lowest mortality class instead
          writeOGR(segments, results.dir, sprintf("%s_pred_map_%s", fire.out, best.colname), driver="ESRI Shapefile", overwrite_layer=TRUE)   ## write prediction map shapefile for the left-out fire
          
          ## Join the predicted labels for the segments to the corresponding segID in the complete vector of all images
          setkey(y.pred.segID.dt, segID) ## set key as the segment IDs column
          segID.allpixels.out.dt <- data.table(segID=allpixels.out.dt[,segID], key = "segID") #
          Y.predicted.object.single[idx.pix.out] <- segID.allpixels.out.dt[y.pred.segID.dt, ypred]  ## join to data.table based on a common key with this command allpixels.out.segID.dt[y.pred.segID.dt], then select only ypred as a column
          
        } ## end if best parameters
        
        ## if segmentation scale is among the set specified by params$multisc, store for every pixel the summarized predictor values at the object level
        if (ranger %in% params$multisc$ranger[multisc.idx] & 
            spatialr %in% params$multisc$spatialr[multisc.idx] & 
            minsize %in% params$multisc$minsize[multisc.idx]) {
          
          single.scale.obj.dt <- summarize.all(allpixels.out.dt, c("NAME", "segID"), c("X", "Y"))  ## overwrite single.scale.obj.dt bc it is not needed
          
          multi.scale.pix.out.dt[, segID:=allpixels.out.dt[,segID]]
          setkey(multi.scale.pix.out.dt, NAME, segID)
          cols <- c("NAME", "segID", params$obj.predictors)  ## colnames of columns to be joined
          multi.scale.pix.out.dt <- multi.scale.pix.out.dt[single.scale.obj.dt[, cols, with=F]]   ## merge data.tables by the defined keys (NAME and segID)
          setnames(multi.scale.pix.out.dt, old=params$obj.predictors, new=sprintf("Sc%s_%s", multisc.idx, params$obj.predictors))
          multisc.idx <- multisc.idx + 1  
          
        } ## end if multisc parameters
        
      }
    }
  } ## end for on MeanShift parameters
  
  multi.scale.pix.dt <- rbind(multi.scale.pix.in.dt, multi.scale.pix.out.dt)  ## stack together the in and out pixels (unique segments are identified by fire NAME and segID) 
  
  for (appr in c("PixelBased", "ObjectMultiScale")) {
    
    ## depending on the type of approach set on which predictors the RF will run
    if (appr == "PixelBased") {
      predictors <- params$base.predictors    ## base pixel-based predictors
    } else if (appr == "ObjectMultiScale") {
      predictors <- multisc.predictors    ## multiscale predictors (base predictors + different scales set in params$multisc)
    }

    nr.vars <- length(predictors)
    if (params$mtry == 'sqrt_nr_var') {
      mtries <- floor(sqrt(nr.vars))
    } else if (params$mtry == 'nr_var_div_3') {
      mtries <- floor(nr.vars/3)
    }
    
    ## train RF on in pixels and predict on out pixels, always with a given set of "predictors"
    if (params$parallel.RF) {  ## in parallel with foreach()
      set.seed(params$seed)   ## set the same seed every time so that we obtain the same results
      nr.clusters <- min(params$ntree, detectCores())  ## set as minimum between the nr of trees and the nr of cores
      tot.nrtrees <- params$ntree      ## total nr of trees to be shared across cores (clusters)
      nrtrees.clust <- tot.nrtrees%/%(nr.clusters-1)   ## nr of trees per cluster (computed over the total nr of cluster minus 1)
      remaind <- tot.nrtrees%%(nr.clusters-1)   ## remainder of trees for last cluster
      cl <- makeCluster(nr.clusters)    ## initialize cores
      registerDoParallel(cl)   ## regirster them
      ## actual loop that will aggregate the results (output of all trees in the RF object rf.RF, the output of foreach()), use .multicombine=T for increased speed
      RF <- foreach (ntrees=c(rep(nrtrees.clust, nr.clusters-1), remaind), .combine=combine, .multicombine=T, .packages=c('randomForest', 'data.table')) %dopar% {
        randomForest(x=multi.scale.pix.dt[idx.pix.in, predictors, with=FALSE], y=multi.scale.pix.dt[[params$targ]][idx.pix.in], 
                     ntree=ntrees, mtry=mtries, nodesize=params$nodesize)
      }
      stopCluster(cl)
    } else {   ## or classical (sequential)
      set.seed(params$seed)
      RF <- randomForest(x=multi.scale.pix.dt[idx.pix.in, predictors, with=FALSE], y=multi.scale.pix.dt[[params$targ]][idx.pix.in],   ## y has to be a vector and the syntax for data.table is first getting the vector with [[]] then subsetting it from outside by adding [segments.in] 
                         ntree=params$ntree, mtry=mtries, nodesize=params$nodesize, importance=params$plot.importance)
    }
    
    ## predict on dt with all the pixel based variables (at each iteration with different predictors)
    Y.predicted <- predict(RF, multi.scale.pix.dt[idx.pix.out, predictors, with=FALSE], type="response", predict.all=F, nodes=F)
    
    if (appr == "PixelBased") {
      Y.predicted.pixel.raw[idx.pix.out] <- Y.predicted  
    } else if (appr == "ObjectMultiScale") {
      Y.predicted.object.multi[idx.pix.out] <- Y.predicted  
    }

  }  ## end for on params$approach

} ## end for on fire.out

#### ASSESSMENT -----------------------------------------------------------

RES <- list()
for (appr in params$approach) {
  
  if (appr == "PixelBased") {
    Y.predicted <- Y.predicted.pixel.raw  
  } else if (appr == "ObjectSingleScale") {
    Y.predicted <- Y.predicted.object.single  
  } else if (appr == "ObjectMultiScale") {
    Y.predicted <- Y.predicted.object.multi  
  }

  ## overall assessment
  metrics.overall <- fires.classif.metrics(Y.predicted, allpixels.dt$CLASS, params$critical.class.label)
  
  ## by fire assessment
  temp.df <- data.frame(predicted=Y.predicted, observed=allpixels.dt$CLASS, NAME=allpixels.dt$NAME)
  metrics.by.fire <- temp.df %>%
                   group_by(NAME) %>% 
                   do(fires.classif.metrics(.$predicted, .$observed, params$critical.class.label))  ## summarize with custom function returning 2 outputs
  
  metrics.by.fire <- as.data.frame(metrics.by.fire[match(as.character(toupper(params$fires)), as.character(metrics.by.fire$NAME)), ]) ## resort fires to respect initial order
  
  metrics <- rbind(cbind(data.frame(NAME="OVERALL"), metrics.overall), metrics.by.fire)   ## stack overall and by fire metrics together
  
  ## save confusion matrices (and whole output by caret's confusionMatrix()) under the name of each fire
  conf.mat.by.fire <- list()
  for (fire in params$fires) {
    cmd <- sprintf('conf.mat.by.fire <- list.append(conf.mat.by.fire, %s=confusionMatrix(temp.df$predicted[temp.df$NAME==toupper(fire)], temp.df$observed[temp.df$NAME==toupper(fire)]))', fire)
    eval(parse(text=cmd))
  }
  
  cmd <- sprintf('RES$%s$metrics <- metrics; RES$%s$conf.mat.by.fire <- conf.mat.by.fire', appr, appr)
  eval(parse(text=cmd))
  
}  ## end for on params$approach

RES.file = file.path(results.dir, 'RESULTS.Rdata', fsep = .Platform$file.sep) 
save(RES, file = RES.file)

#### PRINT LOGS ---------------------------------------------------------

## clock global time
toc <- proc.time()-tic[3]
end.message <- sprintf("Total elapsed time: %s, finished running on %s", seconds_to_period(toc[3]), Sys.time())
print(end.message)




