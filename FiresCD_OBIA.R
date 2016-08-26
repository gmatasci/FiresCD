#### CODE INFOS --------------------------------------------------------------

## Project Name: Fire severity classification 
## Authors: Giona Matasci (giona.matasci@gmail.com), Ignacio San Miguel (ignaciosanmiguel86@gmail.com)   
## File Name:                            
## Objective: 

#### TO DO -------------------------------------------------------------------

# PRIOR TO ACTUAL RUN:
# - 

## STILL TO DO:
# - integrate multiscale approach using seg.ids.df as a starting point
# - introduce LOOCV within meanshift loop instead of OOB? (think if need to redo segmentation or just RF in 2nd LOOCV)
# - check error: Error in `colnames<-`(`*tmp*`, value = c(NA, NA, NA, "elev_sd")) : length of 'dimnames' [2] not equal to array extent
# - uniformize fire names (caps, etc.)
# - split aspect in 2 components sin/cos
# - new pixel level GT with majority class by surface within pixel
# - save these raster files for visualization
# - run as it is on complete set of new fire images (no masks): big fires will impact more the process than the current test on Tanghe and Liege
# - base parameter selection and assessment on a mixed metric: 0.5*Kappa + 0.5*Fmeasure_PartialMortality

## SOLVED:
# -V issue with images/projection -- OTB does not accept Lambert Conical when producing a shp, we have to save a sqlite file and then convert it to shp with ogr2ogr
# -V fill ypred df based in indices and not on rbind -- done
# -V read about OTB meanshift minsize -- the 2 other parameters control the size and shape of the segments: minsize is just used as a post-processing to merge segments smaller than a threshold

#### INIT --------------------------------------------------------------------

print('OBIA for fires CD') 

rm(list=ls())

# param_file <- "D:/Research/ANALYSES/NationalImputationForestAttributes/BAP_Imputation_working/wkg/AllUTMzones_paramsGL.Rdata"
# load(param_file)

# source("D:/Research/ANALYSES/NationalImputationForestAttributes/BAP_Imputation_working/scripts_NationalImputationForestAttributes/Functions_NatImp.R")

#### PARAMETERS ---------------------------------------------

params <- list()

## General
params$approach <- list("PixelBased", "ObjectSingleScale", "ObjectMultiScale")  ## Type of approach (so far only ObjectSingleScale is implemented):
                                                      ## PixelBased --> classic pixel-based approach 
                                                      ## ObjectSingleScale --> optimizes the parameters of a Meanshift segmentation with an object-based approach  
                                                      ## ObjectMultiScale --> features computes for different segmentations are stacked in the same df with a pixel-based approach
# params$fires <- list("Liege")
params$fires <- list("Tanghe", "Liege")  ## fires to consider
params$subsetting <- T    ## to subset data for development purposes (still to implement)
params$mort.classes <- c("0-5%", "6-25%", "95-100%")
params$critical.class <- "Class: 6-25%"  ## to compute F-measure as an alternative, more detailed measure to Kappa

## Meanshift
params$meanshift.testmode <- F    ## if set to T saves each segmentation result (shp) with parameter values to visually inspect the segments
# params$ranger.vect <- c(10, 50, 200) ## On SF imagebest with 50
# params$spatialr.vect <- c(5, 20 , 50)  ## best with 5  
# params$minsize.vect <- c(10, 100, 1000) ## best with 100
# params$ranger.vect <- c(10, 20, 200, 500) ## On SF imagebest with 50
# params$spatialr.vect <- c(3, 10, 20, 50)  ## best with 5
# params$minsize.vect <- c(10) ## best with 100
params$ranger.vect <- c(5, 20) ## Range radius: on SF imagebest with 50, Liege with 25
params$spatialr.vect <- c(5, 20)  ## Spatial radius: on SF best with 5, Liege with 5
params$minsize.vect <- c(5) ## Minimum object size: on SF best with 100, for fireswith 5, setting it to 0 just produces 1-pixel segments

## RF 
params$base.predictors <- c("db1", "db2", "db3", "db4", "db5", "db7", 
                               "dNBR", "dTCG", "dTCB", "dTCW", "dNDVI", "dNDWI", 
                               "elev", "slope", "aspect", 
                               "EOSD")   ## list of starting predictors (for the classification part, as the segmentation is run on 6-band difference image only)
params$obj.predictors <- c("db1_mean", "db1_sd", "db2_mean", "db2_sd", "db3_mean", "db3_sd", "db4_mean", "db4_sd", "db5_mean", "db5_sd", "db7_mean", "db7_sd", 
                              "dNBR_mean", "dNBR_sd", "dTCG_mean", "dTCG_sd", "dTCB_mean", "dTCB_sd", "dTCW_mean", "dTCW_sd", "dNDVI_mean", "dNDVI_sd", "dNDWI_mean", "dNDWI_sd",
                              "elev_mean", "elev_sd", "slope_mean", "slope_sd", "aspect_mean", "aspect_sd",
                              "EOSD", "EOSD_majpct",
                              "nrpixseg")   ## list of final predictors computed at the object level
# params$metric <- "Kappa"   ## not implemented yet
# params$metric <- "MixKappaFm"  
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
data.dir <- file.path(base.dir, "Data", fsep = .Platform$file.sep)   ## data directory (nothing should be written here)
results.dir <- file.path(base.dir, "Results", fsep = .Platform$file.sep)   ## results directory (outputs go here)
figures.dir <- file.path(base.dir, "Figures", fsep = .Platform$file.sep)   ## figures directory (figures go here)
temp.dir <- file.path(data.dir, "temp")  ## directory for temporary files like the segmentation shps (overwritten each time)
if (!file.exists(temp.dir)) {dir.create(temp.dir, showWarnings=F, recursive=T)}  ## create it

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

## Mean Shift with OTB 
meanShiftOTB <- function(ranger, spatialr, minsize, OTB.dir, OGR.dir, data.dir, temp.dir, img.file.name, seg.file.name, test.mode) {
  
  img.file <- file.path(data.dir, img.file.name)   ## image to segment
  segment.file.sqlt <- file.path(temp.dir, sprintf("%s.sqlite", seg.file.name))  ## initial segmentation result as a sqlite file (only way to work with Lambert conformal conic projection)
  segment.file <- file.path(temp.dir, sprintf("%s.shp", seg.file.name))  ## segmentation result converted to shp via ogr2ogr
  
  ## delete existing segmentation files
  if (file.exists(segment.file)){ 
    unlink(segment.file, recursive = T, force = T) 
    unlink(segment.file.sqlt, recursive = T, force = T)
  } 
  
  ## create text for command to run OTB (equivalent to saving this text in a batch file .bat and running it with system2() )
  command.text <- sprintf("%s -in %s -filter meanshift -filter.meanshift.ranger %d -filter.meanshift.spatialr %d -filter.meanshift.minsize %d -mode vector -mode.vector.out %s -mode.vector.layername %s & %s %s %s", 
                      file.path(OTB.dir, "otbcli_Segmentation"),   ## complete path to OTB bat file doing the segmentation
                      img.file,
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

## summarize numerical variables
summarize.custom.num <- function(df, gvar1, gvar2, uniq.var, oper) {
  df %>%
    group_by_(gvar1, gvar2) %>%
    summarise_(npixel = interp(~get(oper)(number), number=as.name(uniq.var)))
}

## summarize factor variables
summarize.custom.fac <- function(df, gvar1, gvar2, gvar3, uniq.var, oper) {
  df %>%
    group_by_(gvar1, gvar2, gvar3) %>%
    summarise_(npixel = interp(~get(oper)(number), number=as.name(uniq.var))) %>%  
    slice(which.max(npixel))
}

fires.classif.metrics <- function(predicted, observed, critical.class) {
  RES <- confusionMatrix(predicted, observed)
  sens <- RES$byClass[critical.class, "Sensitivity"]
  spec <- RES$byClass[critical.class, "Specificity"]
  return(data.frame(FmeasCrit=(2*sens*spec)/(sens+spec), Kappa=as.vector(RES$overall[2])))
}

#### START --------------------------------------------------------------

tic <- proc.time() ## start clocking global time

load(file.path(data.dir, "allpixels.Rdata"))

#### LOOCV OVER FIRES ---------------------------------------------------

## initialize empty factor vector with appropriate levels to store final class predictions at each round of the Leave-one-out cross-validation loop
Y.predicted <- factor(rep(NA, nrow(allpixels)), levels=levels(allpixels$CLASS))  
for (fire.out in params$fires) {  ## LOO-CV loop over the fires to leave out

  fires.in <- params$fires[!params$fires %in% fire.out]   ## fires to keep in at this round of the loop 
  idx.pix.in <- allpixels$NAME %in% toupper(fires.in)    ## logical indices of the pixels to keep in
  idx.pix.out <- !allpixels$NAME %in% toupper(fires.in)   ## logical indices of the pixels to leave out
  
  allpixels.in <- allpixels[idx.pix.in,]   ## df with pixels in
  allpixels.out <- allpixels[idx.pix.out,]  ## df with pixels out

#### SEGMENTATION PARAMETERS OPTIMIZATION --------------------------------

  multi.scale.pix.df <- data.frame(matrix(nrow=0, ncol=0))   ## stil to implement
  OAs.OOB.meanshift <- data.frame(matrix(nrow=0, ncol=4))   ## initialize empty matrix to store Out-of-Bag Overall Accuracies and associated parameters for each segmentation
  colnames(OAs.OOB.meanshift) <- c("ranger", "spatialr", "minsize", "OA")
  ## df storing segment IDs for all the pixels kept in (as many columns as combinations of parameters)
  seg.ids.df <- data.frame(matrix(nrow=nrow(allpixels.in), ncol=length(params$ranger.vect)*length(params$spatialr.vect)*length(params$minsize.vect)))  
  ms.par.col <- 1  ## start at 1 the index over columns of dataframe of segment IDs
  
  ## fill object level df with per segment mean and std dev
  multi.scale.pix.dt <- data.table()  ## initialize dt to store object-level features for all the kept-in pixels
  
  ## 3 nested loops over the Meanshift parameters
  for (ranger in params$ranger.vect) {
    for (spatialr in params$spatialr.vect) {
      for (minsize in params$minsize.vect) {
        
        ## loop to segment all the fires kept in
        for (fire.in in fires.in) {
          
          img.file.name <- sprintf("%s_6bands.tif", fire.in)
          if (params$meanshift.testmode) {  ## if TRUE save segmentation result with a telling name 
            seg.file.name <- sprintf("Segment_%s_ranger%d_spatial%d_minsize%d", fire.in, ranger, spatialr, minsize)
          } else {   ## otherwise just overwrite it in a file called Segment
            seg.file.name <- "Segment"
          }
          segments <- meanShiftOTB(ranger, spatialr, minsize, OTB.dir, OGR.dir, data.dir, temp.dir, img.file.name, seg.file.name, params$meanshift.testmode)  ## main Meanshift segmentation

          fire.image.proj <- CRS(proj4string(segments))   ## get projection string from segments
          
          ## create spatial points object with coordinates of all the pixels in current fire.in to extract segments
          coords <- SpatialPoints(allpixels.in[allpixels.in$NAME==toupper(fire.in), c("X", "Y")], proj4string=fire.image.proj)  
          id.segment <- over(coords, segments)$dn  ## get segment IDs contained in dn column of the output of function over()
          
          seg.ids.df[allpixels.in$NAME==toupper(fire.in), ms.par.col] <- id.segment  ## dynamically fill the part of the column number ms.par.col of the df of segment IDs corresponding to fire.in
          
        }
        
        colnames(seg.ids.df)[ms.par.col] <- sprintf("rr%ssr%sms%s", ranger, spatialr, minsize)  ## name that same column with a string containing Meanshift parameter values
        
        allpixels.in$segID <- seg.ids.df[, ms.par.col]  ## temporarily add segment IDs to allpixel.in df
        allpixels.in$npixel <- 1   ## add column of ones to be used to count the nr of pixels per segment
        
        ## fill object level df with per segment mean and std dev
        single.scale.obj.df <- data.frame()  ## initialize df to store object-level features
        count <- 0  ## initialize iterator
        
        ## loop over the columns to summarize only
        for (col in colnames(allpixels.in)[!(colnames(allpixels.in) %in% c("NAME", "segID", "X", "Y"))]){
          if (col %in% c("CLASS", "EOSD")){ ## if the column is a factor apply the summarizing function for factors with majority as summarizer
            target <- summarize.custom.fac(allpixels.in, "NAME", "segID", col, "npixel", "sum")
            colnames(target)[4] <- paste(col, "_npixel", sep="")
          } else if (col %in% c("npixel")){   ## else if it is the columns of ones apply the summarzing function for continuous variables with sum as summarizer
            target <- summarize.custom.num(allpixels.in, "NAME", "segID", col, "sum")
          } else {  ## else if it is a continuous column apply the summarzing function for continuous variables with mean (final column ending with "_mean") and standard deviation ("_sd") as summarizers
            target1 <- summarize.custom.num(allpixels.in, "NAME", "segID", col, "mean")
            colnames(target1)[3] <- paste(col, "_mean", sep="") 
            target2 <- summarize.custom.num(allpixels.in, "NAME", "segID", col, "sd")
            colnames(target2)[3] <- paste(col, "_sd", sep="") 
            target <- data.frame(target1[,1:3], target2[,3])
          }
          if (count == 0){  ## if first round of the loop single.scale.obj.df takes the values of the df target
            single.scale.obj.df <- rbind(single.scale.obj.df, as.data.frame(target))
          } else {   ## then, just append the newly computed columns columnwise (cbind)
            newcols <- !(colnames(as.data.frame(target)) %in% colnames(single.scale.obj.df) )
            single.scale.obj.df <- cbind(single.scale.obj.df, target[,newcols])
          }
          count <- count + 1
        }
        colnames(single.scale.obj.df)[ncol(single.scale.obj.df)] <- "nrpixseg"   ## last column is the number of pixels per segment
        single.scale.obj.df$CLASS_npixel <- single.scale.obj.df$CLASS_npixel / single.scale.obj.df$nrpixseg  ## compute percentages of the majority class for both... 
        colnames(single.scale.obj.df)[colnames(single.scale.obj.df) == "CLASS_npixel"] <- "CLASS_majpct"  ## ...the ground truth...
        single.scale.obj.df$EOSD_npixel <- single.scale.obj.df$EOSD_npixel / single.scale.obj.df$nrpixseg
        colnames(single.scale.obj.df)[colnames(single.scale.obj.df) == "EOSD_npixel"] <- "EOSD_majpct"  ## ...and the EOSD layer
        
        ## if segmentation level is among a specified set, store for every pixel the summarized predictor values at the object level
        
        data.table(single.scale.obj.df"NAME", "segID"
        
        ## set mtry parameter according to params$mtry
        nr.vars <- length(params$obj.predictors) 
        if (params$mtry == 'sqrt_nr_var') {
          mtries <- floor(sqrt(nr.vars))
        } else if (params$mtry == 'nr_var_div_3') {
          mtries <- floor(nr.vars/3)
        }
        set.seed(params$seed)   
        ## apply RF on df with object-level values using as predictors the columns listed in params$obj.predictors and with response variable the column specified in params$targ
        RF <- randomForest(x=single.scale.obj.df[,params$obj.predictors], y=single.scale.obj.df[,params$targ], ntree=params$ntree, mtry=mtries, nodesize=params$nodesize, importance=params$plot.importance)
        OAs.OOB.meanshift <- rbind( OAs.OOB.meanshift, data.frame(ranger, spatialr, minsize, OA=1-RF$err.rate[params$ntree, "OOB"]) )   ## grow df with results
      
        ms.par.col <- ms.par.col+1  ## increment index over columns of dataframe of IDs
        
      }  ## end for on ranger
    }  ## end for on spatialr
  }  ## end for on minsize
  
  OAs.OOB.meanshift.sorted <- arrange(OAs.OOB.meanshift, desc(OA))   ## sort df with results by decreasing OA
  
#### SEGMENTATION WITH BEST PARAMETERS ---------------------------------------------
  
  ## segment left out fire with best parameters
  img.file.name <- sprintf("%s_6bands.tif", fire.out)
  if (params$meanshift.testmode) {
    seg.file.name <- sprintf("Segment_%s_ranger%d_spatial%d_minsize%d", fire.out, OAs.OOB.meanshift.sorted$ranger[1], OAs.OOB.meanshift.sorted$spatialr[1], OAs.OOB.meanshift.sorted$minsize[1])
  } else {
    seg.file.name <- "Segment"
  }
  segments <- meanShiftOTB(OAs.OOB.meanshift.sorted$ranger[1], OAs.OOB.meanshift.sorted$spatialr[1], OAs.OOB.meanshift.sorted$minsize[1], OTB.dir, OGR.dir, data.dir, temp.dir, img.file.name, seg.file.name, params$meanshift.testmode)
  fire.image.proj <- CRS(proj4string(segments)) 
  coords <- SpatialPoints(allpixels.out[, c("X", "Y")], proj4string=fire.image.proj)  ## create spatial points for extrating segments from OTB
  id.segment.out <- over(coords, segments)$dn
  allpixels.out$segID <- id.segment.out  ## add segment IDs to the df with the left-out pixels
  allpixels.out$npixel <- 1   ## add column of ones to count nr of pixels per segment
  
  ## retrieve segment IDs for kept-in fires associated with the best combination of parameters
  best.colname <- sprintf("rr%ssr%sms%s", OAs.OOB.meanshift.sorted$ranger[1], OAs.OOB.meanshift.sorted$spatialr[1], OAs.OOB.meanshift.sorted$minsize[1])
  id.segment.in <- seg.ids.df[, best.colname]   
  allpixels.in$segID <- id.segment.in  ## add IDs to the kept-in pixels df
  
  allpixels.final <- rbind(allpixels.in, allpixels.out)  ## stack together the in and out pixels (unique segments are identified by fire NAME and segID) 
  
  ## fill object-level df with per segment mean and std dev (same as above) for both in and out fires 
  count <- 0
  single.scale.obj.df <- data.frame()
  for (col in colnames(allpixels.final)[!(colnames(allpixels.final) %in% c("NAME", "segID", "X", "Y"))]){
    if (col %in% c("CLASS", "EOSD")){
      target <- summarize.custom.fac(allpixels.final, "NAME", "segID", col, "npixel", "sum")
      colnames(target)[4] <- paste(col, "_npixel", sep="")
    } else if (col %in% c("npixel")){
      target <- summarize.custom.num(allpixels.final, "NAME", "segID", col, "sum")
    } else {
      target1 <- summarize.custom.num(allpixels.final, "NAME", "segID", col, "mean")
      colnames(target1)[3] <- paste(col, "_mean", sep="") 
      target2 <- summarize.custom.num(allpixels.final, "NAME", "segID", col, "sd")
      colnames(target2)[3] <- paste(col, "_sd", sep="") 
      target <- data.frame(target1[,1:3], target2[,3])
    }
    if (count == 0){
      single.scale.obj.df <- rbind(single.scale.obj.df, as.data.frame(target))
    } else {
      newcols <- !(colnames(as.data.frame(target)) %in% colnames(single.scale.obj.df) )
      single.scale.obj.df <- cbind(single.scale.obj.df, target[,newcols])
    }
    count <- count + 1
  }
  colnames(single.scale.obj.df)[ncol(single.scale.obj.df)] <- "nrpixseg"
  single.scale.obj.df$CLASS_npixel <- single.scale.obj.df$CLASS_npixel / single.scale.obj.df$nrpixseg
  colnames(single.scale.obj.df)[colnames(single.scale.obj.df) == "CLASS_npixel"] <- "CLASS_majpct"
  single.scale.obj.df$EOSD_npixel <- single.scale.obj.df$EOSD_npixel / single.scale.obj.df$nrpixseg
  colnames(single.scale.obj.df)[colnames(single.scale.obj.df) == "EOSD_npixel"] <- "EOSD_majpct"
  
  

#### CLASSIFICATION OF LEFT-OUT FIRE ----------------------------------------
  
  nr.vars <- length(params$obj.predictors)
  if (params$mtry == 'sqrt_nr_var') {
    mtries <- floor(sqrt(nr.vars))
  } else if (params$mtry == 'nr_var_div_3') {
    mtries <- floor(nr.vars/3)
  }
  segments.in <- single.scale.obj.df$NAME %in% toupper(fires.in)   ## retrieve indices of kept-in segments in the object-level df
  segments.out <- single.scale.obj.df$NAME == toupper(fire.out)    ## retrieve indices of left-out segments in the object-level df
  
  ## train RF on in segments and predict on out segments, always with params$obj.predictors
  set.seed(params$seed)
  RF <- randomForest(x=single.scale.obj.df[segments.in, params$obj.predictors], y=single.scale.obj.df[segments.in, params$targ], ntree=params$ntree, mtry=mtries, nodesize=params$nodesize, importance=params$plot.importance)
  Y.predicted.segments.out <- predict(RF, single.scale.obj.df[segments.out, params$obj.predictors], type="response", predict.all=F, nodes=F)
  
  ## build a df with the predicted class for each segment (Y.predicted.segments.out) and the associated segment ID (single.scale.obj.df[segments.out, "segID"])
  y.pred.segID.df <- data.frame(segID=single.scale.obj.df[segments.out, "segID"], ypred=Y.predicted.segments.out)
  
  ## assign the predicted class to the "data" df of the segments shp by merging by segment ID ("dn" or "segID"), all.x=T is used to keep the segments for which there is no prediction (outside of fire)
  segments@data <- merge(segments@data, y.pred.segID.df, by.x="dn", by.y="segID", all.x=T)  
  segments@data$ypred[is.na(segments@data$ypred)] <- params$mort.classes[1]    ## if NA are assigned to polygons outside of fire, assign the lowest mortality class instead
  writeOGR(segments, results.dir, sprintf("%s_pred_map_%s", fire.out, best.colname), driver="ESRI Shapefile", overwrite_layer=TRUE)   ## write prediction map shapefile for the left-out fire
  
  y.pred.segID.dt <- data.table(y.pred.segID.df, key = "segID") ## create data.table with key being the segment IDs column
  allpixels.out.segID.dt <- data.table(segID=allpixels.out$segID, key = "segID") #
  
  Y.predicted[idx.pix.out] <- allpixels.out.segID.dt[y.pred.segID.dt, ypred]  ## join to data.table based on a common key with this command allpixels.out.segID.dt[y.pred.segID.dt], then select only ypred as a column
  
} ## end for on fire.out

#### ASSESSMENT -----------------------------------------------------------

## overall assessment
RES.LOOCV <- fires.classif.metrics(Y.predicted, allpixels$CLASS, params$critical.class)

## by fire assessment
temp.df <- data.frame(predicted=Y.predicted, observed=allpixels$CLASS, NAME=allpixels$NAME)
metrics.by.fire <- temp.df %>%
                 group_by(NAME) %>% 
                 do(fires.classif.metrics(.$predicted, .$observed, params$critical.class))  ## summarize with custom function returning 2 outputs

metrics.by.fire <- as.data.frame(metrics.by.fire[match(as.character(toupper(params$fires)), as.character(metrics.by.fire$NAME)), ]) ## resort fires to respect initial order

conf.mat.by.fire <- list()
for (fire in params$fires) {
  conf.mat.by.fire <- list.append(conf.mat.by.fire, confusionMatrix(temp.df$predicted[temp.df$NAME==toupper(fire)], temp.df$observed[temp.df$NAME==toupper(fire)]))
}

RES.FINAL <- rbind(cbind(data.frame(NAME="OVERALL"), RES.LOOCV), metrics.by.fire)   ## stack overall and by fire metrics together

RES.file = file.path(results.dir, 'RESULTS.Rdata', fsep = .Platform$file.sep) 
save(RES.FINAL, file = RES.file)

# ## implement mixed measure to get final ranking?
# if (params$metric == "Kappa") {
#   metric <- as.vector(RES.LOOCV$Kappa)  ## to remove the column name
# } else if (params$metric == "MixKappaFm") {
#   metric <- as.vector(0.5*RES.LOOCV$Kappa + 0.5*RES.LOOCV$Fmeas)
# }

#### PRINT LOGS ---------------------------------------------------------

## clock global time
toc <- proc.time()-tic[3]
end.message <- sprintf("Total elapsed time: %s, finished running on %s", seconds_to_period(toc[3]), Sys.time())
print(end.message)




