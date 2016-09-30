#########################################################
#
#  by 
#################################################################   

#%% IMPORT MODULES

import arcpy, os, sys

sys.path.append(r'D:/Research/MyPythonModules')
from mybasemodule import*

#%% PARAMETERS

params = {}

params['fires'] = ["Abraham", "Cone", "Flett", "Keane", "Leggo", "Levellers", "Liege", "Mcarther", "Overflow", "Perry", "Rail", "Rainbow", "Steephill", "Tanghe"]  ## list of 14 fires to map
#params['fires'] = ["Abraham", "Leggo", "Flett", "Tanghe"] 

root_dir = r'D:\Research\ANALYSES\FiresCD'   ## base directory

fire_img_dir = os.path.join(root_dir, 'Data', r'3_bands')   ## directory of tiff images 
gt_dir = os.path.join(root_dir, 'Data', r'rasterized_API_polygons')  ## directory of GT layers 
maps_dir = os.path.join(root_dir, 'Results', r'Maps')  ## directory of predicted maps 

fire_img_lyr = os.path.join(root_dir, 'MXDs', 'LYRs', r'Abraham_3bands.tif.lyr')  ## path to layer file providing the symbology for the RGB composite
maps_lyr = os.path.join(root_dir, 'MXDs', 'LYRs', r'Abraham_polygon_rasterized.tif.lyr')  ## path to layer file providing the symbology for the categorical rasters (unique values)

new_empty_gr_lyr = os.path.join(root_dir, 'MXDs', 'LYRs', r'NewGroupLayer.lyr')  ## path to empty layer file to be loaded to initialize each group

input_mxd = os.path.join(root_dir, 'MXDs', r'empty.mxd')  ## path to empty mxd to load
output_mxd =  os.path.join(root_dir, 'MXDs', r'Fires.mxd')  ## path to final mxd to save results


#%% START

print('Start')

tic = tic()

file_dirs = [fire_img_dir, gt_dir, maps_dir]  ## create list of directories to iterate on

## Initialize MXD and dataframe to contain the layers
mxd = arcpy.mapping.MapDocument(input_mxd)
df = arcpy.mapping.ListDataFrames(mxd)[0]
mxd.activeView = df.name
mxd.title = df.name


#%% ADD LAYERS AND DEFINE SYMBOLOGY

print('Applying symbology')

for file_dir in file_dirs:  ## iterate over the directories
    arcpy.env.workspace = file_dir   ## set this directory as workspace on which arcpy.ListFiles() will be called on 
    for file_to_add in arcpy.ListFiles('*.tif'):   ## iterate over the files in the each directory
        layer = arcpy.mapping.Layer(os.path.join(file_dir, file_to_add))   ## set layer as an arcpy object 
        if any(x in layer.name for x in params['fires']):   ## do all the rest only if the layer name is in the list of fires
            arcpy.mapping.AddLayer(df, layer, "BOTTOM")   ## add it to the current dataframe df
            if layer.isRasterLayer:   ## do all the symbology operation only if it is a raster layer
                layer_to_update = arcpy.mapping.ListLayers(mxd, layer.name, df)[0]   ## redefine layer_to_update as an arcpy object (1st element, i.e. index 0, of the list of layers with that layer name)
                if file_dir == maps_dir or file_dir == gt_dir:   ## if it is a categorical layer
                    arcpy.ApplySymbologyFromLayer_management(layer_to_update, maps_lyr)  ## apply the symbology of the layer located at the maps_lyr path
                elif file_dir == fire_img_dir:   ## else if it is a rgb layer
                    arcpy.ApplySymbologyFromLayer_management(layer_to_update, fire_img_lyr)  ## first apply the symbology with ApplySymbologyFromLayer_management (needs a path to the layer providing the symbology)...
                    source_layer = arcpy.mapping.Layer(fire_img_lyr)   ## ...then set as an arpy object (source_layer) the layer providing the symbology (path set as fire_img_lyr) because...
                    arcpy.mapping.UpdateLayer(df, layer_to_update, source_layer, True)  ## ...UpdateLayer needs a source layer that should be an arcpy object, not a path!
            arcpy.RefreshTOC()
            arcpy.RefreshActiveView()


#%% GROUP LAYERS    
   
print('Grouping layers')
   
## Group layers by fire
for fire in params['fires']:
    group_layer = arcpy.mapping.Layer(new_empty_gr_lyr)   ## add empty group layer created beforehand
    group_layer.name = fire
    arcpy.mapping.AddLayer(df, group_layer, "BOTTOM")  ## add group layer to dataframe
    group_layer = arcpy.mapping.ListLayers(mxd, fire, df)[0]    ## redefine group_layer as being the group layer named by fire (to avoid error later on)
    list_layers_fire = arcpy.mapping.ListLayers(mxd, fire+'*', df)   ## list layers starting with the name of the fire
    for layer in list_layers_fire:
        if layer.isRasterLayer:    ## avoids taking layers other than raster, including group layers
            if "ClassMap" in layer.name:
                placement = "TOP"
            else:
                placement = "BOTTOM"
            arcpy.mapping.AddLayerToGroup(df, group_layer, layer, placement)
            arcpy.mapping.RemoveLayer(df, layer)    ## remove layer to avoid duplicates
                     
mxd.saveACopy(output_mxd)

toc(tic)

