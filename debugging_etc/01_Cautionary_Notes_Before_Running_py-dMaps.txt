
# ----------------------------------------------
#
# CAUTIONARY NOTES BEFORE RUNNING py-dMaps
# Last Update March,10th 2022 
# ----------------------------------------------

# It is strongly recommended to: 

# (a) Mask your land with 1e+20 instead of NaNs;
# (b) Have you netcdf attribute coordinates in the order: lon,lat,time;
# (c) Squeeze any 1-D depth dimension, i.e. to have a truly 2D*time field (without ghost dimensions) 
#
# Below there is a list of commands useful to preprocess your data for the purposes above.  



# 0. If you want to know: the order of your coordinates, the miss.values and if there's any ghost dimension, then type: 
ncdump -h yourfile.nc 

# 1. Remove NaNs from missing values, and substitute them with a mask number meaningful in py-dMaps;

#For example, set Missing Values to 1e+20 instead of NaN
cdo setmissval,1e+20 input.nc output.nc


# 2. Re-order the dimensions in the NetCDF attributes, because py-dMaps runs with: lon,lat,time (and only like this.)
# type: 

ncks -A -v lon input.nc reorder.nc
ncks -A -v lat input.nc reorder.nc
ncks -A -v time input.nc reorder.nc
#ncks -A -v depth input.nc reorder.nc    # Uncomment this if you also have a ghost depth dimension to squeeze later 
ncks -A -v yourvarname input.nc reorder.nc

ncdump -h reorder.nc # This only checks if the order is now correct in your output file "reorder.nc"


# Note: yourvarname is the name of the variable in your netcdf file (for example, see it with: cdo sinfon yourfile.nc)
# Note: reorder.nc should be always the same file-name, it's over written. input.nc is the netcdf file you want to re-order.  
# Note: check whether it's "lon" or "longitude", "lat" or "latitude" in the attributes and adjust the commands accordingly. 


# 3. Squeeze the depth dimension, in case it is of dim=1. (i.e. you want a truly 2D*time field)
ncwa -a depth input.nc output.nc

#-----------------------------
# *** IMPORTANT NOTES *** 
#-----------------------------

# py-dMaps works with anomalies. 
# If you want to detrend and deseasonalized your data with cdo too, you can do as follows. (as stated in the README on GitHub)
 
# Remove seasonality (i.e. for monthly data): 
cdo -L -ymonsub input.nc -ymonmean input.nc file_deseas.nc

# Remove linear trend too: 
cdo detrend file_deseas.nc file_deseas-and-detrended.nc

# If you want to deseasonalize and detrend your data with cdo, you must do that ***BEFORE*** point 2.
# Otherwise cdo will not be able to work with NetCDF data with coordinates in that order, and you'll have troubles. 
  

# If you work on a remote server, it's likely you'll have to load modules before running the commands above.
# For example:  
module load cdo 
module load nco 



