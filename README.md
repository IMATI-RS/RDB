# RDB
MATLAB scripts for computing Remote Derived Bathymetry (RDB) from multispectral optical imagery, specifically designed for shallow coastal and riverine environments.

# Overview
The proposed methodology builds upon the band ratio technique of *Stumpf et al. (2003)* (https://doi.org/10.4319/lo.2003.48.1_part_2.0547) and enables bathymetry retrieval from multispectral optical imagery in shallow water areas, such as coastal or riverine environments. This approach leverages the correlation between water depth and the ratio of log-transformed reflectance values from appropriate spectral bands, using in situ ground truth bathymetric data for model calibration and validation.

-	The script **RDB_BandRatio_Lin.m** aims to faithfully replicate the method proposed by *Stumpf et al. (2003)* and utilizes a linear expression to correlate water depth with the value of the spectral band ratio.
-	In contrast, the script **RDB_BandRatio_Exp.m** uses an exponential regression model, which can prove more flexible in fitting the observed data trends depending on the application context.

Both scripts also include a preliminary preprocessing phase, which enables the resampling of high-density ground truth bathymetric points to align with the spatial resolution of the multispectral imagery. This prevents multiple depth values from being assigned to a single pixel.

# Requirements
MATLAB (Curve Fitting, Map and Statistics Toolboxes) is required in order to run the scripts.

# User-specified inputs
All user-specified inputs are located in the “INPUT DATA AND PARAMETERS” section of the codes. These include:
-	The name of the orthophoto (.tif file) which must be preprocessed beforehand by applying a land mask to each spectral band (masked pixels must be assigned a value of zero).
-	The name of the vector file (.shp file) containing ground truth bathymetric points to be used for model calibration and validation.

Additionally, users can also specify:
-	A series of parameters to filter the ground truth bathymetric points upstream (e.g., minimum and maximum water depth values).
-	The method to be used for resampling the ground truth bathymetric points (options include "median" and "mean").
-	The percentage of available ground truth bathymetric points to be randomly selected and used for model calibration; the remaining points will automatically be used for model validation.
-	The spectral band ratio to be considered (options include "blue/red" and "blue/green").

The scripts include additional comments to help users understand each step in the process.
