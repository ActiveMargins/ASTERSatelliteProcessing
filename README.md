## ASTER Handling Functions

This project was an attempt to work on functional scoping, the details of the ASTER sensor, and geospatial techniques in R. 


### THE ASTER SENSOR
The [Advanced Spaceborne Thermal Emissions and Reflection Radiometer (ASTER)](https://asterweb.jpl.nasa.gov/) onboard the TERRA satellite produces visible near infrared (VNIR), shortwave infrared (SWIR), and thermal infrared (TIR) that, among other uses, is commonly used for geological mapping. This repository contains functions that I use to process and interpret ASTER data (freely available from NASA). I've tried to write the functions to use a few packages as possible, but I incorporate terra and ggplot2. Most importantly I use the [terra package](https://rspatial.org/terra/pkg/index.html) to handle rasters and perform some calculations. 

The ASTER sensor is quite important for geologic mapping as the range of wavelengths from NIR to TIR are very helpful for the characterization of the geologic units. Additionally, the sensor has quite high spatial resolution (NIR - 15m, SWIR - 30m, TIR - 90m), which allows for reasonable interpreations to be made. 

### THE ASTER LIBRARY
The ASTER Library is one of the strengths of the ASTER/TERRA platform. A duplicate sensor has beeen used by  NASA to collect reflectance spectra profiles samples of vegetation, soil, man-made materials, and minerals. This allows for simple classfication methods to be used for mapping ASTER scenes. The ASTER library comes in the form of numerous .txt files. In the case of minerals, multiple example of each mineral (e.g., albite) have been measured. Each record is an individual .txt file that with one header line that contains the record number, the mineral name, and sample name. The reflectance in the first 9 bands (SWIR and NIR) are individual rows. Thermal bands (10-14) are excluded in the library.

<p align="center">
  <img width="1000" height="525" src="https://github.com/ActiveMargins/ASTERSatelliteProcessing/blob/master/images/ASTER_workflow.png">
</p>

### FUNCTIONS

### 1) L1B_to_TOA(scene_loc,save=FALSE,save_loc=NULL)

**Description**

This function takes a folder of L1B data in geotiff format and returns a raster stack object (from Terra package) that has been converted to top of atmosphere (TOA) reflectance. The function arguments "save" and "save_loc" allow for the raster stack to be saved to a specific location once processed. Due to lazy evaluation of arguments , "save_loc" does not need to be specified if save=FALSE.  

**Inputs**

scene_loc  -  A folder location for input scene .tiff files and associated .met files.
save  -  A boolean argument that specifies if the processed scene bands should be saved in .tiff format in a output folder (save_loc)
save_loc  -  A folder location for processed images if save=TRUE

L1B ASTER scenes downloaded from NASA's EarthData portal contain two files:
1. A zip file containing the geotiffs along with other metadata files (e.g., .met and .txt files)
2. A scene metadata file (.met)

Unzip the zipped file to a working folder and drop the scene metadata file within the folder with the geotiff files..

*NOTE: THERE SHOULD ONLY BE ONE SCENE PER FOLDER*

**Outputs**

The output of this function is a single raster stack object (object class outlined in the terra package) consisting of a raster for each band. The data will be processed to TOA reflectance.

### 2) folder_to_stack(scene_loc)

**Description**

This function loads procssed .tifs files (TOA reflectance) and compiles them into a raster stack (as per terra). 

**Inputs**

scene_loc  -  A folder location for processed ASTER .tifs. It's best if these have the same naming as when they are saved by L1B_to_TOA(). The minimum requirement is that the final characters in the file names are the band numbers separated from the scene number by an underscore.

**Outputs**

The output is a raster stack.

### 3) create_cover_mask(stack,calc_NDVI=FALSE,NDVI_value=NULL,calc_NSDI=FALSE,NSDI_value=NULL)

**Description**

This function creates a binary raster that will be used to mask other rasters. The mask is used mask vegetation and snow.

**Inputs**

stack  -  a processed raster stack
calc_NDVI=FALSE  -  boolean to specify if NDVI should be included in the mask.
NDVI_value=NULL  -  NDVI values used to threshold vegetation/non-vegetation. 0.3 is a traditional value.
calc_NDSI=FALSE  -  boolean to specify if NDSI should be included in the mask.
NDSI_value=NULL  -  NDSI value used to threshold snow/non-snow. 0.4 is a traditional value. 

**Outputs**

A single raster that can be used to mask indices or ASTER bands. Areas of ground cover (vegetation and snow) have values of 0, while areas of exposed rock have a value of 1.

### 4) calculate_indices(stack, outliers.rm=TRUE)

**Description**

This function calculates numerous indicies that can be used to map exposed bedrock. The mineral/alteration indicies are well explained in the code.

**Inputs**

stack  -  a processed raster stack (TOA reflectance)
outliers.rm=TRUE  -  boolean to discern whether outliers should be removed. In some cases there are index values that are calculated that are not representative of the geology. Largely occuring in the margins of the scene.  

**Outputs**

This function creates numerous rasters in the global environment and may challenge the ram of your workstation...

### 5) mask_raster(raster, cover_mask)

**Description**

This function applies a mask an input raster.

**Inputs**

raster  -  input single raster (e.g., index) or satck of multiple rasters
cover_mask  -  cover mask generated through the create_cover_mask() function.

**Outputs**

A raster that has been masked by the binary mask.

### 6) raster_PCA(stack)

**Description**

This function performs PCA on a raster stack. Only the first three components are computed.

**Inputs**

stack  -  a raster stack

**Outputs**

An PCA object. Individual principle components can be accessed through "object$PC1" through "object$PC3". Individual principle components can be plotted by "plot(object$PC1)". RGB can be plotted through "plotRGB(object)".


### 7) assemble_ASTER_lib(ASTER_lib_folder)

**Description**

This function assembles an ASTER library object from a folder containing ASTER library .txt. The function takes individual .txt and creates a tidy dataframe from them that can be used in plotting or classification.

**Inputs**

ASTER_lib_folder  -  a folder location that contains ASTER library .txt files.

**Outputs**

The output of this function is a single tidy dataframe the lists the reflectance of each band each .txt. file input into the folder. The Ref_offset is used for plotting spectral profiles if the user desires them to be offset and not overlapping. This dataframe is referred to later as a *ASTER library object*

| Reflectance | Mineral    | Band | Wavelength | Ref_offset |
| ----------- | ---------- | ---- | ---------- | ---------- |
| 0.53672278  | Actinolite | 1    | 0.560      | 0.5367228  |
| 0.51193970  | Actinolite | 2    | 0.660      | 0.5119397  |
| ...         | ...        | ...  | ...        | ...        |

### 8) plot_ASTER_lib(lib_obj, offset=FALSE)

**Description**

Thie fucntion creates a plot of desired spectra from an ASTER library object (described in the previous function).

**Inputs**

lib_obj  -  An ASTER library object (a dataframe formated similar to the above)

offset=FALSE  -  A boolen argument that specifies if the spectra should be overlapping/not offset. It is common for spectra to be offset. If offset = FALSE, the plot will be reflectance vs. wavelength. If offset=TRUE, then it will be Ref_offset vs. wavelength. 

**Outputs**

This function creates a single ggplot for all the mineral species that are within the ASTER library object. Lines will be coloured by mineral. 

### 9) SAM(lib_obj, input_scene_spectra)

**Description**

This function performs a spectral angle mapper classification for input spectra. The input spectra will be compared to the 

**Inputs**

lib_obj  - An ASTER library object (a dataframe formated similar to the above)

input_scene_spectra - a 1D vector of reflectance values (ordered 1-9).

**Outputs**

A list containing two parts:
1. A dataframe that descibes the spectral angle between the minerals within the library object and the input spectra
2. A a classification (i.e., the mineral with the smallest angle)
