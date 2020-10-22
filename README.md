# ASTER Data Processing - Script by D.Coutts

The Advanced Spaceborne Thermal Emissions and Reflection Radiometer (ASTER) onboard the TERRA satellite produces VNIR, SWIR, and TIR that is commonly used for geological mapping. This script is used to process L1B data in GeoTiff format into top of atmostphere reflectance (for bands 1-9) and surface temperature (for bands 10-14). 

### Inputs
L1B ASTER scenes downloaded from NASA's EarthData portal contain two files:
1. A zip file containing the geotiffs along with other metadata files (e.g., .met and .txt files)
2. A scene metadata file (.met)

Unzip the zipped file to a working folder and drop the scene metadata file within this folder.
*NOTE: THERE SHOULD ONLY BE ONE SCENE PER FOLDER*

To process this scene, point the script at the working folder.
