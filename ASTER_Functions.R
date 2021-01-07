library(stringr)
library(dplyr)
library(tidyr)
library(lubridate)
library(terra)
library(ggplot2)

L1B_to_TOA <- function(scene_loc, save=FALSE, save_loc=NULL){
    #1.1 Load constants/tables
    #1.1.1 Calculated unit conversion coefficient table
    df_concoeff_table <- data.frame(Band=c(1,2,3,4,5,6,7,8,9), 
                                    HGH=c(0.676,0.708,0.423, 0.1087, 0.0348, 0.0313, 0.0299, 0.0209,0.0159), 
                                    NOR=c(1.688, 1.415, 0.862, 0.2174, 0.0696, 0.0625, 0.0597, 0.0417,0.0318), 
                                    LOW=c(2.25,1.89,1.15,0.290,0.0925,0.0830,0.0795, 0.0556,0.0424))
    
    #1.1.2 ESUN table (Smith)
    df_esun <- data.frame(Band=c(1,2,3,4,5,6,7,8,9),
                          ESUN=c(1845.99,1555.74,1119.47,231.25,79.81,74.99,68.66,59.74,56.92))
    
    #1.1.3 ASTER thermal bands K1, K2, and effective wavelength
    df_thermal_table <- data.frame(Band=c(10,11,12,13,14),
                                   Gain=c(0.006822, 0.006780, 0.006590, 0.005693, 0.005225),
                                   EffWave=c(8.291, 8.634, 9.075, 10.657, 11.318),
                                   K1=c(3040.136402, 2482.375199, 1935.060183, 866.468575, 641.326517),
                                   K2=c(1735.337945, 1666.398761, 1585.420044, 1350.069147, 1271.221673))
    
    #1.1.4 Create dummy band to help in resizing (will be one of bands 1-3)
    ras_dummy <- NULL
    
    #1.2 Create filelist to aid in navigation
    #1.2.1 Rename .met files to .txt so they are readable
    bandmetadata_filelist <- list.files(scene_loc,pattern=".met$", full.names = TRUE)
    newfiles <- gsub(".met$", ".txt", bandmetadata_filelist)
    file.rename(bandmetadata_filelist, newfiles)
    bandmetadata_filelist <- list.files(scene_loc,pattern="tif.txt$", full.names = TRUE)
    
    #1.2.2 Get the name of all the tif files
    tif_filelist <- list.files(scene_loc,pattern=".tif$", full.names = TRUE)
    
    #1.2.3 Remove band 3B from the file list (tif_filelist).
    band_file_list <- str_sub(tif_filelist, start=-6,end=-5)
    band_3B_file <- which(band_file_list=="3B")
    tif_filelist <- tif_filelist[-band_3B_file]
    
    #1.2.3 Organize the file list (tif_filelist) so that band 1 is first.
    band_file_list <- str_sub(tif_filelist, start=-7)
    band_file_list <- regmatches(band_file_list,gregexpr("[[:digit:]]+", band_file_list))
    band_file_list <- as.numeric(unlist(band_file_list))
    band_1_file <- which(band_file_list==1)
    tif_filelist <- c(tif_filelist[band_1_file],tif_filelist[-band_1_file])
    
    #1.2.3 Get the satellite meta file
    scene_meta <- list.files(scene_loc,pattern=".zip.txt$", full.names = TRUE)
    scene_meta <- readLines(con = scene_meta)
    scene_id <- substr(scene_meta[grep("VALUE",scene_meta)[1]],43,48)
    
    #1.3 Extract constants used in calculation of radiance at sensor, reflectance, and surface temp
    #1.3.1 Extract conversion ceofficients from the aster scene metadata (scene_meta)
    concoeff_rows <- grep("ASTERGains",scene_meta)
    concoeff_filter <- scene_meta[concoeff_rows:(concoeff_rows+10)]
    concoeff_filter <- concoeff_filter[grep("3B",concoeff_filter)]
    concoeff_filter <-substr(concoeff_filter,nchar(concoeff_filter)-81,nchar(concoeff_filter)) 
    concoeff <- concoeff_filter %>% str_replace_all("[^A-Z1-9,]","" ) %>% str_split(pattern=",") #also there maybe a second low (low2)
    df_concoeff <- data.frame(bands=concoeff, stringsAsFactors = FALSE) 
    colnames(df_concoeff) <- "bands"
    df_concoeff <- df_concoeff %>% separate(bands, into=c("bands","gain"), sep=1)
    df_concoeff <- df_concoeff[-4,]
    df_concoeff$bands <- as.numeric(df_concoeff$bands)
    
    if(df_concoeff[3,2]=="NLOW"){ #also there maybe a second low (low2)
        df_concoeff[3,2] <- "LOW"
    } else if (df_concoeff[3,2]=="NNOR"){
        df_concoeff[3,2] <- "NOR"
    } else {
        df_concoeff[3,2]=="HGH"
    }
    
    #1.3.2 Extract satellite data from the first band_meta file
    bandmeta_file <- bandmetadata_filelist[1]
    band_meta <- readLines(con = bandmeta_file) #read in ugly band metadatafile
    
    #1.3.3 Extract the Julian Date by a regrex filter and clculate solar distance (d)
    JD_rows <- grep("CALENDARDATE",band_meta)
    JD_filter <- band_meta[JD_rows[1]:JD_rows[2]]
    JD_filter <- JD_filter[grep("VALUE",JD_filter)]
    JD <- str_extract_all(JD_filter,"[-\\.\\d]", simplify=TRUE)
    JD <- paste(JD, collapse="")
    JD <- yday(ymd(JD))
    d <- (1-0.01672*cos((0.986*(JD-1))*pi/180))
    
    #1.3.4 Extract Solar Angle by a regex filter and calculate ____ (z) 
    SA_rows <- grep("SOLAR_ELEVATION_ANGLE",band_meta)
    SA_filter <- band_meta[SA_rows[1]:SA_rows[2]]
    SA_filter <- SA_filter[grep("VALUE",SA_filter)]
    SA <- str_extract_all(SA_filter,"[-\\.\\d]", simplify=TRUE)
    SA <- as.numeric(paste(SA, collapse=""))
    z <- (90-SA)*pi/180
    
    #2.0 Process reflectance (Bands 1-10) or surface temperature (Bands 11-13) for every tif file in the file list
    for (i in 1:length(tif_filelist)){
        #2.1 Load band .tif and extract band number 
        #2.1.1 Read  and project the raster that contains the digital number
        tif_file <- tif_filelist[i]
        ras_DN <- rast(tif_file)
        ras_DN <- project(ras_DN, crs(ras_DN))
        ras_DN[ras_DN==0] <- NA
        
        #2.1.2 Extract band number from the .tif file name
        band <- substr(tif_file, nchar(tif_file)-10, nchar(tif_file))
        band <- str_extract_all(band,"[0-9]", simplify=TRUE)
        band <- as.numeric(paste(band, collapse=""))
        
        #2.2 For SWIR bands (1-9) calculate the radiance at sensor then reflectance. For TIR bands (10-14) calculate radiance at sensor then surface temperature 
        if(band>=1 & band<=9){
            print(paste("Processing band:", band))
            #2.2.1 Select the conversion coefficient/gain
            con_coeff_gain <- df_concoeff[band,2]
            con_coeff <- df_concoeff_table %>% filter(Band==band) %>% dplyr::select(con_coeff_gain)
            con_coeff <- con_coeff[1,1]
            
            #2.2.2 Calculate radiance at sensor
            ras_RAD <- (ras_DN-1)*con_coeff
            
            #2.2.3 Calculate TOA reflectance!
            ESUN <- df_esun[band,2]
            ras_REF <- (ras_RAD*pi*(d^2))/(ESUN*z)
            
            #2.2.4 Set "rast_dummy" to bands 1-3, we will use rast_dummy beneath for resizing.
            if(band %in% c(1,2,3) & is.null(ras_dummy)){
                ras_dummy <- ras_REF
            }
            
            #2.2.5 Resize SWIR bands to to 15m
            if(band %in% c(4,5,6,7,8,9)){
                ras_REF <- resample(x=ras_REF,y=ras_dummy,method="near")  
            }
            
            #2.2.6 Assign to a global variable
            assign(paste("band",band,"_rast",sep=""), ras_REF, envir=.GlobalEnv)
            
            #2.2.7 Assign to a variable
            if(save==TRUE){
                writeRaster(ras_REF,filename=paste(save_loc,"/",scene_id,"_band",band,".tif",sep=""))
            }
            
        } else {
            print(paste("Processing band:", band))
            #2.3.1 Select the conversion coefficients/gain
            #extract conversion coefficient
            thermal_filter <- df_thermal_table %>% filter(Band==band)
            con_coeff <- thermal_filter[[2]]
            #extract L
            L <- thermal_filter[[3]]
            #extract k1
            K1 <- thermal_filter[[4]]
            #extract k2
            K2 <- thermal_filter[[5]]
            
            #2.3.2 Calculate radiance at sensor
            ras_RAD <- (ras_DN-1)*con_coeff
            
            #2.3.3 Calculate TOA brightness temperature!
            ras_T <-  K2/(log((K1/ras_RAD)+1))
            ras_Tc <- ras_T - 273.15
            
            #2.3.4 Resize resolution from 90m to 15m
            ras_Tc <- resample(x=ras_Tc,y=ras_dummy,method="near")
            
            #2.3.5 Assign to a varaible/object in rstudio
            assign(paste("band",band,"_rast",sep=""), ras_Tc, envir=.GlobalEnv)
            
            #2.3.6 Write the file to the save location if desired.
            if(save==TRUE){
                writeRaster(ras_Tc,filename=paste(save_loc,"/",scene_id,"_band",band,".tif",sep=""))
            }
        } #end of processing thermal bands (bands 10-14)
    } #End of band loop
    
    #3.0 Create final raster stack (r)
    r <- c(band1_rast,band2_rast,band3_rast,band4_rast,band5_rast,band6_rast,band7_rast,band8_rast,band9_rast,band10_rast,band11_rast,band12_rast,band13_rast,band14_rast)
    names(r) <- c("band1","band2","band3","band4", "band5", "band6", "band7", "band8", "band9", "band10", "band11", "band12", "band13", "band14")
    r
}



folder_to_stack <- function(scene_loc){
    #read .tif files in the folder
    tif_filelist <- list.files(scene_loc,pattern=".tif$", full.names = TRUE)
    
    #isolate band numbers and order files
    band_file_list <- str_sub(tif_filelist, start=-6,end=-5)
    band_file_list <- as.numeric(unlist(regmatches(band_file_list,gregexpr("[[:digit:]]+", band_file_list))))
    tif_filelist <- tif_filelist[order(band_file_list)]
    
    #Create stack and name the bands
    stack <- rast(tif_filelist)
    names(stack) <- c("b1","b2","b3","b4", "b5", "b6", "b7", "b8", "b9", "b10", "b11", "b12", "b13", "b14")
    stack
}



create_cover_mask <- function(stack,calc_NDVI=FALSE,NDVI_value=NULL,calc_NSDI=FALSE,NSDI_value=NULL){
    #Set up dummy variables
    veg_mask <- 1
    snow_mask <- 1
    
    #1. Create individual masks
    #1.1 Vegetation mask through NDVI
    if(calc_NDVI==TRUE | calc_NSDI==TRUE){ #Sanity check for calc_NDVI=FALSE, calc_NSDI=FALSE
        if(calc_NDVI==TRUE & !is.null(NDVI_value)){ #Genearlly, >0.3 vegetation, <0.3 rocks
            NDVI <- (stack$b3-stack$b2)/(stack$b3+stack$b2) #(band3-band2)/(band3+band2) #Calculate NDVI
            vclmat <- matrix(c(-Inf, NDVI_value, 1,NDVI_value, Inf, 0), ncol=3, byrow=TRUE) #Create classification matrix
            veg_mask <- classify(NDVI, vclmat, include.lowest=TRUE) #Classify NDVI layer by matrix NDVI<cut_off=1, NDVI>cutoffvalue=0
        }
        
        #1.2 Snow mask through NSDI
        if(calc_NSDI==TRUE & !is.null(NSDI_value)){ #Genearlly, >0.4 snow, <0.4 rocks (very low to negative) sunlit rock
            NSDI <- (stack$b1-stack$b4)/(stack$b1+stack$b4) #calculate NSDI
            sclmat <- matrix(c(-Inf,NSDI_value,1, NSDI_value,Inf,0), ncol=3, byrow=TRUE) #Create classification matrix 
            snow_mask <- classify(NSDI, sclmat, include.lowest=TRUE) #Classify NSDI layer based on the index
        }
        
        #2. Combine masks
        veg_mask * snow_mask #Masks are multiplied. If one is not calculated, (e.g., calc_NDVI==TRUE, calc_NSDI==FALSE) it will the calculated one will simply be multiplied by 1
    } else {
        print("Nothing calculated. calc_NDVI or calc_NSDI must equal TRUE")
    }
}



calculate_indices <- function(stack,outliers.rm=TRUE){
    
    remove_index_outliers <- function(index){
        if(outliers.rm==TRUE){
            
            clmat <- matrix(c(-Inf,0,0), ncol=3, byrow=TRUE) #classification matrix for -Inf<x<0 == 0
            index <- classify(index, clmat, include.lowest=TRUE) #Classify index layer based on the matrix
            
            if(global(index, mean,na.rm=TRUE)[1,1]==Inf){ #if we have an infinite mean we need to get rid the few infinite values
                clmat <- matrix(c(-Inf,0,0,100,Inf,100), ncol=3, byrow=TRUE) #Values above 100 (not a real index value) will be changed to 100
                index <- classify(index, clmat, include.lowest=TRUE) #Classify index layer based on the matrix
            }
            
            index_max <- global(index, mean,na.rm=TRUE)[1,1] + global(index, sd,na.rm=TRUE)[1,1]
            clmat <- matrix(c(-Inf,0,0, index_max,Inf,index_max), ncol=3, byrow=TRUE) #Create classification matrix 
            index <- classify(index, clmat, include.lowest=TRUE) #Classify index layer based on the matrix
        }
        index
    }
    
    #VEGETATION
    #NDVI
    NDVI <<- (stack$b3-stack$b2)/(stack$b3+stack$b2)#(band3-band2)/(band3+band2)
    print("Completed NDVI")
    
    #SILICA
    #Quartz rich rocks
    index_quartz <<- stack$b14/stack$b12 #same as b13/b12
    index_quartz <<- remove_index_outliers(index_quartz)
    print("Completed index_quartz")
    
    #Silica index
    index_silica <<- (stack$b11*stack$b11)/stack$b10/stack$b12 
    index_silica <<- remove_index_outliers(index_silica)
    print("Completed index_silica")
    
    #Siliceous rocks index
    index_siliceous <<- (stack$b11*stack$b11)/(stack$b10*stack$b12)
    index_siliceous <<- remove_index_outliers(index_siliceous)
    print("Completed index_siliceous")
    
    #Opaline Silica 
    index_opal <<- (stack$b5+stack$b8)/(stack$b6+stack$b7)
    index_opal <<- remove_index_outliers(index_opal)
    print("Completed index_opal")
    
    ##CARBONATES AND MAFIC MINERALS
    #Carbonate Index
    index_carbonate <<- (stack$b13/stack$b14) #(D13/D14)
    index_carbonate <<- remove_index_outliers(index_carbonate)
    print("Completed index_carbonate")
    
    #Mafic Index
    index_mafic <<- (stack$b12/stack$b13) #(D12/D13)
    index_mafic<<- remove_index_outliers(index_mafic)
    print("Completed index_mafic")
    
    #Chlorite index (chlorite, epidote, amphibole)
    index_chlor <<- (stack$b6+stack$b9)/(stack$b7+stack$b8)
    index_chlor<<- remove_index_outliers(index_chlor)
    print("Completed index_chlor")
    
    #MgOH index (and amphiboles)
    index_MgOH <<- (stack$b6+stack$b9)/(stack$b8)
    index_MgOH <<- remove_index_outliers(index_MgOH)
    print("Completed index_MgOH")
    
    #Amphibole Index
    index_amphibole <<- stack$b6/stack$b8
    index_amphibole <<- remove_index_outliers(index_amphibole)
    print("Completed index_amphibole")
    
    
    #SILICATES
    #AlOH (Clay Idex)
    index_aloh <<- (stack$b4*stack$b7)/(stack$b6*stack$b6)
    index_aloh <<- remove_index_outliers(index_aloh)
    print("Completed index_aloh")
    
    #Kaolinite
    index_koal <<- stack$b7/stack$b5 #is poor in most cases
    index_koal <<- remove_index_outliers(index_koal)
    print("Completed index_koal")
    
    #Phyllic alteration (Sericite, muscovite, illite, smectite)
    index_phyllic <<- (stack$b5+stack$b7)/stack$b6
    index_phyllic <<- remove_index_outliers(index_phyllic)
    print("Completed index_phyllic")
    
    #Alunite
    index_alunite <<- (stack$b7*stack$b7)/(stack$b5*stack$b8)
    index_alunite <<- remove_index_outliers(index_alunite)
    print("Completed index_alunite")
    
    #Phengitic
    index_phengitic <<- stack$b5/stack$b6
    index_phengitic <<- remove_index_outliers(index_phengitic)
    print("Completed index_phengitic")
    
    #Muscovite
    index_muscovite <<- stack$b7/stack$b6
    index_muscovite <<- remove_index_outliers(index_muscovite)
    print("Completed index_muscovite")
    
    #Alteration
    index_alteration <<- stack$b4/stack$b5
    index_alteration <<- remove_index_outliers(index_alteration)
    print("Completed index_alteration")
    
    #Host Rock
    index_host_rock <<- stack$b5/stack$b6
    index_host_rock <<- remove_index_outliers(index_host_rock)
    print("Completed index_host_rock")
    
    
    #IRON BARING MINERALS
    #Ferrous silicates
    index_ferrous_silicates <<- stack$b5/stack$b4
    index_ferrous_silicates <<- remove_index_outliers(index_ferrous_silicates)
    print("Completed index_ferrous_silicates")
    
    #Gossan
    index_gossan <<- stack$b4/stack$b2
    index_gossan <<- remove_index_outliers(index_gossan)
    print("Completed index_gossan")
    
    #Ferric Iron (Fe 3+)
    index_ferric <<- stack$b2/stack$b1
    index_ferric <<- remove_index_outliers(index_ferric)
    print("Completed index_ferric")
    
    #Ferrous Iron (Fe 2+)
    index_ferrous <<- (stack$b5+stack$b3)/(stack$b1/stack$b2)
    index_ferrous <<- remove_index_outliers(index_ferrous)
    print("Completed index_ferrous")
    
    #Laterite
    index_laterite <<- stack$b4/stack$b5
    index_laterite <<- remove_index_outliers(index_laterite)
    print("Completed index_laterite")
}



mask_raster <- function(raster,cover_mask){
    raster <- raster * cover_mask # multiply the input raster by the binary mask
    clmat <- matrix(c(0,NA), ncol=2, byrow=TRUE) #Create classification matrix 
    raster <- classify(raster, clmat, include.lowest=TRUE) #Classify input layer based on the classify matrix
    raster
}



raster_pca <- function(stack){

    #Sample 1000000 points from it because it would be too much to do it on the whole stack
    sr_sample <- spatSample(stack, 100000)
    
    #Make into dataframe, give better names, and remove NA's
    sr_sample <- as.data.frame(sr_sample)
    sr_sample <- sr_sample[complete.cases(sr_sample),]
    
    #Compute PCA
    pca <<- prcomp(x=sr_sample, retx=TRUE, center=TRUE, scale. = TRUE, na.action=na.omit)
    
    #little function that only computes pca on the first three columns
    pca_predict3 <- function(model, data, ...) {
        predict(model, data, na.rm=TRUE, ...)[,1:3]
    }
    
    predict(stack, pca, fun=pca_predict3, na.rm=FALSE)
}



assemble_ASTER_lib <- function(ASTER_lib_folder){ #This function assembles an ASTER Library dataframe that can be plotted or used.
    lib_spectra_files <- list.files(ASTER_lib_folder,full.names=TRUE)
    lib_spectra <- NULL
    for(i in lib_spectra_files){
        individual_spectra <- read.table(i, sep=",", stringsAsFactors=FALSE, header=TRUE)
        min_name <- colnames(individual_spectra)
        min_name_dots <- gregexpr("\\.",min_name)
        name_start <- min_name_dots[[1]][4]+1
        name_stop <- min_name_dots[[1]][5]-1
        individual_spectra$Mineral <- substr(min_name,name_start,name_stop)
        individual_spectra$Band <- seq_len(nrow(individual_spectra))
        individual_spectra$Wavelength <- c(0.56,0.66,0.82,1.65,2.165,2.205,2.26,2.33,2.395)
        colnames(individual_spectra)[1] <- "Reflectance"
        individual_spectra$Reflectance[which(individual_spectra$Reflectance==-1.23e+34)] <- 0.00001 #Fix null values
        
        if(i==lib_spectra_files[1]){
            individual_spectra$Ref_Offset <- individual_spectra$Reflectance
        } else {
            max_reflect_prev <- max(lib_spectra$Ref_Offset)
            min_reflect_current <- min(individual_spectra$Reflectance)
            individual_spectra$Ref_Offset <- individual_spectra$Reflectance - (min_reflect_current-max_reflect_prev) 
        }
        lib_spectra <- rbind(lib_spectra,individual_spectra)
    }
    return(lib_spectra)
}



plot_ASTER_lib <- function(lib_obj,offset=FALSE){
    if(offset==TRUE){
        ggplot(data=lib_obj, mapping=aes(x=Wavelength,y=Ref_Offset,color=Mineral)) + 
            geom_line() + theme_minimal() + 
            xlab("Wavelength (um)") + 
            ylab("Relative Reflectance")
    } else {
        ggplot(data=lib_obj, mapping=aes(x=Wavelength,y=Reflectance,color=Mineral)) + 
            geom_line() + 
            theme_minimal() + 
            xlab("Wavelength (um)") +
            ylab("Reflectance")
    }
}



SAM <- function(lib_obj, input_scene_spectra){
    SAM_class <- NULL
    for (r in unique(lib_obj$Mineral)){
        lib_temp <- lib_obj[which(lib_obj$Mineral==r),]
        lib_temp <- cbind(lib_temp,input_scene_spectra)
        lib_temp$t2 <- lib_temp$Reflectance^2
        lib_temp$r2 <- lib_temp$input_scene_spectra^2
        lib_temp$tr <- lib_temp$Reflectance * lib_temp$input_scene_spectra
        
        a <- acos((sum(lib_temp$tr))/(sqrt(sum(lib_temp$t2))*sqrt(sum(lib_temp$r2))))
        
        SAM_class <- rbind(SAM_class, data.frame(Mineral=unique(lib_temp$Mineral), Contribution=a))
        print(SAM_class)
    }
    list(Contributions=SAM_class, Classification=SAM_class[which.min(SAM_class$Contribution),1])
}