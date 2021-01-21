#----------------------------------------------------
#  R SCRIPT FOR MAPPING NATIONAL SOIL SALINITY 
#----------------------------------------------------
#
# Author: Omuto, CT (cthine@gmail.com)
#
# ---------------------------------------------------------------------------
# Introduction: 
# This script is intended to support development of country-level
# soil salinity map using the following input data:
# (1) georeferenced EC,pH,ESP as input soil factors (between #0-100 cm deep).
# (2) climate raster map like mean annual rainfall amounts and/or temperature
# (3) parent material raster map such as geology map or soil map
# (4) reflectance from remote sensing image bands such as MODIS (7 bands of MOD9A1)
# (4) altitude/relief parameters such as DEM, slope, length of slope (ls factor),
#     longitudinal curvature, valley depth, and channel network to basin level
# (5) land cover/use raster map
#-----------------------------------------------------------------------------
#
# IMPORTANT CONSIDERATIONS FOR THE INPUT DATA
# (1) The soil data should have the following columns: Profile number, Longitudes,
#     Latitudes, upper soil depth and lower soil depth for each profile,EC, pH, ESP
# (2) Climate, geology, relief parameters, land cover, and remote sensing images MUST
#     have uniform projection, pixel sizes, coordinate system, and spatial extent (and should overlay and confirmed)
#     in software such as ILWIS
# -----------------------------------------------------------------------------

###### PART ONE: IMPORTING DATA AND DATA PREPARATION #####
#Step 1: Load the libraries and set your working directory
#The libraries may need to be installed if they are missing
#Otherwise load the following packages if R is started afresh

{
install.packages(c("raster", "sp", "rgdal", "car", "carData", "dplyr", "spacetime", "gstat", "automap", "randomForest", "fitdistrplus", "e1071", "caret", "soilassessment", "soiltexture", "GSIF", "aqp", "plyr", "Hmisc", "corrplot", "factoextra", "spup", "purrr", "lattice", "ncf", "qrnn", "nnet", "mda", "RColorBrewer", "vcd", "readxls","maptools","neuralnet","psych","ggrepel", "plotly")) 
}#Use this script to install the packages if using R for the first time. Internet connectivity is required

### step 1: Load the libraries #######--------------------------------------------------

library(sp); library(foreign); library(rgdal); library(car);library(carData); library(maptools)
library(spacetime); library(gstat); library(automap);library(randomForest);library(fitdistrplus);
library(e1071); library(caret); library(raster); library(soilassessment); library(soiltexture);
library(GSIF); library(aqp); library(plyr); library(Hmisc); library(corrplot); library(factoextra)
library(spup); library(purrr); library(lattice);library(ncf);library(npsurvSS); 
library(nnet); library(class); library(mda); library(RColorBrewer); library(vcd); library(grid); 
library(neuralnet);library(readxl); library(psych);library(qrnn); library(dplyr)

#### Load test data: This is for practice with test sample

### For test data, you do not need Step 2 to Step 3 (data import). Start from Step 4

#####Step 2: Import the soil data and predictors ####
## Import soil data and inspect soil depth ranges

soil=read.csv("soildata.csv",header = T)
#soil = read_xlsx ("soildata.xlsx")

### Import predictors

{
  predictors=readGDAL("dems.tif");
  predictors$slope=readGDAL("SLOPEs.mpr")$band1;
  predictors$loncurve=readGDAL("loncurves.mpr")$band1;
  predictors$ls=readGDAL("lss.mpr")$band1;
  predictors$cnbl=readGDAL("cnbls.mpr")$band1;
  predictors$valley=readGDAL("valleys.mpr")$band1;
  predictors$lcover=readGDAL("lcover.mpr")$band1;
  predictors$geology=readGDAL("mintemp.mpr")$band1;
  predictors$pgeology=readGDAL("maxtemp.mpr")$band1;
  predictors$rain=readGDAL("rain.mpr")$band1;
  predictors$swir1=readGDAL("BBand6.mpr")$band1;
  predictors$swir2=readGDAL("BBand7.mpr")$band1;
  predictors$BBlue=readGDAL("BBand3.mpr")$band1;
  predictors$BGreen=readGDAL("BBand4.mpr")$band1;
  predictors$BRed=readGDAL("BBand1.mpr")$band1;
  predictors$BIRed=readGDAL("BBand2.mpr")$band1;
  predictors$soilmap=readGDAL("soilmap.mpr")$band1;
  
}
predictors$dem=predictors$band1
predictors$band1=NULL

###### PREDICTOR DATA PREPARATION #####
#Step 4: Assess the distribution of the predictors and transform where necessary 
#Check and remove NA if available in any layer

summary(predictors)

#Step 5: Derive the remote sensing indices of salinity
predictors$SI1=imageIndices(predictors$BBlue,predictors$BGreen,predictors$BRed,predictors$BIRed,predictors$swir1,predictors$swir2,"SI1");summary(predictors$SI1)
predictors$SI2=imageIndices(predictors$BBlue,predictors$BGreen,predictors$BRed,predictors$BIRed,predictors$swir1,predictors$swir2,"SI2");summary(predictors$SI2)
predictors$SI3=imageIndices(predictors$BBlue,predictors$BGreen,predictors$BRed,predictors$BIRed,predictors$swir1,predictors$swir2,"SI3");summary(predictors$SI3)
predictors$SI4=imageIndices(predictors$BBlue,predictors$BGreen,predictors$BRed,predictors$BIRed,predictors$swir1,predictors$swir2,"SI4");summary(predictors$SI4)
predictors$SI5=imageIndices(predictors$BBlue,predictors$BGreen,predictors$BRed,predictors$BIRed,predictors$swir1,predictors$swir2,"SI5");summary(predictors$SI5)
predictors$SI6=imageIndices(predictors$BBlue,predictors$BGreen,predictors$BRed,predictors$BIRed,predictors$swir1,predictors$swir2,"SI6");summary(predictors$SI6)
predictors$SAVI=imageIndices(predictors$BBlue,predictors$BGreen,predictors$BRed,predictors$BIRed,predictors$swir1,predictors$swir2,"SAVI");summary(predictors$SAVI)
predictors$VSSI=imageIndices(predictors$BBlue,predictors$BGreen,predictors$BRed,predictors$BIRed,predictors$swir1,predictors$swir2,"VSSI");summary(predictors$VSSI)
predictors$NDSI=imageIndices(predictors$BBlue,predictors$BGreen,predictors$BRed,predictors$BIRed,predictors$swir1,predictors$swir2,"NDSI");summary(predictors$NDSI)
predictors$NDVI=imageIndices(predictors$BBlue,predictors$BGreen,predictors$BRed,predictors$BIRed,predictors$swir1,predictors$swir2,"NDVI");summary(predictors$NDVI)
predictors$SR=imageIndices(predictors$BBlue,predictors$BGreen,predictors$BRed,predictors$BIRed,predictors$swir1,predictors$swir2,"SR");summary(predictors$SR)
predictors$CRSI=imageIndices(predictors$BBlue,predictors$BGreen,predictors$BRed,predictors$BIRed,predictors$swir1,predictors$swir2,"CRSI");summary(predictors$CRSI)
predictors$BI=imageIndices(predictors$BBlue,predictors$BGreen,predictors$BRed,predictors$BIRed,predictors$swir1,predictors$swir2,"BI");summary(predictors$BI)

##Check the distribution of the image indices
{#: SI1, SI2, SI3, SI4, SI5, SI6, SAVI, VSSI, NDSI, NDVI, SR, CRSI, BI
  
hist(predictors@data[,"SI2"],cex=1)
summary(predictors$SI1)
predictors$SI1=sqrt(predictors$SI1)
predictors$SI4=ifelse((predictors$SI4)<0,-1*sqrt(abs(predictors$SI4)),sqrt(predictors$SI4))
hist(predictors$SI1)
summary(predictors$SI1)

}

##Step 6: Convert the predictors to dataframe and perform PCA,
predicters=predictors@data[,c("SI1","SI2","SI3","SI4","SI5","SI6","SAVI","VSSI","NDSI","NDVI","SR","CRSI", "BI")]

##Step 7: Check and plot correlation between predictors and perform PCA analysis 
soil.cor=cor(predicters)
corrplot(soil.cor,method="number",number.cex = 0.7)

pca<-prcomp(predicters[], scale=TRUE)
fviz_eig(pca) #Plot the predictor importance to know the number of PCs to include in the list of predictors
summary(pca)
fviz_pca_var(pca, col.var = "contrib",gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE)

#Check for variable contribution
{
  var_coord_func <- function(loadings, comp.sdev){loadings*comp.sdev}
  loadings <- pca$rotation;sdev <- pca$sdev
  var.coord <- t(apply(loadings, 2, var_coord_func, sdev)) 
  var.cos2 <- var.coord^2;  comp.cos2 <- apply(var.cos2, 2, sum)
  contrib <- function(var.cos2, comp.cos2){var.cos2*100/comp.cos2}
  var.contrib <- t(apply(var.cos2,1, contrib, comp.cos2))
  head(var.contrib[, 1:11])
}
var.pca=as.data.frame(var.contrib[, 1:11])
write.csv(var.pca,file="pca_variable_contribution.csv")

#Step 8: Convert the PCs to grid maps and attach their corresponding values to the soil dataset
Pred.pcs<-predict(pca,predicters[])
predictors@data$PCA1=Pred.pcs[,1]
predictors@data$PCA2=Pred.pcs[,2]
predictors@data$PCA3=Pred.pcs[,3]
predictors@data$PCA4=Pred.pcs[,4]
predictors@data$PCA5=Pred.pcs[,5]
predictors@data$PCA6=Pred.pcs[,6]
predictors@data$PCA7=Pred.pcs[,7]
predictors@data$PCA8=Pred.pcs[,8]

#Step 9: The following data can now be removed to free up computer space
rm(pca);rm(Pred.pcs);rm(predicters); rm(soil.cor); rm(var.pca)
rm(loadings); rm(var.contrib); rm(var.coord); rm(var.cos2); rm(comp.cos2);rm(sdev);rm(contrib);rm(var_coord_func)
predictors@data[,c("swir1","swir2","BBlue","BGreen","BRed","BIRed","NSI","SI1","SI2","SI3","SI4","SI5","SI6","SR","NDSI","NDVI","SAVI","VSSI","CRSI","BI")]=NULL

#Save the harmonized predictors for future use

save(predictors, file="predictors.RData")

#### PART TWO: DIGITAL MAPPING OF INDICATORS ######
# Part A: Mapping EC   #######
###Step 1a: Harmonize the depths for EC
## If EC was not determined from saturated paste extract, then, it needs harmonization
## Harmonization may require soil texture data
## If EC is from saturated paste extract, then jump to Step 1a-3 

## Step 1a-1: First harmonize the texture classes
str(soil) # Check and note the columns containing the texture components

{
variable.names(soil[11]); variable.names(soil[12]); variable.names(soil[13]); 
soil$dummy= rowSums(soil[, 11:13])
soil1=subset(soil,!is.na(soil$dummy))
soil1$dummy=NULL # remove the dummy
soil0=data.frame(soil1)
SSCP=soil0[,c("Sand","Silt","Clay")]
names(SSCP) <- c('SAND', 'SILT', 'CLAY')
SSCP <- round(SSCP, 2)
SSCP_norm <- TT.normalise.sum(tri.data = SSCP[,1:3], residuals = T)
colnames(SSCP_norm)[1:3] <- paste0(colnames(SSCP_norm)[1:3],"_nm")
SSCP <- cbind(SSCP, round(SSCP_norm, 2))
SSCP$CLAY=SSCP$CLAY_nm;SSCP$SILT=SSCP$SILT_nm;SSCP$SAND=SSCP$SAND_nm
rm(SSCP_norm)
soil0=cbind(soil0,"TEXCLASS" =TT.points.in.classes(tri.data =SSCP[, c('SAND', 'SILT',  'CLAY')],class.sys = "USDA.TT", PiC.type  = "t",collapse  = ', '))
soil0$TEXCLASS=as.factor(soil0$TEXCLASS)
soil0$TEXCLASS1=as.numeric(soil0$TEXCLASS)
summary(soil0$TEXCLASS)
rm(SSCP)
#View(classnames("texture")) # In order to get the classess and reclassify using numeric values
soil0$TEXCLASS=car::recode(soil0$TEXCLASS,"'Cl, ClLo'='ClLo'") #Here, the double classes are change one at a time
summary(soil0$TEXCLASS)
soil0$TEXCLASS1=dplyr::recode(soil0$TEXCLASS,Cl=1, ClLo =7,Lo=11,LoSa=10,Sa=12,SaCl=8,SaClLo=9,SaLo=5,SiCl=2,SiClLo=3,SiLo=4,Si=6,CS=13,MS=14,HCL=16,FS=15)
summary(soil0$TEXCLASS1)
soil1=soil0
soil1=subset(soil1,!is.na(soil1$TEXCLASS1))

#This part may be used to develop pedotransfer model for a target variable from other variables
#For example, to predict EC from cations: Mg, Na, K, Ca  
soil33=subset(soil1,!is.na(soil1$EC))
soil33=soil33[,c("EC","Mg","Na","K","Ca")]
soil33=subset(soil33,!is.na(soil33$Mg))
EC.pf=pedoTrasnfer("randomforest",soil33,EC,Mg, Na,K,Ca)
soil33$EC1=predict(EC.pf,soil33)
cor(soil33$EC,soil33$EC1)^2
soil1$EC1=ifelse(is.na(soil1$EC),predict(EC.pf,soil1),soil1$EC)
summary(soil1$EC1)

### Step 1a-2: Harmonize the EC values
soil1$ECse=ECconversion1(soil1$EC1,soil1$OC,soil1$CLAY,soil1$TEXCLASS1,"1:5","FAO") 

summary(soil1$ECse); summary(soil1$EC)
soil1=subset(soil1,!is.na(soil1$ECse))

}

### Step 1a-3: Now harmonize the depths
soil1 =subset(soil,!is.na(soil$EC))

lon=soil1$Longitude
lat=soil1$Latitude
id=soil1$Pits
top=soil1$Upper
bottom=soil1$Lower
horizon=soil1$Horizon
ECdp=soil1$EC #(or use ECdp=soil1$EC if data is from saturated paste extract)
prof1=join(data.frame(id,top,bottom, ECdp, horizon),data.frame(id,lon,lat),type="inner")
depths(prof1)=id~top+bottom
site(prof1)=~lon+lat
coordinates(prof1) = ~lon+lat
proj4string(prof1)=CRS("+proj=longlat +datum=WGS84 +no_defs")
depth.s = mpspline(prof1, var.name= "ECdp", lam=0.21,d = t(c(0,30,100,150)))
plot(prof1, color= "ECdp", name="horizon",color.palette = rev(brewer.pal(8, 'Accent')),par=c(cex.lab=2.0)) #Figure 4.5

#Step 2a: Convert the harmonized depths into spatial dataframe and harmonize the projection
soilhrmdepths=data.frame(depth.s$idcol, depth.s$var.std, check.names = TRUE)
soil2a=merge(soil1,soilhrmdepths,by=intersect(names(soil1),names(soilhrmdepths)),  by.x="Pits",by.y="depth.s.idcol",all=TRUE)
coordinates(soil2a)=~Longitude+Latitude
proj4string(soil2a)=CRS("+proj=longlat +datum=WGS84 +no_defs")

#Harmonize CRS. #Make sure to use correct +proj and +zone
soil1=spTransform(soil2a,CRS("+proj=utm +zone=36 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")) #Make sureto use correct CRS

#Step 3a: Further harmonize EC values to ECse and its statistical distribution
hist(soil1$X0.30.cm)
summary(soil1$X0.30.cm)
soil1=subset(soil1,!is.na(soil1$X0.30.cm))
bubble(soil1,"X0.30.cm", main="Harmonized Electrical Conductivity (0-30 cm)")

soil1$dummy=(soil1$X0.30.cm)+0.001# add "+0.001" if minimum in summary(soil1$X0.30.cm) is zero
hist(soil1$dummy, main="Frequency distribution (before transformation)", xlab="Harmonized EC (dS/m)")
soil1$Tran=(soil1$dummy^(as.numeric(car::powerTransform(soil1$dummy, family ="bcPower")["lambda"]))-1)/(as.numeric(car::powerTransform(soil1$dummy, family ="bcPower")["lambda"]))
hist(soil1$Tran, main="Frequency distribution (after transformation)",xlab="Harmonized EC (dS/m)")
crs(predictors); crs(soil1); #soil1$Tran=soil1$X0.30.cm

#Step 4a: Attach the spatial predictors
{predictors.ov=over(soil1, predictors)
  soil1$dem=predictors.ov$dem
  soil1$slope=predictors.ov$slope
  soil1$cnbl=predictors.ov$cnbl
  soil1$ls=predictors.ov$ls
  soil1$valley=predictors.ov$valley
  soil1$loncurve=predictors.ov$loncurve
  soil1$lcover=predictors.ov$lcover
  soil1$rain=predictors.ov$rain
  soil1$pgeology=predictors.ov$pgeology
  soil1$geology=predictors.ov$geology
  soil1$PCA1=predictors.ov$PCA1
  soil1$PCA2=predictors.ov$PCA2
  soil1$PCA3=predictors.ov$PCA3
  soil1$PCA4=predictors.ov$PCA4
  soil1$PCA5=predictors.ov$PCA5
}

#Step 5a: Choosing the appropriate model for predictive mapping of EC 
summary(soil1); 
soil1=subset(soil1,!is.na(soil1$dem))
soil11a=soil1@data[,c("Tran","dem","slope","ls","loncurve","cnbl","valley","rain","lcover","geology","PCA1","PCA2","PCA3","PCA4","PCA5")]

{
#The next line will search for suitable DSM model for the data. It may take time to complete
regmodelSuit(soil11a,Tran,dem,geology,loncurve,slope,rain,cnbl,valley,lcover,ls,PCA1,PCA2,PCA3, PCA4,PCA5)

# Note the best model produced by regmodelSuit above

{soil4=as.data.frame(soil1)
  bound <- floor((nrow(soil4)/4)*3)
  soil3 <- soil4[sample(nrow(soil4)), ]
  df.traina <- soil3[1:bound, ]
  df.testa <- soil3[(bound+1):nrow(soil3), ]}

# Insert the model xxxx noted from regmodelSuit into the section for method = "xxxx" in the next line for rf.ec
rf.ec=train(Tran~(slope+ls+geology+cnbl+rain+loncurve+valley+lcover+dem+PCA1+PCA2+PCA3+PCA4+PCA5),  data = df.traina,  method = "qrf", trControl=trainControl( method = "cv",number=5,returnResamp = "all",savePredictions = TRUE, search = "random",verboseIter = FALSE))

print(rf.ec)

# Show the prediction interval
df.testa$Strain=predict(rf.ec,newdata=df.testa)
hist(df.testa$Strain,xlab="Box-Cox Transformed ECse (0-30cm)", main=NULL)
abline(v = quantile(df.testa$Strain, probs = c(0.05, 0.95)),lty = 5, col = "red")

lmbda=(as.numeric(powerTransform(soil1$dummy, family ="bcPower")["lambda"]))
df.testa$Strain=predict(rf.ec,newdata=df.testa)
df.testa$Strain=(df.testa$Strain*lmbda+1)^(1/lmbda)
summary(df.testa$Strain)
summary(df.testa$dummy)

cor(df.testa$Strain,df.testa$dummy)^2
{plot(df.testa$Strain~df.testa$dummy, cex.lab=1.5, cex.axis=1.5,cex=1.5,xlim=c(0,117), ylim=c(0,117), xlab="Harmonized Electrical Conductivity",ylab="Modelled Electrical Conductivity", main="Accuracy assessment on hold-out samples")
  abline(a=0,b=1,lty=20, col="blue")}

Bias=mean(df.testa$Strain-df.testa$dummy,na.rm=TRUE)
RMSE=sqrt(sum((df.testa$Strain-df.testa$dummy)^2,na.rm=TRUE)/length((df.testa$Strain-df.testa$dummy)))
Rsquared=cor(df.testa$Strain,df.testa$dummy)^2
NSE=1-sum((df.testa$Strain-df.testa$dummy)^2,na.rm=TRUE)/sum((df.testa$Strain-mean(df.testa$dummy,na.rm=TRUE))^2,na.rm=TRUE)
statia=data.frame(Bias,RMSE,Rsquared,NSE);View(statia)
write.csv(statia,file = "EC_30_validmodel_stats.csv")
}

##Step 6a: This part is optional. Plough back all the data and make aggregate model for mapping EC
{
soil5=as.data.frame(soil1)
rf.ec=train(Tran~(dem+slope+ls+cnbl+rain+loncurve+valley+lcover+geology+PCA1+PCA2+PCA3+PCA4+PCA5),  data = soil5,  method = "qrf",trControl=trainControl( method = "cv",number=5,returnResamp = "all",savePredictions = TRUE, search = "random",verboseIter = FALSE))

}

#Step 7a: Determine the final EC map
lmbda1=(as.numeric(powerTransform(soil1$dummy, family ="bcPower")["lambda"]))
predictors$ECte=predict(rf.ec,predictors)
predictors$ECse=(predictors$ECte*lmbda1+1)^(1/lmbda1)

#Step 8a: Compare the predicted map and validation dataset
featureRep(predictors["ECse"],soil1)
summary(predictors$ECse);summary(df.testa$dummy)


#Step 9a: Export the final map
plot(predictors["ECse"]) #This part may be run to displaye the map
writeRaster(raster(predictors["ECse"]), drivername = "GTiff", "Top0_30ECs.tif", overwrite=F)
writeRaster(raster(predictors["ECte"]), drivername = "GTiff", "Top0_30ECte.tif",overwrite=F)

#Step10a: Determine uncertainty for EC
soil6a=soil1[,c("Tran")]
predictors6a=predictors[c("dem","slope","rain","cnbl","lcover","loncurve","pgeology","geology","ls","valley","PCA1","PCA2","PCA3","PCA4","PCA5")]

#This part runs simulations and can take time
pred_uncerta=predUncertain(soil6a,predictors6a,4,95,"qrandomforest")
summary(pred_uncerta$pred_width)
spplot(pred_uncerta, "pred_width")

EC0_30_uncertain=(pred_uncerta$pred_width)
summary(EC0_30_uncertain)
writeRaster(EC0_30_uncertain, filename="EC0_30_uncertain.tif",format="GTiff",overwrite=F)
writeRaster(pred_uncerta$pred_sd, filename="EC0_30_sd.tif",format="GTiff",overwrite=F)

#########Part B: PH Mapping######
###Step 1b: Harmonize the depths for pH
soil11=subset(soil,!is.na(soil$pH))
summary(soil11$pH)

lon2=soil11$Longitude
lat2=soil11$Latitude
id2=soil11$Pits
top2=soil11$Upper
bottom2=soil11$Lower
horizon2=soil11$Horizon
PHdp=soil11$pH
prof2=join(data.frame(id2,top2,bottom2, PHdp, horizon2),data.frame(id2,lon2,lat2),type="inner")
depths(prof2)=id2~top2+bottom2
site(prof2)=~lon2+lat2
coordinates(prof2) = ~lon2+lat2
proj4string(prof2)=CRS("+proj=longlat +datum=WGS84 +no_defs")
depth.s2 <- mpspline(prof2, var.name= "PHdp", lam=0.21, d = t(c(0,30,100,150)))
plot(prof2, color= "PHdp", name="horizon",color.palette = rev(brewer.pal(20, 'Spectral'))) #Figure 4.5

#Step 2b: Convert the harmonized depths into spatial dataframe and harmonize the projection
soilhrmdepths2=data.frame(depth.s2$idcol, depth.s2$var.std, check.names = TRUE)
soil21=merge(soil11,soilhrmdepths2,by=intersect(names(soil11),names(soilhrmdepths2)),  by.x="Pits",by.y="depth.s2.idcol",all=TRUE)
coordinates(soil21)=~Longitude+Latitude
proj4string(soil21)=CRS("+proj=longlat +datum=WGS84")
soil11=spTransform(soil21,CRS("+proj=utm +zone=36 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

#Step 3b: Further harmonize pH statistical distribution
soil11=subset(soil11,!is.na(soil11$X0.30.cm))
hist(soil11$X0.30.cm)
summary(soil11$X0.30.cm);summary(soil11$pH)
bubble(soil11,"X0.30.cm",main="Harmonized pH (0-30 cm)")

soil11$dummy=(soil11$X0.30.cm) 
soil11$Tran=(soil11$dummy^(as.numeric(powerTransform(soil11$dummy, family ="bcPower")["lambda"]))-1)/(as.numeric(powerTransform(soil11$dummy, family ="bcPower")["lambda"]))
hist(soil11$Tran)


#Step 4b: Extract and attach pixel values of the spatial predictors
{predictors.ov1=over(soil11, predictors)
  soil11$dem=predictors.ov1$dem
  soil11$slope=predictors.ov1$slope
  soil11$cnbl=predictors.ov1$cnbl
  soil11$ls=predictors.ov1$ls
  soil11$loncurve=predictors.ov1$loncurve
  soil11$valley=predictors.ov1$valley
  soil11$lcover=predictors.ov1$lcover
  soil11$rain=predictors.ov1$rain
  soil11$pgeology=predictors.ov1$pgeology
  soil11$geology=predictors.ov1$geology
  soil11$PCA1=predictors.ov1$PCA1
  soil11$PCA2=predictors.ov1$PCA2
  soil11$PCA3=predictors.ov1$PCA3
  soil11$PCA4=predictors.ov1$PCA4
  soil11$PCA5=predictors.ov1$PCA5
}

#Step 5b: Choosing the appropriate model for predictive mapping of pH
summary(soil11)
soil11=subset(soil11,!is.na(soil11$dem))
soil11b=soil11@data[,c("Tran","dem","slope","ls","cnbl","loncurve","valley","rain","lcover","pgeology","geology","PCA1","PCA2","PCA3","PCA4","PCA5")] 

#This part searches for suitable DSM model. It can take a long time
regmodelSuit(soil11b,Tran,dem,geology,pgeology,slope,rain,loncurve,cnbl,valley,lcover,ls,PCA1,PCA2,PCA3,PCA4,PCA5)

{soil4b=as.data.frame(soil11)
  bound <- floor((nrow(soil4b)/4)*3)
  soil3b <- soil4b[sample(nrow(soil4b)), ]
  df.trainb <- soil3b[1:bound, ]
  df.testb <- soil3b[(bound+1):nrow(soil3b), ]}

rf.ph=train(Tran~(slope+loncurve+rain+cnbl+ls+valley+geology+pgeology+lcover+dem+PCA1+PCA2+PCA3+PCA4+PCA5),  data = df.trainb,  method = "qrf",trControl=trainControl( method = "cv",number=5,returnResamp = "all",savePredictions = TRUE, search = "random",verboseIter = FALSE))

print(rf.ph)

# Show the prediction interval
df.testb$Strain=predict(rf.ph,newdata=df.testb)
hist(df.testb$Strain,xlab="Box-Cox Transformed PH (0-30cm)", main=NULL)
abline(v = quantile(df.testb$Strain, probs = c(0.05, 0.95)),lty = 5, col = "red")

lmbda2=(as.numeric(powerTransform(soil11$dummy, family ="bcPower")["lambda"]))
df.testb$Strain=predict(rf.ph,newdata=df.testb)
df.testb$Strain=(df.testb$Strain*lmbda2+1)^(1/lmbda2)
summary(df.testb$Strain)
summary(df.testb$dummy)

cor(df.testb$Strain,df.testb$dummy)^2
{plot(df.testb$Strain~df.testb$dummy,xlim=c(2,8),ylim=c(2,8) ,xlab="Measured PH",ylab="Modelled PH", main="Accuracy assessment for pH")
  abline(a=0,b=1,lty=20, col="blue")}

Biasb=mean(df.testb$Strain-df.testb$dummy,na.rm=TRUE)
RMSEb=sqrt(sum((df.testb$Strain-df.testb$dummy)^2,na.rm=TRUE)/length((df.testb$Strain-df.testb$dummy)))
Rsquaredb=cor(df.testb$Strain,df.testb$dummy)^2
NSEb=1-sum((df.testb$Strain-df.testb$dummy)^2,na.rm=TRUE)/sum((df.testb$Strain-mean(df.testb$dummy,na.rm=TRUE))^2,na.rm=TRUE)
statib=data.frame(Biasb,RMSEb,Rsquaredb,NSEb);View(statib)
write.csv(statib,file = "PH0_30_validmodel_stats.csv")

##Step 6b: This part is optional: Plough back all the data and make aggregate model for mapping pH
{
soil5b=as.data.frame(soil11)
rf.ph=train(Tran~(slope+loncurve+cnbl+rain+ls+geology+pgeology+valley+lcover+dem+PCA1+PCA2+PCA3+PCA4+PCA5),  data = soil5b,  method = "qrf",trControl=trainControl( method = "cv",number=5,returnResamp = "all",savePredictions = TRUE, search = "random",verboseIter = FALSE))

}


##Step 7b: Determine the final pH map
lmbda2=(as.numeric(powerTransform(soil11$dummy, family ="bcPower")["lambda"]))
predictors$PHt=predict(rf.ph,predictors)
predictors$PH=(predictors$PHt*lmbda2+1)^(1/lmbda2)
summary(predictors$PH)

#Step 8b: Compare the predicted map and validation dataset
featureRep(predictors["PH"],soil11)
summary(predictors$PH);summary(df.testb$dummy)


#Step 9b: Export the final map
hist(predictors$PHt)
writeRaster(raster(predictors["PH"]), drivername = "GTiff", "Top0_30PH.tif",overwrite=F)
writeRaster(raster(predictors["PHt"]), drivername = "GTiff", "Top0_30PHt.tif",overwrite=F)


##Step10b: Determine location uncertainty for pH 
soil6b=soil11[,c("Tran")]
predictors6b=predictors[c("dem","rain","cnbl","slope","lcover","loncurve","pgeology","geology","ls","valley","PCA1","PCA2","PCA3","PCA4","PCA5")]

#This part uses simulations and can take time
pred_uncertb=predUncertain(soil6b,predictors6b,4,95,"qrandomforest")
spplot(pred_uncertb, "pred_width", scales = list(draw = TRUE), sp.layout = list(list("sp.points", soil6b, pch = 3,cex=0.4,col="cyan")))


ph0_30_uncertain=(pred_uncertb$pred_width)
summary(ph0_30_uncertain)
writeRaster(ph0_30_uncertain, filename="pH0_30_uncertain.tif",format="GTiff",overwrite=F)
writeRaster(pred_uncertb$pred_sd, filename="pH0_30_sd.tif",format="GTiff",overwrite=F)


#########Part C: ESP Mapping####
#########Step 1c: Harmonize the depths
soil12=subset(soil,!is.na(soil$ESP)); #soil12$ESP=ifelse(soil12$ESP>100,0.01*soil12$ESP,soil12$ESP)

lon3=soil12$Longitude
lat3=soil12$Latitude
id3=soil12$Pits
top3=soil12$Upper
bottom3=soil12$Lower
horizon3=soil12$Horizon
ESPdp=soil12$ESP
prof3=join(data.frame(id3,top3,bottom3, ESPdp, horizon3),data.frame(id3,lon3,lat3),type="inner")
depths(prof3)=id3~top3+bottom3
site(prof3)=~lon3+lat3
coordinates(prof3) = ~lon3+lat3
proj4string(prof3)=CRS("+proj=longlat +datum=WGS84 +no_defs")
depth.s3 <- mpspline(prof3, var.name= "ESPdp", lam=0.21, d = t(c(0,30,100,150)))
plot(prof3, color= "ESPdp", name="horizon",color.palette = rev(brewer.pal(20, 'Spectral'))) 

#Step 2c: Convert the harmonized depths into spatial dataframe and harmonize the projection
soilhrmdepths3=data.frame(depth.s3$idcol, depth.s3$var.std, check.names = TRUE)
soil21=merge(soil12,soilhrmdepths3,by=intersect(names(soil12),names(soilhrmdepths3)),  by.x="Pits",by.y="depth.s3.idcol",all=TRUE)
coordinates(soil21)=~Longitude+Latitude
proj4string(soil21)=CRS("+proj=longlat +datum=WGS84")
soil13=spTransform(soil21,CRS("+proj=utm +zone=36 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

#Step 3c: Further harmonize ESP statistical distribution
soil13=subset(soil13,!is.na(soil13$X0.30.cm))
hist(soil13$X0.30.cm)
summary(soil13$X0.30.cm);summary(soil13$ESP)
bubble(soil13,"X0.30.cm",main="Harmonized Exchangeable Sodium Percent (0-30 cm)")

soil13$dummy=(soil13$X0.30.cm)+0.001 #A small value +0.001 is added to remove zero discontinuity

soil13$Tran=(soil13$dummy^(as.numeric(powerTransform(soil13$dummy, family ="bcPower")["lambda"]))-1)/(as.numeric(powerTransform(soil13$dummy, family ="bcPower")["lambda"]))
hist(soil13$Tran)


#Step 4c: Attach the spatial predictors
{predictors.ov3=over(soil13, predictors)
  soil13$dem=predictors.ov3$dem
  soil13$slope=predictors.ov3$slope
  soil13$cnbl=predictors.ov3$cnbl
  soil13$ls=predictors.ov3$ls
  soil13$valley=predictors.ov3$valley
  soil13$loncurve=predictors.ov3$loncurve
  soil13$lcover=predictors.ov3$lcover
  soil13$rain=predictors.ov3$rain
  soil13$pgeology=predictors.ov3$pgeology
  soil13$geology=predictors.ov3$geology
  soil13$PCA1=predictors.ov3$PCA1
  soil13$PCA2=predictors.ov3$PCA2
  soil13$PCA3=predictors.ov3$PCA3
  soil13$PCA4=predictors.ov3$PCA4
  soil13$PCA5=predictors.ov3$PCA5

  }


#Step 5c: Choosing the appropriate model for predictive mapping of ESP and validating it
summary(soil13)
soil13=subset(soil13,!is.na(soil13$dem))
soil11c=soil13@data[,c("Tran","dem","slope","ls","loncurve","valley","lcover","pgeology","geology","PCA1","PCA2","PCA3","PCA4","PCA5")]
regmodelSuit(soil11c,Tran,dem,geology,pgeology,slope,loncurve,valley,lcover,ls,PCA1,PCA2,PCA3, PCA4,PCA5)

{soil4c=as.data.frame(soil13)
  bound <- floor((nrow(soil4c)/4)*3)
  soil3c <- soil4c[sample(nrow(soil4c)), ]
  df.trainc <- soil3c[1:bound, ]
  df.testc <- soil3c[(bound+1):nrow(soil3c), ]}

rf.esp=train(Tran~(slope+loncurve+ls+geology+pgeology+valley+lcover+dem+PCA1+PCA2+PCA3+PCA4+PCA5),  data = df.trainc,  method = "qrf",trControl=trainControl( method = "cv",number=5,returnResamp = "all",savePredictions = TRUE, search = "random",verboseIter = FALSE))

print(rf.esp)

# Show the prediction interval
df.testc$Strain=predict(rf.esp,newdata=df.testc)
hist(df.testc$Strain,xlab="Box-Cox Transformed ESP (0-30cm)", main=NULL)
abline(v = quantile(df.testc$Strain, probs = c(0.05, 0.95)),lty = 5, col = "red")


lmbda3=(as.numeric(powerTransform(soil13$dummy, family ="bcPower")["lambda"]))
df.testc$Strain=predict(rf.esp,newdata=df.testc)
df.testc$Strain=(df.testc$Strain*lmbda3+1)^(1/lmbda3)
summary(df.testc$Strain)
summary(df.testc$dummy)

cor(df.testc$Strain,df.testc$dummy)^2
{plot(df.testc$Strain~df.testc$dummy, xlim=c(0,78),ylim=c(0,78),xlab="Harmonized ESP",ylab="Modelled ESP", main="Accuracy assessment on holdout samples for ESP")
  abline(a=0,b=1,lty=20, col="blue")}

Biasc=mean(df.testc$Strain-df.testc$dummy,na.rm=TRUE)
RMSEc=sqrt(sum((df.testc$Strain-df.testc$dummy)^2,na.rm=TRUE)/length((df.testc$Strain-df.testc$dummy)))
Rsquaredc=cor(df.testc$Strain,df.testc$dummy)^2
NSEc=1-sum((df.testc$Strain-df.testc$dummy)^2,na.rm=TRUE)/sum((df.testc$Strain-mean(df.testc$dummy,na.rm=TRUE))^2,na.rm=TRUE)
static=data.frame(Biasc,RMSEc,Rsquaredc,NSEc);View(static)
write.csv(static,file = "ESP0_30_validmodel_stats.csv")

##Step 6c: This part is optional: Plough back all the data and make the aggregate model for mapping EC
{
soil5c=as.data.frame(soil13)
rf.esp=train(Tran~(slope+loncurve+cnbl+rain+ls+geology+pgeology+valley+lcover+dem+PCA1+PCA2+PCA3+PCA4+PCA5),  data = soil5c,  method = "qrf",trControl=trainControl( method = "cv",number=5,returnResamp = "all",savePredictions = TRUE, search = "random",verboseIter = FALSE))

}

#Step 7c: Determine the final ESP map
lmbda3=(as.numeric(powerTransform(soil13$dummy, family ="bcPower")["lambda"]))
predictors$ESPt=predict(rf.esp,predictors)
predictors$ESP=(predictors$ESPt*lmbda3+1)^(1/lmbda3)
summary(predictors$ESP)

#Step 8c: Compare the predicted map and validation dataset
featureRep(predictors["ESP"],soil13)
summary(predictors$ESP);summary(df.testc$dummy)

#Step9c: Export the final map
hist(predictors$ESP)
writeRaster(raster(predictors["ESP"]), drivername = "GTiff", "Top0_30ESP.tif",overwrite=F)
writeRaster(raster(predictors["ESPt"]), drivername = "GTiff", "Top0_30ESPt.tif",overwrite=F)

#Step 10c: Determine location uncertainty for ESP
soil6c=soil13[,c("Tran")]
predictors6c=predictors[c("dem","slope","lcover","loncurve","pgeology","geology","ls","valley","PCA1","PCA2","PCA3","PCA4","PCA5")]

#This part runs simulations and can take a long time
pred_uncertc=predUncertain(soil6c,predictors6c,4,95,"qrandomforest")
spplot(pred_uncertc, "pred_width", scales = list(draw = TRUE),sp.layout = list(list("sp.points", soil12, pch = 3,cex=0.4,col="green")),legend=list())


ESP0_30_uncertain=(pred_uncertc$pred_width)
writeRaster(ESP0_30_uncertain, filename="ESP0_30_uncertain.tif",format="GTiff",overwrite=F)
writeRaster(pred_uncertc$pred_sd, filename="ESP0_30_sd.tif",format="GTiff",overwrite=F)

############PART THREE: MAPPING SALT-AFFECTED SOILS####
##### Part D: Classifying salt-affected soils #####
###Step 1d: Develop the salt-affected map
#Use this part if you do not have the raster maps in the predictors 
{
predictors=readGDAL("Top0_31ECse.tif")
predictors$PH=readGDAL("Top0_31PH.tif")$band1
predictors$ESP=readGDAL("Top0_31ESP.tif")$band1
predictors$ECse=predictors$band1
}
predictors$band1=NULL

#Else continue from this line
predictors$salty=saltClass(predictors$ECse,predictors$PH,predictors$ESP,"FAO")
summary(predictors$salty)
predictors$saltiness=classCode(predictors$salty,"saltclass")
spplot(predictors["salty"])
spplot(predictors["saltiness"],col.regions=heat.colors(10,rev = TRUE))

predictors$Salt_affected=saltSeverity(predictors$ECse,predictors$PH,predictors$ESP,"FAO")
predictors$saltaffectedness=classCode(predictors$Salt_affected,"saltseverity")
summary(predictors$saltaffectedness)
spplot(predictors["Salt_affected"])
spplot(predictors["saltaffectedness"],col.regions=heat.colors(20,rev = TRUE))

predictors$Saltclass2=as.numeric(predictors$saltaffectedness)
salinity_LUT2=classLUT(predictors["saltaffectedness"],"saltseverity")

writeRaster(raster(predictors["Saltclass2"]), drivername = "GTiff", "Saltaffected0_31.tif",overwrite=FALSE)
write.table(salinity_LUT2,file = "saltaffected_LUT0_30.txt",row.names = FALSE)


###Step 2d: Import the validation dataset and classify the salt classes
#Step 2d-1: Import the validation dataset 

{  soilv=readOGR(".","soilvalid")
  #or soilv=read.csv("validation_harmonized.csv",header=T)
  summary(soilv)}
  
#Step ii: convert the validation dataset into spatial dataframe
  coordinates(soilv) = ~Longitude+Latitude
  proj4string(soilv)=CRS("+proj=utm +zone=36 +datum=WGS84 +units=m +no_defs")
  #soilv=spTransform(soilv,CRS("+proj=utm +zone=36 +datum=WGS84 +units=m +no_defs"))

  

#Step 2d-2: CLassify the salt in the validation dataset
soilv=subset(soilv,soilv$Horizon==1)
soilv$salt_affected1=saltSeverity(soilv$EC,soilv$pH,soilv$ESP,"FAO")
summary(soilv$salt_affected1)
soilv$saltaffectedness1=classCode(soilv$salt_affected1,"saltseverity")
summary(soilv$saltaffectedness1)

#Step 3d: Extract the salt classes from the map using the validation samples 
soilv=subset(soilv,!is.na(soilv$saltaffectedness1))
predictors.ovv=over(soilv, predictors)
soilv$salt_affected=predictors.ovv$Salt_affected
soilv$saltaffectedness=predictors.ovv$saltaffectedness
summary(soilv$salt_affected)
summary(soilv$saltaffectedness)
soilv=subset(soilv,!is.na(soilv$saltaffectedness))
summary(soilv$saltaffectedness)

#Step 4d: Accuracy assessment for the map

confusion(soilv$salt_affected1, soilv$salt_affected) 
confusion(soilv$salt_affected1, soilv$salt_affected)
agreementplot(confusion(soilv$salt_affected1, soilv$salt_affected),margins = (mar=c(3,3,3,3)) ,xlab = "Holdout SAS classes", ylab = "Modelled SAS classes")$Bangdiwala

{
#Step 5d: Use this part for accuracy assessment in case of unequal cases in holdout and calibration samples

soild=soilv
soild$salt_affected=as.factor(soild$salt_affected);summary(soild$salt_affected)
levels(soild$salt_affected)=c(levels(soild$salt_affected),"15")
summary(soild$salt_affected)

soild$salt_affected1=as.factor(soild$salt_affected1);summary(soild$salt_affected1)
levels(soild$salt_affected1)=c(levels(soild$salt_affected1),"11")
summary(soild$salt_affected1)
confusion(soild$salt_affected1, soild$salt_affected)
agreementplot(confusion(soild$salt_affected1, soild$salt_affected),margins = (mar=c(3,3,3,3)) ,xlab = "Holdout SAS classes", ylab = "Modelled SAS classes")$Bangdiwala
}

######Part E: Uncertainty assessment #####
###Step 1e: Load the raster maps

#Step 1e-1: Import the data if notalready in R 
# This part if for importing the maps in R
# Use this part if the data is not yet in R environment
{
ECte=raster("Top0_31ECte.tif"); ECsd=raster("ECse0_31_sd.tif")
PHde=raster("Top0_31PHt.tif"); PHsd=raster("pH0_31_sd.tif")
ESPt=raster("Top0_31ESPt.tif"); ESPsd=raster("ESP0_31_sd.tif")

EC=raster("Top0_31ECse.tif");names(EC)=c(" EC")
PH=raster("Top0_31PH.tif"); names(PH)=c("PH")
ESP=raster("Top0_31ESP.tif"); names(ESP)=c("ESP")

EC1=as(EC,"SpatialPixelsDataFrame")
PH1=as(PH,"SpatialPixelsDataFrame")
ESP1=as(ESP,"SpatialPixelsDataFrame")
}

#Step 1e-2: COnvert R Grid data into Raster data
#Convert R grid objects into raster objects
#Use this part if the data is already in R as Grid data type
{
EC=raster(predictors["ECse"]); names(EC)=c("EC"); EC1=as(EC,"SpatialPixelsDataFrame")
PH=raster(predictors["PH"]); names(PH)=c("PH"); PH1=as(PH,"SpatialPixelsDataFrame")
ESP=raster(predictors["ESP"]); names(EC)=c("ESP");ESP1=as(ESP,"SpatialPixelsDataFrame") 

ECte=raster(predictors["ECte"]);ECsd=pred_uncerta$pred_sd; names(ECsd)=c("ECsd")
PHde=raster(predictors["PHt"]);PHsd=pred_uncertb$pred_sd; names(PHsd)=c("PHsd")
ESPt=raster(predictors["ESPt"]);ESPsd=pred_uncertc$pred_sd; names(ESPsd)=c("ESPsd")
}

###Step 2e: Estimate the spatial correlation models
slot(slot(EC1, "grid"), "cellsize")=rep(mean(slot(slot(EC1, "grid"), "cellsize")), 2)
b=nrow(EC1)
c=trunc(0.01*b)
jj=EC1[sample(b,c),] 

vrm=autofitVariogram(EC~1,jj)
plot(vrm)
acf((EC1$EC)) 

EC_crm <- makeCRM(acf0 = 0.85, range = 20000, model = "Sph")
plot(EC_crm, main = "EC correlogram")

slot(slot(PH1, "grid"), "cellsize")=rep(mean(slot(slot(PH1, "grid"), "cellsize")), 2)
b=nrow(PH1)
c=trunc(0.01*b)
jj=PH1[sample(b,c),] 

vrm=autofitVariogram(PH~1,jj)
plot(vrm)
acf((PH1$PH)) 

PH_crm <- makeCRM(acf0 = 0.85, range = 20000, model = "Sph")
plot(PH_crm, main = "PH correlogram")

slot(slot(ESP1, "grid"), "cellsize")=rep(mean(slot(slot(ESP1, "grid"), "cellsize")), 2)
b=nrow(ESP1)
c=trunc(0.01*b)
jj=ESP1[sample(b,c),] 

vrm=autofitVariogram(ESP~1,jj)
plot(vrm)
acf((PH1$PH)) 

ESP_crm <- makeCRM(acf0 = 0.85, range = 20000, model = "Sph")
plot(ESP_crm, main = "ESP correlogram")

##Step 3e: Define uncertainty model for the EC, PH and ESP
EC_UM <- defineUM(distribution = "norm", distr_param = c(ECte, ECsd),crm = EC_crm, id = "EC")
PH_UM <- defineUM(distribution = "norm", distr_param = c(PHde, PHsd),crm = PH_crm, id = "PH")
ESP_UM <- defineUM(distribution = "norm", distr_param = c(ESPt, ESPsd), crm = ESP_crm, id = "ESP")
class(EC_UM);class(PH_UM);class(ESP_UM)

cor(values(ECte),values(PHde)); cor(values(ECte),values(ESPt)); cor(values(PHde),values(ESPt)) #Get the correlation values and use them in the matrix in the next line for EC:PH:ESP
salinityMUM = defineMUM(UMlist = list(EC_UM, PH_UM, ESP_UM), cormatrix = matrix(c(1, cor(values(ECte),values(PHde)), cor(values(ECte),values(ESPt)),
                                                                                  cor(values(ECte),values(PHde)),1,cor(values(PHde),values(ESPt)),
                                                                                  cor(values(ECte),values(ESPt)),cor(values(PHde),values(ESPt)),1), nrow = 3, ncol = 3))

class(salinityMUM)

##Step 4e: create possible realizations from the joint distribution of OC and TN
MC <- 100
input_sample = genSample(UMobject = salinityMUM, n = MC, samplemethod = "ugs",  nmax = 20, asList = FALSE)

##Step 5e: compute and plot input sample statistics such as mean and standard deviation
EC_sample = input_sample[[1:MC]]
PH_sample = input_sample[[(MC+1):(2*MC)]]
ESP_sample = input_sample[[(2*MC+1):(3*MC)]]

EC_sample_mean <- mean(EC_sample)
PH_sample_mean <- mean(PH_sample)
ESP_sample_mean <- mean(ESP_sample)

EC_sample_sd <- calc(EC_sample, fun = sd)
PH_sample_sd <- calc(PH_sample, fun = sd)
ESP_sample_sd <- calc(ESP_sample, fun = sd)

par(mfrow=c(2,2),mar = c(1, 1, 2, 2), mgp = c(1.7, 0.5, 0), oma = c(0, 0, 0, 1),
    las = 1, cex.main = 1, tcl = -0.2, cex.axis = 0.8, cex.lab = 0.8)
plot(EC_sample_mean, main = "Mean of ECt realizations", xaxt = "n", yaxt = "n")
plot(PH_sample_mean, main = "Mean of PHt realizations", xaxt = "n", yaxt = "n")
plot(ESP_sample_mean, main = "Mean of ESPt realizations", xaxt = "n", yaxt = "n")

###Step 6e: Define the salinity model
Salinity_model_raster <- function (EC1,PH1,ESP1){
  ww=EC1
  ww=raster(ww)
  ww$salt=saltSeverity(values(EC1),values(PH1),values(ESP1),"FAO")
  ww=ww$salt; names(ww)=c("salt")
  ww
}

l <- list()
l[[1]] <- map(1:50, function(x){input_sample[[x]]})
l[[2]] <- map(51:100, function(x){input_sample[[x]]})
l[[3]] <- map(101:50, function(x){input_sample[[x]]})
input_sample=l

salinity_sample <- propagate(realizations = input_sample, model = Salinity_model_raster, n = MC)

# coerce salinity list to a raster stack and compute and plot the slope sample statistics
salinity_sample <- raster::stack(salinity_sample)
names(salinity_sample) <- paste("salt.", c(1:nlayers(salinity_sample)), sep = "")
salinity_freq = modal(salinity_sample, freq=TRUE)
salinity_prop = salinity_freq/100
salinity_SErr = sqrt(salinity_prop*(1-salinity_prop)/100)
CL=0.95
z_star=round(qnorm((1-CL)/2,lower.tail=F),digits = 2)
salinity_MErr=z_star*salinity_SErr

plot(salinity_MErr, main = "95% Prediction width of error margin on proportions", xaxt = "n", yaxt = "n")

##Step 7e: Export the final modelling uncertainty to output
writeRaster(salinity_MErr,filename="Salt0_30cm_uncertain1.tif",format="GTiff", overwrite=F)

###### PART FOUR: MAP UPDATE #####
##Step 1f: Import the soil indicators and classify
# use this part to import the data if they are not already in R
{
  predictors=readGDAL("Top0_31ECse.tif")#
  predictors$PH=readGDAL("Top0_31PH.tif")$band1
  predictors$ESP=readGDAL("Top0_31ESP.tif")$band1
  predictors$EC=predictors$band1
  predictors$band1=NULL

predictors$Salt_affected=saltSeverity(predictors$EC,predictors$PH,predictors$ESP,"FAO")
predictors$saltaffectedness=classCode(predictors$Salt_affected,"saltseverity")

}

##Step 2f: Identify the class to monitor or update
summary(predictors$saltaffectedness)
hist(predictors$Salt_affected, main = "Topsoil salt-affecetd classes ", xlab="Codes for salt-affected classes")
salts=predictors["saltaffectedness"]; salts=as(salts,"SpatialPixelsDataFrame")
salty=as.data.frame(salts)
salty1=data.frame(plyr::count(salty$saltaffectedness))
colnames(salty1)=c("Saltclass","Cases")
salty1$Props=round(salty1$Cases/sum(salty1$Cases)*100,1)
barchart(Saltclass~Props, data=salty1, xlab="Proportion (%)")

#Step 3f: Determine the number of sampling points to update the class
predictors$Saltclass=as.numeric(predictors$saltaffectedness)
salinity_LUT30=classLUT(predictors["saltaffectedness"],"saltseverity")
View(salinity_LUT30)

salt_affected_class=5
clorp_factors=4
soilsample=predictors["Saltclass"]
survey=surveyPoints(soilsample,clorp_factors,salt_affected_class,10)
length(survey$new)
spplot(soilsample, scales=list(draw=TRUE),sp.layout=list("sp.points",survey,pch=8,col="cyan"))

#Step 4f: Export the survey points to ESRI shapefile
writeOGR(survey,".","SurveyPointsClass5",driver = "ESRI Shapefile",overwrite_layer=FALSE)
