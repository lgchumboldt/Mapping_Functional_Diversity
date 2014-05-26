rm(list=ls())
########################################################################################################################
###############################################I-Load packages##########################################################
#########################################################################################################################

library(ade4)
library(ape)
library (FD)
library(geometry)
library(maptools)
library(raster)
library(rgdal)

########################################################################################################################
#########################II-Read data necessary for script and set user determined properties############################
########################################################################################################################


{# set working directory
print("Choose your working directory:")
working_directory<-choose.dir()
setwd(working_directory)

#Select folder containing distribution maps
print("Choose de folder containing distribution maps")
maps_folder<-choose.dir()

#Select folder to store functional diversity and taxonomic diversity maps
print("Select the folder where both the functional diversity and species richness maps will be stored: ")
save_files<-choose.dir()

#Trait data: select file containing trait data
print("Select the file with the functional characters:")
trait=read.csv(file.choose())

#Select the mask to be used in the analysis. The mask must be in a shapefile format, *.shp (Colombia, Paramos , cordilleras etc.)
mask1<-file.choose()

#Chose the names for the output files
functional_diversity_map<-readline("Name of the functional diversity map: NO spaces and INCLUDE extension, .asc (Ej: Functional_diversity.asc): ")

species_richness<-readline("Name of the species richness map: NO spaces and INCLUDE extension (Ej: Species_richness.asc): ")
#Determinar resolucion
res<-readline("Raster resolution: ")
}


###############################################################################################################
########################################III-Definition of functions to be used#####################################
###############################################################################################################


### 1- Function: Functional diversity (FRic, FEve, FDiv, FSpe)  ###


### Author: Sebastien Vill�ger, adapted by Claire Fortunel (Please acknowledge as appropriate)  

#  Notations corresponds with Vill�ger et al. (2008) Ecology, 89: 2290-2301 for FRic, FEve, FDiv; and Bellwood et al. (2006) Proc. R. Soc. B., 273: 101-107 for FSpe

# Function to calculate the four Functional diversity indices (So far we only use the FRic index since no abundance data is available from distribution maps)

FDind=function(trait,abund) {
  # T = number of traits
  T=dim(trait)[2]
  # c = number of communities
  C=dim(abund)[1]
  # check coherence of number of species in 'traits' and 'abundances'
  if (dim(abund)[2]!=dim(trait)[1]) stop(" Error : different number of species in 'trait' and 'abund' matrices ")
  # check format of traits values
  if (ncol(trait)<2) stop ("'Trait' must have at least 2 columns")
  if (is.numeric(trait)==F) stop ("Traits values must be numeric")
  # check absence of NA in 'traits'
  if (length(which(is.na(trait)==T))!=0) stop(" Error : NA in 'trait' matrix ")
  # replacement of NA in 'abund' by '0'
  abund[which(is.na(abund))]=0
  # definition of vector for results, with communities'names as given in 'abund'
  Nbsp=rep(NA,C) ; names(Nbsp)=row.names(abund)
  FRic=rep(NA,C) ; names(FRic)=row.names(abund)
  FEve=rep(NA,C) ; names(FEve)=row.names(abund)
  FDiv=rep(NA,C) ; names(FDiv)=row.names(abund)
  FSpe=rep(NA,C) ; names(FSpe)=row.names(abund)
  # scaling and centering of each trait according to all species values
  traitCS=scale(trait, center=TRUE, scale=TRUE)
  # functional specialization of each species (distance to point 0,0 in the standardized functional space)
  FSpeS=(apply(traitCS, 1, function(x) {x%*%x}))^0.5
  # loop to compute FRic, FEve, FDiv and FSpe on each community
  for (i in 1:C){
    # selection of species present in the community
    esppres=which(abund[i,]>0)
    #  number of species in the community
    S=length(esppres) ; Nbsp[i]=S
    # check if more species than traits
   # if (S<=T) stop(paste("Number of species must be higher than number of traits in community:",row.names(abund)[i]))
    # filter on 'trait' and 'abund' to keep only values of species present in the community
    tr=traitCS[esppres,] ; ab=as.matrix(abund[i,esppres])
    # scaling of abundances
    abondrel=ab/sum(ab)
    # Functional Diversity Indices
      # FRic
        # Using convhulln function
              # volume
              FRic[i]=round(convhulln(tr,"FA")$vol,6)
              # identity of vertices
              vert0=convhulln(tr,"Fx TO 'vert.txt'")
              vert1=scan("vert.txt",quiet=T)
			        vert2=vert1+1
              vertices=vert2[-1]
      # FEve
        # computation of inter-species euclidian distances
              distT=dist(tr, method="euclidian")
        # computation of Minimum Spanning Tree and conversion of the 'mst' matrix into 'dist' class
              linkmst=mst(distT) ; mstvect=as.dist(linkmst)
        # computation of the pairwise cumulative relative abundances and conversion into 'dist' class
              abond2=matrix(0,nrow=S,ncol=S)
              for (q in 1:S)
              for (r in 1:S)
              abond2[q,r]=abondrel[q]+abondrel[r]
              abond2vect=as.dist(abond2)  # end of q,r
        # computation of weighted evenness (EW) for the (S-1) branches to link S species
              EW=rep(0,S-1)
              flag=1
              for (m in 1:((S-1)*S/2)){if (mstvect[m]!=0) {EW[flag]=distT[m]/(abond2vect[m]) ; flag=flag+1}}  # end of m
        # computation of the partial weighted evenness (PEW) and comparison with 1/S-1, and computation of FEve
              minPEW=-rep(0,S-1) ; OdSmO=1/(S-1)
              for (l in 1:(S-1))
                minPEW[l]=min((EW[l]/sum(EW)), OdSmO)  # end of l
              FEve[i]=round(((sum(minPEW))- OdSmO)/(1-OdSmO),6)
      # FDiv
        # traits values of vertices of the convex hull
              trvertices=tr[vertices,]
        # coordinates of the center of gravity of the vertices (Gv) of the convex hull
              baryv=apply(trvertices,2,mean)
        # euclidian distances to Gv (dB) of each of S species (centro de gravedad)
              distbaryv=rep(0,S)
              for (j in 1:S)
                distbaryv[j]=(sum((tr[j,]-baryv)^2) )^0.5  # end of j
        # mean euclidian distance to the center of gravity of the S species (i.e. mean of dB values)  (mean centro gravedad)
# bigger effect if you are very abundant (andrea)
              meandB=mean(distbaryv)
        # deviation of each species dB from mean dB
              devdB=distbaryv-meandB
        # relative abundances-weighted mean deviation
              abdev=abondrel*devdB
        # relative abundances-weighted mean of absolute deviations
              ababsdev=abondrel*abs(devdB)
        # computation of FDiv
              FDiv[i]=round((sum(abdev)+meandB)/(sum(ababsdev)+meandB),6)              
      # FSpe
        # mean functional specialization in the communities
              FSpe[i]=(abund[i,]/sum(abund[i,]))%*%FSpeS
  } # end of i
  # result storage
  res=data.frame(Nbsp=Nbsp, FRic=FRic, FEve=FEve, FDiv=FDiv, FSpe=FSpe) ; row.names(res)=row.names(abund)
  invisible(res)
}# end of function


### 2- Function: Rasterize distribution maps and cut to mask extent  ###

rasterize_species= function (x,mask=mascara) {
  r<-raster(ncol=1462,nrow=624)
  res(r)<-resolucion #resolution
  r<-crop(r, extent(mascara))
  values(r)<-0
  map<-readOGR(dsn=distribution_maps_folder,layer=x)
  r<-rasterize(map,r,1,update=T,background=0)
  r<-mask(r,mascara)
  valor<-unique(getValues(r))
  
  if(length(valor)==1&&is.na(valor)==TRUE){
    
    
  }
  else {
    x
    writeRaster(r,paste(x,".asc"))
    return (raster(paste(x,".asc")))
  }
}#end of function

####################################################################################################################
########################################IV- Genarate a map of Functional Diversity######################################
####################################################################################################################


#######################################
### Script begins ###
######################################
#log transform all variables in trait
trait$L_Culmen_expuesto<-log(trait$L_Culmen_expuesto)
trait$Ancho_Pico<-log(trait$Ancho_Pico)
trait$L_Ala<-log(trait$L_Ala)
trait$L_Cola<-log(trait$L_Cola)
trait$L_Hallux<-log(trait$L_Hallux)
# Trait data: select file containing trait data

str(trait)
head(trait)
summary(trait)
traitnb=dim(trait)[2]-1
spnb=dim(trait)[1]

####Comunidades####

#set working directory to folder containing distribution maps
setwd(maps_folder)

#Determine folder containing distribution maps and generate list of species
distribution_maps_folder<-maps_folder
distribution_files<-list.files(path=distribution_maps_folder, pattern= "*.shp$")
species_names<-sub(".shp","",distribution_files)
tabla<-as.data.frame(species_names)
colnames(tabla)<-"Grilla"

#Read polygon with mask
mask<-readShapePoly(mask1)

#Determine working resolution
resolution<-as.numeric(res)

#Load all distribution maps and rasterize
r<-raster(ncol=1462,nrow=624)
res(r)<-resolution #resolution
r<-crop(r,extent(mask))
grilla=r
names(grilla)="grilla"
grilla[1:ncell(grilla)]<-1:ncell(grilla)
lista_grilla<-list(grilla)
names(lista_grilla)<-"Grilla"
layers<-lapply(species_names,rasterize_species)
names(layers)<-as.vector(tabla$Grilla)
layers[sapply(layers,is.null)]<-NULL
lista_completa<-c(lista_grilla,layers)
Stack<-stack(lista_completa)

#Turn maps into dataframe for computation of FD
marco<-as.data.frame(Stack)
marco=na.omit(marco)


# Select traits

traitc=trait[,2:6]
traitcn=apply(traitc,2,as.numeric)

# Calculate FRic, FEve, FDiv and FSpe

FD=as.numeric()
Community=unique(marco$Grilla)

for (i in 1:length(Community)) {
  abund.i=marco[marco$Grilla==Community[i],1:250]
  abundm.i=as.matrix(abund.i)
  FD.i=FDind(traitcn,abundm.i)
  FD=rbind(FD,FD.i)
}


marco$fd<-FD[,2]

r<-raster(ncol=1462,nrow=624)
res(r)<-resolution
r<-crop(r,extent(mask))
values(r)<-0
map<-readOGR(dsn=distribution_maps_folder,layer=species_names[[1]])
r<-rasterize(map,r,map$PRESENCE,update=T)
r<-mask(r,mask)
fd_ras<-r
values(fd_ras)<-NA #se eliminan todos los valores del modelo de distribuci�n


#Asignar al raster los valores de PD que corresponden a cada pixel
fd_ras[marco$Grilla]<-marco$fd


#Compute taxonomic diversity as number of species
aves_col<-stack(layers)
TD<-calc(aves_col,sum)


#Plot both maps in R to the left functional diversity and to the right taxonomic diversity
par(mfrow=c(1,2))

plot(fd_ras, main="Diversidad Funcional")
plot(TD, main="Riqueza de especies")

#Write rasters (FD and TD) to file
setwd(save_files)
writeRaster(fd_ras,functional_diversity_map)
writeRaster(TD,species_richness)






