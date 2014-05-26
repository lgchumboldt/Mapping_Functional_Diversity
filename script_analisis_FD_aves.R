##################################
#####Load required packages######
#################################

library(maptools)
library(raster)
library(rgdal)
###############################
####Load areas for analysis####
###############################

    #Cargar paramos
    ruta_paramos<-"C:\\Users\\GIC 40\\Documents\\Andrea\\SIG\\Colombia\\paramos"
    paramos<-readOGR(dsn=ruta_paramos,layer="paramos_project_2012")
    #Cargar Colombia
    colombia_shape<-"C:\\Users\\GIC 66\\Documents\\Andrea\\SIG\\Colombia\\adm"
    colombia<-readOGR(dsn=colombia_shape,layer="COL_adm0")
    #Cargar cordilleras (cortado por encima de 2600m)
    cordillera<-"C:\\Users\\GIC 66\\Documents\\Andrea\\SIG\\Colombia\\cordillera\\cordillera.shp"
    cordillera<-readShapePoly(cordillera)
    
#Load rasters of Functional Diversity and Species richness
    #AVES
Fd_aves<-raster("C:\\Users\\GIC 66\\Documents\\Andrea\\Proyectos\\Diversidad_funcional\\riqueza_17marzo\\por_complejos\\sumas\\riqueza_funcional_por_complejo.asc")
Td_aves<-raster("C:\\Users\\GIC 66\\Documents\\Andrea\\Proyectos\\Diversidad_funcional\\riqueza_17marzo\\por_complejos\\sumas\\riqueza_especies_por_complejo.asc")

#Crear un marco de datos con valores de Pd y Td
  #AVES
marco_Fd_aves<-as.data.frame(Fd_aves)
marco_Td_aves<-as.data.frame(Td_aves)
marco_aves<-data.frame(marco_Fd_aves$riqueza_funcional_por_complejo,marco_Td_aves$riqueza_especies_por_complejo)


#Generar modelos de regresi�n 
  #LINEAL
  modelo_lineal<-lm(marco_aves)
  #LOESS
  modelo_loess<-loess(marco_aves)

#Calcular y guardar los residuales de las regresiones
  #LINEAL
  residuos_lineal<-resid(modelo_lineal)
  #LOESS
  residuos_loess<-resid(modelo_loess)
#Agregar los resiudales al marco de datos
marco_aves$residuos_Fd_lineal<-rep(NA,length(marco_aves[,1]))
marco_aves$residuos_Fd_lineal[as.numeric(names(residuos_lineal))]<-residuos_lineal[names(residuos_lineal)]

marco_aves$residuos_Fd_loess<-rep(NA,length(marco_aves[,1]))
marco_aves$residuos_Fd_loess[as.numeric(names(residuos_loess))]<-residuos_loess[names(residuos_loess)]

#Generar rasters con los datos de residuos

Fd_aves<-raster("C:\\Users\\GIC 66\\Documents\\Andrea\\Proyectos\\Diversidad_funcional\\riqueza_17marzo\\por_complejos\\sumas\\riqueza_funcional_por_complejo.asc")

residuos_lineal_raster<-Fd_aves
values(residuos_lineal_raster)<-NA
values(residuos_lineal_raster)<-marco_aves$residuos_Fd_lineal

Fd_aves<-raster("C:\\Users\\GIC 66\\Documents\\Andrea\\Proyectos\\Diversidad_funcional\\riqueza_17marzo\\por_complejos\\sumas\\riqueza_funcional_por_complejo.asc")

residuos_loess_raster<-Fd_aves
values(residuos_loess_raster)<-NA
values(residuos_loess_raster)<-marco_aves$residuos_Fd_loess
writeRaster(residuos_lineal_raster,"C:\\Users\\GIC 66\\Documents\\Andrea\\Proyectos\\Diversidad_funcional\\riqueza_17marzo\\por_complejos\\sumas\\residuals_lineal_FD.asc","ascii")
writeRaster(residuos_loess_raster,"C:\\Users\\GIC 66\\Documents\\Andrea\\Proyectos\\Diversidad_funcional\\riqueza_17marzo\\por_complejos\\sumas\\residuals_loess_FD.asc","ascii")




#Graphs and maps    
  #Regressions
  par(mfrow=c(2,2),mar=c(4,4,3,3))

  plot(marco_aves$marco_Td_aves.riqueza_especies_por_complejo,marco_aves$marco_Fd_aves.riqueza_funcional_por_complejo,xlab="N�mero de especies",ylab="Diversidad funcional", main="Regresi�n lineal")
  abline(modelo_lineal,col="red")
  plot(marco_aves$marco_Td_aves.riqueza_especies_por_complejo,marco_aves$marco_Fd_aves.riqueza_funcional_por_complejo,xlab="N�mero de especies",ylab="Diversidad funcional", main="Regresi�n loess")
  Td_vals<-seq(0,600,1)
  loess_vals<-predict(modelo_loess,Td_vals)
  lines(loess_vals ~ Td_vals,col="red")

#Residuals

plot(marco_aves$marco_Td_aves.riqueza_especies_por_complejo,marco_aves$residuos_Fd_lineal,xlab="Number of species",ylab="Residuals of functional diversity")
abline(0,0)
plot(marco_aves$marco_Td_aves.riqueza_especies_por_complejo,marco_aves$residuos_Fd_loess,xlab="Number of species",ylab="Residuals of functional diversity")
abline(0,0)

#Maps
par(mfrow=c(1,1))
plot(colombia,main="Linear regression")
plot(residuos_lineal_raster,box=F,axes=F,add=T)
plot(colombia, main="Loess regression")
plot(residuos_loess_raster,box=F,axes=F,add=T)
