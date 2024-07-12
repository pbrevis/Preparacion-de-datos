##Cargando paquetes de R
library(polymapR)
library(hexbin)
library(dplyr)

rm(list = ls()) # Limpiando el ambiente global

# Guardando la ubicación del directorio de trabajo,
# una vez que fuese definida manualmente
directorio <- getwd()
head(directorio)

##Leyendo archivo *.csv con datos de genotipo (SNPs)
##en población segregante bi-parental F1
##Rosa L. ‘Radbrite’ (Brite Eyes) x Rosa L. ‘BAIgirl’ (Easy Elegance My Girl)
##(BExMG) n=157

bexmg_df <- read.csv("BExMG_dosage.csv",
                        stringsAsFactors = FALSE,
                        row.names = 1)
                        #primera columna tiene nombres de marcadores

#Vista rápida al df
glimpse(bexmg_df)

#Confirmando que los datos fueron ingresados a un 'data frame'
class(bexmg_df)

#polymapR requiere los datos en forma de matriz
bexmg_pop <- as.matrix(bexmg_df)

#Confirmando que los datos ahora están en forma de matriz
class(bexmg_pop)

#Vista rápida a las primeras columnas (1-a-6) de la matriz de datos
head(bexmg_pop[,1:6])


#Dimensiones de la matriz: filas (marcadores) x columnas (individuos)
#m=35565 marcadores
#n=161 individuos F1 + 2 padres
dim(bexmg_pop)

#Chequeando que proporciones de segregación corresponden a lo esperado
#según dosis de los padres y herencia tetrasómica
f1checked <- checkF1(dosage_matrix = bexmg_pop, parent1 = "P1", parent2 = "P2",
                     F1 = colnames(bexmg_pop)[3:ncol(bexmg_pop)],
                     polysomic = TRUE, disomic = FALSE, mixed = FALSE,
                     ploidy = 4)
                     #Rosa de jardín se considera autopoliploide


#Revisando resultados de desviación de la distribución de probabilidad.
#Altos valores (~1) de 'qall_mult' y 'qall_weights' indican consistencia
#entre dosis de los padres y segregación de marcadores en la población F1
glimpse(f1checked$checked_F1)


#check class
class(f1checked$checked_F1)

#check dimensions
dim(f1checked$checked_F1)

#Varios pasos a continuación para filtrar marcadores, dejando afuera
#aquellos que se desvían de proporciones de segregación esperadas
all_markers <- f1checked$checked_F1
nonskewed_markers <- filter(all_markers, qall_weights > 0.9)
                     # sólo conservar marcadores cuyo valor qall_weights > 0.9

#Valores promedio y mínimo de qall_mult entre los marcadores seleccionados
mean(nonskewed_markers$qall_mult)
min(nonskewed_markers$qall_mult)

nonskewed_marker_names <- nonskewed_markers$MarkerName

#Vista previa: nombre de marcadores confiables
head(nonskewed_marker_names)

#Número de marcadores confiables
length(nonskewed_marker_names)

#Porcentaje de marcadores seleccionados
100*(length(nonskewed_marker_names) / length(all_markers$MarkerName))

#Nuevo dataframe, excluyendo marcadores no confiables (skewed)
#Conservando sólo 21254 de 35565 marcadores (~60%)
bexmg_df2 <- bexmg_df %>% filter(row.names(bexmg_df) %in% nonskewed_marker_names)

#polymapR requiere los datos en forma de matriz
#nueva matriz bexmg_pop2 ahora sólo con marcadores confiables
bexmg_pop2 <- as.matrix(bexmg_df2)

#Confirmando que los datos ahora están en forma de matriz
class(bexmg_pop2)

#PCA de la dosis de marcadores, para detectar posibles individuos problemáticos
PCA_progeny(dosage_matrix = bexmg_pop2, 
            highlight = list(c("P1", "P2")), 
            colors = "red")

#Opcional: eliminar individuo extremo por su nombre
#bexmg_pop2 <- subset(bexmg_pop, select=-c(X16401_N019))


#Creando resumen sobre calidad de datos
markerdatasum <- marker_data_summary(dosage_matrix = bexmg_pop2,
                           ploidy = 4,
                           pairing = "random",
                           parent1 = "P1",
                           parent2 = "P2",
                           progeny_incompat_cutoff = 0.05)

#Histograma de dosis en los padres (todos los marcadores, n=35,565)
pq_before_convert <- parental_quantities(dosage_matrix = bexmg_pop2, 
                                         las = 2)

#Convirtiendo tipos de segregación a su dosis más simple, para facilitar computación
segregating_data <- convert_marker_dosages(dosage_matrix = bexmg_pop2, ploidy = 4)

#Histograma de dosis en los padres, luego de convertir tipos de segregación
pq_after_convert <- parental_quantities(dosage_matrix = segregating_data)

#Control de calidad de datos: eliminando marcadores (filas) con problemas
screened_data <- screen_for_NA_values(dosage_matrix = segregating_data, 
                                      margin = 1, # opción 1 indica marcadores
                                      cutoff =  0.1, # tolerancia de 10% NA
                                      print.removed = FALSE) 

#Control de calidad de datos: eliminando individuos (columnas) con problemas
screened_data2 <- screen_for_NA_values(dosage_matrix = screened_data, 
                                       margin = 2, # opción 2 indica individuos                                       
                                       cutoff = 0.1, # tolerancia de 10% NA

                                       print.removed = FALSE)

#Revisando individuos duplicados
screened_data3 <- screen_for_duplicate_individuals(dosage_matrix = screened_data2, 
                                                   cutoff = 0.95, 
                                                   plot_cor = TRUE)

#Revisando marcadores duplicados
screened_data4 <- screen_for_duplicate_markers(dosage_matrix = screened_data3)

#Opcional: identificar marcadores confiables en base al tamaño del bin
#Formula M/(4NLy)
#reliable.markers <- names(which(sapply(screened_data4$bin_list,length) >= 6))
#reliable_data <- screened_data4$filtered_dosage_matrix[reliable.markers,]

#Extrayendo la matriz de dosis, elemento al interior de 'screened_data4', que
#es producto de la función 'screen_for_duplicate_markers'
filtered_data <- screened_data4$filtered_dosage_matrix

#Histograma de dosis en los padres, luego de eliminar marcadores duplicados
#Esta gráfica muestra el número final de marcadores a usar en el desarrollo
#del mapa genético
pq_screened_data <- parental_quantities(dosage_matrix = filtered_data)



# Fin de la preparación de datos para desarrollar mapas de ligamiento