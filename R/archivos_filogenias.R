library(tidyverse)
library(treeio)
library(rlang)

setwd('/home/biocomp/HIP1/HIP_FIGURES/')

# Las tablas de los palindromos significativos las obtengo con:
# significativos.R

# CREAMOS df_tippoint.csv

# primero agrego la clasificación taxonomica a cada especie
# esta clasificación la hice con:
# perl taxonomy_id_to_tax.pl pico_2022/gbff_files/ pico_2022
# ADICIONALMENTE EDITE EL TIPO DE AGUA A MANO CON LA TABLA DEL ARTICULO

# Genome size
genomesize <- c()
genomesize$specie <-(tabla[,1])
genomesize$genomesize <-(tabla[,9])
write.table(genomesize, file=paste0("pico_df_genome_size.txt"), sep="\t", row.names = F)
GenomeSize = read.table(file=paste0("pico_df_genome_size.txt"), header=TRUE, sep="\t")
genomesize =unique(GenomeSize)

# cargamos tabla de clasificación
taxonomy<- read.table(file="filogenia_pico/pico_2022_taxonomy.txt", header=TRUE, sep="\t")
# pico_32_df_tippoint
for (i in 1:length(taxonomy[,1])){
  for (j in 1:length(genomesize[,1])){
    if ((isTRUE(taxonomy[i,1]==genomesize[j,1])==TRUE)){
      taxonomy[i,4] = genomesize[j,2]
    }
  }
}
write.table(taxonomy, file=paste0("pico_df_tippoint.txt"), sep="\t", row.names = F)

# CARGAMOS SIGNIFICTIVOS
tabla = read.table(file="Markov_count_all_pico_2022_fna_2022-9-7_7hrs54mins_.txt", header=TRUE, sep="\t")
sel32 = read.table(file="Marcov_pico_2022_significatives1e-32.txt", header=TRUE, sep="\t")
sel64 = read.table(file="Marcov_pico_2022_significatives1e-64.txt", header=TRUE, sep="\t")
sel128 = read.table(file="Marcov_pico_2022_significatives1e-128.txt", header=TRUE, sep="\t")
sel256 = read.table(file="Marcov_pico_2022_significatives1e-256.txt", header=TRUE, sep="\t")
## Funcion para obtener nombres
SigNames <- function(x) {
  d <- substitute(x)
  n <- sapply(d[-1],deparse)
  return(n)
}

df_sel <- list(sel32,sel64,sel128,sel256)
df_names <- SigNames(c(sel32,sel64,sel128,sel256))
names(df_sel)<-lapply(df_names, function(x) paste0("df_",x))

for (l in 1:length(df_sel)){
  df_l <- df_sel[[l]]
  ### *** ###
  
  # HACEMOS UNA LISTA DE ESPECIES 
  spp_sel <- c()
  for(i in 1:nrow(df_l)) {
    spp_sel <- append(spp_sel,df_l[i,1])
  }
  # QUITAMOS LAS REPETICIONES
  spp_sel <-unique(spp_sel)
  spp_sel
  
  # CREAMOS UNA LISTA VACIA QUE CONTENDRÁ AL PALINDROMO MAS ABUNDANTE POR ESPECIE
  highest_counts <- c()
  
  # BUSCAMOS EL PALINDROMO MAS ABUNDANTE:
  
  for(j in 1:length(spp_sel)) {         # PARA CADA ESPECIE:
    spp_pals_counts <- c()              # CREAMOS UNA LISTA DE LOS PALINDROMOS EXISTENTES
    print(spp_sel[j])
    for(i in 1:length(df_l[,1])) {     # PARA CADA RENGLON DE LA TABLA DE SIGNIFICATIVOS
      if (isTRUE(spp_sel[j] == df_l[i,1])==TRUE){ # SI LA ESPECIE COINCIDE CON LA QUE ESTAMOS CHECANDO
        spp_pals_counts <- append(spp_pals_counts, c(df_l[i,18])) # ENTONCES GUARDAMOS EL OE DEL PALINDROMO
      }                                                            # EN LA LISTA DE PALINDROMOS EXISTENTES
    }
    for (k in 1:length(df_l[,1])) {    # PARA CADA REGLON DE LA TABLA DE SIGNIFICATIVOS
      if (isTRUE((df_l[k,18])==spp_pals_counts[which.max(abs(spp_pals_counts))])==TRUE){ # SI EL OE EN CUESTION ES EL OE MAS GRANDE
        highest_counts <- rbind(highest_counts, c(df_l[k,])) # ENTONCES GUARDAMOS LA LINEA EN  HIGHEST COUNTS
        highest_counts <- unique(highest_counts)  # QUITAMOS LOS REPETIDOS
      }
    }
  }
  
  # CREAMOS LA TABLA CON LOS PALINDROMOS MAS SOBREREPRESENTADOS
  write.table(highest_counts, file=paste0("Highest_counts_",names(df_sel)[l],".txt"), sep="\t", row.names = F)
  
  # ABRIMOS LA TABLA PARA CREAR EL ARCHIVO
  tabla_l <- read.csv(file=paste0("Highest_counts_",names(df_sel)[l],".txt"), header=TRUE, sep="\t")
  tabla_l <- unique(tabla_l)
  
  # TABLA DE ESPECIES COMPLETAS CON 0's
  # En elcaso de las especies sin palindromos les agrego HIP con valor de 0
  # este no se cuenta en el grafico
  df_barplot_attr_zeros <- c()
  for (i in 1:length(tabla[,1])){
    vec_zeros <- c(tabla[i,1], "GCGATCGC", 0)
    df_barplot_attr_zeros <- rbind(df_barplot_attr_zeros, vec_zeros)
  }
  df_barplot_attr_zeros <- unique(df_barplot_attr_zeros)
  
  pico_df_barplot_attr <- c()
  # CREAMOS df_barplot_attr.txt
  for (i in 1:length(tabla_l[,1])){
    vec <- c(tabla_l[i,1], tabla_l[i,3], tabla_l[i,18])
    pico_df_barplot_attr <- rbind(pico_df_barplot_attr, vec)
  }
  # Combino las tablas con especies que tienen palindromos y las que no
  for (i in 1:length(df_barplot_attr_zeros[,1])){
    for (j in 1:length(pico_df_barplot_attr[,1])){
      if ((isTRUE(df_barplot_attr_zeros[i,1]==pico_df_barplot_attr[j,1])==TRUE)){
        df_barplot_attr_zeros[i,2] = pico_df_barplot_attr[j,2]
        df_barplot_attr_zeros[i,3] = pico_df_barplot_attr[j,3]
      }
    }
  }
  colnames(df_barplot_attr_zeros) = c("ID", "Sites", "HigherAbundance")
  write.table(df_barplot_attr_zeros, file=paste0("pico_",names(df_sel)[l],"_barplot_attr.txt"), sep="\t", row.names = F)
  
  # CREAMOS df_ring_heatmap
  # Primero hago una lista con todos los palindromos significativos
  pals_sel <- c()
  for (i in 1:nrow(df_l)) {
    pals_sel <- append(pals_sel,df_l[i,3])
  }
  pals_sel <- unique(pals_sel)
  # Para cada especie agrego el conteo para cada uno de los palindromos significativos (pals_sel)
  df_ring_heatmap_counts <- c()
  for (i in 1:length(df_l[,1])){
    vec_counts <- c(df_l[i,1], df_l[i,3], df_l[i,18])
    df_ring_heatmap_counts <- rbind(df_ring_heatmap_counts, vec_counts)
    df_ring_heatmap_counts <- unique(df_ring_heatmap_counts)
  }
  # Para cada especie agrego un CERO al conteo para cada uno de los palindromos significativos (pals_sel)
  df_ring_heatmap_zeros <- c()
  for (i in 1:length(tabla[,1])){
    for (j in 1:length(pals_sel)){
      vec_zeros <- c(tabla[i,1], pals_sel[j], 0)
      df_ring_heatmap_zeros <- rbind(df_ring_heatmap_zeros, vec_zeros)
    }
    df_ring_heatmap_zeros <- unique(df_ring_heatmap_zeros)
  }
  
  # Combino la tabla de CEROS y la de conteo de cada palindromo significativo
  for (i in 1:length(df_ring_heatmap_zeros[,1])){
    for (j in 1:length(df_ring_heatmap_counts[,1])){
      if ((isTRUE(df_ring_heatmap_zeros[i,1]==df_ring_heatmap_counts[j,1])==TRUE)&(isTRUE(df_ring_heatmap_zeros[i,2]==df_ring_heatmap_counts[j,2])==TRUE)){
        df_ring_heatmap_zeros[i,3] = df_ring_heatmap_counts[j,3]
      }
    }
  }
  colnames(df_ring_heatmap_zeros) = c("ID", "Sites", "Abundance")
  write.table(df_ring_heatmap_zeros, file=paste0("pico_",names(df_sel)[l],"_ring_heatmap.txt"), sep="\t", row.names = F)

  ### *** ###
  }
######
# HACEMOS UNA LISTA DE ESPECIES 
spp_sel <- c()
for(i in 1:nrow(sel32)) {
  spp_sel <- append(spp_sel,sel32[i,1])
}
# QUITAMOS LAS REPETICIONES
spp_sel <-unique(spp_sel)
spp_sel

# CREAMOS UNA LISTA VACIA QUE CONTENDRÁ AL PALINDROMO MAS ABUNDANTE POR ESPECIE
highest_counts <- c()

# BUSCAMOS EL PALINDROMO MAS ABUNDANTE:

for(j in 1:length(spp_sel)) {         # PARA CADA ESPECIE:
  spp_pals_counts <- c()              # CREAMOS UNA LISTA DE LOS PALINDROMOS EXISTENTES
  print(spp_sel[j])
  for(i in 1:length(sel32[,1])) {     # PARA CADA RENGLON DE LA TABLA DE SIGNIFICATIVOS
    if (isTRUE(spp_sel[j] == sel32[i,1])==TRUE){ # SI LA ESPECIE COINCIDE CON LA QUE ESTAMOS CHECANDO
      spp_pals_counts <- append(spp_pals_counts, c(sel32[i,18])) # ENTONCES GUARDAMOS EL OE DEL PALINDROMO
    }                                                            # EN LA LISTA DE PALINDROMOS EXISTENTES
  }
  for (k in 1:length(sel32[,1])) {    # PARA CADA REGLON DE LA TABLA DE SIGNIFICATIVOS
    if (isTRUE((sel32[k,18])==spp_pals_counts[which.max(abs(spp_pals_counts))])==TRUE){ # SI EL OE EN CUESTION ES EL OE MAS GRANDE
      highest_counts <- rbind(highest_counts, c(sel32[k,])) # ENTONCES GUARDAMOS LA LINEA EN  HIGHEST COUNTS
      highest_counts <- unique(highest_counts)  # QUITAMOS LOS REPETIDOS
    }
  }
}

# CREAMOS LA TABLA CON LOS PALINDROMOS MAS SOBREREPRESENTADOS
write.table(highest_counts, file="Highest_counts_sel32.txt", sep="\t", row.names = F)

# ABRIMOS LA TABLA PARA CREAR EL ARCHIVO
tabla_32 <- read.csv(file="Highest_counts_sel32.txt", header=TRUE, sep="\t")
tabla_32 <- unique(tabla_32)
pico_df_barplot_attr <- c()

# TABLA DE ESPECIES COMPLETAS CON 0's
# En elcaso de las especies sin palindromos les agrego HIP con valor de 0
# este no se cuenta en el grafico
df_barplot_attr_zeros <- c()
for (i in 1:length(tabla[,1])){
  vec_zeros <- c(tabla[i,1], "GCGATCGC", 0)
  df_barplot_attr_zeros <- rbind(df_barplot_attr_zeros, vec_zeros)
}
df_barplot_attr_zeros <- unique(df_barplot_attr_zeros)

# CREAMOS df_barplot_attr.txt
for (i in 1:length(tabla_32[,1])){
  vec <- c(tabla_32[i,1], tabla_32[i,3], tabla_32[i,18])
  pico_df_barplot_attr <- rbind(pico_df_barplot_attr, vec)
}
# Combino las tablas con especies que tienen palindromos y las que no
for (i in 1:length(df_barplot_attr_zeros[,1])){
  for (j in 1:length(pico_df_barplot_attr[,1])){
    if ((isTRUE(df_barplot_attr_zeros[i,1]==pico_df_barplot_attr[j,1])==TRUE)){
      df_barplot_attr_zeros[i,2] = pico_df_barplot_attr[j,2]
      df_barplot_attr_zeros[i,3] = pico_df_barplot_attr[j,3]
    }
  }
}
colnames(df_barplot_attr_zeros) = c("ID", "Sites", "HigherAbundance") 
write.table(df_barplot_attr_zeros, file="pico_32_df_barplot_attr.txt", sep="\t", row.names = F)

# CREAMOS df_ring_heatmap
# Primero hago una lista con todos los palindromos significativos
pals_sel <- c()
for (i in 1:nrow(sel32)) {
  pals_sel <- append(pals_sel,sel32[i,3])
}
pals_sel <- unique(pals_sel)
# Para cada especie agrego el conteo para cada uno de los palindromos significativos (pals_sel)
df_ring_heatmap_counts <- c()
for (i in 1:length(sel32[,1])){
  vec_counts <- c(sel32[i,1], sel32[i,3], sel32[i,18])
  df_ring_heatmap_counts <- rbind(df_ring_heatmap_counts, vec_counts)
  df_ring_heatmap_counts <- unique(df_ring_heatmap_counts)
}
# Para cada especie agrego un CERO al conteo para cada uno de los palindromos significativos (pals_sel)
df_ring_heatmap_zeros <- c()
for (i in 1:length(tabla[,1])){
  for (j in 1:length(pals_sel)){
    vec_zeros <- c(tabla[i,1], pals_sel[j], 0)
    df_ring_heatmap_zeros <- rbind(df_ring_heatmap_zeros, vec_zeros)
  }
  df_ring_heatmap_zeros <- unique(df_ring_heatmap_zeros)
}

# Combino la tabla de CEROS y la de conteo de cada palindromo significativo
for (i in 1:length(df_ring_heatmap_zeros[,1])){
  for (j in 1:length(df_ring_heatmap_counts[,1])){
    if ((isTRUE(df_ring_heatmap_zeros[i,1]==df_ring_heatmap_counts[j,1])==TRUE)&(isTRUE(df_ring_heatmap_zeros[i,2]==df_ring_heatmap_counts[j,2])==TRUE)){
      df_ring_heatmap_zeros[i,3] = df_ring_heatmap_counts[j,3]
    }
  }
}
colnames(df_ring_heatmap_zeros) = c("ID", "Sites", "Abundance") 
write.table(df_ring_heatmap_zeros, file="pico_32_df_ring_heatmap.txt", sep="\t", row.names = F)

# CREAMOS df_tippoint.csv

# primero agrego la clasificación taxonomica a cada especie
# esta clasificación la hice con:
# perl taxonomy_id_to_tax.pl pico_2022/gbff_files/ pico_2022
# ADICIONALMENTE EDITE EL TIPO DE AGUA A MANO CON LA TABLA DEL ARTICULO

# Genome size
genomesize <- c()
genomesize$specie <-(tabla[,1])
genomesize$genomesize <-(tabla[,9])
write.table(genomesize, file="genome_size_test.txt", sep="\t", row.names = F)
GenomeSize = read.table(file="genome_size_test.txt", header=TRUE, sep="\t")
genomesize =unique(GenomeSize)

# cargamos tabla de clasificación
taxonomy<- read.table(file="filogenia_pico/pico_2022_taxonomy.txt", header=TRUE, sep="\t")
# pico_32_df_tippoint
for (i in 1:length(taxonomy[,1])){
  for (j in 1:length(genomesize[,1])){
    if ((isTRUE(taxonomy[i,1]==genomesize[j,1])==TRUE)){
      taxonomy[i,4] = genomesize[j,2]
    }
  }
}
write.table(taxonomy, file="pico_32_df_tippoint.txt", sep="\t", row.names = F)
######

## Estas listas de especies las saco para ver si las
## etiquetas de todos los archivos estan bien y evitar problemas al graficas:
##  ***
##  ***
spp <- tabla[,1]
spp <- unique(spp)

tabla_heat = read.table(file="pico_32_df_ring_heatmap.txt", header=TRUE, sep="\t")
spp_heat <- tabla_heat[,1] 
spp_heat <- unique(spp_heat)

tabla_bar = read.table(file="pico_32_df_barplot_attr.txt", header=TRUE, sep="\t")
spp_bar <-tabla_bar [,1]
spp_bar <- unique(spp_bar)

tax = read.table(file="pico_2022_taxonomy.txt", header=TRUE, sep="\t")
spp_tax <-tax[,1]
spp_tax <- unique(spp_tax)

tree <- read.tree("filogenia_pico/SpeciesTree_rooted.txt")
spp_tree <- tree[["tip.label"]]
spp_tree <- unique(spp_tree)

etiquetas <- c()
etiquetas$Table <- spp
etiquetas$HeatMap <- spp_heat
etiquetas$BarPlot <- spp_bar
etiquetas$Taxonomy <- spp_tax
etiquetas$Tree <- spp_tree

write.table(etiquetas, file=paste0("Etiquetas_pico.txt"), sep="\t", row.names = F)

##  ***
##  ***