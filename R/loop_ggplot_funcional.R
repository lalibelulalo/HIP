library(rlang)
library(tidyverse)

# CARGAMOS DATOS
setwd('/home/biocomp/HIP1/HIP_FIGURES/')
tabla = read.table(file="Markov_count_all_pico_2022_fna_2022-9-7_7hrs54mins_.txt", header=TRUE, sep="\t")

# CALCULAMOS PVALUES
tabla$pval = pbinom((tabla$obs-1),
                    (tabla$genomesize -8+1),
                    (tabla$markov3/(tabla$genomesize-8+1)),lower.tail = FALSE) # log.p = "TRUE"
# CALCULAMOS FDRS, FRECUENCIAS OBSERVADAS, OE y CONTENIDO GC
tabla$fdrs <- p.adjust(tabla$pval, method="fdr")
tabla$frecObs <- (1000*(tabla$obs)/(tabla$genomesize))
tabla$OE <- ((tabla$obs)/(tabla$markov3))
tabla$GC <-(((tabla$G)+(tabla$C))/((tabla$A)+(tabla$Th)+(tabla$G)+(tabla$C)))

# IDENTIFICAMOS LOS CONJUNTOS DE PALINDROMOS SIGNIFICATIVOS
# Y CREAMOS UNA TABLA PARA CADA PVALUE 
prefijo = "Marcov_pico_2022_significatives"
pvalues <- c(1e-32,1e-64,1e-128,1e-256)
filenames =c()
for (pvalue in pvalues) {
  FileName <- paste(prefijo,pvalue,".txt",sep="")
  filenames <- c(filenames,FileName)
  significatives.pvalue <- c()
  
  for(i in 1:nrow(tabla)) {
    if(tabla[i,15] < pvalue){
      significatives.pvalue <- rbind(significatives.pvalue,c(tabla[i,]))
      write.table(significatives.pvalue, file=FileName, sep="\t", row.names = F)
    }
  }
}

# CARGAMOS TABLAS PARAQ CADA PVALUE DE CORTE 
sel32=read.table(file="Marcov_pico_2022_significatives1e-32.txt", header=TRUE, sep="\t")
sel64=read.table(file="Marcov_pico_2022_significatives1e-64.txt", header=TRUE, sep="\t")
sel128=read.table(file="Marcov_pico_2022_significatives1e-128.txt", header=TRUE, sep="\t")
sel256=read.table(file="Marcov_pico_2022_significatives1e-256.txt", header=TRUE, sep="\t")
df2 <- list(sel32,sel64,sel128,sel256)

# GRAFICAMOS CADA TABLA
library(tidyverse)
#Create dataframes(In this example n = 3)

df_1 <- sel32  
df_2 <- sel64
df_3 <- sel128
df_4 <- sel256
SigNames <- function(x) {
  d <- substitute(x)
  n <- sapply(d[-1],deparse)
  return(n)
}

# GUARDAMOS DATAFRAMES EN UNA LISTA
all_data <- list(sel32,sel64,sel128,sel256)

df_names <- SigNames(c(sel32,sel64,sel128,sel256))
names(all_data)<-lapply(df_names, function(x) paste0("df_",x))

#Graph and save for each dataframe

for (i in 1:length(all_data)){
  df_i <- all_data[[i]]
  benp <-  
    #df_i %>%
    ggplot(data = tabla,
           aes(x = log10(frecObs),
               y = log10(obs/markov3))) +
    geom_point(color = "gray", size = 0.5, alpha = 1/5) +
    geom_point(data=df_i, aes(x = log10(frecObs), y = log10(obs/markov3), col = palindrome)) +
    labs(title=paste0("Palindromos ", names(all_data)[i]), subtitle="") 

    ggsave(benp, file=paste0(names(all_data)[i],"_significativepalindromes.png"),width=7, height=6, units="in", scale=1.5)
}

# BORRAMOS LAS TABLAS
sapply(paste0("Marcov_pico_2022_significatives1e-", 32:256, ".txt"), unlink)

