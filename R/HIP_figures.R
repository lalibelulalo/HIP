library(ggplot2)

setwd('/home/biocomp/HIP1/HIP_FIGURES/')
tabla = read.table(file="Markov_count_all_pico_2022_fna_2022-9-21_11hrs24mins_HIP_LIKE.txt", header=TRUE, sep="\t")
tabla$pval = pbinom((tabla$obs-1),
                    (tabla$genomesize -8+1),
                    (tabla$markov3/(tabla$genomesize-8+1)),lower.tail = FALSE) # log.p = "TRUE"
tabla$fdrs <- p.adjust(tabla$pval, method="fdr")
tabla$frecObs <- (1000*(tabla$obs)/(tabla$genomesize))
tabla$OE <- ((tabla$obs)/(tabla$markov3))
tabla$GC <-(((tabla$G)+(tabla$C))/((tabla$A)+(tabla$Th)+(tabla$G)+(tabla$C)))

write.table(tabla, file="Markov_count_all_pico_2022_fna_2022-9-21_11hrs24mins_HIP_LIKE_e3.txt", sep="\t", row.names = F)

# "spp"
# "spp.1" tabla_s[,2]
# "palindrome" tabla_s[,3]
# "obs"
# "markov0"
# "markov1"
# "markov2"
# "markov3"
# "genomesize"
# "A"
# "T"
# "C"
# "G"
# "N"
# "pval"  tabla_s[,15]
# "fdrs"  tabla_s[,16]
# "frecObs" tabla_s[,17]
# "OE" tabla_s[,18]
# "GC"tabla_s[,19]

# calothrix_df_ring_heatmap.csv
# ID  Sites Abundance
# Calothrix_sp._336_3,ACGCGCGT,0

# calothrix_df_barplot_attr.csv
# ID,Sites,HigherAbundance
# Calothrix_sp._336_3,GCGATCGC,65.3963064436089

# calothrix_df_tippoint.csb
# ID,Phylum,Type,Size
# Calothrix_sp._336_3,Nostocales,freshwater lake,6.798186

#----------------------
## TABLAS SIGNIFICANTES
prefijo = "Marcov_pico_2022_significatives_HIP_like"
pvalues <- c(1e-32,1e-64,1e-128,1e-256)
filenames = c()
plotnames = c()
for (pvalue in pvalues) {
  FileName <- paste(prefijo,pvalue,".txt",sep="")
  PlotName <- paste(prefijo,pvalue,".png",sep="")
  filenames <- c(filenames,FileName)
  plotnames <- c(plotnames,PlotName)
  significatives.pvalue <- c()
  
  for(i in 1:nrow(tabla)) {
    if(tabla[i,15] < pvalue){
      significatives.pvalue <- rbind(significatives.pvalue,c(tabla[i,]))
      write.table(significatives.pvalue, file=FileName, sep="\t", row.names = F)
    }
  }
}

colores = c("blueviolet","blue","green","magenta")
tablas = filenames
temp_tab_2 = read.table(file=tablas[i], header=TRUE, sep="\t")
sel32=read.table(file="Marcov_pico_2022_significatives1e-32.txt", header=TRUE, sep="\t")
sel64=read.table(file="Marcov_pico_2022_significatives1e-64.txt", header=TRUE, sep="\t")
sel128=read.table(file="Marcov_pico_2022_significatives1e-128.txt", header=TRUE, sep="\t")
sel256=read.table(file="Marcov_pico_2022_significatives1e-256.txt", header=TRUE, sep="\t")

selections = c(sel32, sel64, sel128, sel256)

## Figura significancias
ggplot(data = tabla,
       aes(x = log10(frecObs),
           y = log10(obs/markov3))) +
  geom_point(color = "gray", size = 3, alpha = 1/5) +
  geom_point(data = sel32, color = "green") +
  geom_point(data = sel64, color = "yellow") +
  geom_point(data = sel128, color = "orange") +
  geom_point(data = sel256, color = "red")

## Figura palindromos significantes

ggplot(data = tabla,
       aes(x = log10(frecObs),
           y = log10(obs/markov3))) +
  geom_point(color = "gray", size = 3, alpha = 1/5) +
  geom_point(data=sel32, aes(x = log10(frecObs), y = log10(obs/markov3), col = palindrome))
ggsave("palindromos_significantes.pdf",p3, width=6, height=6, units="in", scale=3)
