library(ggtreeExtra)
library(ggtree)
library(treeio)
library(tidytree)
library(ggstar)
library(ggplot2)
library(ggnewscale)
library(TDbook)

setwd('/home/biocomp/HIP1/HIP_FIGURES/filogenia_pico/')
tree <- read.tree("arbol_raiz.txt")
dat1 <- read.table(file="pico_df_tippoint.txt", header=TRUE, sep="\t", stringsAsFactors = FALSE)
pwidths <- c(0.8,0.4,0.2,0.2)
contador = 1
## phylums
phylum <- c()
for (i in 1:nrow(dat1)) {
  phylum <- append(phylum,dat1[i,2])
}
phylum <- unique(phylum)
phylum
## colores para el phylum
cols_phylum <- rainbow(length(phylum))

dataframes <- c("sel32","sel64","sel128","sel256")
for (l in dataframes){
  dat2_l<-read.table(file=paste0("pico_df_",l,"_ring_heatmap.txt"), header=TRUE, sep="\t")
  dat3_l<-read.table(file=paste0("pico_df_",l,"_barplot_attr.txt"), header=TRUE, sep="\t")
  
  etiquetas_l <- c()
  etiquetas_l$tip_point <- sort(dat1[,1]) 
  etiquetas_l$ring_heatmap <- sort(unique(dat2_l[,1]))
  etiquetas_l$barplot <- sort(dat3_l[,1])
  etiquetas_l$arbol <- sort(tree[["tip.label"]])
  write.table(etiquetas_l, file=paste0("pico_",l,"_etiquetas.tsv"), sep="\t", row.names = F)
  
  ###
  ## palindromos
  pals_l <- c()
  for (i in 1:nrow(dat2_l)) {
    pals_l <- append(pals_l,dat2_l[i,2])
  }
  pals_l <- unique(pals_l)
  pals_l
  ## colores para los palindromos
  cols_l <- rainbow(length(pals_l))
  ###
  
  nodeids <- nodeid(tree, tree$node.label[nchar(tree$node.label)>4])
  nodedf <- data.frame(node=nodeids)
  nodelab <- gsub("[\\.0-9]", "", tree$node.label[nchar(tree$node.label)>4])
  
  poslist <- c(0.1)
  labdf <- data.frame(node=nodeids, label=nodelab, pos=poslist)
  
  p <- ggtree(tree, size=0.1, open.angle=5)+ # layout="fan" #, open.angle=5**
    geom_tiplab(size=2, align=TRUE, offset=0.01) +
    geom_hilight(data=nodedf, mapping=aes(node=node),
                 extendto=0.6, alpha=0.3, fill="grey", color="grey50",
                 size=0.05) +
    geom_cladelab(data=labdf, 
                  mapping=aes(node=node, 
                              label=label,
                              offset.text=pos),
                  hjust=0.5,
                  barsize=NA,
                  horizontal=FALSE, 
                  fontsize=1.4,
                  fontface="italic"
    ) 
  #p
  p2 <- p
  p2 <- p2 + theme(plot.margin = unit(c(14,8,14,8), "mm"))
  #p2
  
  p2 <- p2 %<+% dat1 + geom_star(
    mapping=aes(fill=Phylum, starshape=Type, size=Size),
    position="identity",starstroke=0.1)+
    scale_fill_manual(values=cols_phylum,
                      guide=guide_legend(keywidth = 0.5, 
                                         keyheight = 0.5, order=1,
                                         override.aes=list(starshape=15)),
                      na.translate=FALSE)+
    scale_starshape_manual(values=c(12, 8, 4, 1),
                           guide=guide_legend(keywidth = 0.5, 
                                              keyheight = 0.5, order=4),
                           na.translate=FALSE) +
    scale_size_continuous(range = c(1, 2.5),
                          guide = guide_legend(keywidth = 0.5, 
                                               keyheight = 0.5, order=3,
                                               override.aes=list(starshape=15)))
  
  p2
  p1 <-p2
  
  p1 <- p1 + new_scale_fill() +
    geom_fruit(data=dat2_l,
               geom=geom_tile,
               mapping=aes(y=ID,
                           x=Sites, 
                           alpha=Abundance,
                           fill=Sites),
               color = "grey50", offset = 0.28, size = 0.09,## Offset QUE TAN LEJOS ESTA EL HEATMAP 0.4 0.02
               axis.params = list(axis=unique(dat2_l[i,2]), ## etiquetas pal_sel
                                  text.angle = 0, text.size = 1.85, hjust=0.5), pwidth = pwidths[contador])+ # pwidth : ancho de celda
    scale_alpha_continuous(range=c(0, 1),                                                   # hjust mueve las etiquetas del heatmap
                           guide=guide_legend(keywidth = 0.3, 
                                              keyheight = 0.3, order=5)) +
    geom_fruit(data=dat3_l, geom=geom_bar,
               mapping=aes(y=ID, x=HigherAbundance, fill=Sites),
               pwidth=0.38, 
               orientation="y", 
               stat="identity",
               axis.params=list(
                 axis       = "x",
                 text.size  = 1.8,
                 hjust      = 1,
                 vjust      = 0.5,
                 nbreak     = 3,
               ),
               grid.params=list(),offset = 0.05) +
    scale_fill_manual(values=cols_l,
                      guide=guide_legend(keywidth = 0.3, 
                                         keyheight = 0.3, order=4))+
    geom_treescale(fontsize=2, linesize=0.3, x=0.001, y=0.001) + # X es que tan largas son las ramas
    theme(legend.position=c(0.93, 0.5),
          legend.background=element_rect(fill=NA),
          legend.title=element_text(size=6.5),
          legend.text=element_text(size=4.5),
          legend.spacing.y = unit(0.02, "cm"),
    )
  p1
  p4=p1
  p5=p4 + theme_tree()
  ggsave(paste0("pico_",l,"_filogenia.pdf"),p5, width=6, height=6, units="in", scale=3)
  contador = contador+1
}
####################
## SEL32
p1

dat2_32<-read.table(file="pico_df_sel32_ring_heatmap.txt", header=TRUE, sep="\t")
dat3_32<-read.table(file="pico_df_sel32_barplot_attr.txt", header=TRUE, sep="\t")

etiquetas_32 <- c()
etiquetas_32$tip_point <- sort(dat1[,1]) 
etiquetas_32$ring_heatmap <- sort(unique(dat2_32[,1]))
etiquetas_32$barplot <- sort(dat3_32[,1])
etiquetas_32$arbol <- sort(tree[["tip.label"]])
write.table(etiquetas_32, file="test_etiquetas.tsv", sep="\t", row.names = F)

###
## palindromos
pals_32 <- c()
for (i in 1:nrow(dat2_32)) {
  pals_32 <- append(pals_32,dat2_32[i,2])
}
pals_32 <- unique(pals_32)
pals_32
## colores para los palindromos
cols_32 <- rainbow(length(pals_32))
###

# adjust the order
dat2_32$Sites <- factor(dat2_32$Sites, 
                       levels=pals_32)
dat3_32$Sites <- factor(dat3_32$Sites, 
                       levels=pals_32)
# extract the clade label information. Because some nodes of tree are
# annotated to genera, which can be displayed with high light using ggtree.
# Extraer la información de la etiqueta del clado.
# Debido a que algunos nodos del árbol están anotados en los géneros,
# que se pueden mostrar con luz alta usando ggtree.
nodeids <- nodeid(tree, tree$node.label[nchar(tree$node.label)>4])
nodedf <- data.frame(node=nodeids)
nodelab <- gsub("[\\.0-9]", "", tree$node.label[nchar(tree$node.label)>4])

# The layers of clade and hightlight
poslist <- c(0.1)
labdf <- data.frame(node=nodeids, label=nodelab, pos=poslist)

# The circular layout tree.
p <- ggtree(tree, size=0.1, open.angle=5)+ # layout="fan" #, open.angle=5**
  geom_tiplab(size=2, align=TRUE, offset=0.01) +
  geom_hilight(data=nodedf, mapping=aes(node=node),
               extendto=0.6, alpha=0.3, fill="grey", color="grey50",
               size=0.05) +
  geom_cladelab(data=labdf, 
                mapping=aes(node=node, 
                            label=label,
                            offset.text=pos),
                hjust=0.5,
                barsize=NA,
                horizontal=FALSE, 
                fontsize=1.4,
                fontface="italic"
  ) 
p
p2 <- p
p2 <- p2 + theme(plot.margin = unit(c(14,8,14,8), "mm"))
p2
#----
p2 <- p2 %<+% dat1 + geom_star(
  mapping=aes(fill=Phylum, starshape=Type, size=Size),
  position="identity",starstroke=0.1)+
  scale_fill_manual(values=cols_phylum_32,
                    guide=guide_legend(keywidth = 0.5, 
                                       keyheight = 0.5, order=1,
                                       override.aes=list(starshape=15)),
                    na.translate=FALSE)+
  scale_starshape_manual(values=c(12, 8, 4, 1),
                         guide=guide_legend(keywidth = 0.5, 
                                            keyheight = 0.5, order=4),
                         na.translate=FALSE) +
  scale_size_continuous(range = c(1, 2.5),
                        guide = guide_legend(keywidth = 0.5, 
                                             keyheight = 0.5, order=3,
                                             override.aes=list(starshape=15)))

p2
p1 <-p2

p1 <- p1 + new_scale_fill() +
  geom_fruit(data=dat2_32,
             geom=geom_tile,
             mapping=aes(y=ID,
                         x=Sites, 
                         alpha=Abundance,
                         fill=Sites),
             color = "grey50", offset = 0.28, size = 0.09,## Offset QUE TAN LEJOS ESTA EL HEATMAP 0.4 0.02
             axis.params = list(axis=unique(dat2_32[i,2]), ## etiquetas pal_sel
                                text.angle = 0, text.size = 1.85, hjust=0.5), pwidth = 0.8)+ # pwidth : ancho de celda
  scale_alpha_continuous(range=c(0, 1),                                                   # hjust mueve las etiquetas del heatmap
                         guide=guide_legend(keywidth = 0.3, 
                                            keyheight = 0.3, order=5)) +
  geom_fruit(data=dat3_32, geom=geom_bar,
             mapping=aes(y=ID, x=HigherAbundance, fill=Sites),
             pwidth=0.38, 
             orientation="y", 
             stat="identity",
             axis.params=list(
               axis       = "x",
               text.size  = 1.8,
               hjust      = 1,
               vjust      = 0.5,
               nbreak     = 3,
             ),
             grid.params=list(),offset = 0.05) +
  scale_fill_manual(values=cols_32,
                    guide=guide_legend(keywidth = 0.3, 
                                       keyheight = 0.3, order=4))+
  geom_treescale(fontsize=2, linesize=0.3, x=0.001, y=0.001) + # X es que tan largas son las ramas
  theme(legend.position=c(0.93, 0.5),
        legend.background=element_rect(fill=NA),
        legend.title=element_text(size=6.5),
        legend.text=element_text(size=4.5),
        legend.spacing.y = unit(0.02, "cm"),
  )
p1

p4=p1
p5=p4 + theme_tree() 
ggsave("sel_32.pdf",p5, width=6, height=6, units="in", scale=3)

