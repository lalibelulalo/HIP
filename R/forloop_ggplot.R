library(tidyverse)

#Create dataframes(In this example n = 3)
df_1 <- data.frame(a1 = 1:5,
                   b1 = 1:5)  
df_2 <- data.frame(a1 = 1:10,
                   b1 = 1:10)
df_3 <- data.frame(a1 = 1:15,
                   b1 = 1:15)
df_4 <- data.frame(a1 = 1:20,
                   b1 = 1:20)

##Store dataframes in list
example.list<-lapply(1:4, function(x) eval(parse(text=paste0("df_", x)))) #In order to store all datasets in one list using their name
names(example.list)<-lapply(1:4, function(x) paste0("df_", x))

#Graph and save for each dataframe

for (i in 1:length(example.list)){
  df_i <- example.list[[i]]
  benp <-  
    df_i %>%
    ggplot(aes(x=b1)) + 
    geom_histogram(fill="steelblue", aes(y=..density.., alpha=..count..), bins=60) + 
    labs(title="Beneficios", subtitle="") + ylab("Densidad") + 
    xlab("Beneficios ($millones)") + 
    geom_vline(aes(xintercept=mean(b1)), color="red4",linetype="dashed") +
    theme(legend.position = "none") + 
    annotate("text", x= mean(df_i$b1), y=0, label=round(mean(df_i$b1), digits = 2), 
             colour="red4", size=3.5, vjust=-1.5, hjust=-0.5) 
  ggsave(benp, file=paste0(names(example.list)[i],"_histogram.png"))
}


