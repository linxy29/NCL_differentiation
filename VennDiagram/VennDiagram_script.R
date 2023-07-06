install.packages("ggVennDiagram")
library(ggVennDiagram)
library(ggplot2)

#upload data
vivo<-read.delim("genelist_vivo.txt",header=TRUE)
vitro<-read.delim("genelist_vitro.txt",header=TRUE)
mouse<-read.delim("genelist_mouse.txt",header=TRUE)

#make list
vivo<-as.list(vivo)
vitro<-as.list(vitro)
mouse <- as.list(mouse)
whole<-c(vivo,vitro,mouse)

ggVennDiagram(whole,label_alpha = 0)
t<-process_region_data(Venn(whole))
write.csv(t$item[[1]],"only_vivo.csv")
write.csv(t$item[[2]],"only_vitro.csv")
write.csv(t$item[[3]],"only_mouse.csv")
write.csv(t$item[[4]],"vivo_intersect_vitro.csv")
write.csv(t$item[[5]],"vivo_intersect_mouse.csv")
write.csv(t$item[[6]],"vitro_intersect_mouse.csv")
write.csv(t$item[[7]],"common_to_all.csv")
