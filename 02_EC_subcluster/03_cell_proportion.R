library(monocle3)
library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(ggplot2)
library(dittoSeq)
library(dplyr)
library(tidyverse)
library(harmony)
library(cluster)
library(SeuratWrappers)
CAV_snRNA<- readRDS("./02_EC/CAV/EC_CAV_snRNA_add_gene_symbol.rds")
CAV_snATAC<- readRDS("./02_EC/CAV/EC_CAV_snATAC_add_peaks_fragment.rds")
# the proportion in differnent condition 

# treatment: CTRL,IZ,BZ,RZ,FZ

# for CAV_subtype
CAV_snRNA$major_labl<- factor(CAV_snRNA$major_labl,levels=c("CTRL","IZ","BZ","RZ","FZ"))
CAV_snATAC$major_labl<- factor(CAV_snATAC$major_labl,levels=c("CTRL","IZ","BZ","RZ","FZ"))
str(df)
pdf("/md01/nieyg/project/EC/human_MI/02_EC/CAV/CAV_detail_vs_major_labl_proportion.pdf",width=10,height=5)
df <- as.data.frame(CAV_snRNA@meta.data)
# major_labl
df_ct <- df %>%
    group_by(major_labl, CAV_detail) %>%
    summarise(counts = n()) %>%
    mutate(cell_proportion = counts / sum(counts))
df_ct$CAV_detail<- factor(df_ct$CAV_detail,levels=levels(CAV_snRNA$CAV_detail))
cols <- myUmapcolors
p1 <- ggplot(df_ct, aes(major_labl, cell_proportion, fill=CAV_detail)) + 
    geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = cols) +
    theme_cowplot() +
    xlab("CAV_snRNA") + ylab("") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
df <- as.data.frame(CAV_snATAC@meta.data)
# major_labl
df_ct <- df %>%
    group_by(major_labl, CAV_detail) %>%
    summarise(counts = n()) %>%
    mutate(cell_proportion = counts / sum(counts))
df_ct$CAV_detail<- factor(df_ct$CAV_detail,levels=levels(CAV_snATAC$CAV_detail))
cols <- myUmapcolors
p2 <- ggplot(df_ct, aes(major_labl, cell_proportion, fill=CAV_detail)) + 
    geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = cols) +
    theme_cowplot() +
    xlab("CAV_snATAC") + ylab("") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
p1|p2
dev.off()

pdf("/md01/nieyg/project/EC/human_MI/02_EC/CAV/CAV_vs_major_labl_proportion.pdf",width=10,height=5)
df <- as.data.frame(CAV_snRNA@meta.data)
# major_labl
df_ct <- df %>%
    group_by(major_labl, CAV) %>%
    summarise(counts = n()) %>%
    mutate(cell_proportion = counts / sum(counts))
cols <- myUmapcolors
df_ct$CAV<- factor(df_ct$CAV,levels=levels(CAV_snRNA$CAV))
p1 <- ggplot(df_ct, aes(major_labl, cell_proportion, fill=CAV)) + 
    geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = cols) +
    theme_cowplot() +
    xlab("CAV_snRNA") + ylab("") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
df <- as.data.frame(CAV_snATAC@meta.data)
# major_labl
df_ct <- df %>%
    group_by(major_labl, CAV) %>%
    summarise(counts = n()) %>%
    mutate(cell_proportion = counts / sum(counts))
cols <- myUmapcolors
df_ct$CAV<- factor(df_ct$CAV,levels=levels(CAV_snATAC$CAV))
p2 <- ggplot(df_ct, aes(major_labl, cell_proportion, fill=CAV)) + 
    geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = cols) +
    theme_cowplot() +
    xlab("CAV_snATAC") + ylab("") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
p1|p2
dev.off()


library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(ggplot2)
library(dittoSeq)
library(dplyr)
library(tidyverse)
library(harmony)
library(cluster)
CAV_snRNA<- readRDS("./02_EC/CAV/EC_CAV_snRNA_add_gene_symbol.rds")
CAV_snATAC<- readRDS("./02_EC/CAV/EC_CAV_snATAC_add_peaks_fragment.rds")
# add the metadata from clinical data 
clinical_data<- read.csv("/data/R02/nieyg/project/EC/human_MI/00_data/clinical_data.csv")
CAV_snRNA$age<- clinical_data[match(CAV_snRNA$donor_id,clinical_data$patient),5]
CAV_snRNA$days.after.infarction<- clinical_data[match(CAV_snRNA$donor_id,clinical_data$patient),]$days.after.infarction
CAV_snRNA$disease<- clinical_data[match(CAV_snRNA$donor_id,clinical_data$patient),]$disease
table(CAV_snRNA$days.after.infarction)
CAV_snRNA$days.after.infarction<- factor(CAV_snRNA$days.after.infarction,levels=c("control","2","3","4","5","6","11","31","40","45","62","88","101","153","166"))


CAV_snATAC$age<- clinical_data[match(CAV_snATAC$donor_id,clinical_data$patient),5]
CAV_snATAC$days.after.infarction<- clinical_data[match(CAV_snATAC$donor_id,clinical_data$patient),]$days.after.infarction
CAV_snATAC$disease<- clinical_data[match(CAV_snATAC$donor_id,clinical_data$patient),]$disease
table(CAV_snATAC$days.after.infarction)
CAV_snATAC$days.after.infarction<- factor(CAV_snATAC$days.after.infarction,levels=c("control","2","3","4","5","6","11","31","40","45","62","88","101","153","166"))

saveRDS(CAV_snATAC,"./02_EC/CAV/EC_CAV_snATAC_add_peaks_fragment.rds")
saveRDS(CAV_snRNA,"./02_EC/CAV/EC_CAV_snRNA_add_gene_symbol.rds")

library(cowplot)
myUmapcolors <- c(  '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', '#AA9A59', '#E63863', '#E39A35', 
         '#C1E6F3', '#6778AE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#625D9E', 
         '#68A180', '#3A6963', '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', 
         '#116530', '#678983', '#A19882', '#FFBCBC', '#24A19C', '#FF9A76', "#8DD3C7",
         "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", 
         "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#E41A1C", "#377EB8", "#4DAF4A", 
         "#FF7F00", "#FFFF33", "#A65628", "#F781BF")

# days.after.infarction of all region vs major CAV in snRNA and snATAC
Idents(CAV_snRNA)<- CAV_snRNA$major_labl
CAV_snRNA$region_dai<- paste(CAV_snRNA$major_labl,CAV_snRNA$days.after.infarction,sep="_dai")
CAV_snRNA$region_dai[which(CAV_snRNA$region_dai=="CTRL_daicontrol")]<-"CTRL"
CAV_snRNA$region_dai<- factor(CAV_snRNA$region_dai,levels=c("CTRL","IZ_dai2","IZ_dai4","IZ_dai5","IZ_dai6","IZ_dai11","IZ_dai45","BZ_dai2","BZ_dai5","BZ_dai31","RZ_dai2" ,"RZ_dai3" ,"RZ_dai5","RZ_dai31","FZ_dai40","FZ_dai62","FZ_dai88","FZ_dai101","FZ_dai153","FZ_dai166"))

Idents(CAV_snATAC)<- CAV_snATAC$major_labl
CAV_snATAC$region_dai<- paste(CAV_snATAC$major_labl,CAV_snATAC$days.after.infarction,sep="_dai")
CAV_snATAC$region_dai[which(CAV_snATAC$region_dai=="CTRL_daicontrol")]<-"CTRL"
CAV_snATAC$region_dai<- factor(CAV_snATAC$region_dai,levels=c("CTRL","IZ_dai2","IZ_dai4","IZ_dai5","IZ_dai6","IZ_dai11","IZ_dai45","BZ_dai2","BZ_dai5","BZ_dai31","RZ_dai2" ,"RZ_dai3" ,"RZ_dai5","RZ_dai31","FZ_dai40","FZ_dai62","FZ_dai88","FZ_dai101","FZ_dai153","FZ_dai166"))
  pdf("/md01/nieyg/project/EC/human_MI/02_EC/CAV/CAV_vs_dpi_proportion.pdf",width=10,height=5)
  df <- as.data.frame(CAV_snRNA@meta.data)
  # region_dai
  df_ct <- df %>%
      group_by(region_dai, CAV) %>%
      summarise(counts = n()) %>%
      mutate(cell_proportion = counts / sum(counts))
  cols <- myUmapcolors
  df_ct$CAV<- factor(df_ct$CAV,levels=levels(CAV_snRNA$CAV))
  p1 <- ggplot(df_ct, aes(region_dai, cell_proportion, fill=CAV)) + 
      geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
      scale_fill_manual(values = cols) +
      theme_cowplot() +
      xlab("snRNA") + ylab("") +
      theme(axis.text.x = element_text(angle = 60, hjust = 1))
  df <- as.data.frame(CAV_snATAC@meta.data)
  # region_dai
  df_ct <- df %>%
      group_by(region_dai, CAV) %>%
      summarise(counts = n()) %>%
      mutate(cell_proportion = counts / sum(counts))
  cols <- myUmapcolors
  df_ct$CAV<- factor(df_ct$CAV,levels=levels(CAV_snATAC$CAV))
  p2 <- ggplot(df_ct, aes(region_dai, cell_proportion, fill=CAV)) + 
      geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
      scale_fill_manual(values = cols) +
      theme_cowplot() +
      xlab("snATAC") + ylab("") +
      theme(axis.text.x = element_text(angle = 60, hjust = 1))
  p1
  p2
  dev.off()

Idents(CAV_snRNA)<- CAV_snRNA$major_labl
CAV_snRNA$region_dai<- paste(CAV_snRNA$major_labl,CAV_snRNA$days.after.infarction,sep="_dai")
CAV_snRNA$region_dai[which(CAV_snRNA$region_dai=="CTRL_daicontrol")]<-"CTRL"
CAV_snRNA$region_dai<- factor(CAV_snRNA$region_dai,levels=c("CTRL","IZ_dai2","IZ_dai4","IZ_dai5","IZ_dai6","IZ_dai11","IZ_dai45","BZ_dai2","BZ_dai5","BZ_dai31","RZ_dai2" ,"RZ_dai3" ,"RZ_dai5","RZ_dai31","FZ_dai40","FZ_dai62","FZ_dai88","FZ_dai101","FZ_dai153","FZ_dai166"))
df <- as.data.frame(CAV_snRNA@meta.data)
df_ct <- df %>%
      group_by(region_dai, CAV) %>%
      summarise(counts = n()) %>%
      mutate(cell_proportion = counts / sum(counts))
library(ggpubr)
library(ggplot2)
subtype<- levels(CAV_snRNA$CAV)
pdf("/md01/nieyg/project/EC/human_MI/02_EC/CAV/snRNA_CAV_vs_dpi_proportion_dot_line_plot.pdf",width=5,height=5)
for (i in subtype){
    print(i);
    df_tmp<- df_ct[df_ct$CAV==i,]
    dai0<- df_tmp[df_tmp$region_dai=="CTRL",] 
    df_tmp<- df_tmp[-which(df_tmp$region_dai=="CTRL"),];
    dai0_df<- data.frame(region_dai=c("IZ_dai0","BZ_dai0","FZ_dai0"),CAV=rep(dai0$CAV,3),counts=rep(dai0$counts,3),cell_proportion=rep(dai0$cell_proportion,3))
    df_tmp<- rbind(df_tmp,dai0_df)
    df_tmp$region<- gsub("_dai.*","",df_tmp$region_dai)
    df_tmp$days<- gsub(".Z_dai","",df_tmp$region_dai)
    df_tmp$days<- as.numeric(df_tmp$days);
    df_tmp<- df_tmp[-which(df_tmp$region=="FZ"),]
    p<- ggplot(data=df_tmp,aes(x=days,y=cell_proportion))+
           geom_line(aes(color=region))+
           geom_point(aes(color=region),size=5)+
           scale_color_manual(values = c("IZ"="#BE017C",
                                         "BZ"="#FF8A4E",
                                         "RZ"="#FFC324"))+
           labs(y="% cell proportion",
                x="days after injury")+
           theme_classic()+ggtitle(i)
    print(p)
}
dev.off()

pdf("/md01/nieyg/project/EC/human_MI/02_EC/CAV/snRNA_IZ_CAV_vs_dpi_proportion_dot_line_plot.pdf",width=8,height=4)

    df_tmp<-df_ct
    dai0<- df_tmp[df_tmp$region_dai=="CTRL",] 
    df_tmp<- df_tmp[-which(df_tmp$region_dai=="CTRL"),];
    dai0_df<- data.frame(region_dai=c("IZ_dai0","BZ_dai0","FZ_dai0"),CAV=rep(dai0$CAV,3),counts=rep(dai0$counts,3),cell_proportion=rep(dai0$cell_proportion,3))
    df_tmp<- rbind(df_tmp,dai0_df)
    df_tmp$region<- gsub("_dai.*","",df_tmp$region_dai)
    df_tmp$days<- gsub(".Z_dai","",df_tmp$region_dai)
    df_tmp$days<- as.numeric(df_tmp$days);
    df_tmp<- df_tmp[which(df_tmp$region=="IZ"),]
    ggplot(data=df_tmp,aes(x=days,y=cell_proportion))+
           geom_line(aes(color=CAV))+
           geom_point(aes(color=CAV),size=5)+
           scale_color_manual(values = c("Artery"=myUmapcolors[1],
                                         "C_A"=myUmapcolors[2],
                                         "Capillary"=myUmapcolors[3],
                                         "C_V"=myUmapcolors[4],
                                         "Venous"=myUmapcolors[5]))+
           labs(y="% cell proportion",
                x="days after injury")+
           theme_classic()+ggtitle("IZ region")

dev.off()


Idents(CAV_snATAC)<- CAV_snATAC$major_labl
CAV_snATAC$region_dai<- paste(CAV_snATAC$major_labl,CAV_snATAC$days.after.infarction,sep="_dai")
CAV_snATAC$region_dai[which(CAV_snATAC$region_dai=="CTRL_daicontrol")]<-"CTRL"
CAV_snATAC$region_dai<- factor(CAV_snATAC$region_dai,levels=c("CTRL","IZ_dai2","IZ_dai4","IZ_dai5","IZ_dai6","IZ_dai11","IZ_dai45","BZ_dai2","BZ_dai5","BZ_dai31","RZ_dai2" ,"RZ_dai3" ,"RZ_dai5","RZ_dai31","FZ_dai40","FZ_dai62","FZ_dai88","FZ_dai101","FZ_dai153","FZ_dai166"))
df <- as.data.frame(CAV_snATAC@meta.data)
df_ct <- df %>%
      group_by(region_dai, CAV) %>%
      summarise(counts = n()) %>%
      mutate(cell_proportion = counts / sum(counts))
library(ggpubr)
library(ggplot2)
subtype<- levels(CAV_snATAC$CAV)
pdf("/md01/nieyg/project/EC/human_MI/02_EC/CAV/snATAC_CAV_vs_dpi_proportion_dot_line_plot.pdf",width=5,height=5)
for (i in subtype){
    print(i);
    df_tmp<- df_ct[df_ct$CAV==i,]
    dai0<- df_tmp[df_tmp$region_dai=="CTRL",] 
    df_tmp<- df_tmp[-which(df_tmp$region_dai=="CTRL"),];
    dai0_df<- data.frame(region_dai=c("IZ_dai0","BZ_dai0","FZ_dai0"),CAV=rep(dai0$CAV,3),counts=rep(dai0$counts,3),cell_proportion=rep(dai0$cell_proportion,3))
    df_tmp<- rbind(df_tmp,dai0_df)
    df_tmp$region<- gsub("_dai.*","",df_tmp$region_dai)
    df_tmp$days<- gsub(".Z_dai","",df_tmp$region_dai)
    df_tmp$days<- as.numeric(df_tmp$days);
    df_tmp<- df_tmp[-which(df_tmp$region=="FZ"),]
    p<- ggplot(data=df_tmp,aes(x=days,y=cell_proportion))+
           geom_line(aes(color=region))+
           geom_point(aes(color=region),size=5)+
           scale_color_manual(values = c("IZ"="#BE017C",
                                         "BZ"="#FF8A4E",
                                         "RZ"="#FFC324"))+
           labs(y="% cell proportion",
                x="days after injury")+
           theme_classic()+ggtitle(i)
    print(p)
}
dev.off()

pdf("/md01/nieyg/project/EC/human_MI/02_EC/CAV/snATAC_IZ_CAV_vs_dpi_proportion_dot_line_plot.pdf",width=8,height=4)
    df_tmp<-df_ct
    dai0<- df_tmp[df_tmp$region_dai=="CTRL",] 
    df_tmp<- df_tmp[-which(df_tmp$region_dai=="CTRL"),];
    dai0_df<- data.frame(region_dai=c("IZ_dai0","BZ_dai0","FZ_dai0"),CAV=rep(dai0$CAV,3),counts=rep(dai0$counts,3),cell_proportion=rep(dai0$cell_proportion,3))
    df_tmp<- rbind(df_tmp,dai0_df)
    df_tmp$region<- gsub("_dai.*","",df_tmp$region_dai)
    df_tmp$days<- gsub(".Z_dai","",df_tmp$region_dai)
    df_tmp$days<- as.numeric(df_tmp$days);
    df_tmp<- df_tmp[which(df_tmp$region=="IZ"),]
    ggplot(data=df_tmp,aes(x=days,y=cell_proportion))+
           geom_line(aes(color=CAV))+
           geom_point(aes(color=CAV),size=5)+
           scale_color_manual(values = c("Artery"=myUmapcolors[1],
                                         "C_A"=myUmapcolors[2],
                                         "Capillary"=myUmapcolors[3],
                                         "C_V"=myUmapcolors[4],
                                         "Venous"=myUmapcolors[5]))+
           labs(y="% cell proportion",
                x="days after injury")+
           theme_classic()+ggtitle("IZ region")

dev.off()

  pdf("/md01/nieyg/project/EC/human_MI/02_EC/CAV/CAV_detail_vs_dpi_proportion.pdf",width=10,height=5)
  df <- as.data.frame(CAV_snRNA@meta.data)
  # region_dai
  df_ct <- df %>%
      group_by(region_dai, CAV_detail) %>%
      summarise(counts = n()) %>%
      mutate(cell_proportion = counts / sum(counts))
  cols <- myUmapcolors
  df_ct$CAV_detail<- factor(df_ct$CAV_detail,levels=levels(CAV_snRNA$CAV_detail))
  p1 <- ggplot(df_ct, aes(region_dai, cell_proportion, fill=CAV_detail)) + 
      geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
      scale_fill_manual(values = cols) +
      theme_cowplot() +
      xlab("snRNA") + ylab("") +
      theme(axis.text.x = element_text(angle = 60, hjust = 1))
  df <- as.data.frame(CAV_snATAC@meta.data)
  # region_dai
  df_ct <- df %>%
      group_by(region_dai, CAV_detail) %>%
      summarise(counts = n()) %>%
      mutate(cell_proportion = counts / sum(counts))
  cols <- myUmapcolors
  df_ct$CAV_detail<- factor(df_ct$CAV_detail,levels=levels(CAV_snATAC$CAV_detail))
  p2 <- ggplot(df_ct, aes(region_dai, cell_proportion, fill=CAV_detail)) + 
      geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
      scale_fill_manual(values = cols) +
      theme_cowplot() +
      xlab("snATAC") + ylab("") +
      theme(axis.text.x = element_text(angle = 60, hjust = 1))
  p1
  p2
  dev.off()
Idents(CAV_snRNA)<- CAV_snRNA$major_labl
CAV_snRNA$region_dai<- paste(CAV_snRNA$major_labl,CAV_snRNA$days.after.infarction,sep="_dai")
CAV_snRNA$region_dai[which(CAV_snRNA$region_dai=="CTRL_daicontrol")]<-"CTRL"
CAV_snRNA$region_dai[which(CAV_snRNA$region_dai%in%c("BZ_dai2","BZ_dai5","BZ_dai31"))]<-"BZ"
CAV_snRNA$region_dai[which(CAV_snRNA$region_dai%in%c("RZ_dai2" ,"RZ_dai3" ,"RZ_dai5","RZ_dai31"))]<-"RZ"
CAV_snRNA$region_dai[which(CAV_snRNA$region_dai%in%c("FZ_dai40","FZ_dai62","FZ_dai88","FZ_dai101","FZ_dai153","FZ_dai166"))]<-"FZ"
CAV_snRNA$region_dai<- factor(CAV_snRNA$region_dai,levels=c("CTRL","IZ_dai2","IZ_dai4","IZ_dai5","IZ_dai6","IZ_dai11","IZ_dai45","BZ","RZ","FZ"))

Idents(CAV_snATAC)<- CAV_snATAC$major_labl
CAV_snATAC$region_dai<- paste(CAV_snATAC$major_labl,CAV_snATAC$days.after.infarction,sep="_dai")
CAV_snATAC$region_dai[which(CAV_snATAC$region_dai=="CTRL_daicontrol")]<-"CTRL"
CAV_snATAC$region_dai[which(CAV_snATAC$region_dai%in%c("BZ_dai2","BZ_dai5","BZ_dai31"))]<-"BZ"
CAV_snATAC$region_dai[which(CAV_snATAC$region_dai%in%c("RZ_dai2" ,"RZ_dai3" ,"RZ_dai5","RZ_dai31"))]<-"RZ"
CAV_snATAC$region_dai[which(CAV_snATAC$region_dai%in%c("FZ_dai40","FZ_dai62","FZ_dai88","FZ_dai101","FZ_dai153","FZ_dai166"))]<-"FZ"
CAV_snATAC$region_dai<- factor(CAV_snATAC$region_dai,levels=c("CTRL","IZ_dai2","IZ_dai4","IZ_dai5","IZ_dai6","IZ_dai11","IZ_dai45","BZ","RZ","FZ"))
  pdf("/md01/nieyg/project/EC/human_MI/02_EC/CAV/CAV_vs_IZ_dpi_proportion.pdf",width=10,height=5)
  df <- as.data.frame(CAV_snRNA@meta.data)
  # region_dai
  df_ct <- df %>%
      group_by(region_dai, CAV) %>%
      summarise(counts = n()) %>%
      mutate(cell_proportion = counts / sum(counts))
  cols <- myUmapcolors
  df_ct$CAV<- factor(df_ct$CAV,levels=levels(CAV_snRNA$CAV))
  p1 <- ggplot(df_ct, aes(region_dai, cell_proportion, fill=CAV)) + 
      geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
      scale_fill_manual(values = cols) +
      theme_cowplot() +
      xlab("snRNA") + ylab("") +
      theme(axis.text.x = element_text(angle = 60, hjust = 1))
  df <- as.data.frame(CAV_snATAC@meta.data)
  # region_dai
  df_ct <- df %>%
      group_by(region_dai, CAV) %>%
      summarise(counts = n()) %>%
      mutate(cell_proportion = counts / sum(counts))
  cols <- myUmapcolors
  df_ct$CAV<- factor(df_ct$CAV,levels=levels(CAV_snATAC$CAV))
  p2 <- ggplot(df_ct, aes(region_dai, cell_proportion, fill=CAV)) + 
      geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
      scale_fill_manual(values = cols) +
      theme_cowplot() +
      xlab("snATAC") + ylab("") +
      theme(axis.text.x = element_text(angle = 60, hjust = 1))
  p1|p2
  dev.off()

  pdf("/md01/nieyg/project/EC/human_MI/02_EC/CAV/CAV_detail_vs_IZ_dpi_proportion.pdf",width=10,height=5)
  df <- as.data.frame(CAV_snRNA@meta.data)
  # region_dai
  df_ct <- df %>%
      group_by(region_dai, CAV_detail) %>%
      summarise(counts = n()) %>%
      mutate(cell_proportion = counts / sum(counts))
  cols <- myUmapcolors
  df_ct$CAV_detail<- factor(df_ct$CAV_detail,levels=levels(CAV_snRNA$CAV_detail))
  p1 <- ggplot(df_ct, aes(region_dai, cell_proportion, fill=CAV_detail)) + 
      geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
      scale_fill_manual(values = cols) +
      theme_cowplot() +
      xlab("snRNA") + ylab("") +
      theme(axis.text.x = element_text(angle = 60, hjust = 1))
  df <- as.data.frame(CAV_snATAC@meta.data)
  # region_dai
  df_ct <- df %>%
      group_by(region_dai, CAV_detail) %>%
      summarise(counts = n()) %>%
      mutate(cell_proportion = counts / sum(counts))
  cols <- myUmapcolors
  df_ct$CAV_detail<- factor(df_ct$CAV_detail,levels=levels(CAV_snATAC$CAV_detail))
  p2 <- ggplot(df_ct, aes(region_dai, cell_proportion, fill=CAV_detail)) + 
      geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
      scale_fill_manual(values = cols) +
      theme_cowplot() +
      xlab("snATAC") + ylab("") +
      theme(axis.text.x = element_text(angle = 60, hjust = 1))
  p1|p2
  dev.off()


# only the IZ region 
# dpi vs major CAV in snRNA and snATAC
IZ_snRNA<- subset(CAV_snRNA_rmCTRL,idents="IZ")
IZ_snATAC<- subset(CAV_snATAC_rmCTRL,idents="IZ")

pdf("/md01/nieyg/project/EC/human_MI/02_EC/CAV/IZ_CAV_vs_dpi_proportion.pdf",width=10,height=5)
df <- as.data.frame(IZ_snRNA@meta.data)
# dpi
df_ct <- df %>%
    group_by(dpi, CAV) %>%
    summarise(counts = n()) %>%
    mutate(cell_proportion = counts / sum(counts))
cols <- myUmapcolors
df_ct$CAV<- factor(df_ct$CAV,levels=levels(IZ_snRNA$CAV))
p1 <- ggplot(df_ct, aes(dpi, cell_proportion, fill=CAV)) + 
    geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = cols) +
    theme_cowplot() +
    xlab("IZ snRNA") + ylab("") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
df <- as.data.frame(IZ_snATAC@meta.data)
# dpi
df_ct <- df %>%
    group_by(dpi, CAV) %>%
    summarise(counts = n()) %>%
    mutate(cell_proportion = counts / sum(counts))
cols <- myUmapcolors
df_ct$CAV<- factor(df_ct$CAV,levels=levels(IZ_snATAC$CAV))
p2 <- ggplot(df_ct, aes(dpi, cell_proportion, fill=CAV)) + 
    geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = cols) +
    theme_cowplot() +
    xlab("IZ snATAC") + ylab("") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
p1|p2
dev.off()

# check the C_A_1 annotation and the P16 sample in snATAC
# the P16 sample distribution in Umap 
# snRNA 

Idents(CAV_snATAC)<- CAV_snATAC$donor_id;
pdf("./02_EC/CAV/P16_snATAC_CAV_Umap.pdf",width=6,height=5)
DimPlot(object = CAV_snATAC,cols.highlight=c("darkblue","darkred"),cells.highlight=WhichCells(CAV_snATAC,idents="P16"),reduction = "umap_harmony",cols=myUmapcolors,pt.size=0.7, group.by = "donor_id")
DimPlot(object = CAV_snATAC,reduction = "umap_harmony",cols=myUmapcolors,pt.size=0.7, group.by = "CAV")
DimPlot(object = CAV_snATAC,reduction = "umap_harmony",cols=myUmapcolors,pt.size=0.7, group.by = "CAV_detail")
DimPlot(object = CAV_snATAC,reduction = "umap_harmony",cols=myUmapcolors,pt.size=0.7, group.by = "donor_id")
dev.off()

# plot the sample cell number 
# snRNA 
CAV_snRNA$sample_name<- paste(CAV_snRNA$region_dai,CAV_snRNA$donor_id,sep="-")
cell_number=as.data.frame(table(CAV_snRNA$sample_name))
colnames(cell_number)<-c("sample","number")
level<- c("CTRL-P1","CTRL-P7","CTRL-P8","CTRL-P17",
"IZ_dai2-P3"  ,"IZ_dai2-P9","IZ_dai4-P16"  ,"IZ_dai5-P2"  ,
"IZ_dai6-P10","IZ_dai11-P15" ,"IZ_dai45-P13" ,
"BZ_dai2-P3"   ,"BZ_dai5-P2"  ,"BZ_dai31-P12" ,
"RZ_dai2-P3"   ,"RZ_dai2-P9"   ,"RZ_dai3-P6"  ,
"RZ_dai5-P2","RZ_dai31-P11","FZ_dai40-P19","FZ_dai62-P14",
"FZ_dai88-P4"  ,"FZ_dai101-P5","FZ_dai153-P18","FZ_dai166-P20")
cell_number$sample<- factor(cell_number$sample,levels=level)
cell_number$color<-c(myUmapcolors,myUmapcolors)[1:length(level)]
pdf("./02_EC/CAV/snRNA_sample_cellnumber.pdf",width=8,height=4)
p<-ggplot(data = cell_number, aes_string(x = "sample", y = "number", 
        fill = "sample")) +  xlab(" ") + ylab("# of cells") + 
        scale_fill_manual(values = cell_number$color) + 
        geom_bar( stat = "identity", width = 0.6) +
        theme_classic()+ggtitle("snRNA_sample_cellnumber")+
        theme(axis.text.x = element_text(angle = 270,hjust=0));
p
p+geom_text(aes(label = number), size = 3) 
dev.off();

# snATAC
CAV_snATAC$sample_name<- paste(CAV_snATAC$region_dai,CAV_snATAC$donor_id,sep="-")
cell_number=as.data.frame(table(CAV_snATAC$sample_name))
colnames(cell_number)<-c("sample","number")
level<- c("CTRL-P1","CTRL-P7","CTRL-P8","CTRL-P17",
"IZ_dai2-P3"  ,"IZ_dai2-P9","IZ_dai4-P16"  ,"IZ_dai5-P2"  ,
"IZ_dai6-P10","IZ_dai11-P15" ,"IZ_dai45-P13" ,
"BZ_dai2-P3"   ,"BZ_dai5-P2"  ,"BZ_dai31-P12" ,
"RZ_dai2-P3"   ,"RZ_dai2-P9"   ,"RZ_dai3-P6"  ,
"RZ_dai5-P2","RZ_dai31-P11","FZ_dai40-P19","FZ_dai62-P14",
"FZ_dai88-P4"  ,"FZ_dai101-P5","FZ_dai153-P18","FZ_dai166-P20")
cell_number$sample<- factor(cell_number$sample,levels=level)
cell_number$color<-c(myUmapcolors,myUmapcolors)[1:length(level)]
pdf("./02_EC/CAV/snATAC_sample_cellnumber.pdf",width=8,height=4)
p<-ggplot(data = cell_number, aes_string(x = "sample", y = "number", 
        fill = "sample")) +  xlab(" ") + ylab("# of cells") + 
        scale_fill_manual(values = cell_number$color) + 
        geom_bar( stat = "identity", width = 0.6) +
        theme_classic()+ggtitle("snATAC_sample_cellnumber")+
        theme(axis.text.x = element_text(angle =270,hjust=0));
p
p+geom_text(aes(label = number), size = 3) 
dev.off();



Idents(CAV_snATAC)<- CAV_snATAC$major_labl
CAV_snATAC$region_dai<- paste(CAV_snATAC$major_labl,CAV_snATAC$days.after.infarction,sep="_dai")
CAV_snATAC$region_dai[which(CAV_snATAC$region_dai=="CTRL_daicontrol")]<-"CTRL"
CAV_snATAC$region_dai[which(CAV_snATAC$region_dai%in%c("BZ_dai2","BZ_dai5","BZ_dai31"))]<-"BZ"
CAV_snATAC$region_dai[which(CAV_snATAC$region_dai%in%c("RZ_dai2" ,"RZ_dai3" ,"RZ_dai5","RZ_dai31"))]<-"RZ"
CAV_snATAC$region_dai[which(CAV_snATAC$region_dai%in%c("FZ_dai40","FZ_dai62","FZ_dai88","FZ_dai101","FZ_dai153","FZ_dai166"))]<-"FZ"
CAV_snATAC$region_dai[which(CAV_snATAC$region_dai%in%c("IZ_dai2"))]<-"IZ_early"
CAV_snATAC$region_dai[which(CAV_snATAC$region_dai%in%c("IZ_dai4","IZ_dai5","IZ_dai6"))]<-"IZ_middle"
CAV_snATAC$region_dai[which(CAV_snATAC$region_dai%in%c("IZ_dai11","IZ_dai45"))]<-"IZ_late"
CAV_snATAC$region_dai<- factor(CAV_snATAC$region_dai,levels=c("CTRL","IZ_early","IZ_middle","IZ_late","BZ","RZ","FZ"))

Idents(CAV_snRNA)<- CAV_snRNA$major_labl
CAV_snRNA$region_dai<- paste(CAV_snRNA$major_labl,CAV_snRNA$days.after.infarction,sep="_dai")
CAV_snRNA$region_dai[which(CAV_snRNA$region_dai=="CTRL_daicontrol")]<-"CTRL"
CAV_snRNA$region_dai[which(CAV_snRNA$region_dai%in%c("BZ_dai2","BZ_dai5","BZ_dai31"))]<-"BZ"
CAV_snRNA$region_dai[which(CAV_snRNA$region_dai%in%c("RZ_dai2" ,"RZ_dai3" ,"RZ_dai5","RZ_dai31"))]<-"RZ"
CAV_snRNA$region_dai[which(CAV_snRNA$region_dai%in%c("FZ_dai40","FZ_dai62","FZ_dai88","FZ_dai101","FZ_dai153","FZ_dai166"))]<-"FZ"
CAV_snRNA$region_dai[which(CAV_snRNA$region_dai%in%c("IZ_dai2"))]<-"IZ_early"
CAV_snRNA$region_dai[which(CAV_snRNA$region_dai%in%c("IZ_dai4","IZ_dai5","IZ_dai6"))]<-"IZ_middle"
CAV_snRNA$region_dai[which(CAV_snRNA$region_dai%in%c("IZ_dai11","IZ_dai45"))]<-"IZ_late"
CAV_snRNA$region_dai<- factor(CAV_snRNA$region_dai,levels=c("CTRL","IZ_early","IZ_middle","IZ_late","BZ","RZ","FZ"))

  pdf("/md01/nieyg/project/EC/human_MI/02_EC/CAV/CAV_vs_IZ_stages_proportion.pdf",width=10,height=5)
  df <- as.data.frame(CAV_snRNA@meta.data)
  # region_dai
  df_ct <- df %>%
      group_by(region_dai, CAV) %>%
      summarise(counts = n()) %>%
      mutate(cell_proportion = counts / sum(counts))
  cols <- myUmapcolors
  df_ct$CAV<- factor(df_ct$CAV,levels=levels(CAV_snRNA$CAV))
  p1 <- ggplot(df_ct, aes(region_dai, cell_proportion, fill=CAV)) + 
      geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
      scale_fill_manual(values = cols) +
      theme_cowplot() +
      xlab("snRNA") + ylab("") +
      theme(axis.text.x = element_text(angle = 60, hjust = 1))
  df <- as.data.frame(CAV_snATAC@meta.data)
  # region_dai
  df_ct <- df %>%
      group_by(region_dai, CAV) %>%
      summarise(counts = n()) %>%
      mutate(cell_proportion = counts / sum(counts))
  cols <- myUmapcolors
  df_ct$CAV<- factor(df_ct$CAV,levels=levels(CAV_snATAC$CAV))
  p2 <- ggplot(df_ct, aes(region_dai, cell_proportion, fill=CAV)) + 
      geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
      scale_fill_manual(values = cols) +
      theme_cowplot() +
      xlab("snATAC") + ylab("") +
      theme(axis.text.x = element_text(angle = 60, hjust = 1))
  p1|p2
  dev.off()


# cell number distribution of all celltype 
# IZ early middle late stage 

library(ggplot2)
library(reshape2)
library(cowplot)
snRNA<- readRDS("./01_all_celltype/all_celltype_snRNA_add_gene_symbol_change_celltype.rds")
snATAC<- readRDS("./01_all_celltype/all_celltype_snATAC_add_peaks_fragment.rds")
clinical_data<- read.csv("/data/R02/nieyg/project/EC/human_MI/00_data/clinical_data.csv")
snRNA$age<- clinical_data[match(snRNA$donor_id,clinical_data$patient),5]
snRNA$days.after.infarction<- clinical_data[match(snRNA$donor_id,clinical_data$patient),]$days.after.infarction
snRNA$disease<- clinical_data[match(snRNA$donor_id,clinical_data$patient),]$disease
table(snRNA$days.after.infarction)
snRNA$days.after.infarction<- factor(snRNA$days.after.infarction,levels=c("control","2","3","4","5","6","11","31","40","45","62","88","101","153","166"))

snRNA$region_dai<- paste(snRNA$major_labl,snRNA$days.after.infarction,sep="_dai")
snRNA$region_dai[which(snRNA$region_dai=="CTRL_daicontrol")]<-"CTRL"
snRNA$region_dai[which(snRNA$region_dai%in%c("BZ_dai2","BZ_dai5","BZ_dai31"))]<-"BZ"
snRNA$region_dai[which(snRNA$region_dai%in%c("RZ_dai2" ,"RZ_dai3" ,"RZ_dai5","RZ_dai31"))]<-"RZ"
snRNA$region_dai[which(snRNA$region_dai%in%c("FZ_dai40","FZ_dai62","FZ_dai88","FZ_dai101","FZ_dai153","FZ_dai166"))]<-"FZ"
snRNA$region_dai[which(snRNA$region_dai%in%c("IZ_dai2"))]<-"IZ_early"
snRNA$region_dai[which(snRNA$region_dai%in%c("IZ_dai11","IZ_dai4","IZ_dai5","IZ_dai6"))]<-"IZ_middle"
snRNA$region_dai[which(snRNA$region_dai%in%c("IZ_dai45"))]<-"IZ_late"
snRNA$region_dai<- factor(snRNA$region_dai,levels=c("CTRL","IZ_early","IZ_middle","IZ_late","BZ","RZ","FZ"))

snATAC$age<- clinical_data[match(snATAC$donor_id,clinical_data$patient),5]
snATAC$days.after.infarction<- clinical_data[match(snATAC$donor_id,clinical_data$patient),]$days.after.infarction
snATAC$disease<- clinical_data[match(snATAC$donor_id,clinical_data$patient),]$disease
table(snATAC$days.after.infarction)
snATAC$days.after.infarction<- factor(snATAC$days.after.infarction,levels=c("control","2","3","4","5","6","11","31","40","45","62","88","101","153","166"))

snATAC$region_dai<- paste(snATAC$major_labl,snATAC$days.after.infarction,sep="_dai")
snATAC$region_dai[which(snATAC$region_dai=="CTRL_daicontrol")]<-"CTRL"
snATAC$region_dai[which(snATAC$region_dai%in%c("BZ_dai2","BZ_dai5","BZ_dai31"))]<-"BZ"
snATAC$region_dai[which(snATAC$region_dai%in%c("RZ_dai2" ,"RZ_dai3" ,"RZ_dai5","RZ_dai31"))]<-"RZ"
snATAC$region_dai[which(snATAC$region_dai%in%c("FZ_dai40","FZ_dai62","FZ_dai88","FZ_dai101","FZ_dai153","FZ_dai166"))]<-"FZ"
snATAC$region_dai[which(snATAC$region_dai%in%c("IZ_dai2"))]<-"IZ_early"
snATAC$region_dai[which(snATAC$region_dai%in%c("IZ_dai11","IZ_dai4","IZ_dai5","IZ_dai6"))]<-"IZ_middle"
snATAC$region_dai[which(snATAC$region_dai%in%c("IZ_dai45"))]<-"IZ_late"
snATAC$region_dai<- factor(snATAC$region_dai,levels=c("CTRL","IZ_early","IZ_middle","IZ_late","BZ","RZ","FZ"))

  df <- as.data.frame(snRNA@meta.data)
  # region_dai
  df_ct <- df %>%
      group_by(region_dai, celltype) %>%
      summarise(counts = n()) %>%
      mutate(cell_proportion = counts / sum(counts))


pdf("/md01/nieyg/project/EC/human_MI/01_all_celltype/Allcelltype_cellnumber_vs_IZ_stages.pdf",width=10,height=5)
ggplot(df_ct, aes(x=region_dai, y=counts, fill=celltype)) +
  geom_bar(stat="identity", position=position_dodge(),
           color="black", width=.8) +
  scale_fill_manual(values = cols) +
  #scale_y_continuous(expand=c(0,0))+
  #coord_cartesian(ylim = c(0, 8))+
  theme(axis.text.x = element_text(size = 14, color = "black"))+##设置x轴字体大小
  theme(axis.text.y = element_text(size = 14, color = "black"))+##设置y轴字体大小
  theme(title=element_text(size=13))+#设置标题字体大小
  theme_bw()
  dev.off()

snATAC <- RenameIdents(
  object = snATAC,
  'neuronal receptor cell' = 'Neuronal',
  'cardiac muscle myoblast' = 'vCMs',
  'smooth muscle myoblast' = 'vSMCs',
  'pericyte' = 'Pericytes',
  'lymphoid lineage restricted progenitor cell' = 'Lymphoid',
  'immature innate lymphoid cell' = 'Myeloid',
  'fibroblast of cardiac tissue' = 'Fibroblasts',
  'cardiac endothelial cell' = 'Endothelial',
  'mast cell' = 'Mast',
  'adipocyte of epicardial fat of left ventricle' = 'Adipocytes',
  'unknown' = 'Cycling cells'
  )
snATAC@meta.data$celltype<-Idents(snATAC)
table(snATAC$celltype,snATAC$major_labl)

  pdf("/md01/nieyg/project/EC/human_MI/01_all_celltype/Allcelltype_vs_IZ_stages_proportion.pdf",width=10,height=5)
  df <- as.data.frame(snRNA@meta.data)
  # region_dai
  df_ct <- df %>%
      group_by(region_dai, celltype) %>%
      summarise(counts = n()) %>%
      mutate(cell_proportion = counts / sum(counts))
  cols <- myUmapcolors
  df_ct$celltype<- factor(df_ct$celltype,levels=levels(snRNA$celltype))
  p1 <- ggplot(df_ct, aes(region_dai, cell_proportion, fill= celltype)) + 
      geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
      scale_fill_manual(values = cols) +
      theme_cowplot() +
      xlab("snRNA") + ylab("") +
      theme(axis.text.x = element_text(angle = 60, hjust = 1))
  df <- as.data.frame(snATAC@meta.data)
  # region_dai
  df_ct <- df %>%
      group_by(region_dai, celltype) %>%
      summarise(counts = n()) %>%
      mutate(cell_proportion = counts / sum(counts))
  cols <- myUmapcolors
  df_ct$celltype<- factor(df_ct$celltype,levels=levels(snATAC$celltype))
  p2 <- ggplot(df_ct, aes(region_dai, cell_proportion, fill=celltype)) + 
      geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
      scale_fill_manual(values = cols) +
      theme_cowplot() +
      xlab("snATAC") + ylab("") +
      theme(axis.text.x = element_text(angle = 60, hjust = 1))
  p1|p2
  dev.off()




# cell number distribution of all celltype 
# IZ days

snRNA$region_dai<- paste(snRNA$major_labl,snRNA$days.after.infarction,sep="_dai")
snRNA$region_dai[which(snRNA$region_dai=="CTRL_daicontrol")]<-"CTRL"
snRNA$region_dai[which(snRNA$region_dai%in%c("BZ_dai2","BZ_dai5","BZ_dai31"))]<-"BZ"
snRNA$region_dai[which(snRNA$region_dai%in%c("RZ_dai2" ,"RZ_dai3" ,"RZ_dai5","RZ_dai31"))]<-"RZ"
snRNA$region_dai[which(snRNA$region_dai%in%c("FZ_dai40","FZ_dai62","FZ_dai88","FZ_dai101","FZ_dai153","FZ_dai166"))]<-"FZ"
snRNA$region_dai<- factor(snRNA$region_dai,levels=c("CTRL","IZ_dai2","IZ_dai4","IZ_dai5","IZ_dai6","IZ_dai11","IZ_dai45","BZ","RZ","FZ"))

snATAC$region_dai<- paste(snATAC$major_labl,snATAC$days.after.infarction,sep="_dai")
snATAC$region_dai[which(snATAC$region_dai=="CTRL_daicontrol")]<-"CTRL"
snATAC$region_dai[which(snATAC$region_dai%in%c("BZ_dai2","BZ_dai5","BZ_dai31"))]<-"BZ"
snATAC$region_dai[which(snATAC$region_dai%in%c("RZ_dai2" ,"RZ_dai3" ,"RZ_dai5","RZ_dai31"))]<-"RZ"
snATAC$region_dai[which(snATAC$region_dai%in%c("FZ_dai40","FZ_dai62","FZ_dai88","FZ_dai101","FZ_dai153","FZ_dai166"))]<-"FZ"
snATAC$region_dai<- factor(snATAC$region_dai,levels=c("CTRL","IZ_dai2","IZ_dai4","IZ_dai5","IZ_dai6","IZ_dai11","IZ_dai45","BZ","RZ","FZ"))

  df <- as.data.frame(snRNA@meta.data)
  # region_dai
  df_ct <- df %>%
      group_by(region_dai, celltype) %>%
      summarise(counts = n()) %>%
      mutate(cell_proportion = counts / sum(counts))


pdf("/md01/nieyg/project/EC/human_MI/01_all_celltype/Allcelltype_cellnumber_vs_IZ_days.pdf",width=10,height=5)
ggplot(df_ct, aes(x=region_dai, y=counts, fill=celltype)) +
  geom_bar(stat="identity", position=position_dodge(),
           color="black", width=.8) +
  scale_fill_manual(values = cols) +
  #scale_y_continuous(expand=c(0,0))+
  #coord_cartesian(ylim = c(0, 8))+
  theme(axis.text.x = element_text(size = 14, color = "black"))+##设置x轴字体大小
  theme(axis.text.y = element_text(size = 14, color = "black"))+##设置y轴字体大小
  theme(title=element_text(size=13))+#设置标题字体大小
  theme_bw()
  dev.off()

  pdf("/md01/nieyg/project/EC/human_MI/01_all_celltype/Allcelltype_vs_IZ_days_proportion.pdf",width=10,height=5)
  df <- as.data.frame(snRNA@meta.data)
  # region_dai
  df_ct <- df %>%
      group_by(region_dai, celltype) %>%
      summarise(counts = n()) %>%
      mutate(cell_proportion = counts / sum(counts))
  cols <- myUmapcolors
  df_ct$celltype<- factor(df_ct$celltype,levels=levels(snRNA$celltype))
  p1 <- ggplot(df_ct, aes(region_dai, cell_proportion, fill= celltype)) + 
      geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
      scale_fill_manual(values = cols) +
      theme_cowplot() +
      xlab("snRNA") + ylab("") +
      theme(axis.text.x = element_text(angle = 60, hjust = 1))
  df <- as.data.frame(snATAC@meta.data)
  # region_dai
  df_ct <- df %>%
      group_by(region_dai, celltype) %>%
      summarise(counts = n()) %>%
      mutate(cell_proportion = counts / sum(counts))
  cols <- myUmapcolors
  df_ct$celltype<- factor(df_ct$celltype,levels=levels(snATAC$celltype))
  p2 <- ggplot(df_ct, aes(region_dai, cell_proportion, fill=celltype)) + 
      geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
      scale_fill_manual(values = cols) +
      theme_cowplot() +
      xlab("snATAC") + ylab("") +
      theme(axis.text.x = element_text(angle = 60, hjust = 1))
  p1|p2
  dev.off()
