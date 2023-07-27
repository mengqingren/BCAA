setwd("K:/TCGA/Anlysis/BCAA")
library(tidyverse)
######
######## Confirm genes => 18 #########
library(GSVA)
library(GSEABase)
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(limma)
library(tidyverse)
Pathway <- read.gmt("../../msigdb/c2.all.v7.5.1.symbols.gmt")
BCAA.SYMBOLs <- Pathway[str_detect(Pathway$term,"BRANCHED") | str_detect(Pathway$term,"VALINE"),]
#### REACTOME 21 GENES ####
BCAA.SYMBOLs[str_detect(BCAA.SYMBOLs$term,"REACTOME"),] %>% dim()
#### KEGG 55 GENES ####
BCAA.SYMBOLs[str_detect(BCAA.SYMBOLs$term,"KEGG"),] %>% dim()
#### Total 56 GENES ####
#### Venn ####
dat<-c("REACTOME" = 21, 
       "KEGG" = 55, 
       "REACTOME&KEGG" = 18)

library(Vennerable)
x<-(BCAA.SYMBOLs[str_detect(BCAA.SYMBOLs$term,"REACTOME"),])$gene
y<-(BCAA.SYMBOLs[str_detect(BCAA.SYMBOLs$term,"KEGG"),])$gene
Figure.D<-Venn(list("REACTOME"=x,"KEGG"=y))
plot(Figure.D,doWeight=T)
#### 18 GENES ####
Symbols <- intersect(x,y)
GeneType <- read.table("../../../../Metastasis/TissueMicrobiome/Human.GRC38.GeneType.txt",header = T,sep = ",")
GeneType.Unique <- GeneType %>% dplyr::select(Gene.stable.ID,Gene.name,NCBI.gene..formerly.Entrezgene..ID) %>% unique()
data <- GeneType.Unique %>% filter(Gene.name %in% Symbols) %>%
  magrittr::set_colnames(c("ENSEMBL","SYMBOL","ENTREZ"))
write.csv(data,"BCAA.ID.csv")

######## CNV #########
BCAA.ID <- read.csv("BCAA.ID.csv")
library(readxl)
CNV.Data <- read_xlsx("GSCA/CnvSummaryTable.xlsx")
CNV.Data <- CNV.Data %>%
  mutate(cancertype=factor(cancertype,levels = rev(sort(unique(.$cancertype)))))
Amp.Figure <- ggplot(CNV.Data,aes(symbol,cancertype,fill=a_total))+
  geom_tile(color="black")+
  geom_text(aes(label=round(a_total,1)),size=2)+
  ggthemes::theme_few()+
  theme(#legend.position = "bottom",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9,angle = 90,vjust = 0.5,hjust = 1),
        axis.text.y = element_text(size=9),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="",y="",title = "Amplication")+
  scale_fill_gradientn("Percentage(%)",colors = brewer.pal(11,"RdBu")[6:1])+
  theme(plot.title = element_text(hjust = 0.5))

ggsave(Amp.Figure,filename = "CNV.Amp.Figure.pdf",width = 5.5,height = 7)

Del.Figure <- ggplot(CNV.Data,aes(symbol,cancertype,fill=d_total))+
  geom_tile(color="black")+
  geom_text(aes(label=round(d_total,1)),size=2)+
  ggthemes::theme_few()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9,angle = 90,vjust = 0.5,hjust = 1),
    axis.text.y = element_text(size=9),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="",y="",title = "Deletion")+
  scale_fill_gradientn("Percentage(%)",colors = brewer.pal(11,"RdBu")[6:11])+
  theme(plot.title = element_text(hjust = 0.5))

ggsave(Del.Figure,filename = "CNV.Del.Figure.pdf",width = 5.5,height = 7)

######## CNV + Expression #########
library(tidyverse)
library(readxl)
library(RColorBrewer)
C.Data <- read_xlsx("GSCA/CnvAndExpressionTable.xlsx")
C.Data2 <- C.Data %>%
  mutate(symbol=factor(symbol,levels = rev(sort(unique(.$symbol)))))

p1031.1 <- ggplot()+
  geom_tile(data=C.Data2,
            aes(cancertype,symbol,fill=spm),
            color="white")+
  geom_point(data=C.Data2 %>% filter(fdr > 0.05),
             aes(cancertype,symbol,size=-log10(fdr)),
             shape=4,color="black")+
  geom_point(data=C.Data2 %>% filter(fdr <= 0.05),
             aes(cancertype,symbol,color=spm,size=-log10(fdr)),
             shape=1,color="black")+
  scale_size_continuous("-log10(FDR)",range = c(0.5,4))+
  ggthemes::theme_few()+
  theme(legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9,angle = 90,vjust = 0.5,hjust = 1),
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="",y="")+
  scale_fill_gradient2("Spearman'CC",
                       high = "#E64B35FF", #brewer.pal(11,"RdBu")[1],
                       low = "#053061", #brewer.pal(11,"RdBu")[11],
                       mid = "white",
                       breaks = c(-0.4,0,0.4,0.8))

data3 <- C.Data2 %>% 
  mutate(Label=if_else(fdr>0.05,"P>0.05",if_else(spm>0,"Positive","Negative"))) %>%
  group_by(symbol,Label) %>%
  dplyr::summarise(Count=n()) %>%
  #mutate(symbol=factor(symbol,levels = Gene.Symble)) %>%
  ungroup() %>%
  mutate(Label=factor(Label,levels = c("P>0.05","Positive","Negative")))

p1031.top <- ggplot(data3,aes(Count,symbol,fill=Label))+
  geom_col()+
  ggthemes::theme_few()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size=9),
    axis.ticks.y = element_blank(),
    axis.title = element_text(size=10),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="Count",y="")+
  scale_fill_manual("Spearman'CC",values = c("grey50",brewer.pal(11,"RdBu")[c(2,10)]))

p1031.F <- p1031.1 %>% aplot::insert_right(p1031.top,width = 0.1)
ggsave(p1031.F,filename = "CNV.Expression.Spearman.Figure.pdf",width = 10,height = 3)

######## SNV #########
BCAA.ID <- read.csv("BCAA.ID.csv")
library(maftools)
library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(data.table)
#### Mutation Freq ####
SNP.Data <- read_xlsx("GSCA/SnvSummaryTable.xlsx")
SNP.Data <- SNP.Data %>%
  mutate(Count=EffectiveMut+NonEffectiveMut) %>%
  mutate(cancertype=paste(cancertype,"(n=",sample_size,")",sep = "")) %>%
  dplyr::select(cancertype,symbol,Count) %>%
  spread(symbol,Count,fill=0) %>%
  gather(symbol,Count,-cancertype) %>%
  mutate(SampleSize=str_remove_all(cancertype,".*=") %>% str_remove_all("\\)")) %>%
  mutate(SampleSize=as.numeric(as.character(SampleSize))) %>%
  mutate(Percentage=round(Count/SampleSize*100,2)) %>%
  mutate(cancertype=factor(cancertype,levels = rev(sort(unique(.$cancertype)))))
write.csv(SNP.Data,file = "SNP.Frequency.csv",row.names = F)

SNP.Figure <- ggplot(SNP.Data,aes(symbol,cancertype,fill=Percentage))+
  geom_tile(color="black")+
  geom_text(aes(label=round(Count,1)),size=2)+
  ggthemes::theme_few()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9,angle = 90,vjust = 0.5,hjust = 1),
    axis.text.y = element_text(size=9),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="",y="",title = "Mutation")+
  scale_fill_gradientn("Mutation frequency(%)",colors = brewer.pal(11,"BrBG")[6:11])+
  theme(plot.title = element_text(hjust = 0.5))

ggsave(SNP.Figure,filename = "SNP.Frequency.Figure.pdf",width = 6,height = 7)

#### Oncoplot => UCEC ####
SNP.Data <- read.csv("../AnimoAcidSensing/TCGA.PanCancer.Symbol.MUTECT2.SNP.MAF.csv")
SNP.Clinical <- read.csv("../AnimoAcidSensing/TCGA.PanCancer.Symbol.MUTECT2.SNP.MAF.Clinical.csv")

UCEC.Clinical <- SNP.Clinical %>% filter(Cancer=="UCEC")
UCEC.SNP <- SNP.Data %>% filter(Tumor_Sample_Barcode %in% UCEC.Clinical$Tumor_Sample_Barcode)

UCEC.SNP2 <- UCEC.SNP %>% mutate(Variant_Classification = 
                                   if_else(Variant_Classification == "inframe_insertion","In_Frame_Ins",
                                           if_else(Variant_Classification == "inframe_deletion","In_Frame_Del",
                                                   if_else(Variant_Classification == "frameshift_variant","Frame_Shift_Var",
                                                           if_else(Variant_Classification == "missense_variant","Missense_Mutation",
                                                                   if_else(Variant_Classification == "synonymous_variant","Nonsense_Mutation",
                                                                           if_else(Variant_Classification == "3_prime_UTR_variant" | 
                                                                                     Variant_Classification == "5_prime_UTR_variant","Prime_UTR_Var",
                                                                                   if_else(Variant_Classification == "intron_variant","Intron_Var",
                                                                                           if_else(str_detect(Variant_Classification,";"),"Multi_Hit",
                                                                                                   if_else(str_detect(Variant_Classification,"splice"),"Splice_Site","Multi_Hit"))))))))))
library(maftools)
library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(data.table)
tcga = read.maf(maf = UCEC.SNP2,clinicalData=UCEC.Clinical,
                vc_nonSyn=c("In_Frame_Ins",
                            "In_Frame_Del",
                            "Frame_Shift_Var",
                            "Missense_Mutation",
                            "Nonsense_Mutation",
                            "Prime_UTR_Var",
                            "Intron_Var",
                            "Splice_Site"
                )) # 
library(RColorBrewer)
vc_cols = brewer.pal(9,"Paired")
names(vc_cols) <- c("In_Frame_Ins",
                    "Missense_Mutation",
                    "In_Frame_Del",
                    "Frame_Shift_Var",
                    "Nonsense_Mutation",
                    "Prime_UTR_Var",
                    "Intron_Var",
                    "Splice_Site",
                    "Multi_Hit"
)

SYMBOLs <- (BCAA.ID[BCAA.ID$SYMBOL %in% UCEC.SNP2$Hugo_Symbol,])$SYMBOL

maftools::oncoplot(maf = tcga,genes = SYMBOLs,bgCol = "white",
                   colors = vc_cols,
                   additionalFeatureCol = "white",
                   showTitle =F,
                   annoBorderCol = "white",
                   clinicalFeatures=c("OS"))

#### Oncoplot => SKCM ####
SKCM.Clinical <- SNP.Clinical %>% filter(Cancer=="SKCM")
SKCM.SNP <- SNP.Data %>% filter(Tumor_Sample_Barcode %in% SKCM.Clinical$Tumor_Sample_Barcode)

SKCM.SNP2 <- SKCM.SNP %>% mutate(Variant_Classification = 
                                   if_else(Variant_Classification == "inframe_insertion","In_Frame_Ins",
                                           if_else(Variant_Classification == "inframe_deletion","In_Frame_Del",
                                                   if_else(Variant_Classification == "frameshift_variant","Frame_Shift_Var",
                                                           if_else(Variant_Classification == "missense_variant","Missense_Mutation",
                                                                   if_else(Variant_Classification == "synonymous_variant","Nonsense_Mutation",
                                                                           if_else(Variant_Classification == "3_prime_UTR_variant" | 
                                                                                     Variant_Classification == "5_prime_UTR_variant","Prime_UTR_Var",
                                                                                   if_else(Variant_Classification == "intron_variant","Intron_Var",
                                                                                           if_else(str_detect(Variant_Classification,";"),"Multi_Hit",
                                                                                                   if_else(str_detect(Variant_Classification,"splice"),"Splice_Site","Multi_Hit"))))))))))
library(maftools)
library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(data.table)
tcga = read.maf(maf = SKCM.SNP2,clinicalData=SKCM.Clinical,
                vc_nonSyn=c("In_Frame_Ins",
                            "In_Frame_Del",
                            "Frame_Shift_Var",
                            "Missense_Mutation",
                            "Nonsense_Mutation",
                            "Prime_UTR_Var",
                            "Intron_Var",
                            "Splice_Site"
                )) # 
library(RColorBrewer)
vc_cols = brewer.pal(9,"Paired")
names(vc_cols) <- c("In_Frame_Ins",
                    "Missense_Mutation",
                    "In_Frame_Del",
                    "Frame_Shift_Var",
                    "Nonsense_Mutation",
                    "Prime_UTR_Var",
                    "Intron_Var",
                    "Splice_Site",
                    "Multi_Hit"
)

SYMBOLs <- (BCAA.ID[BCAA.ID$SYMBOL %in% UCEC.SNP2$Hugo_Symbol,])$SYMBOL

maftools::oncoplot(maf = tcga,genes = SYMBOLs,bgCol = "white",
                   colors = vc_cols,
                   additionalFeatureCol = "white",
                   showTitle =F,
                   annoBorderCol = "white",
                   clinicalFeatures=c("OS"))

######## Diagnosis #######
#### Logistic model => ≥5 Normal ####
library(pROC)
library(ROCR)
library(caret)
library(tidyverse)
BCAA.ID <- read.csv("BCAA.ID.csv")
PhenoType <- read.csv("../../mRNA.Phenotype.csv")
PhenoType <- PhenoType %>% filter(Sample.Type == "Primary Tumor" | 
                                    Sample.Type == "Solid Tissue Normal" | 
                                    Sample.Type == "Primary Blood Derived Cancer - Peripheral Blood") %>% 
  dplyr::select(Project.ID,Sample.ID,Sample.Type) %>%
  mutate(Type=ifelse(Sample.Type=="Solid Tissue Normal","Normal","Tumor")) %>%
  mutate(Cancer=str_remove(Project.ID,"TCGA-"))

#files <- list.files("../../mRNA.PanCancer.Exp/")
files <- list.files("../Pancancer/",pattern = "[DEG]",full.names = T)
files <- files[str_detect(files,"DEG.csv")]

dir = "../../Clinical.XENA/"
Dia.AUC <- data.frame()

for (file in files) {
  project <- str_remove(file,".*/") %>% str_remove_all("_.*") %>% str_remove_all(".*-")
  phenotype <- read.csv(paste(dir,project,"/TCGA-",project,".GDC_phenotype.tsv",sep = ""),sep = "\t",header = T) %>% 
    dplyr::filter(sample_type.samples %in% c("Primary Tumor","Primary Blood Derived Cancer - Peripheral Blood","Solid Tissue Normal"))
  TPM <- read.csv(paste("../../mRNA.PanCancer.Exp/TCGA-",project,".mRNA.TPM.csv",sep = ""),check.names = F,row.names = 1)
  TPM <- log2(TPM+1)
  
  Inter.Samples <- intersect(phenotype$submitter_id.samples,colnames(TPM))
  Pheno.mid <- phenotype %>% filter(submitter_id.samples %in% Inter.Samples) %>%
    mutate(Label = ifelse(sample_type.samples %in% c("Primary Tumor","Primary Blood Derived Cancer - Peripheral Blood"),1,0)) %>%
    dplyr::select(submitter_id.samples,Label)
  TPM.Target <- TPM[BCAA.ID$ENSEMBL,Pheno.mid$submitter_id.samples] %>% t() %>% data.frame(check.names = F)
  colnames(TPM.Target) <- BCAA.ID$SYMBOL
  TPM.Target$Label = Pheno.mid$Label
  if (sum(TPM.Target$Label == 0) >= 10) {
    Input.Formula <- as.formula(paste("Label ~ ",paste0(BCAA.ID$SYMBOL,collapse = "+"),sep = ""))
    TPM.Target$Label <- factor(TPM.Target$Label,levels = c(0,1),labels = c("Normal","Tumor"))
    
    for (i in 1:500) {
      print(project)
      print(i)
      index <- createDataPartition(TPM.Target$Label, p=0.7, list = F)
      train <- TPM.Target[index, ]
      test <- TPM.Target[-index, ]
      
      control = trainControl(method = "cv",number = 2,classProbs = TRUE)
      
      glm.model = train(Input.Formula,
                        data= train,
                        method = "glm",
                        metric = "ROC",
                        trControl = control)
      glm.probs = predict(glm.model,test,type = "prob")
      glm.ROC = roc(response = test[,c("Label")],
                    predictor = glm.probs$Tumor,
                    levels = levels(test[,c("Label")]))
      AUC <- glm.ROC$auc %>% str_remove(".* ") %>% as.numeric() %>% round(.,4)
      Dia.AUC <- data.frame(AUC = AUC, Repeat = i) %>% 
        mutate(Cancer = project) %>%
        rbind.data.frame(Dia.AUC)
    }
  }else{
    print(project)
  }
}

write.csv(Dia.AUC,file = "Diagnosis.BCAA.csv",row.names = T)

setwd("K:/TCGA/Anlysis/BCAA")
library(tidyverse)
data <- read.csv("Diagnosis.BCAA.csv",row.names = 1) %>%
  group_by(Cancer) %>%
  mutate(AUC.Mean=mean(AUC),
            AUC.Median=median(AUC),
            AUC.SD=sd(AUC)) %>%
  dplyr::select(-Repeat,-AUC) %>%
  unique()
read.csv("Diagnosis.BCAA.csv",row.names = 1) %>% 
  write.csv(file="TableS1.csv")
write.csv(data,file="TableS1.Sup.csv")

#### Figure ####
library(tidyverse)
Dia.AUC <- read.csv("Diagnosis.BCAA.csv",row.names = 1)
#Dia.AUC <- Dia.AUC %>% 
#  mutate(Cancer=factor(Cancer,levels = rev(sort(unique(.$Cancer)))))

library(RColorBrewer)
Cancer.Col <- c(brewer.pal(9,"Set1"),
                brewer.pal(8,"Set2")[1:7],
                brewer.pal(12,"Set3")[1:3],
                brewer.pal(8,"Dark2")[6],
                brewer.pal(12,"Set3")[4:12],
                brewer.pal(8,"Dark2")[1:4])

p1104 <- ggplot(Dia.AUC,aes(x=AUC,y=reorder(Cancer,AUC),fill=Cancer))+
  #geom_violin(cex=1.2)+           
  geom_boxplot()+
  #geom_jitter()
  ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9),
        axis.text.y = element_text(size=9),
        #axis.ticks.x = element_blank(),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="AUROC",y="")+
  scale_fill_manual(values = Cancer.Col)

ggsave(p1104,filename = "Diagnosis.BCAA.pdf",width = 2,height = 4)

######## PCA + UMAP #########
library(tidyverse)
BCAA.ID <- read.csv("BCAA.ID.csv")
PhenoType <- read.csv("../../mRNA.Phenotype.csv")
PhenoType <- PhenoType %>% filter(Sample.Type == "Primary Tumor" | 
                                    Sample.Type == "Solid Tissue Normal" | 
                                    Sample.Type == "Primary Blood Derived Cancer - Peripheral Blood") %>% 
  dplyr::select(Project.ID,Sample.ID,Sample.Type) %>%
  mutate(Type=ifelse(Sample.Type=="Solid Tissue Normal","Normal","Tumor")) %>%
  mutate(Cancer=str_remove(Project.ID,"TCGA-"))

#files <- list.files("../../mRNA.PanCancer.Exp/")
files <- list.files("../Pancancer/",pattern = "[DEG]",full.names = T)
files <- files[str_detect(files,"DEG.csv")]

dir = "../../Clinical.XENA/"
BCAA.Data <- data.frame()

for (file in files) {
  project <- str_remove(file,".*/") %>% str_remove_all("_.*") %>% str_remove_all(".*-")
  phenotype <- read.csv(paste(dir,project,"/TCGA-",project,".GDC_phenotype.tsv",sep = ""),sep = "\t",header = T) %>% 
    dplyr::filter(sample_type.samples %in% c("Primary Tumor","Primary Blood Derived Cancer - Peripheral Blood","Solid Tissue Normal"))
  TPM <- read.csv(paste("../../mRNA.PanCancer.Exp/TCGA-",project,".mRNA.TPM.csv",sep = ""),check.names = F,row.names = 1)
  TPM <- log2(TPM+1)
  
  Inter.Samples <- intersect(phenotype$submitter_id.samples,colnames(TPM))
  Pheno.mid <- phenotype %>% filter(submitter_id.samples %in% Inter.Samples) %>%
    mutate(Label = ifelse(sample_type.samples %in% c("Primary Tumor","Primary Blood Derived Cancer - Peripheral Blood"),1,0)) %>%
    dplyr::select(submitter_id.samples,Label)
  TPM.Target <- TPM[BCAA.ID$ENSEMBL,Pheno.mid$submitter_id.samples] %>% 
    t() %>% data.frame(check.names = F)
  colnames(TPM.Target) <- BCAA.ID$SYMBOL
  TPM.Target$Label = Pheno.mid$Label
  BCAA.Data <- TPM.Target %>% 
    rownames_to_column("SampleID") %>%
    mutate(Cancer=project) %>%
    rbind.data.frame(BCAA.Data)
}
BCAA.Data2 <- BCAA.Data %>% mutate(Type=if_else(Label==1,"Tumor","Normal"))
write.csv(BCAA.Data2,file = "BCAA.PCA.csv",row.names = F)
#### PCA ####
BCAA.Data2 <- read.csv(file = "BCAA.PCA.csv")
PCA.Data <- BCAA.Data2 %>%
  dplyr::select(-Label,-Cancer,-Type) %>%
  remove_rownames() %>%
  column_to_rownames("SampleID")
PCA.Condition <- BCAA.Data2 %>%
  dplyr::select(SampleID,Label,Cancer,Type)

library(ggbiplot)
#data(wine)
BCAA.PCA <- prcomp(PCA.Data, scale. = TRUE)
PCA.Data <- BCAA.PCA$x[,1:5] %>% data.frame(check.names = F) %>%
  cbind.data.frame(PCA.Condition)
write.csv(PCA.Data,file = "BCAA.PCA.Figure.Data.csv")
Immune.Col <- c("#293462","#a64942","#fe5f55","#fff1c1",
                "#5bd1d7","#348498","#004d61","#ff502f",
                "#e41749","#f5587b","#ff8a5c","#fff591",
                "#001871","#ff585d","#ffb549","#41b6e6",
                "#515bd4","#8134af","#dd2a7b","#feda77",
                "#96ceb4","#ffeead","#d9534f","#ffad60",
                "#05445c","#f2317f","#5c4f74","#040000",
                "#de4307","#f29c2b","#f6d04d","#8bc24c",
                "#fef4a9","#3b9a9c","#4bc2c5","#78fee0",
                "#fad3cf","#a696c8","#2470a0","#c50d66",
                "#f07810","#eec60a","#dd7777","#1687a7",
                "#014955","#dd0a35","#ED5485","#FFE869",
                "#bc8420","#F0B775","#D25565","#2E94B9")
##Immune.Col[c(1:4,6:12,15,17,18,20,30:38)]
p0105 <- ggplot(PCA.Data,aes(x=PC1,y=PC2,color=Cancer,shape=Type))+
  geom_point(alpha=0.9)+
  #geom_violin(cex=1.2)+           
  #geom_boxplot()+
  #geom_jitter()
  ggthemes::theme_few()+
  theme(#legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9),
        axis.text.y = element_text(size=9),
        #axis.ticks.x = element_blank(),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="PC1(40.6%)",y="PC2(9.5%)")+
  scale_color_manual(values = Immune.Col[c(1:4,6:12,15,17,18,20,30:38)])
ggsave(p0105,filename = "BCAA.PCA.Legend.pdf",width = 6,height = 4)
ggsave(p0105,filename = "BCAA.PCA.pdf",width = 4,height = 2)
#### UMAP ####
BCAA.UMAP <- uwot::umap(PCA.Data, n_neighbors = 15, min_dist = 0.001, target_weight = 0.5)
UMAP.Data <- BCAA.UMAP %>% data.frame() %>% 
  cbind.data.frame(PCA.Condition)
p0106 <- ggplot(UMAP.Data,aes(x=X1,y=X2,color=Cancer,shape=Type))+
  geom_point(alpha=0.9)+
  #geom_violin(cex=1.2)+           
  #geom_boxplot()+
  #geom_jitter()
  ggthemes::theme_few()+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9),
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="UMAP_1",y="UMAP_2")+
  scale_color_manual(values = Immune.Col[c(1:4,6:12,15,17,18,20,30:38)])
ggsave(p0106,filename = "BCAA.UMAP.Legend.pdf",width = 6,height = 4)
ggsave(p0106,filename = "BCAA.UMAP.pdf",width = 4,height = 2)

######## DEGs => TCGA => All Samples ###########
BCAA.ID <- read.csv("BCAA.ID.csv")

PhenoType <- read.csv("../../mRNA.Phenotype.csv")
PhenoType <- PhenoType %>% filter(Sample.Type == "Primary Tumor" | 
                                    Sample.Type == "Solid Tissue Normal" | 
                                    Sample.Type == "Primary Blood Derived Cancer - Peripheral Blood") %>% 
  dplyr::select(Project.ID,Sample.ID,Sample.Type) %>%
  mutate(Type=ifelse(Sample.Type=="Solid Tissue Normal","Normal","Tumor")) %>%
  mutate(Cancer=str_remove(Project.ID,"TCGA-"))

#files <- list.files("../../mRNA.PanCancer.Exp/")
files <- list.files("../Pancancer/",pattern = "[DEG]",full.names = T)
files <- files[str_detect(files,"DEG.csv")]
BCAA.DEG <- data.frame()
for (file in files) {
  project = file %>% str_remove_all(".*\\/") %>%
    str_remove_all("_.*") %>% str_remove_all("TCGA-")
  data <- read.csv(file.path("../Pancancer/",file))
  BCAA.DEG <- data %>% filter(X %in% BCAA.ID$ENSEMBL) %>%
    mutate(Cancer=project) %>%
    rbind.data.frame(BCAA.DEG)
}
BCAA.DEG <- merge(BCAA.DEG,BCAA.ID,by.x="X",by.y = "ENSEMBL")
write.csv(BCAA.DEG,file="BCAA.DEGs.csv")
#### Figure ####
BCAA.DEG <- read.csv("BCAA.DEGs.csv",row.names = 1)
BCAA.DEG2 <- BCAA.DEG %>% 
  mutate(Threshold=if_else(padj>0.5,"No",
                           if_else(log2FoldChange>=1,"Up",
                                   if_else(log2FoldChange<=-1,"Down","No")))) %>%
  mutate(Label=cut(log2FoldChange,breaks = c(-Inf,-2,-1,0,1,2,Inf),
                   labels = c("≤-2","-2~-1","-1~0","0~1","1~2","≥2"))) %>%
  mutate(Cancer=factor(Cancer,levels = rev(sort(unique(.$Cancer))))) #%>%
  #mutate(SYMBOL=factor(SYMBOL,levels = ))

BCAA.DEG2$Sig <- if_else(BCAA.DEG2$padj > 0.05,"P>0.05",as.character(BCAA.DEG2$Label))
BCAA.DEG2$Sig <- factor(BCAA.DEG2$Sig,levels = c("≤-2","-2~-1","-1~0","0~1","1~2","≥2","P>0.05"))

library(RColorBrewer)
p1030.2 <- ggplot(BCAA.DEG2,aes(SYMBOL,Cancer,fill=Sig))+
  geom_tile(color="black")+
  ggthemes::theme_few()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9,angle = 90,vjust = 0.5,hjust = 1),
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="",y="")+
  scale_fill_manual("log2FoldChange",values = brewer.pal(11,"RdBu")[c(11,10,8,4,2,1,6)])

ggsave(p1030.2,filename = "BCAA.DEGs.Heatmap.pdf",height = 4,width = 6)

######## DEGs => TCGA => Paired Samples #######
library(readxl)
Data <- read_xlsx("GSCA/DifferentialExpressionTable.xlsx")
Data2 <- Data %>% 
  mutate(log2FoldChange=log2(fc)) %>%
  mutate(Threshold=if_else(fdr>0.5,"No",
                           if_else(log2FoldChange>=1,"Up",
                                   if_else(log2FoldChange<=-1,"Down","No")))) %>%
  mutate(Label=cut(log2FoldChange,breaks = c(-Inf,-2,-1,0,1,2,Inf),
                   labels = c("≤-2","-2~-1","-1~0","0~1","1~2","≥2"))) %>%
  mutate(cancertype=factor(cancertype,levels = rev(sort(unique(.$cancertype)))))
  
Data2$Sig <- if_else(Data2$fdr > 0.05,"P>0.05",as.character(Data2$Label))
Data2$Sig <- factor(Data2$Sig,levels = c("≤-2","-2~-1","-1~0","0~1","1~2","≥2","P>0.05"))

library(RColorBrewer)
p1030.3 <- ggplot(Data2,aes(symbol,cancertype,fill=Sig))+
  geom_tile(color="black")+
  ggthemes::theme_few()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9,angle = 90,vjust = 0.5,hjust = 1),
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="",y="")+
  scale_fill_manual("log2FoldChange",values = brewer.pal(11,"RdBu")[c(11,10,8,4,2,1,6)])

ggsave(p1030.3,filename = "BCAA.DEGs.Heatmap.PairedSamples.pdf",height = 4,width = 6)

######## Survival => OS ##########
library(tidyverse)
library(readxl)
library(survival)
library(survminer)
dir.create("Survival.KM")
BCAA.ID <- read.csv("BCAA.ID.csv")
dirs <- list.dirs("../../Clinical.XENA/",full.names = T)
dirs <- dirs[str_detect(dirs,"//")]
BCAA.KM.Pancancer <- data.frame()
for (dir in dirs) {
  print(dir)
  project <- str_remove(dir,".*/")
  #dir.create(paste("AASS.KM/",project,sep = ""))
  ##submitter_id.samples sample_type.samples
  phenotype <- read.csv(paste(dir,"/TCGA-",project,".GDC_phenotype.tsv",sep = ""),sep = "\t",header = T) %>% 
    dplyr::filter(sample_type.samples == "Primary Tumor" | 
                    sample_type.samples == "Primary Blood Derived Cancer - Peripheral Blood")
  TPM <- read.csv(paste("../../mRNA.PanCancer.Exp/TCGA-",project,".mRNA.TPM.csv",sep = ""),check.names = F,row.names = 1)
  TPM.log <- log2(TPM+1)
  #### Calculate AASS ssGSEA 
  library(GSVA)
  BCAA=list()
  BCAA[["BCAA"]]<-BCAA.ID$ENSEMBL
  res.ssgsea <- gsva(as.matrix(TPM.log), 
                     gset.idx.list=BCAA, 
                     method = "ssgsea", 
                     kcdf = "Gaussian")
  
  #### survival datasets
  survival.data <- read.csv(paste(dir,"/TCGA-",project,".survival.tsv",sep = ""),sep = "\t",header = T,check.names = F)
  
  Inter.Samples <- intersect(intersect(phenotype$submitter_id.samples,colnames(TPM)),survival.data$sample)
  #### symbol expression
  TPM.ENSEMBL.Target <- TPM.log[BCAA.ID$ENSEMBL,Inter.Samples] %>% 
    t() %>% data.frame(check.names = F) %>% 
    mutate(SampleID = rownames(.)) %>% 
    merge(.,survival.data,by.x="SampleID",by.y="sample")
  #### BCAA survival
  TPM.log.Target <- res.ssgsea %>% t() %>% data.frame(check.names = F) %>%
    rownames_to_column("SampleID") %>%
    merge(TPM.ENSEMBL.Target,by="SampleID") %>%
    mutate(Cancer=project)
  
  #write.csv(TPM.log.Target,file = paste("AASS.KM/TCGA-",project,".AASS.SurvivalData.csv",sep = ""),row.names = F)
  
  res.cut <- surv_cutpoint(TPM.log.Target,
                           time = paste("OS",".time",sep = ""),
                           event = "OS",
                           variables = "BCAA")
  #res.cat <- surv_categorize(res.cut)
  
  theme.sur <- ggthemes::theme_few()+
    theme(legend.position = "none",
          axis.text = element_text(color = "black"),
          axis.text.x = element_text(size=15),
          axis.text.y = element_text(size=15),
          axis.title = element_text(size=18),
          panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))
  
  dat <- surv_categorize(res.cut) #%>% data.frame() %>% mutate(Sample=SNHGs.ICI.data$Sample.ID)
  fit <- survfit(Surv(OS.time, OS) ~ BCAA,
                 data = dat)
  
  index="OS"
  diff=survdiff(as.formula(paste("Surv(",index,".time,",index,") ~ ","BCAA",sep = "")),data = dat)
  pValue=1-pchisq(diff$chisq,df=1)
  
  cox <- coxph(as.formula(paste("Surv(OS.time, OS) ~ ","BCAA",sep = "")), data = TPM.log.Target)
  coxSummary = summary(cox)
  
  BCAA.KM.Pancancer=rbind(BCAA.KM.Pancancer,
                          data.frame(ENSEMBL="BCAA",
                                     KM.pvalue=pValue,
                                     HR=coxSummary$conf.int[,"exp(coef)"],
                                     HR.95L=coxSummary$conf.int[,"lower .95"],
                                     HR.95H=coxSummary$conf.int[,"upper .95"],
                                     HR.pvalue=coxSummary$coefficients[,"Pr(>|z|)"],
                                     Cancer=project))
  
  #### GENE SURVIVAL
  #OS.Tab=data.frame()
  for(i in BCAA.ID$ENSEMBL){
    symbol <- BCAA.ID[BCAA.ID$ENSEMBL == i,]$SYMBOL %>% as.character()
    cox <- coxph(as.formula(paste("Surv(OS.time, OS) ~ ",i,sep = "")), data = TPM.log.Target)
    coxSummary = summary(cox)
    #group=ifelse(TPM.log.Target[[i]] > median(TPM.log.Target[[i]]),"high","low")
    #if(length(table(group))==1) return(NULL) #
    res.cut <- surv_cutpoint(TPM.log.Target,
                             time = paste("OS",".time",sep = ""),
                             event = "OS",
                             variables = i)
    #res.cat <- surv_categorize(res.cut)
    
    theme.sur <- ggthemes::theme_few()+
      theme(legend.position = "none",
            axis.text = element_text(color = "black"),
            axis.text.x = element_text(size=15),
            axis.text.y = element_text(size=15),
            axis.title = element_text(size=18),
            panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))
    
    dat <- surv_categorize(res.cut) #%>% data.frame() %>% mutate(Sample=SNHGs.ICI.data$Sample.ID)
    fit <- survfit(as.formula(paste("Surv(",index,".time,",index,") ~ ",i,sep = "")),
                   data = dat)
    
    index="OS"
    diff=survdiff(as.formula(paste("Surv(",index,".time,",index,") ~ ",i,sep = "")),data = dat)
    pValue=1-pchisq(diff$chisq,df=1)
    
    #cox <- coxph(as.formula(paste("Surv(OS.time, OS) ~ ",i,sep = "")), data = TPM.log.Target)
    #coxSummary = summary(cox)
    
    theme.sur <- ggthemes::theme_few()+
      theme(legend.position = "none",
            axis.text = element_text(color = "black"),
            axis.text.x = element_text(size=15),
            axis.text.y = element_text(size=15),
            axis.title = element_text(size=18),
            panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))
    library(ggthemes)
    p<-ggsurvplot(fit,pval =TRUE, data = dat, 
                  surv.median.line = "hv",
                  legend.title = symbol,
                  conf.int.style = "step",
                  xlab = "Time in days",
                  #break.time.by = 500,
                  risk.table = "abs_pct",
                  #risk.table.y.text.col = T,
                  #risk.table.y.text = FALSE,
                  legend.labs = c("High", "Low"),
                  #pval = TRUE,
                  conf.int = TRUE,
                  palette = "Set1",
                  #ggtheme = theme.sur,
                  risk.table.y.text.col = T,
                  risk.table.y.text = F,
                  ggtheme = theme.sur)
    print(p)
    #dev.off()
    graph2pdf(file=paste("Survival.KM/TCGA-",project,".OS.BestCutoff.KM.",symbol,".pdf",sep = ''),height = 5,width = 4)
    
    BCAA.KM.Pancancer=rbind(BCAA.KM.Pancancer,
                            data.frame(ENSEMBL=symbol,
                                       KM.pvalue=pValue,
                                       HR=coxSummary$conf.int[,"exp(coef)"],
                                       HR.95L=coxSummary$conf.int[,"lower .95"],
                                       HR.95H=coxSummary$conf.int[,"upper .95"],
                                       HR.pvalue=coxSummary$coefficients[,"Pr(>|z|)"],
                                       Cancer=project))
  }
}
write.csv(BCAA.KM.Pancancer,file = "Survival.OS.BCAA.csv",row.names = F)

if (F) {
  result.cox <- coxph(Surv(time, status) ~ sex, data =  lung)
  summary(result.cox)
  
  ## Method1
  data <- lung
  covariates <- c("age", "sex", "ph.karno", "ph.ecog", "wt.loss")
  
  results <- data.frame(covariate = character(0),
                        beta = numeric(0),
                        HR.confint.upper=numeric(0),
                        HR.confint.upper=numeric(0),
                        wald.test = numeric(0),
                        p.value = numeric(0))
  
  for (covariate in covariates) {
    formula <- as.formula(paste("Surv(time, status) ~", covariate))
    model <- coxph(formula, data = data)
    summary <- summary(model)
    p.value <- signif(summary$wald["pvalue"], digits = 2)
    wald.test <- signif(summary$wald["test"], digits = 2)
    beta <- signif(summary$coef[1], digits = 2)
    HR.confint.lower <- signif(summary$conf.int[,"lower .95"], 2)
    HR.confint.upper <- signif(summary$conf.int[,"upper .95"], 2)
    results <- rbind(results, data.frame(covariate = covariate,
                                         beta = beta,
                                         HR.confint.lower=HR.confint.lower,
                                         HR.confint.upper=HR.confint.upper,
                                         wald.test= wald.test,
                                         p.value = p.value))
  }
  
  #results
  
  ### Method2
  data <- lung
  
  covariates <- c("age", "sex", "ph.karno", "ph.ecog", "wt.loss")
  
  univ_models <- map(covariates, function(covariate) {
    formula <- as.formula(paste("Surv(time, status) ~", covariate))
    model <- coxph(formula, data = data)
    summary <- summary(model)
    p.value <- signif(summary$wald["pvalue"], digits = 2)
    wald.test <- signif(summary$wald["test"], digits = 2)
    beta <- signif(summary$coef[1], digits = 2)
    HR.confint.lower <- signif(summary$conf.int[,"lower .95"], 2)
    HR.confint.upper <- signif(summary$conf.int[,"upper .95"], 2)
    HR <- paste0(signif(exp(summary$coef[2]), digits = 2), " (", HR.confint.lower, "-", HR.confint.upper, ")")
    
    return(list(covariate = covariate,
                beta = beta,
                HR = HR,
                wald.test = wald.test,
                p.value = p.value))
  })
  results <- reduce(univ_models, rbind)
  results <- do.call(rbind, univ_models)
}
#### Figure ####
library(tidyverse)
library(RColorBrewer)
BCAA.KM.Pancancer <- read.csv("Survival.OS.BCAA.csv")

BCAA.KM.Pancancer2 <- BCAA.KM.Pancancer %>%
  mutate(Label = if_else(KM.pvalue > 0.05,"p>0.05",
                         if_else(HR < 1,"Protective","Risky"))) %>%
  #mutate(ENSEMBL=factor(ENSEMBL,levels = rev(c("AASS",Gene.Symble)))) %>%
  mutate(Label=factor(Label,levels = c("Protective","Risky","p>0.05"))) %>%
  filter(ENSEMBL!="BCAA") %>%
  mutate(ENSEMBL=factor(ENSEMBL,levels = rev(sort(unique(.$ENSEMBL)))))

p1 <- ggplot()+
  geom_tile(data = BCAA.KM.Pancancer2,
            aes(y=ENSEMBL,x=Cancer,fill=Label),
            color="white")+
  geom_point(data = BCAA.KM.Pancancer2 %>% filter(KM.pvalue>0.05),
             aes(y=ENSEMBL,x=Cancer,size=-log10(KM.pvalue)),
             shape=4,color="black")+
  geom_point(data = BCAA.KM.Pancancer2 %>% filter(KM.pvalue<= 0.05),
             aes(y=ENSEMBL,x=Cancer,size=-log10(KM.pvalue)),
             shape=1,color="black")+
  scale_size_continuous(range = c(1,4))+
  ggthemes::theme_few()+
  theme(legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9,angle = 90,vjust = 0.5,hjust = 1),
    axis.text.y = element_text(size=9),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="",y="",size="-log10(KM p)",fill="Type")+
  scale_fill_manual(values = c("#4DBBD5FF","#E5086A","grey90"))

ggsave(p1,filename = "Survival.BCAA.OS.KM.pdf",height = 4,width = 8)

BCAA.KM.Pancancer3 <- BCAA.KM.Pancancer %>%
  mutate(Label = if_else(HR.pvalue > 0.05,"p>0.05",
                         if_else(HR < 1,"Protective","Risky"))) %>%
  #mutate(ENSEMBL=factor(ENSEMBL,levels = rev(c("AASS",Gene.Symble)))) %>%
  mutate(Label=factor(Label,levels = c("Protective","Risky","p>0.05"))) %>%
  filter(ENSEMBL!="BCAA") %>%
  mutate(Cancer=factor(Cancer,levels = rev(sort(unique(.$Cancer)))))

p2 <- ggplot()+
  geom_point(data = BCAA.KM.Pancancer3,
             aes(y=Cancer,x=ENSEMBL,fill=Label,size=-log10(HR.pvalue)),
             shape=21,color="white",stroke=1)+
  geom_point(data = BCAA.KM.Pancancer3 %>% filter(HR.pvalue<= 0.05),
             aes(y=Cancer,x=ENSEMBL,fill=Label,size=-log10(HR.pvalue)),
             shape=21,color="black",stroke=1.2)+
  scale_size_continuous(range = c(1,4))+
  ggthemes::theme_few()+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9,angle = 90,vjust = 0.5,hjust = 1),
    axis.text.y = element_text(size=9),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="",y="",size="-log10(COX p)",fill="Type")+
  scale_fill_manual(values = c("#469990","Orange","grey"))

ggsave(p2,filename = "Survival.BCAA.OS.COX.pdf",height = 5,width = 5)

######## Survival => DSS PFI #########
library(tidyverse)
library(readxl)
library(survival)
library(survminer)
#dir.create("Survival.KM")
BCAA.ID <- read.csv("BCAA.ID.csv")
dirs <- list.dirs("../../Clinical.XENA/",full.names = T)
dirs <- dirs[str_detect(dirs,"//")]
BCAA.KM.Pancancer <- data.frame()
for (dir in dirs[16:33]) {
  print(dir)
  project <- str_remove(dir,".*/")
  #dir.create(paste("AASS.KM/",project,sep = ""))
  ##submitter_id.samples sample_type.samples
  phenotype <- read.csv(paste(dir,"/TCGA-",project,".GDC_phenotype.tsv",sep = ""),sep = "\t",header = T) %>% 
    dplyr::filter(sample_type.samples == "Primary Tumor" | 
                    sample_type.samples == "Primary Blood Derived Cancer - Peripheral Blood") %>%
    mutate(Sample.ID=substr(submitter_id.samples,1,15))
  TPM <- read.csv(paste("../../mRNA.PanCancer.Exp/TCGA-",project,".mRNA.TPM.csv",sep = ""),check.names = F,row.names = 1)
  TPM.log <- log2(TPM+1)
  #### Calculate AASS ssGSEA 
  library(GSVA)
  BCAA=list()
  BCAA[["BCAA"]]<-BCAA.ID$ENSEMBL
  res.ssgsea <- gsva(as.matrix(TPM.log), 
                     gset.idx.list=BCAA, 
                     method = "ssgsea", 
                     kcdf = "Gaussian")
  
  #### survival datasets
  survival.data <- read.csv(paste(dir,"/../XENA.TCGA.",project,"_survival.txt",sep = ""),
                            sep = "\t",header = T,check.names = F,row.names = 1) %>%
    rownames_to_column("Sample.ID")
  
  #Inter.Samples <- intersect(intersect(phenotype$submitter_id.samples,colnames(TPM)),
  #                           survival.data$sample)
  #### symbol expression
  TPM.ENSEMBL.Target <- TPM.log[BCAA.ID$ENSEMBL,] %>% 
    t() %>% data.frame(check.names = F) %>% 
    mutate(SampleID = rownames(.)) %>% 
    mutate(Sample.ID=substr(SampleID,1,15)) %>%
    merge(.,survival.data,by.x="Sample.ID") %>%
    filter(Sample.ID %in% phenotype$Sample.ID)
  #### BCAA survival
  TPM.log.Target <- res.ssgsea %>% t() %>% data.frame(check.names = F) %>%
    rownames_to_column("SampleID") %>%
    merge(TPM.ENSEMBL.Target,by="SampleID") %>%
    mutate(Cancer=project)
  
  #write.csv(TPM.log.Target,file = paste("AASS.KM/TCGA-",project,".AASS.SurvivalData.csv",sep = ""),row.names = F)
  #index="PFI"
  for (index in c("DSS","PFI")) {
    middata <- TPM.log.Target %>% dplyr::select(c("BCAA",BCAA.ID$ENSEMBL,"DSS","DSS.time","PFI","PFI.time"))
    TPM.log.Target2 <- middata %>% 
      dplyr::select(c("BCAA",index,paste(index,".time",sep = ""))) %>%
      na.omit()
    if (nrow(TPM.log.Target2) >= 10) {
      res.cut <- surv_cutpoint(TPM.log.Target,
                               time = paste(index,".time",sep = ""),
                               event = index,
                               variables = "BCAA")
      #res.cat <- surv_categorize(res.cut)
      
      dat <- surv_categorize(res.cut) #%>% data.frame() %>% mutate(Sample=SNHGs.ICI.data$Sample.ID)
      #fit <- survfit(Surv(OS.time, OS) ~ BCAA,data = dat)
      
      #index="OS"
      diff=survdiff(as.formula(paste("Surv(",index,".time,",index,") ~ ","BCAA",sep = "")),data = dat)
      pValue=1-pchisq(diff$chisq,df=1)
      
      cox <- coxph(as.formula(paste("Surv(",index,".time,",index,") ~ ","BCAA",sep = "")), data = TPM.log.Target)
      coxSummary = summary(cox)
      
      BCAA.KM.Pancancer=rbind(BCAA.KM.Pancancer,
                              data.frame(ENSEMBL="BCAA",
                                         KM.pvalue=pValue,
                                         HR=coxSummary$conf.int[,"exp(coef)"],
                                         HR.95L=coxSummary$conf.int[,"lower .95"],
                                         HR.95H=coxSummary$conf.int[,"upper .95"],
                                         HR.pvalue=coxSummary$coefficients[,"Pr(>|z|)"],
                                         Cancer=project,
                                         Type=index))
      
      #### GENE SURVIVAL
      #OS.Tab=data.frame()
      for(i in BCAA.ID$ENSEMBL){
        symbol <- BCAA.ID[BCAA.ID$ENSEMBL == i,]$SYMBOL %>% as.character()
        cox <- coxph(as.formula(paste("Surv(OS.time, OS) ~ ",i,sep = "")), data = TPM.log.Target)
        coxSummary = summary(cox)
        #group=ifelse(TPM.log.Target[[i]] > median(TPM.log.Target[[i]]),"high","low")
        #if(length(table(group))==1) return(NULL) #
        res.cut <- surv_cutpoint(TPM.log.Target,
                                 time = paste(index,".time",sep = ""),
                                 event = index,
                                 variables = i)
        #res.cat <- surv_categorize(res.cut)
        
        dat <- surv_categorize(res.cut) #%>% data.frame() %>% mutate(Sample=SNHGs.ICI.data$Sample.ID)
        fit <- survfit(as.formula(paste("Surv(",index,".time,",index,") ~ ",i,sep = "")),
                       data = dat)
        
        #index="OS"
        diff=survdiff(as.formula(paste("Surv(",index,".time,",index,") ~ ",i,sep = "")),data = dat)
        pValue=1-pchisq(diff$chisq,df=1)
        
        #cox <- coxph(as.formula(paste("Surv(OS.time, OS) ~ ",i,sep = "")), data = TPM.log.Target)
        #coxSummary = summary(cox)
        
        theme.sur <- ggthemes::theme_few()+
          theme(legend.position = "none",
                axis.text = element_text(color = "black"),
                axis.text.x = element_text(size=15),
                axis.text.y = element_text(size=15),
                axis.title = element_text(size=18),
                panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))
        library(ggthemes)
        library(export)
        p<-ggsurvplot(fit,pval =TRUE, data = dat, 
                      surv.median.line = "hv",
                      legend.title = symbol,
                      conf.int.style = "step",
                      xlab = "Time in days",
                      #break.time.by = 500,
                      risk.table = "abs_pct",
                      #risk.table.y.text.col = T,
                      #risk.table.y.text = FALSE,
                      legend.labs = c("High", "Low"),
                      #pval = TRUE,
                      conf.int = TRUE,
                      palette = "Set1",
                      #ggtheme = theme.sur,
                      risk.table.y.text.col = T,
                      risk.table.y.text = F,
                      ggtheme = theme.sur)
        print(p)
        #dev.off()
        graph2pdf(file=paste("Survival.KM/TCGA-",project,".",index,".BestCutoff.KM.",symbol,".pdf",sep = ''),height = 5,width = 4)
        
        BCAA.KM.Pancancer=rbind(BCAA.KM.Pancancer,
                                data.frame(ENSEMBL=symbol,
                                           KM.pvalue=pValue,
                                           HR=coxSummary$conf.int[,"exp(coef)"],
                                           HR.95L=coxSummary$conf.int[,"lower .95"],
                                           HR.95H=coxSummary$conf.int[,"upper .95"],
                                           HR.pvalue=coxSummary$coefficients[,"Pr(>|z|)"],
                                           Cancer=project,
                                           Type=index))
      }
    }
    #write.csv(BCAA.KM.Pancancer,file = paste("Survival.",index,".BCAA.csv",sep = ""),row.names = F)
  }
}
BCAA.KM.Pancancer <- BCAA.KM.Pancancer %>% unique()
write.csv(BCAA.KM.Pancancer,file = "Survival.DSS.PFI.BCAA.csv")
#### Figure => DSS ####
library(tidyverse)
library(RColorBrewer)
BCAA.KM.Pancancer <- read.csv("Survival.DSS.PFI.BCAA.csv")

BCAA.KM.Pancancer2 <- BCAA.KM.Pancancer %>%
  filter(Type=="DSS") %>%
  mutate(Label = if_else(KM.pvalue > 0.05,"p>0.05",
                         if_else(HR < 1,"Protective","Risky"))) %>%
  #mutate(ENSEMBL=factor(ENSEMBL,levels = rev(c("AASS",Gene.Symble)))) %>%
  mutate(Label=factor(Label,levels = c("Protective","Risky","p>0.05"))) %>%
  filter(ENSEMBL!="BCAA") %>%
  mutate(ENSEMBL=factor(ENSEMBL,levels = rev(sort(unique(.$ENSEMBL)))))
  #mutate(Cancer=factor(Cancer,levels = rev(sort(unique(.$Cancer)))))

p1 <- ggplot()+
  geom_tile(data = BCAA.KM.Pancancer2,
            aes(y=ENSEMBL,x=Cancer,fill=Label),
            color="white")+
  geom_point(data = BCAA.KM.Pancancer2 %>% filter(KM.pvalue>0.05),
             aes(y=ENSEMBL,x=Cancer,size=-log10(KM.pvalue)),
             shape=4,color="black")+
  geom_point(data = BCAA.KM.Pancancer2 %>% filter(KM.pvalue<= 0.05),
             aes(y=ENSEMBL,x=Cancer,size=-log10(KM.pvalue)),
             shape=1,color="black")+
  scale_size_continuous(range = c(1,4))+
  ggthemes::theme_few()+
  theme(legend.position = "bottom",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9,angle = 90,vjust = 0.5,hjust = 1),
        axis.text.y = element_text(size=9),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="",y="",size="-log10(KM p)",fill="Type")+
  scale_fill_manual(values = c("#4DBBD5FF","#E5086A","grey90"))

ggsave(p1,filename = "Survival.BCAA.DSS.KM.pdf",height = 4,width = 8)



p1 <- ggplot()+
  geom_point(data = BCAA.KM.Pancancer2,
             aes(y=Cancer,x=ENSEMBL,fill=Label,size=-log10(KM.pvalue)),
             shape=21,color="white",stroke=1)+
  geom_point(data = BCAA.KM.Pancancer2 %>% filter(KM.pvalue<= 0.05),
             aes(y=Cancer,x=ENSEMBL,fill=Label,size=-log10(KM.pvalue)),
             shape=21,color="black",stroke=1.2)+
  scale_size_continuous(range = c(1,4))+
  ggthemes::theme_few()+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9,angle = 90,vjust = 0.5,hjust = 1),
    axis.text.y = element_text(size=9),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="",y="",size="-log10(KM p)",fill="Type")+
  scale_fill_manual(values = c("#0094FF","#E5086A","#CFB99E"))

ggsave(p1,filename = "Survival.BCAA.DSS.KM.pdf",height = 5,width = 4)

BCAA.KM.Pancancer3 <- BCAA.KM.Pancancer %>%
  filter(Type=="DSS") %>%
  mutate(Label = if_else(HR.pvalue > 0.05,"p>0.05",
                         if_else(HR < 1,"Protective","Risky"))) %>%
  #mutate(ENSEMBL=factor(ENSEMBL,levels = rev(c("AASS",Gene.Symble)))) %>%
  mutate(Label=factor(Label,levels = c("Protective","Risky","p>0.05"))) %>%
  filter(ENSEMBL!="BCAA") %>%
  mutate(Cancer=factor(Cancer,levels = rev(sort(unique(.$Cancer)))))

p2 <- ggplot()+
  geom_point(data = BCAA.KM.Pancancer3,
             aes(y=Cancer,x=ENSEMBL,fill=Label,size=-log10(HR.pvalue)),
             shape=21,color="white",stroke=1)+
  geom_point(data = BCAA.KM.Pancancer3 %>% filter(HR.pvalue<= 0.05),
             aes(y=Cancer,x=ENSEMBL,fill=Label,size=-log10(HR.pvalue)),
             shape=21,color="black",stroke=1.2)+
  scale_size_continuous(range = c(1,4))+
  ggthemes::theme_few()+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9,angle = 90,vjust = 0.5,hjust = 1),
    axis.text.y = element_text(size=9),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="",y="",size="-log10(COX p)",fill="Type")+
  scale_fill_manual(values = c("#469990","Orange","grey"))

ggsave(p2,filename = "Survival.BCAA.DSS.COX.pdf",height = 5,width = 5)
#### Figure => PFI ####
library(tidyverse)
library(RColorBrewer)
BCAA.KM.Pancancer <- read.csv("Survival.DSS.PFI.BCAA.csv")

BCAA.KM.Pancancer2 <- BCAA.KM.Pancancer %>%
  filter(Type=="PFI") %>%
  mutate(Label = if_else(KM.pvalue > 0.05,"p>0.05",
                         if_else(HR < 1,"Protective","Risky"))) %>%
  #mutate(ENSEMBL=factor(ENSEMBL,levels = rev(c("AASS",Gene.Symble)))) %>%
  mutate(Label=factor(Label,levels = c("Protective","Risky","p>0.05"))) %>%
  filter(ENSEMBL!="BCAA") %>%
  mutate(ENSEMBL=factor(ENSEMBL,levels = rev(sort(unique(.$ENSEMBL)))))
#mutate(Cancer=factor(Cancer,levels = rev(sort(unique(.$Cancer)))))

p1 <- ggplot()+
  geom_tile(data = BCAA.KM.Pancancer2,
            aes(y=ENSEMBL,x=Cancer,fill=Label),
            color="white")+
  geom_point(data = BCAA.KM.Pancancer2 %>% filter(KM.pvalue>0.05),
             aes(y=ENSEMBL,x=Cancer,size=-log10(KM.pvalue)),
             shape=4,color="black")+
  geom_point(data = BCAA.KM.Pancancer2 %>% filter(KM.pvalue<= 0.05),
             aes(y=ENSEMBL,x=Cancer,size=-log10(KM.pvalue)),
             shape=1,color="black")+
  scale_size_continuous(range = c(1,4))+
  ggthemes::theme_few()+
  theme(legend.position = "bottom",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9,angle = 90,vjust = 0.5,hjust = 1),
        axis.text.y = element_text(size=9),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="",y="",size="-log10(KM p)",fill="Type")+
  scale_fill_manual(values = c("#4DBBD5FF","#E5086A","grey90"))

ggsave(p1,filename = "Survival.BCAA.PFI.KM.pdf",height = 4,width = 8)

BCAA.KM.Pancancer3 <- BCAA.KM.Pancancer %>%
  filter(Type=="PFI") %>%
  mutate(Label = if_else(HR.pvalue > 0.05,"p>0.05",
                         if_else(HR < 1,"Protective","Risky"))) %>%
  #mutate(ENSEMBL=factor(ENSEMBL,levels = rev(c("AASS",Gene.Symble)))) %>%
  mutate(Label=factor(Label,levels = c("Protective","Risky","p>0.05"))) %>%
  filter(ENSEMBL!="BCAA") %>%
  mutate(Cancer=factor(Cancer,levels = rev(sort(unique(.$Cancer)))))

p2 <- ggplot()+
  geom_point(data = BCAA.KM.Pancancer3,
             aes(y=Cancer,x=ENSEMBL,fill=Label,size=-log10(HR.pvalue)),
             shape=21,color="white",stroke=1)+
  geom_point(data = BCAA.KM.Pancancer3 %>% filter(HR.pvalue<= 0.05),
             aes(y=Cancer,x=ENSEMBL,fill=Label,size=-log10(HR.pvalue)),
             shape=21,color="black",stroke=1.2)+
  scale_size_continuous(range = c(1,4))+
  ggthemes::theme_few()+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9,angle = 90,vjust = 0.5,hjust = 1),
    axis.text.y = element_text(size=9),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="",y="",size="-log10(COX p)",fill="Type")+
  scale_fill_manual(values = c("#469990","Orange","grey"))

ggsave(p2,filename = "Survival.BCAA.PFI.COX.pdf",height = 5,width = 5)

######## Methylation ############
BCAA.SYMBOLs <- read.csv("BCAA.ID.csv")
library(readxl)
# methylation <- read_xlsx("GSCA/ExpressionAndMethylationTable.xlsx")
library(tidyverse)
library(magrittr)

PhenoData <- read.csv("../../mRNA.Phenotype.csv") %>%
  filter(Sample.Type %in% c("Primary Tumor",
                            "Primary Blood Derived Cancer - Peripheral Blood",
                            "Solid Tissue Normal")) %>%
  mutate(SampleType = ifelse(Sample.Type == "Solid Tissue Normal","Normal","Tumor"))

GPL <- read.csv("../../Clinical.XENA/Methylation450K.GPL13534-11288.txt",sep = "\t",header = T,comment.char = "#")
GPL <- GPL %>% dplyr::select(ID,
                             UCSC_RefGene_Name,
                             UCSC_RefGene_Group,
                             Regulatory_Feature_Group,
                             Relation_to_UCSC_CpG_Island)
BCAA.GPL <- data.frame()
for (i in 1:nrow(GPL)) {
  print(i)
  SYMBOL=unique(unlist(str_split(GPL$UCSC_RefGene_Name[i],";")))
  if (length(SYMBOL) >= 1) {
    for (GENE in SYMBOL) {
      if (GENE %in% BCAA.SYMBOLs$SYMBOL) {
        BCAA.GPL <- data.frame(ID=GPL$ID[i],
                               UCSC_RefGene_Name=GENE,
                               UCSC_RefGene_Group=paste(unique(sort(unlist(str_split(GPL$UCSC_RefGene_Group[i],";")))),collapse = ";"),
                               Regulatory_Feature_Group=GPL$Regulatory_Feature_Group[i],
                               Relation_to_UCSC_CpG_Island=GPL$Relation_to_UCSC_CpG_Island[i]) %>%
          rbind.data.frame(BCAA.GPL)
      }
    }
  }
}
BCAA.GPL2 <- BCAA.GPL %>% filter(str_detect(ID,"cg")) %>%
  arrange(UCSC_RefGene_Name)

write.csv(BCAA.GPL2,file = "Methylation450K.BCAA.GPL.ID.csv",row.names = F)
#### Spearman => Data preparation ####
library(tidyverse)
BCAA.GPL2 <- read.csv(file = "Methylation450K.BCAA.GPL.ID.csv")
#PhenoData <- read.csv("../../mRNA.Phenotype.csv") %>%
PhenoData <- read.csv("../../mRNA.Phenotype.csv") %>%
  filter(Sample.Type %in% c("Primary Tumor",
                            "Primary Blood Derived Cancer - Peripheral Blood",
                            "Solid Tissue Normal")) %>%
  mutate(SampleType = ifelse(Sample.Type == "Solid Tissue Normal","Normal","Tumor"))

PROJECTS <- unique(PhenoData$Project.ID) %>% str_remove_all("TCGA-")

=> Datrfor (project in PROJECTS) {
  print(project)
  Cancer=paste("TCGA-",project,sep = "")
  mid.Pheno <- PhenoData %>% filter(Project.ID == Cancer) %>%
    dplyr::select(Project.ID,Sample.ID,Sample.Type,SampleType)
  File.path = paste("../TotalRNA.Microbiome/Clinical.XENA/",project,"/",sep = "")
  #File.path = paste(project,"/",sep = "")
  File.name = paste("TCGA-",project,".methylation450.tsv",sep = "")
  Data <- read.csv(file = file.path(File.path,File.name),sep = "\t",check.names = F,row.names = 1)
  Data2 <- Data %>% rownames_to_column("ID") %>% 
    #filter(ID %in% BCAA.GPL2$ID) %>%
    gather(SampleID,Methylation,-ID) %>%
    na.omit() %>%
    filter(SampleID %in% mid.Pheno$Sample.ID) %>%
    merge(mid.Pheno,by.x="Sample.ID",by.y="Sample.ID") %>%
    merge(BCAA.GPL2,by="ID")
  
  Output.name=paste("XENA.BCAA.Methylation450K.",project,".csv",sep = "")
  write.csv(Data2,file=Output.name)
}
#### Spearman => Calculate #### 
library(tidyverse)
BCAA.GPL2 <- read.csv(file = "Methylation450K.BCAA.GPL.ID.csv")
PhenoData <- read.csv("../../mRNA.Phenotype.csv") %>%
  filter(Sample.Type %in% c("Primary Tumor",
                            "Primary Blood Derived Cancer - Peripheral Blood",
                            "Solid Tissue Normal")) %>%
  mutate(SampleType = ifelse(Sample.Type == "Solid Tissue Normal","Normal","Tumor")) %>%
  mutate(Project.ID=str_remove(Project.ID,"TCGA-"))

BCAA.ID <- read.csv("BCAA.ID.csv") #%>% filter(ENSEMBL!="ENSG00000285069")

library(survival)
library(survminer)

Methylation.Wilcox <- data.frame()
Methylation.Exp.Spearman <- data.frame()
Methylation.Survival <- data.frame()
files <- list.files(path = "Methylation450K",pattern = "txt$")
for (file in files) {
  project = (str_split(file,"\\.") %>% unlist())[2]
  print(project)
  methylated.data <- read.csv(file.path("Methylation450K",file),sep = "\t",row.names = 1,check.names = F)
  #methylated.data[is.na(methylated.data)] = 0
  methylated.data <- methylated.data %>% t() %>% data.frame(check.names = F)
  NA.count <- sapply(methylated.data, function(x) sum(is.na(x)))
  methylated.data2 <- methylated.data[,NA.count < nrow(methylated.data)] %>%
    data.frame(check.names = F) %>%
    rownames_to_column("Sample.ID")
  if (any(str_detect(methylated.data2$Sample.ID,"-??A$"))) {
    mid.Pheno <- PhenoData %>% filter(Project.ID==project) %>%
      dplyr::select(Sample.ID,SampleType) %>%
      merge(methylated.data2,by="Sample.ID") %>%
      mutate(SampleID=substr(Sample.ID,1,15))
  }else{
    mid.Pheno <- PhenoData %>% filter(Project.ID==project) %>%
      dplyr::select(Sample.ID,SampleType) %>% 
      mutate(SampleID=substr(Sample.ID,1,15)) %>%
      merge(methylated.data2,by.x="SampleID",by.y="Sample.ID")
  }
  
  
  #### Spearman => Methylation + Expression ####
  Exp.filename <- paste("TCGA-",project,".mRNA.TPM.csv",sep = "")
  Exp <- read.csv(file = file.path("../../mRNA.PanCancer.Exp/",Exp.filename),row.names = 1,check.names = F)
  Exp <- log2(Exp+1)
  
  BCAA.Exp <- Exp[BCAA.ID$ENSEMBL,] %>% data.frame(check.names = F) %>%
    rownames_to_column("ENSEMBL") %>%
    merge(BCAA.ID,by="ENSEMBL") %>%
    dplyr::select(-ENSEMBL,-ENTREZ) %>%
    column_to_rownames("SYMBOL") %>%
    t() %>%
    data.frame(check.names = F) %>%
    rownames_to_column("Sample.ID") %>%
    merge(mid.Pheno,by="Sample.ID")
  
  for (symbol in unique(BCAA.GPL2$UCSC_RefGene_Name)) {
    mid.cg <- BCAA.GPL2 %>% filter(UCSC_RefGene_Name==symbol)
    
    if (symbol == "GATSL2") {
      symbol="CASTOR2"
    }else if(symbol == "GATSL3"){
      symbol="CASTOR1"
    }
    
    for (cgid in mid.cg$ID) {
      if (cgid %in% colnames(BCAA.Exp)) {
        test1 <- cor.test(BCAA.Exp[[cgid]],BCAA.Exp[[symbol]],
                          method = "spearman",
                          exact = F)
        Methylation.Exp.Spearman <- data.frame(SYMBOL=symbol,
                                               ID=cgid,
                                               Rho=test1$estimate,
                                               Pvalue=test1$p.value,
                                               Cancer=project,
                                               Count=length(BCAA.Exp[[cgid]])) %>%
          rbind.data.frame(Methylation.Exp.Spearman)
      }
    }
  }
  
  #### Survival => Methylation ####
  TUMOR <- mid.Pheno %>% filter(SampleType=="Tumor") %>%
    mutate(SampleID=substr(Sample.ID,1,15))
  Sur.filename <- paste("XENA.TCGA.",project,"_survival.txt",sep = "")
  Sur.Data <- read.csv(file.path("../../Clinical.XENA/",Sur.filename),row.names = 1,sep = "\t",check.names = F) %>%
    rownames_to_column("SampleID") %>%
    merge(TUMOR,by="SampleID")
  
  for (id in colnames(TUMOR)[colnames(TUMOR) %>% str_detect("cg")]) {
    Sur.Data2 <- Sur.Data %>% dplyr::select(c("OS","OS.time",id)) %>%
      na.omit()
    Sur.Data2$group=ifelse(Sur.Data2[[id]] > median(Sur.Data2[[id]]),"high","low")
    if(length(table(Sur.Data2$group))==1) print("Only one group")
    
    index="OS"
    diff=survdiff(as.formula(paste("Surv(",index,".time,",index,") ~ ","group",sep = "")),
                  data = Sur.Data2)
    pValue=1-pchisq(diff$chisq,df=1)
    
    cox <- coxph(as.formula(paste("Surv(OS.time, OS) ~ ",id,sep = "")), 
                 data = Sur.Data2)
    coxSummary = summary(cox)
    
    Methylation.Survival=rbind(Methylation.Survival,
                               data.frame(ID=id,
                                          KM.pvalue=pValue,
                                          HR=coxSummary$conf.int[,"exp(coef)"],
                                          HR.95L=coxSummary$conf.int[,"lower .95"],
                                          HR.95H=coxSummary$conf.int[,"upper .95"],
                                          HR.pvalue=coxSummary$coefficients[,"Pr(>|z|)"],
                                          Cancer=project))
  }
  
  CTRL <- mid.Pheno %>% filter(SampleType=="Normal")
  #### Wilcox ####
  print(nrow(CTRL))
  if (nrow(CTRL) >= 10 & nrow(TUMOR) >= 10) {
    Tumor.Samples <- mid.Pheno %>% filter(SampleType=="Tumor")
    Ctrl.Samples <- mid.Pheno %>% filter(SampleType=="Normal")
    for (cg.id in colnames(methylated.data2[,2:ncol(methylated.data2)])) {
      test <- wilcox.test(Tumor.Samples[[cg.id]],
                          Ctrl.Samples[[cg.id]],
                          exact = F,
                          alternative = "two.sided")
      Methylation.Wilcox <- data.frame(ID=cg.id,
                                       Wilcox.pvalue=test$p.value,
                                       Mean.Ctrl = mean(na.omit(Ctrl.Samples[[cg.id]])),
                                       Median.Ctrl = median(na.omit(Ctrl.Samples[[cg.id]])),
                                       Mean.Tumor = mean(na.omit(Tumor.Samples[[cg.id]])),
                                       Median.Tumor = median(na.omit(Tumor.Samples[[cg.id]])),
                                       Cancer=project) %>%
        rbind.data.frame(Methylation.Wilcox)
    }
  }
}

#Methylation.BCAA <- Methylation.Exp.Spearman %>% remove_rownames() %>%
#  #merge(Methylation.Wilcox,by=c("ID","Cancer")) %>%
#  merge(Methylation.Survival,by=c("ID","Cancer")) %>%
#  merge(BCAA.GPL2,by="ID")

write.csv(Methylation.Wilcox,file = "Methylation.BCAA.Wilcox.csv",row.names = F)
write.csv(Methylation.Exp.Spearman,file = "Methylation.BCAA.Exp.Spearman.csv",row.names = F)
write.csv(Methylation.Survival,file = "Methylation.BCAA.Survival.csv",row.names = F)

#### Figure => All methylation sites ####
#Methylation.BCAA <- read.csv(file = "Methylation.BCAA.Survival.Wilcox.csv")
Methylation.Exp.Spearman <- read.csv(file = "Methylation.BCAA.Exp.Spearman.csv")

read.csv(file = "Methylation.BCAA.Exp.Spearman.csv") %>%
  filter(Pvalue < 0.05) %>%
  write.csv("TableS3.csv",row.names = F)

Methylation.Survival <- read.csv(file = "Methylation.BCAA.Survival.csv")
#Methylation.Survival2 <- round(Methylation.Survival,1e20)
Methylation.BCAA2 <- Methylation.Exp.Spearman %>% 
  dplyr::select(ID,Cancer,Rho) %>%
  #pivot_longer(-gene)
  spread(Cancer,Rho,fill=0) %>%
  remove_rownames() %>%
  column_to_rownames("ID")

library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(ggtree)
library(aplot)
library(ggh4x)
## heatmap
heatmap <- ggplot(Methylation.Exp.Spearman,aes(x=Cancer,y=ID,fill=Rho))+
  geom_tile()+
  scale_fill_gradient2(mid="white",low="#313695",high="#A50026")+
  guides(fill=guide_colorbar(direction = "vertical",
                             reverse = F,barwidth = unit(.6, "cm"),
                             barheight = unit(3,"cm")))+
  ggthemes::theme_few()+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9,angle = 90,vjust = 0.5),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="",y="",fill="Spearman'CC")

## cancer Cluster
cancerorder <- hclust(dist(t(Methylation.BCAA2)))
phc <- hclust(dist(t(Methylation.BCAA2))) %>% ggtree() + 
  layout_dendrogram()+
  theme_void()
## methylation sites Cluster
phr <- hclust(dist(Methylation.BCAA2)) %>% 
  ggtree(layout="rectangular",branch.length="none")+
  theme_void()
## methylation type
library(RColorBrewer)
Cg.group <- BCAA.GPL2 %>% 
  ggplot(aes(x=1,y=ID,fill=UCSC_RefGene_Group))+
  geom_col(width=1,color="black",size=0.01)+
  theme_void()+
  scale_fill_manual(values = c(brewer.pal(9,"Set1"),brewer.pal(12,"Set3")))

Cg.cpg <- BCAA.GPL2 %>% 
  ggplot(aes(x=1,y=ID,fill=Relation_to_UCSC_CpG_Island))+
  geom_col(width=1,color="black",size=0.01)+
  theme_void()+
  scale_fill_manual(values = c("white",brewer.pal(8,"Dark2")[c(1,2,3,4,6)]))

## methylation genes
Colors17 <- c("#949483","#F47B7B","#9F1F5C","#EF9020",
              "#00AF3E","#85B7E2","#29245C","#FFD616",
              "#E5352B","#E990AB","#0081B4","#96CBB3",
              "#91BE3E","#39A6DD","#EB0973","#DDE2E0",
              "#333C41","#9013FE","#FFC208","#1FADC5",
              "#2B697A")
Cg.genes <- Methylation.Exp.Spearman %>% 
  dplyr::select(ID,SYMBOL) %>%
  unique() %>%
  mutate(SYMBOL=factor(SYMBOL,levels = BCAA.ID$SYMBOL)) %>%
  ggplot(aes(x=1,y=ID,fill=SYMBOL))+
  geom_col(width=1,color="black",size=0.01)+
  theme_void()+
  scale_fill_manual(values = Colors17)

## methylation survival
Cg.survival <- Methylation.Survival %>%
  mutate(Type = if_else(KM.pvalue <= 0.05 | HR.pvalue <= 0.05,
                        if_else(HR > 1,"Risky","Protective"),"p>0.05")) %>%
  mutate(Type=factor(Type,levels = c("Protective","Risky","p>0.05"))) %>%
  mutate(Cancer=factor(Cancer,levels = colnames(Methylation.BCAA2)[cancerorder$order])) %>%
  ggplot(aes(x=Cancer,y=ID,fill=Type))+
  geom_tile()+ #
  ggthemes::theme_few()+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9,angle = 90,vjust = 0.5),
    axis.text.y = element_blank(),
    axis.ticks.y=element_blank(),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="",y="",fill="")+
  scale_fill_manual(values = c(brewer.pal(8,"Set1")[1],brewer.pal(8,"Set1")[2],"white"))

## Final
BCAA.Methylation.All <- heatmap %>%
  #insert_top(cell_type,height=0.05) %>% 
  aplot::insert_top(phc,height=0.05) %>% 
  #insert_left(gene_type,width = 0.06) %>%
  aplot::insert_left(phr,width = 0.1) %>%
  aplot::insert_right(Cg.genes,width = 0.05) %>%
  aplot::insert_right(Cg.group,width = 0.05) %>%
  aplot::insert_right(Cg.cpg,width = 0.05) %>%
  aplot::insert_right(Cg.survival,width = 1)
ggsave(BCAA.Methylation.All,filename = "Methylation.BCAA.Figure.pdf",height = 8,width = 12)
ggsave(BCAA.Methylation.All,filename = "Methylation.BCAA.Figure.Legend.pdf",height = 16,width = 12)

#### Figure => Representative ####
Methylation.Exp.Spearman <- read.csv(file = "Methylation.BCAA.Exp.Spearman.csv") %>%
  mutate(Cancer=paste(Cancer,"(n=",Count,")",sep = "")) %>%
  dplyr::select(-Count)
Methylation.Survival <- read.csv(file = "Methylation.BCAA.Survival.csv")
Methylation.Wilcox <- read.csv(file = "Methylation.BCAA.Wilcox.csv")

BCAA.GPL2 <- read.csv("Methylation450K.BCAA.GPL.ID.csv")

CgID <- BCAA.GPL2 %>%
  filter(!str_detect(UCSC_RefGene_Group,"Body"))

Methylation.Survival2 <- Methylation.Survival %>%
  filter(HR < 1) %>%
  filter(KM.pvalue <= 0.05)

Methylation.Exp <- Methylation.Exp.Spearman %>%
  filter(Pvalue <= 0.05) %>%
  filter(Rho <= 0) %>%
  #filter(ID %in% CgID$ID) %>%
  group_by(ID) %>%
  summarise(Count=n()) %>%
  ungroup() %>%
  arrange(desc(Count)) %>%
  merge(Methylation.Exp.Spearman,by="ID") %>%
  dplyr::select(ID,Count,SYMBOL) %>%
  unique() %>%
  data.frame(check.names = F) %>%
  group_by(SYMBOL) %>%
  dplyr::top_n(2,Count) %>%
  merge(Methylation.Exp.Spearman,by=c("ID","SYMBOL")) %>%
  mutate(SYMBOL=factor(SYMBOL,levels = BCAA.ID$SYMBOL)) %>%
  arrange(SYMBOL,ID) %>%
  mutate(ID=factor(ID,levels = rev(unique(.$ID))))

Figure.Body <- ggplot()+
  geom_point(data=Methylation.Exp,
             aes(x=Cancer,y=ID,fill=Rho,size=-log10(Pvalue)),
             shape=21,color="white",stroke=1) +
  geom_point(data=Methylation.Exp %>% filter(Pvalue<=0.05),
             aes(x=Cancer,y=ID,fill=Rho,size=-log10(Pvalue)),
             shape=21,color="black",stroke=1)+
  scale_size_continuous(range = c(1,5))+
  ggthemes::theme_few()+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9,angle = 90,vjust = 0.5),
    #axis.text.y = element_blank(),
    #axis.ticks.y = element_blank(),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="",y="",fill="Spearman'CC")+
  scale_fill_gradient2(mid="white",high="#A50026",low="#313695")

Figure.R2 <- Methylation.Exp %>%
  mutate(Label=if_else(Pvalue > 0.05,"P>0.05",
                       if_else(Rho > 0,"Positive","Negative","P>0.05"))) %>%
  group_by(ID,Label) %>%
  summarise(Count=n()) %>%
  spread(Label,Count,fill=0) %>% 
  mutate(ID=factor(ID,levels = levels(Methylation.Exp$ID))) %>%
  arrange(ID) %>%
  gather(Type,Count,-ID) %>%
  mutate(ID=factor(ID,levels = levels(Methylation.Exp$ID))) %>%
  mutate(Type=factor(Type,levels = rev(c("Positive","Negative","P>0.05")))) %>%
  ggplot(aes(x=Count,y=ID,fill=Type))+
  geom_col(color="black")+
  ggthemes::theme_few()+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9,angle = 90,vjust = 0.5),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="Cancer count",y="",fill="Correlation")+
  scale_fill_manual(values = rev(c("#E7298A","#1B9E77","grey")))

Colors17 <- c("#949483","#F47B7B","#9F1F5C","#EF9020",
              "#00AF3E","#85B7E2","#29245C","#FFD616",
              "#E5352B","#E990AB","#0081B4","#96CBB3",
              "#91BE3E","#39A6DD","#EB0973","#DDE2E0",
              "#333C41","#9013FE","#FFC208","#1FADC5",
              "#2B697A")

Figure.R1 <- Methylation.Exp %>%
  dplyr::select(ID,SYMBOL) %>%
  unique() %>%
  mutate(SYMBOL=factor(SYMBOL,levels = BCAA.ID$SYMBOL)) %>%
  ggplot(aes(x=1,y=ID,fill=SYMBOL))+
  geom_col(color="black")+
  ggthemes::theme_few()+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="",y="",fill="Symbol")+
  scale_fill_manual(values = Colors17)

library(RColorBrewer)
Cg.group <- BCAA.GPL2 %>% 
  ggplot(aes(x=1,y=ID,fill=UCSC_RefGene_Group))+
  geom_col(width=1,color="black",size=0.01)+
  theme_void()+
  scale_fill_manual(values = c(brewer.pal(9,"Set1"),brewer.pal(12,"Set3")))

Cg.cpg <- BCAA.GPL2 %>% 
  ggplot(aes(x=1,y=ID,fill=Relation_to_UCSC_CpG_Island))+
  geom_col(width=1,color="black",size=0.01)+
  theme_void()+
  scale_fill_manual(values = c("white",brewer.pal(8,"Dark2")[c(1,2,3,4,6)]))

Figure.F <- Figure.Body %>%
  aplot::insert_right(Figure.R1,width = 0.05) %>%
  aplot::insert_right(Cg.group,width = 0.05) %>%
  aplot::insert_right(Cg.cpg,width = 0.05) %>%
  aplot::insert_right(Figure.R2,width = 0.2)

ggsave(Figure.F,filename = "Methylation.BCAA.Figure.Representative.pdf",height = 8,width = 10)
ggsave(Figure.F,filename = "Methylation.BCAA.Figure.Representative.Legend.pdf",height = 16,width = 10)

######## ATAC ########
setwd("K:/TCGA/Anlysis/BCAA")
BCAA.ID <- read.csv("BCAA.ID.csv")
setwd("K:/TCGA/Cancer.ATAC")
Sequencing.Data <- read.csv("Cancer.ATAC.Metadata.csv")
ATAC.BCAA <- data.frame()
for (project in unique(Sequencing.Data$cohort)) {
  Anno <- read.csv(paste("ATAC.",project,".PeakChIPSeeker.csv",sep = ""))
  SNHGs.Peak <- Anno %>% filter(SYMBOL %in% BCAA.ID$SYMBOL) %>%
    dplyr::select(Peak_name,SYMBOL) %>% unique()
  
  FPKM.ATAC <- read.csv(paste("ATAC.",project,".FPKM.csv",sep = ""),row.names = 1)
  ATAC.BCAA <- FPKM.ATAC[SNHGs.Peak$Peak_name,] %>% data.frame() %>% 
    rownames_to_column("Peak_name") %>%
    gather(Sample,FPKM,-Peak_name) %>%
    merge(.,SNHGs.Peak,by="Peak_name") %>%
    mutate(Cancer = project) %>%
    rbind.data.frame(ATAC.BCAA)
}
setwd("K:/TCGA/Anlysis/BCAA")
write.csv(ATAC.BCAA,file = "BCAA.ATAC.FPKM.csv",row.names = F)

#### Data preparation ####
setwd("K:/TCGA/Cancer.ATAC")
Sequencing.Data <- read.csv("Cancer.ATAC.Metadata.csv")
Cancer.ATAC <- data.frame()
for (name in BCAA.ID$SYMBOL) {
  Cancer.ATAC <- rbind(Cancer.ATAC,data.frame(SYMBOL=name,Cancer=unique(Sequencing.Data$cohort),MeanFPKM=0))
}

setwd("K:/TCGA/Anlysis/BCAA")
ATAC.BCAA <- read.csv(file = "BCAA.ATAC.FPKM.csv")
write.csv(ATAC.BCAA,file = "ATAC.BCAA.FPKM.csv",row.names = F)

#### Figure ####
ATAC.BCAA <- read.csv(file = "ATAC.BCAA.FPKM.csv")
ATAC.BCAA.FPKM <- ATAC.BCAA %>% group_by(SYMBOL,Cancer) %>%
  summarise(MeanFPKM = mean(FPKM),
            SDFPKM=sd(FPKM),
            CVFPKM=round(MeanFPKM/SDFPKM,4)*100) %>% ungroup() %>%
  mutate(SYMBOL=factor(SYMBOL,levels = rev(BCAA.ID$SYMBOL)))

write.csv(ATAC.BCAA.FPKM,"ATAC.BCAA.FPKM.MeanCV.csv",row.names = F)


library(RColorBrewer)
p1 <- ggplot(data=ATAC.BCAA.FPKM,aes(x=Cancer,
                                     y=SYMBOL,
                                     fill=MeanFPKM,
                                     size=MeanFPKM))+
  geom_point(shape=22)+
  scale_size_continuous(range = c(1,5))+
  ggthemes::theme_few()+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9,angle = 90,vjust = 0.5),
    axis.text.y = element_text(size=9),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="",y="",size="Mean of FPKM",fill="CV of FPKM")+
  scale_fill_gradientn(colors = brewer.pal(11,"RdYlBu")[4:1])

ggsave(p1,filename = "ATAC.BCAA.FPKM.MeanCV.pdf",height = 3,width = 6)
#### BCAT1 ####
ATAC.BCAA <- read.csv(file = "BCAA.ATAC.FPKM.csv")
BCAT1.ATAC <- ATAC.BCAA %>% filter(SYMBOL=="BCAT1") %>%
  mutate(Cancer=factor(Cancer,levels = rev(unique(sort(.$Cancer)))))

Immune.Col <- c("#293462","#a64942","#fe5f55","#fff1c1",
                "#5bd1d7","#348498","#004d61","#ff502f",
                "#e41749","#f5587b","#ff8a5c","#fff591",
                "#001871","#ff585d","#ffb549","#41b6e6",
                "#515bd4","#8134af","#dd2a7b","#feda77",
                "#96ceb4","#ffeead","#d9534f","#ffad60",
                "#05445c","#f2317f","#5c4f74","#040000",
                "#de4307","#f29c2b","#f6d04d","#8bc24c",
                "#fef4a9","#3b9a9c","#4bc2c5","#78fee0",
                "#fad3cf","#a696c8","#2470a0","#c50d66",
                "#f07810","#eec60a","#dd7777","#1687a7",
                "#014955","#dd0a35","#ED5485","#FFE869",
                "#bc8420","#F0B775","#D25565","#2E94B9")

p3<-ggplot(BCAT1.ATAC ,aes(x=FPKM,y=reorder(Cancer,FPKM),fill=Cancer))+
  #geom_violin()+
  geom_boxplot(color="black")+
  ggthemes::theme_few()+
  theme(legend.position = "none",
        #axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9,color = "black"),
        axis.text.y = element_text(size=9,color = "black"),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  #scale_x_sqrt()+
  scale_x_continuous(breaks = c(0,1,5,10,20,40,60,80),
                     labels = c(0,1,5,10,20,40,60,80),
                     trans = "sqrt")+
  scale_fill_manual(values = Immune.Col)+
  labs(y="",x="FPKM of BCAT1")
ggsave(p3,filename = "ATAC.BCAA.BCAT1.FPKM.pdf",height = 5,width = 2.3) 
#### BCAT2 ####
ATAC.BCAA <- read.csv(file = "BCAA.ATAC.FPKM.csv")
ATAC.BCAA %>% filter(SYMBOL=="BCAT2") %>%
  group_by(Cancer) %>% top_n(2,FPKM) %>% View()

BCAT2.ATAC <- ATAC.BCAA %>% filter(SYMBOL=="BCAT2") %>%
  mutate(Cancer=factor(Cancer,levels = rev(unique(sort(.$Cancer)))))

Immune.Col <- c("#293462","#a64942","#fe5f55","#fff1c1",
                "#5bd1d7","#348498","#004d61","#ff502f",
                "#e41749","#f5587b","#ff8a5c","#fff591",
                "#001871","#ff585d","#ffb549","#41b6e6",
                "#515bd4","#8134af","#dd2a7b","#feda77",
                "#96ceb4","#ffeead","#d9534f","#ffad60",
                "#05445c","#f2317f","#5c4f74","#040000",
                "#de4307","#f29c2b","#f6d04d","#8bc24c",
                "#fef4a9","#3b9a9c","#4bc2c5","#78fee0",
                "#fad3cf","#a696c8","#2470a0","#c50d66",
                "#f07810","#eec60a","#dd7777","#1687a7",
                "#014955","#dd0a35","#ED5485","#FFE869",
                "#bc8420","#F0B775","#D25565","#2E94B9")

p3<-ggplot(BCAT2.ATAC ,aes(x=FPKM,y=reorder(Cancer,FPKM),fill=Cancer))+
  #geom_violin()+
  geom_boxplot(color="black")+
  ggthemes::theme_few()+
  theme(legend.position = "none",
        #axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9,color = "black"),
        axis.text.y = element_text(size=9,color = "black"),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  #scale_x_sqrt()+
  scale_x_continuous(breaks = c(0,1,5,10,20,40,60,80),
                     labels = c(0,1,5,10,20,40,60,80),
                     trans = "sqrt")+
  scale_fill_manual(values = Immune.Col)+
  labs(y="",x="FPKM of BCAT2")
ggsave(p3,filename = "ATAC.BCAA.BCAT2.FPKM.pdf",height = 5,width = 2.3) 
######## NMF ############
#### Cluster => TPM => NMF ####
library(NMF)
library(tidyverse)
dir.create("NMF.BCAA")
BCAA.ID <- read.csv("BCAA.ID.csv")
#PhenoData <- read.csv("../../mRNA.Phenotype.csv") %>%
PhenoData <- read.csv("../mRNA.Phenotype.csv") %>%
  filter(Sample.Type %in% c("Primary Tumor",
                            "Primary Blood Derived Cancer - Peripheral Blood",
                            "Solid Tissue Normal")) %>%
  mutate(SampleType = ifelse(Sample.Type == "Solid Tissue Normal","Normal","Tumor")) %>%
  mutate(Project.ID=str_remove(Project.ID,"TCGA-")) %>%
  filter(SampleType == "Tumor")

for (project in unique(PhenoData$Project.ID)) {
  print(project)
  mid.Pheno <- PhenoData %>% filter(Project.ID == project)
  
  #filepath = "../../mRNA.PanCancer.Exp/"
  filepath = "../mRNA.PanCancer.Exp/"
  filename = paste("TCGA-",project,".mRNA.TPM.csv",sep = "")
  data <- read.csv(file.path(filepath,filename),check.names = F,row.names = 1)
  BCAA.Data <- data[BCAA.ID$ENSEMBL,] %>% t() %>% 
    data.frame(check.names = F) %>%
    rownames_to_column("Sample.ID") %>%
    filter(Sample.ID %in% mid.Pheno$Sample.ID) %>%
    column_to_rownames("Sample.ID") %>%
    t()
  
  ranks <- 2:10
  seed <- 20230206
  figure.name <- paste(project,".NMF.2-10.pdf",sep = "")
  file.name <- paste(project,".NMF.2-10.rds",sep = "")
  result = nmf(BCAA.Data,
               ranks,
               method="brunet", 
               nrun=5,
               seed =seed)
  saveRDS(result,file = file.path("NMF.BCAA",file.name))
  pdf(file = figure.name)
  plot(result)
  dev.off()
}




#### Cluster => Confirm ####
setwd("K:/TCGA/Anlysis/BCAA")
library(tidyverse)
files <- list.files(path = "NMF.BCAA/",pattern = "rds$")
NMF.Rank <- data.frame()
for (file in files) {
  data <- readRDS(file.path("NMF.BCAA/",file))
  NMF.Rank <- data$measures[,c("cophenetic","rank")] %>%
    data.frame(check.names = F) %>%
    mutate(Cancer=file %>% str_remove_all("\\..*")) %>%
    rbind.data.frame(NMF.Rank)
}
write.csv(NMF.Rank,file = "NMF.BCAA/NMF.BCAA.Cophenetic.csv",row.names = F)
ggplot(NMF.Rank,aes(rank,cophenetic,color=Cancer))+
  geom_point(size=2)+
  geom_line()+
  facet_wrap(~Cancer,ncol = 6)+
  ggthemes::theme_few()+
  theme(legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9),
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))

######## miRNA => Network Pearson ########
#### multiMiR ####
library(multiMiR)
library(plyr)
BCAA.ID <- read.csv("BCAA.ID.csv")
multiMiR.Data <- data.frame()
for (symbol in BCAA.ID$SYMBOL) {
  miRNA.Target <- multiMiR::get_multimir(target = symbol,table="all",org = 'hsa',summary = TRUE)
  multiMiR.Data <- miRNA.Target@data %>% rbind.fill(multiMiR.Data)
}
write.csv(multiMiR.Data,file = "miRNA.BCAA.multiMiR.csv")
#### TCGA => miRNA + mRNA => Correlation ####
library(miRBaseVersions.db)
library(vegan)
Inter.miRNA.mRNA <- read.csv(file = "miRNA.BCAA.multiMiR.csv",row.names = 1)
Inter.miRNA.mRNA <- Inter.miRNA.mRNA[,c(2:6)] %>% data.frame() %>%
  unique()
BCAA.ID <- read.csv("BCAA.ID.csv")

dirs <- list.dirs("../../Clinical.XENA/")
dirs <- dirs[str_detect(dirs,"//")]
miRNA.mRNA.Pearson <- data.frame()
for (dir in dirs) {
  project <- dir %>% str_remove(".*\\/")
  print(project)
  mRNA.Filename <- paste("TCGA-",project,".mRNA.TPM.csv",sep = "")
  mRNA.Exp <- read.csv(file.path("../../mRNA.PanCancer.Exp/",mRNA.Filename),row.names = 1,check.names = F)
  mRNA.Exp <- log2(mRNA.Exp+1)
  BCAA.Exp <- mRNA.Exp[BCAA.ID$ENSEMBL,] %>% data.frame(check.names = F) %>%
    rownames_to_column("ENSEMBL") %>%
    merge(BCAA.ID,by="ENSEMBL") %>%
    dplyr::select(-ENSEMBL,-ENTREZ) %>%
    remove_rownames() %>%
    column_to_rownames("SYMBOL") %>%
    t() %>% data.frame(check.names = F) %>%
    rownames_to_column("Sample.ID")
  
  miRNA.Filename <- paste("TCGA-",project,".mature.miRNA.RPM.csv",sep = "")
  miRNA.Exp <- read.csv(file.path("../../miRNA.Mature.PanCancer.Exp/",miRNA.Filename),
                        check.names = F,row.names = 1)
  
  miRNA.mRNA.Exp <- miRNA.Exp %>% t() %>%
    data.frame(check.names = F) %>%
    rownames_to_column("Sample.ID") %>%
    merge(BCAA.Exp,by="Sample.ID")
  
  if (nrow(miRNA.mRNA.Exp) > 0) {
    for (gene in BCAA.ID$SYMBOL) {
      mid.gene.miRNA <- Inter.miRNA.mRNA %>% filter(target_symbol==gene)
      
      if (nrow(mid.gene.miRNA) > 0) {
        print("Samples overlapped")
        for (miRNA.id in mid.gene.miRNA$mature_mirna_id) {
          if (c(miRNA.id) %in% colnames(miRNA.mRNA.Exp)) {
            test <- cor.test(miRNA.mRNA.Exp[[miRNA.id]],
                             miRNA.mRNA.Exp[[gene]],
                             method = "pearson",
                             exact = F,
                             alternative = "two.sided")
            miRNA.mRNA.Pearson <- data.frame(miRNA=miRNA.id,
                                             Target=gene,
                                             Rho=test$estimate,
                                             Pvalue=test$p.value,
                                             Cancer=project) %>%
              rbind.data.frame(miRNA.mRNA.Pearson)
          }
        }
      }else(
        print("No overlap")
      )
    }
  }
}

miRNA.mRNA.Pearson2 <- miRNA.mRNA.Pearson %>%
  na.omit() %>%
  mutate(Label=if_else(Pvalue>0.05,"p>0.05",
                       if_else(Rho>0,"Positive","Negative"))) 

write.csv(miRNA.mRNA.Pearson2,file = "miRNA.TCGA.miRNA.mRNA.Pearson.csv",row.names = F)
#### Figure ####
miRNA.mRNA.Pearson2 <- read.csv(file = "miRNA.TCGA.miRNA.mRNA.Pearson.csv")
library(plyr)
miRNA.mRNA.Pearson3 <- miRNA.mRNA.Pearson2 %>%
  data.frame(check.names = F) %>%
  group_by(Target,miRNA,Label) %>%
  dplyr::summarise(Count=n()) %>%
  ungroup() %>%
  filter(Label=="Negative") %>%
  arrange(desc(Count)) %>%
  dplyr::top_n(30,Count) %>%
  dplyr::select(Target,miRNA) %>%
  arrange(Target,miRNA) %>%
  unique() %>%
  merge(miRNA.mRNA.Pearson2,by=c("Target","miRNA")) %>%
  mutate(Yaxis=paste(miRNA,Target,sep = "->")) %>%
  arrange(Target,miRNA,Cancer) %>% 
  unique() %>%
  mutate(Yaxis=factor(Yaxis,levels = rev(unique(.$Yaxis))))
write.csv(miRNA.mRNA.Pearson3,file = "miRNA.TCGA.miRNA.mRNA.Pearson.Data.csv",row.names = F)

library(RColorBrewer)
Figure.Body <- ggplot()+
  geom_point(data=miRNA.mRNA.Pearson3,
             aes(x=Cancer,y=Yaxis,fill=Rho,size=-log10(Pvalue)),
             shape=22,color="white",stroke=1)+
  geom_point(data=miRNA.mRNA.Pearson3 %>% filter(Pvalue <= 0.05),
             aes(x=Cancer,y=Yaxis,fill=Rho,size=-log10(Pvalue)),
             shape=22,color="black",stroke=1)+
  scale_size_continuous(range = c(2,5))+
  ggthemes::theme_few()+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9,angle = 90,vjust = 0.5),
    axis.text.y = element_text(size=9),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="",y="",size="p value",fill="Pearson'CC")+
  scale_fill_gradient2(high = "#BC3C29FF", mid = "white", low = "#0072B5FF",midpoint = 0)
#"#0072B5FF","#BC3C29FF"  
Figure.R1 <- miRNA.mRNA.Pearson3 %>%
  dplyr::select(Yaxis,Target) %>%
  unique() %>%
  ggplot(aes(x=1,y=Yaxis,fill=Target))+
  geom_col(width = 1,color="black")+
  scale_fill_manual("Symbol",values = brewer.pal(12,"Paired"))+
  theme_void()

Figure.R2 <- miRNA.mRNA.Pearson3 %>%
  group_by(Target,miRNA,Label) %>%
  dplyr::summarise(Count=n()) %>% #arrange(desc(Count))
  ungroup() %>%
  mutate(Yaxis=paste(miRNA,Target,sep = "->")) %>%
  dplyr::select(-Target,-miRNA) %>%
  spread(Label,Count,fill=0) %>%
  gather(Correlation,Count,-Yaxis) %>%
  mutate(Correlation=factor(Correlation,
                            levels = c("p>0.05","Negative","Positive"))) %>%
  ggplot(aes(x=Count,y=Yaxis,fill=Correlation))+
  geom_col(width = 1,color="black")+
  ggthemes::theme_few()+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9,angle = 90,vjust = 0.5),
    axis.text.y = element_blank(),
    axis.ticks.y=element_blank(),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="No. of cancers",y="")+
  scale_fill_manual(values = c("grey90","#66C2A5","#BC3C2999"))

Figure.F <- Figure.Body %>%
  aplot::insert_right(Figure.R1,width = 0.05) %>%
  aplot::insert_right(Figure.R2,width = 0.3)

ggsave(Figure.F,filename = "miRNA.TCGA.miRNA.mRNA.Pearson.pdf",height = 6,width = 9)
ggsave(Figure.F,filename = "miRNA.TCGA.miRNA.mRNA.Pearson.Legend.pdf",height = 12,width = 10)

######## BCAA score ssGSEA #######
library(GSVA)
library(GSEABase)
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
BCAA.ID <- read.csv("BCAA.ID.csv")
BCAA <- list()
BCAA[["BCAA"]] <- BCAA.ID$ENSEMBL

dir.create("BCAA.ssGSEA")
BCAA.ssGSEA <- data.frame()
files <- list.files("../../mRNA.PanCancer.Exp/",pattern = "csv$")
for (file in files) {
  project <- file %>% str_remove_all("\\..*") %>% str_remove_all("TCGA-")
  print(project)
  data <- read.csv(file.path("../../mRNA.PanCancer.Exp/",file),row.names = 1,check.names = F)
  data2 <- log2(data+1)
  
  BCAA.ssGSEA <- gsva(as.matrix(data2), BCAA, 
                      method = "ssgsea", 
                      kcdf = "Gaussian") %>%
    t() %>%
    data.frame(check.names = F) %>%
    rownames_to_column("Sample.ID") %>%
    mutate(Cancer=project) %>%
    rbind.data.frame(BCAA.ssGSEA)
  #write.csv(res.ssgsea,file = paste("HALLMARK.ssGSEA.Pancancer/TCGA-",project,".ssGSEA.csv",sep = ""))
}
write.csv(BCAA.ssGSEA,"BCAA.ssGSEA/BCAA.ssGSEA.csv",row.names = F)
######## BCAA score Survival ####
library(survival)
library(survminer)
library(tidyverse)
library(export)

files <- list.files("../../mRNA.PanCancer.Exp/",pattern = "csv$")
BCAA.Survival <- data.frame()
for (file in files) {
  project <- file %>% str_remove_all("\\..*") %>% str_remove_all("TCGA-")
  print(project)
  #data <- read.csv(file.path("../../mRNA.PanCancer.Exp/",file),row.names = 1,check.names = F)
  #TPM <- log2(data+1)
  TPM <- BCAA.ssGSEA %>% filter(Cancer==project)
  
  phenotype <- read.csv(paste("../../Clinical.XENA/",project,"/TCGA-",project,".GDC_phenotype.tsv",sep = ""),sep = "\t",header = T) %>% 
    dplyr::filter(sample_type.samples == "Primary Tumor" | 
                    sample_type.samples == "Primary Blood Derived Cancer - Peripheral Blood")
  
  #### survival datasets
  survival.data <- read.csv(paste("../../Clinical.XENA/",project,"/TCGA-",project,".survival.tsv",sep = ""),sep = "\t",header = T,row.names = 1,check.names = F) %>%
    rownames_to_column("Sample.ID")
  
  Inter.Samples <- intersect(intersect(phenotype$submitter_id.samples,TPM$Sample.ID),survival.data$Sample.ID)
  #### survival analysis
  TPM.log.Target <- TPM %>% filter(Sample.ID %in% Inter.Samples) %>%
    merge(survival.data,by="Sample.ID")
  res.cut <- surv_cutpoint(TPM.log.Target,
                           time = paste("OS",".time",sep = ""),
                           event = "OS",
                           variables = "BCAA")
  #res.cat <- surv_categorize(res.cut)
  
  theme.sur <- ggthemes::theme_few()+
    theme(legend.position = "none",
          axis.text = element_text(color = "black"),
          axis.text.x = element_text(size=15),
          axis.text.y = element_text(size=15),
          axis.title = element_text(size=18),
          panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))
  
  dat <- surv_categorize(res.cut) #%>% data.frame() %>% mutate(Sample=SNHGs.ICI.data$Sample.ID)
  fit <- survfit(Surv(OS.time, OS) ~ BCAA,
                 data = dat)
  
  index="OS"
  diff=survdiff(as.formula(paste("Surv(",index,".time,",index,") ~ ","BCAA",sep = "")),data = dat)
  pValue=1-pchisq(diff$chisq,df=1)
  
  cox <- coxph(as.formula(paste("Surv(OS.time, OS) ~ ","BCAA",sep = "")), data = TPM.log.Target)
  coxSummary = summary(cox)
  
  library(ggthemes)
  p<-ggsurvplot(fit,pval =TRUE, data = dat, 
                surv.median.line = "hv",
                legend.title = "BCAA",
                conf.int.style = "step",
                xlab = "Time in days",
                #break.time.by = 500,
                risk.table = "abs_pct",
                #risk.table.y.text.col = T,
                #risk.table.y.text = FALSE,
                legend.labs = c("High", "Low"),
                #pval = TRUE,
                conf.int = TRUE,
                palette = "Set1",
                #ggtheme = theme.sur,
                risk.table.y.text.col = T,
                risk.table.y.text = F,
                ggtheme = theme.sur)
  print(p)
  #dev.off()
  graph2pdf(file=paste("BCAA.ssGSEA/TCGA-",project,".OS.BCAA.KM.pdf",sep = ''),height = 5,width = 4)
  
  BCAA.Survival=rbind(BCAA.Survival,
                      data.frame(ENSEMBL="BCAA",
                                 KM.pvalue=pValue,
                                 HR=coxSummary$conf.int[,"exp(coef)"],
                                 HR.95L=coxSummary$conf.int[,"lower .95"],
                                 HR.95H=coxSummary$conf.int[,"upper .95"],
                                 HR.pvalue=coxSummary$coefficients[,"Pr(>|z|)"],
                                 Cancer=project))
  
}
write.csv(BCAA.Survival,file = "BCAA.ssGSEA/BCAA.Survival.csv",row.names = F)
#### Figure ####
BCAA.ssGSEA <- read.csv("BCAA.ssGSEA/BCAA.ssGSEA.csv")

library(survival)
library(survminer)
library(tidyverse)

BCAA.Survival <- read.csv(file = "BCAA.ssGSEA/BCAA.Survival.csv")
BCAA.Survival <- BCAA.Survival %>%
  mutate(Colors=if_else(KM.pvalue >0.05 & HR.pvalue > 0.05,"KM p>0.05 & COX p>0.05",
                        if_else(KM.pvalue >0.05 & HR.pvalue <= 0.05,"KM p>0.05 & COX p≤0.05",
                                if_else(KM.pvalue <=0.05 & HR.pvalue > 0.05,"KM p≤0.05 & COX p>0.05","KM p≤0.05 & COX p≤0.05")))) %>%
  arrange(HR) %>%
  mutate(Color=if_else(HR.pvalue > 0.05,"P>=0.05","P<=0.05")) 
BCAA.Survival$Cancer <- factor(BCAA.Survival$Cancer,levels = unique(BCAA.Survival$Cancer))

library(RColorBrewer)
p1 <- ggplot()+
  geom_segment(data = BCAA.Survival,
               aes(x=0,xend=log10(HR),y=Cancer,yend=Cancer))+
  geom_point(data = BCAA.Survival,
             aes(x=log10(HR),y=Cancer,color=Color),size=3)+
  ggthemes::theme_few()+
  theme(#legend.position = "top",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9,angle = 90,vjust = 0.5),
        axis.text.y = element_text(size=9),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="log10(HR)",y="",color="")+
  scale_color_manual(values = c(brewer.pal(8,"Set1")[1:2],"grey50"))+
  scale_x_continuous(limits = c(-2,3))+
  geom_vline(xintercept = c(0),linetype=2,color="black")+
  coord_flip()

p2 <- ggplot()+
  geom_segment(data = BCAA.Survival,
               aes(x=0,xend=-log10(KM.pvalue),
                   y=Cancer,yend=Cancer))+
  geom_point(data = BCAA.Survival,
             aes(x=-log10(KM.pvalue),
                 y=Cancer,
                 color=-log10(KM.pvalue),size=-log10(KM.pvalue)))+
  ggthemes::theme_few()+
  theme(#legend.position = "top",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9,angle = 90,vjust = 0.5),
    axis.text.y = element_text(size=9),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="-log10(P)",y="",color="")+
  #scale_color_manual(values = c(brewer.pal(8,"Set1")[1:2],"grey50"))+
  scale_color_gradientn(colors = brewer.pal(11,"RdBu")[4:1])+
  scale_x_continuous(limits = c(0,13))+
  geom_vline(xintercept = c(-log10(0.05)),linetype=2,color="black")+
  coord_flip()

ggsave(p2,filename = "BCAA.ssGSEA.Survival2.pdf",height = 2,width = 7)
ggsave(p2,filename = "BCAA.ssGSEA.Survival2.Legend.pdf",height = 12,width = 7)

ggsave(p,filename = "BCAA.ssGSEA.Survival.pdf",height = 2,width = 6) 
#ggsave(p,filename = "BCAA.Survival.Legend.pdf",height = 2,width = 16)

######## BCAA score Normal Tumor ####
setwd("K:/TCGA/Anlysis/BCAA")
library(tidyverse)
BCAA.ssGSEA <- read.csv("BCAA.ssGSEA/BCAA.ssGSEA.csv")
BCAA.ID <- read.csv("BCAA.ID.csv")

PhenoData <- read.csv("../../mRNA.Phenotype.csv") %>%
  filter(Sample.Type %in% c("Primary Tumor",
                            "Primary Blood Derived Cancer - Peripheral Blood",
                            "Solid Tissue Normal")) %>%
  mutate(SampleType = ifelse(Sample.Type == "Solid Tissue Normal","Normal","Tumor")) %>%
  mutate(ProjectID=str_remove_all(Project.ID,"TCGA-")) %>%
  dplyr::select(ProjectID,SampleType,Sample.ID) 
PhenoData2 <- PhenoData %>% filter(SampleType=="Normal")
BCAA.SCORE <- merge(BCAA.ssGSEA,PhenoData %>% filter(ProjectID %in% PhenoData2$ProjectID),by="Sample.ID") %>%
  data.frame(check.names = F) 
library(ggbeeswarm)
library(RColorBrewer)
library(ggpubr)
data2 <- compare_means(formula = BCAA~SampleType,data = BCAA.SCORE,method = "wilcox.test",
                       group.by="Cancer") %>%
  data.frame(check.names = F)
write.csv(data2,file = "BCAA.ssGSEA.Tumor.Normal.Wilcox.csv",row.names = F)

BCAA.SCORE2 <- BCAA.SCORE %>%
  group_by(Cancer,SampleType) %>%
  summarise(Mean=mean(BCAA)) %>%
  ungroup() %>%
  spread(SampleType,Mean,fill=NA) %>%
  mutate(SHAPE=if_else(Normal > Tumor,"Higher in normal","Higher in tumor")) %>%
  mutate(Normal=-Normal) %>%
  gather(Type,BCAA,Normal:Tumor) #%>%
  mutate(Cancer=factor(Cancer,levels = rev(unique(sort(.$Cancer)))))
  
p1 <- ggplot(data=BCAA.SCORE2,
       aes(x=BCAA,y=Cancer))+
  geom_rect(aes(xmin=0, xmax=Inf, ymin=-Inf, ymax=Inf),fill="#E64B35FF",alpha = 0.005)+
  geom_rect(aes(xmin=-Inf, xmax=0, ymin=-Inf, ymax=Inf),fill="#00A087FF",alpha = 0.005)+
  geom_segment(aes(x=0,xend=BCAA,y=Cancer,yend=Cancer,
                   color=Type))+
  geom_point(aes(fill=Type),color="black",
             size=3,shape=21)+
  #geom_text(data=data2,aes(x=10.5,y=Cancer,label=p.signif))+
  #stat_ellipse(aes(fill=Group),type = "norm", geom ="polygon",alpha=0.2,color=NA)+
  #geom_beeswarm(aes(color=SampleType),size=1,dodge.width=0.75,priority = "ascending")+
  labs(x="Mean of BCAA score",y="")+
  #facet_wrap(~Phenotype,scale="free_y",ncol = 6)+
  ggthemes::theme_few()+
  theme(#legend.position = "top",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  #scale_shape_manual(values = c(21,22))+
  scale_color_manual("",values = brewer.pal(9,"Set1")[c(2:1)])+
  scale_fill_manual("",values = brewer.pal(9,"Set1")[c(2:1)])+
  geom_vline(xintercept = c(0),color="black")+
  scale_x_continuous(limits = c(-10,10))
  #coord_cartesian(clip = "off")+
  #geom_text(data = data2,aes(x=11,y=Cancer,label=paste("P=",round(p.format,3),sep = "")))>

data2$label <- paste("P=",data2$p.format,sep = "")
p2 <- ggplot(data2,aes(x=1,y=Cancer))+
  geom_text(aes(label=p.signif))+
  theme_void()
p12 <- p1 %>% aplot::insert_right(p2,width = 0.25)
ggsave(p12,filename = "BCAA.ssGSEA.Tumor.Normal2.pdf",height = 4,width = 3.5)


######## BCAA score Index Spearman ####
setwd("K:/TCGA/Anlysis/BCAA")
library(tidyverse)
library(readxl)
BCAA.ssGSEA <- read.csv("BCAA.ssGSEA/BCAA.ssGSEA.csv") %>%
  mutate(ID=substr(Sample.ID,1,15))
#### ESTIMATE ####
library(readxl)
data <- read_xlsx("../../PanCancer.DiverseIndex.xlsx",sheet = 2)
data2<-data[,-1]
Index.data2 <- merge(BCAA.ssGSEA,data2,by="ID") %>% unique() %>%
  data.frame() %>%
  gather(Attribute,Score,StromalScore:ESTIMATEScore)

Spearman.Data <- data.frame()
for (project in unique(Index.data2$Cancer)) {
  for (score in c("StromalScore","ImmuneScore","ESTIMATEScore")) {
    middata <- Index.data2 %>% filter(Attribute == score) %>%
      filter(Cancer==project)
    for (gene in c("BCAA")) {
      test <- cor.test(middata[[gene]],middata[["Score"]],
                       method = "spearman",
                       alternative = "two.sided",
                       exact = F)
      Spearman.Data <- data.frame(Index=score,
                                  SCORE=gene,
                                  Rho=test$estimate,
                                  Pvalue=test$p.value,
                                  Cancer=project) %>%
        rbind.data.frame(Spearman.Data)
    }
  }
}
write.csv(Spearman.Data,file = "BCAA.ESTIMATE.Spearman.csv",row.names = F)
#### IPS ####
library(readxl)
data <- read_xlsx("../../PanCancer.DiverseIndex.xlsx",sheet = 3)
data2<-data[,-1]
Index.data2 <- merge(BCAA.ssGSEA,data2,by="ID") %>% unique() %>%
  data.frame() %>%
  gather(Attribute,Score,MHC:IPS)

for (project in unique(Index.data2$Cancer)) {
  for (score in c("AZ","CP","EC","MHC","SC","IPS")) {
    middata <- Index.data2 %>% filter(Attribute == score) %>%
      filter(Cancer==project)
    for (gene in c("BCAA")) {
      test <- cor.test(middata[[gene]],middata[["Score"]],
                       method = "spearman",
                       alternative = "two.sided",
                       exact = F)
      Spearman.Data <- data.frame(Index=score,
                                  SCORE=gene,
                                  Rho=test$estimate,
                                  Pvalue=test$p.value,
                                  Cancer=project) %>%
        rbind.data.frame(Spearman.Data)
    }
  }
}
write.csv(Spearman.Data,file = "BCAA.IPS.Spearman.csv",row.names = F)
#### TMB MSI NEO HRD ####
for (index in c(4,6,7,10)) {
  data <- read_xlsx("../../PanCancer.DiverseIndex.xlsx",sheet = index)
  data2<-data[,-1]
  Index.data2 <- merge(BCAA.ssGSEA,data2,by.x="ID",by.y="SampleName") %>% 
    unique() %>% data.frame()
  for (project in unique(Index.data2$Cancer)) {
    middata <- Index.data2 %>% #filter(Attribute == score) %>%
      filter(Cancer==project)
    for (gene in c("BCAA")) {
      test <- cor.test(middata[[gene]],middata[[colnames(data2)[1]]],
                       method = "spearman",
                       alternative = "two.sided",
                       exact = F)
      Spearman.Data <- data.frame(Index=colnames(data2)[1],
                                  SCORE=gene,
                                  Rho=test$estimate,
                                  Pvalue=test$p.value,
                                  Cancer=project) %>%
        rbind.data.frame(Spearman.Data)
    }
  }
}
write.csv(Spearman.Data,file = "BCAA.TMB.MSI.NEO.HRD.Spearman.csv",row.names = F)
#### MATH PURITY PLOIDY ####
for (index in c(5,8,9)) {
  data <- read_xlsx("../../PanCancer.DiverseIndex.xlsx",sheet = index)
  data2<-data[,-1]
  Index.data2 <- merge(BCAA.ssGSEA,data2,by.x="ID",by.y="SampleName") %>% 
    unique() %>% data.frame() 
  for (project in unique(Index.data2$Cancer)) {
    middata <- Index.data2 %>% #filter(Attribute == score) %>%
      filter(Cancer==project)
    for (gene in c("BCAA")) {
      test <- cor.test(middata[[gene]],middata[[colnames(data2)[1]]],
                       method = "spearman",
                       alternative = "two.sided",
                       exact = F)
      Spearman.Data <- data.frame(Index=colnames(data2)[1],
                                  SCORE=gene,
                                  Rho=test$estimate,
                                  Pvalue=test$p.value,
                                  Cancer=project) %>%
        rbind.data.frame(Spearman.Data)
    }
  }
}
write.csv(Spearman.Data,file = "BCAA.MATH.PURITY.PLOIDY.Spearman.csv",row.names = F)
#### RNA stem score ####
for (index in c(13)) {
  data <- read_xlsx("../../PanCancer.DiverseIndex.xlsx",sheet = index)
  data2<-data[,-1]
  Index.data2 <- merge(BCAA.ssGSEA,data2,by.x="ID",by.y="SampleName") %>% 
    unique() %>% data.frame() 
  for (project in unique(Index.data2$Cancer)) {
    middata <- Index.data2 %>% #filter(Attribute == score) %>%
      filter(Cancer==project)
    for (gene in c("BCAA")) {
      test <- cor.test(middata[[gene]],middata[[colnames(data2)[1]]],
                       method = "spearman",
                       alternative = "two.sided",
                       exact = F)
      Spearman.Data <- data.frame(Index=colnames(data2)[1],
                                  SCORE=gene,
                                  Rho=test$estimate,
                                  Pvalue=test$p.value,
                                  Cancer=project) %>%
        rbind.data.frame(Spearman.Data)
    }
  }
}
write.csv(Spearman.Data,file = "BCAA.RNAss.Spearman.csv",row.names = F)


##### Figure ####
Spearman.Data <- read.csv("BCAA.RNAss.Spearman.csv")
write.csv(Spearman.Data,file = "BCAA.ssGSEA.Index.Spearman.csv",row.names = F)

Spearman.Data$Index <- factor(Spearman.Data$Index,levels = unique(Spearman.Data$Index))
Spearman.Data$Cancer <- factor(Spearman.Data$Cancer,levels = sort(unique(Spearman.Data$Cancer)))
library(RColorBrewer)
Cancer.Col <- c(brewer.pal(9,"Set1"),
                brewer.pal(8,"Set2")[1:7],
                brewer.pal(12,"Set3")[1:11],
                brewer.pal(8,"Dark2")[1:6])
p12 <- ggplot()+
  geom_tile(data=Spearman.Data,
            aes(x=Cancer,y=Index,fill=Rho),
            color="white")+
  geom_point(data=Spearman.Data %>% 
               filter(Pvalue <= 0.05),
             aes(x=Cancer,y=Index,size=-log10(Pvalue)),
             shape=1,color="black")+
  geom_point(data=Spearman.Data %>% 
               filter(Pvalue > 0.05),
             aes(x=Cancer,y=Index,size=-log10(Pvalue)),
             shape=4,color="white")+
  scale_size_continuous(range = c(1,4))+
  ggthemes::theme_few()+
  theme(legend.position = "bottom",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9,angle = 90,vjust = 0.5,hjust = 1),
        axis.text.y = element_text(size=9),
        #axis.ticks.x = element_blank(),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="",y="",fill="Spearman'CC")+
  #geom_vline(xintercept = c(0),linetype=2)+
  #theme(strip.text.x = element_text(size = 12))+
  #ggthemes::scale_fill_gradient2_tableau(palette = "Red-Blue Diverging")
  scale_fill_gradient2(low = "#0066CC", mid = "white", high = "#c72e29")

ggsave(p12,filename = "BCAA.ssGSEA.Index.Spearman2.pdf",height = 4,width = 8)

######## BCAA score ImmuneCellAI TIMER ... Spearman ####
setwd("K:/TCGA/Anlysis/BCAA")
library(tidyverse)
BCAA.ssGSEA <- read.csv("BCAA.ssGSEA/BCAA.ssGSEA.csv")
BCAA.ID <- read.csv("BCAA.ID.csv")
TIMER2 <- read.csv("../../TIMER2/infiltration_estimation_for_tcga.csv",check.names = F)

PhenoType <- read.csv("../../mRNA.Phenotype.csv")
PhenoType <- PhenoType %>% filter(Sample.Type == "Primary Tumor" | 
                                    Sample.Type == "Solid Tissue Normal" | 
                                    Sample.Type == "Primary Blood Derived Cancer - Peripheral Blood") %>% 
  dplyr::select(Project.ID,Sample.ID,Sample.Type) %>%
  mutate(Type=ifelse(Sample.Type=="Solid Tissue Normal","Normal","Tumor")) %>%
  mutate(Cancer=str_remove(Project.ID,"TCGA-"))

Immune.Data <- TIMER2 %>% 
  #dplyr::select(c("cell_type",ends_with("CIBERSORT"))) %>%
  remove_rownames() %>%
  column_to_rownames("cell_type") %>%
  #magrittr::set_colnames(str_remove_all(colnames(.),"_CIBERSORT")) %>%
  rownames_to_column("SampleID")

IMMUNECELL.SYMBOL.SPEARMAN <- data.frame()
for (project in unique(PhenoType$Cancer)) {
  print(project)
  
  mid.Pheno <- PhenoType %>% filter(Cancer==project) %>% 
    filter(Type=="Tumor")
  
  Exp.Filename <- paste("TCGA-",project,".mRNA.TPM.csv",sep = "")
  Exp.Data <- read.csv(file.path("../../mRNA.PanCancer.Exp/",Exp.Filename),
                       row.names = 1,check.names = F)
  Exp.Data2 <- log2(Exp.Data+1) %>% data.frame(check.names = F) %>%
    rownames_to_column("ENSEMBL") %>%
    merge(BCAA.ID,by="ENSEMBL") %>%
    dplyr::select(-ENSEMBL,-ENTREZ) %>%
    column_to_rownames("SYMBOL") %>%
    t() %>% data.frame(check.names = F) %>%
    rownames_to_column("Sample.ID") %>%
    filter(Sample.ID %in% mid.Pheno$Sample.ID) %>%
    merge(BCAA.ssGSEA,by = "Sample.ID") %>%
    mutate(SampleID=substr(Sample.ID,1,15)) %>%
    merge(Immune.Data,by="SampleID")
  
  for (gene1 in c("BCAA")) {
    for (gene2 in colnames(Immune.Data)[2:ncol(Immune.Data)]) {
      if (sum(! is.na(Exp.Data2[[gene2]])) >= 5) {
        test <- cor.test(Exp.Data2[[gene1]],Exp.Data2[[gene2]],
                         method = "spearman",
                         exact = F,
                         alternative = "two.sided")
        IMMUNECELL.SYMBOL.SPEARMAN <- data.frame(IMMUNECELL=gene2,
                                                 SYMBOL=gene1,
                                                 Rho=test$estimate,
                                                 Pvalue=test$p.value,
                                                 Cancer=project) %>%
          rbind.data.frame(IMMUNECELL.SYMBOL.SPEARMAN)
      }else{
        IMMUNECELL.SYMBOL.SPEARMAN <- data.frame(IMMUNECELL=gene2,
                                                 SYMBOL=gene1,
                                                 Rho=NA,
                                                 Pvalue=NA,
                                                 Cancer=project) %>%
          rbind.data.frame(IMMUNECELL.SYMBOL.SPEARMAN)
      }
    }
  }
}
write.csv(IMMUNECELL.SYMBOL.SPEARMAN,file = "BCAA.ssGSEA.IMMUNECELL.SPEARMAN.csv",row.names = F)
#### Figure 1 ####
IMMUNECELL.SYMBOL.SPEARMAN <- read.csv("BCAA.ssGSEA.IMMUNECELL.SPEARMAN.csv")
IMMUNECELL.SYMBOL.SPEARMAN2 <- IMMUNECELL.SYMBOL.SPEARMAN %>% 
  filter(str_detect(IMMUNECELL,"score")) %>%
  mutate(Method=str_remove_all(IMMUNECELL,".*_")) %>%
  mutate(ImmuneCell=str_remove_all(IMMUNECELL,"_.*")) %>%
  mutate(ImmuneCell=Hmisc::capitalize(ImmuneCell)) %>%
  mutate(Cancer=factor(Cancer,levels = rev(sort(unique(.$Cancer))))) %>%
  mutate(ImmuneCell=factor(ImmuneCell,levels = c("Microenvironment score",
                                                 "Stroma score",
                                                 "Immune score",
                                                 "Cytotoxicity score")))
Figure.Body <- ggplot()+
  geom_point(data=IMMUNECELL.SYMBOL.SPEARMAN2,
             aes(x=ImmuneCell,y=Cancer,
                 size=-log10(Pvalue),
                 fill=Rho),
             shape=21,color="white")+
  geom_point(data=IMMUNECELL.SYMBOL.SPEARMAN2 %>% 
               filter(Pvalue<=0.05),
             aes(x=ImmuneCell,y=Cancer,
                 size=-log10(Pvalue),
                 fill=Rho),
             shape=21,color="black",stroke=1.2)+
  #facet_wrap(~SCORE)+
  ggthemes::theme_few()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9,angle = 90,vjust = 0.3,hjust = 1),
    axis.text.y = element_text(size=9),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="",y="",fill="Spearman'CC")+
  #geom_vline(xintercept = c(0),linetype=2)+
  #theme(strip.text.x = element_text(size = 12))+
  #ggthemes::scale_fill_gradient2_tableau(palette = "Red-Blue Diverging")
  scale_fill_gradient2(low = "#0066CC", mid = "white", high = "#c72e29")

Figure.Top <- IMMUNECELL.SYMBOL.SPEARMAN2 %>%
  dplyr::select(ImmuneCell,Method) %>%
  unique() %>% 
  mutate(ImmuneCell=factor(ImmuneCell,levels = c("Microenvironment score",
                                                 "Stroma score",
                                                 "Immune score",
                                                 "Cytotoxicity score"))) %>%
  ggplot(aes(ImmuneCell,y=1,fill=Method))+
  geom_col(width = 1,color="black")+
  theme_void()+
  scale_fill_manual(values = brewer.pal(8,"Set1")[2:1])

Figure.F <- Figure.Body %>% aplot::insert_top(Figure.Top,height = 0.03)
ggsave(Figure.F,filename = "BCAA.ssGSEA.ImmuneScore.Spearman.pdf",height = 6,width = 3.5)
#### Figure 2 ####
IMMUNECELL.SYMBOL.SPEARMAN <- read.csv("BCAA.ssGSEA.IMMUNECELL.SPEARMAN.csv")
IMMUNECELL.SYMBOL.SPEARMAN3 <- IMMUNECELL.SYMBOL.SPEARMAN %>% 
  filter(!str_detect(IMMUNECELL,"score")) %>%
  filter(!str_detect(IMMUNECELL,"uncharacter")) %>%
  mutate(Method=str_remove_all(IMMUNECELL,".*_")) %>%
  filter(! str_detect(IMMUNECELL,"EPIC")) %>%
  filter(Method != "CIBERSORT-ABS") %>%
  mutate(ImmuneCell=str_remove_all(IMMUNECELL,"_.*")) %>%
  mutate(ImmuneCell=Hmisc::capitalize(ImmuneCell)) %>%
  mutate(Label=if_else(Pvalue > 0.05,"",if_else(Pvalue > 0.01,"*",if_else(Pvalue > 0.001,"**","***")))) %>%
  mutate(Cancer=factor(Cancer,levels = sort((unique(.$Cancer))))) %>%
  arrange(desc(ImmuneCell),Method) %>%
  mutate(IMMUNECELL=factor(IMMUNECELL,levels = unique(.$IMMUNECELL)))
  #mutate(Cancer=factor(Cancer,levels = rev(sort(unique(.$Cancer)))))

Figure.Body <- ggplot() +
  geom_tile(data=IMMUNECELL.SYMBOL.SPEARMAN3,
            aes(x=Cancer,y=IMMUNECELL,fill=Rho),
            color="white")+
  geom_point(data=IMMUNECELL.SYMBOL.SPEARMAN3 %>% 
               filter(Pvalue <= 0.05),
             aes(x=Cancer,y=IMMUNECELL,size=-log10(Pvalue)),
             shape=1,color="black")+
  geom_point(data=IMMUNECELL.SYMBOL.SPEARMAN3 %>% 
               filter(Pvalue > 0.05),
             aes(x=Cancer,y=IMMUNECELL,size=-log10(Pvalue)),
             shape=4,color="white")+
  scale_size_continuous(range = c(1,3))+
  ggthemes::theme_few()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9,angle = 90,vjust = 0.3,hjust = 1),
    axis.text.y = element_text(size=9),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="",y="",fill="Spearman'CC")+
  scale_fill_gradient2(low = "#0066CC", mid = "white", high = "#c72e29")

#ggsave(Figure.Body,filename = "BCAA.ssGSEA.ImmuneCell.Spearman2.pdf",width = 5.5,height = 6)

Figure.R1 <- IMMUNECELL.SYMBOL.SPEARMAN3 %>%
  dplyr::select(IMMUNECELL,Method) %>%
  unique() %>%
  mutate(IMMUNECELL=factor(IMMUNECELL,levels = unique(.$IMMUNECELL))) %>%
  ggplot(aes(1,y=IMMUNECELL,fill=Method))+
  geom_col(width = 1,color="black")+
  theme_void()+
  #scale_fill_manual(values = brewer.pal(8,"Paired")[1:6])
  ggsci::scale_fill_npg()

Figure.R2 <- IMMUNECELL.SYMBOL.SPEARMAN3 %>%
  mutate(Correlation=if_else(Pvalue > 0.05,"p>0.05",
                             if_else(Rho > 0,"Positive","Negative"))) %>%
  group_by(IMMUNECELL,Correlation) %>%
  dplyr::summarise(Count=n()) %>% #arrange(desc(Count))
  ungroup() %>%
  #mutate(Yaxis=paste(miRNA,Target,sep = "->")) %>%
  #dplyr::select(-Target,-miRNA) %>%
  spread(Correlation,Count,fill=0) %>%
  gather(Correlation,Count,-IMMUNECELL) %>%
  mutate(Correlation=factor(Correlation,
                            levels = c("p>0.05","Negative","Positive"))) %>%
  ggplot(aes(x=Count,y=IMMUNECELL,fill=Correlation))+
  geom_col(width = 1,color="black")+
  ggthemes::theme_few()+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9,angle = 90,vjust = 0.5),
    axis.text.y = element_blank(),
    axis.ticks.y=element_blank(),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="No. of cancers",y="")+
  scale_fill_manual(values = c("grey90","#66C2A5","#BC3C2999"))

Figure.F <- Figure.Body %>% aplot::insert_right(Figure.R1,width = 0.05) %>%
  aplot::insert_right(Figure.R2,width = 0.3)
ggsave(Figure.F,filename = "BCAA.ssGSEA.ImmuneCell.Spearman.pdf",width = 10,height = 12)

######## BCAA score HALLMARKS ####
setwd("K:/TCGA/Anlysis/BCAA")
library(tidyverse)
BCAA.ssGSEA <- read.csv("BCAA.ssGSEA/BCAA.ssGSEA.csv")
BCAA.ID <- read.csv("BCAA.ID.csv")

TERM <- c("HALLMARK_E2F_TARGETS",
          "HALLMARK_G2M_CHECKPOINT",
          "HALLMARK_MYC_TARGETS_V2",
          "HALLMARK_MYC_TARGETS_V1",
          "HALLMARK_MITOTIC_SPINDLE",
          "HALLMARK_NOTCH_SIGNALING",
          "HALLMARK_WNT_BETA_CATENIN_SIGNALING",
          "HALLMARK_MTORC1_SIGNALING",
          "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
          "HALLMARK_TGF_BETA_SIGNALING",
          "HALLMARK_KRAS_SIGNALING_UP",
          "HALLMARK_KRAS_SIGNALING_DN",
          "HALLMARK_IL6_JAK_STAT3_SIGNALING",
          "HALLMARK_INTERFERON_GAMMA_RESPONSE",
          "HALLMARK_INTERFERON_ALPHA_RESPONSE",
          "HALLMARK_INFLAMMATORY_RESPONSE",
          "HALLMARK_HEDGEHOG_SIGNALING",
          "HALLMARK_COMPLEMENT",
          "HALLMARK_DNA_REPAIR",
          "HALLMARK_HYPOXIA",
          "HALLMARK_APOPTOSIS",
          "HALLMARK_P53_PATHWAY",
          "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
          "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")

Spearman.HALLMARK <- data.frame()
for (project in unique(BCAA.ssGSEA$Cancer)) {
  mid.BCAA <- BCAA.ssGSEA %>% filter(Cancer==project)
  filename <- paste("TCGA-",project,".ssGSEA.csv",sep = "")
  ssGSEA.Exp <- read.csv(file.path("../../HALLMARK.ssGSEA.Pancancer/",filename),check.names = F,row.names = 1)
  Spearman.Data <- ssGSEA.Exp %>% rownames_to_column("TERM") %>%
    filter(TERM %in% TERM) %>%
    remove_rownames() %>%
    column_to_rownames("TERM") %>%
    t() %>% data.frame(check.names = F) %>%
    rownames_to_column("Sample.ID") %>%
    merge(mid.BCAA,by="Sample.ID")
  #Spearman.HALLMARK <- rbind.data.frame(Spearman.Data,Spearman.HALLMARK)
  for (HALLMARK in TERM) {
    for (GENE in c("BCAA")) {
      test <- cor.test(Spearman.Data[[HALLMARK]],
                       Spearman.Data[[GENE]],
                       method = "spearman",
                       alternative = "two.sided",
                       exact = F)
      Spearman.HALLMARK <- data.frame(HALLMARK = HALLMARK,
                                      SCORE = GENE,
                                      Rho=test$estimate,
                                      Pvalue=test$p.value,
                                      Cancer=project) %>%
        rbind.data.frame(Spearman.HALLMARK)
    }
  }
  
}
write.csv(Spearman.HALLMARK,"BCAA.ssGSEA.HALLMARK.Spearman.csv",row.names = F)
#### Figure ####
Spearman.HALLMARK <- read.csv("BCAA.ssGSEA.HALLMARK.Spearman.csv")
Spearman.HALLMARK$HALLMARK <- str_remove_all(Spearman.HALLMARK$HALLMARK,"HALLMARK_")
library(tidyverse)
library(RColorBrewer)
Cancer.Col <- c(brewer.pal(9,"Set1"),
                brewer.pal(8,"Set2")[1:7],
                brewer.pal(12,"Set3")[1:11],
                brewer.pal(8,"Dark2")[1:6])
names(Cancer.Col) <- sort(unique(Spearman.HALLMARK$Cancer))
p1 <- ggplot(Spearman.HALLMARK,
             aes(x=Rho,y=reorder(HALLMARK,Rho),
                 size=-log10(Pvalue),
                 color=Cancer))+
  geom_point()+
  #facet_wrap(~SCORE)+
  ggthemes::theme_few()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9),
    axis.text.y = element_text(size=9),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="Spearman'CC",y="")+
  geom_vline(xintercept = c(0),linetype=2)+
  theme(strip.text.x = element_text(size = 12))+
  scale_color_manual(values = Cancer.Col)

Spearman.HALLMARK2 <- Spearman.HALLMARK %>%
  #filter(Pvalue<=0.05) %>%
  mutate(Interaction=if_else(Rho>0,"Positive","Negative")) %>%
  group_by(HALLMARK,Interaction,SCORE) %>%
  summarise(Sig.Count=n()) %>%
  ungroup() %>%
  spread(Interaction,Sig.Count,fill=0) %>%
  magrittr::set_colnames(c("HALLMARK","SCORE","S.N","S.P"))

Spearman.HALLMARK3 <- Spearman.HALLMARK %>%
  filter(Pvalue<=0.05) %>%
  mutate(Interaction=if_else(Rho>0,"Positive","Negative")) %>%
  group_by(HALLMARK,Interaction,SCORE) %>%
  summarise(Sig.Count=n()) %>%
  ungroup() %>%
  spread(Interaction,Sig.Count,fill=0) %>%
  magrittr::set_colnames(c("HALLMARK","SCORE","Sig.N","Sig.P")) %>%
  merge(Spearman.HALLMARK2,by=c("HALLMARK","SCORE")) %>%
  mutate(E.P=S.P-Sig.P) %>%
  mutate(E.N=S.N-Sig.N) %>%
  arrange(S.P) %>%
  mutate(Sig.N=-Sig.N,E.N=-E.N) %>%
  #mutate(HALLMARK=factor(HALLMARK,levels = .$HALLMARK)) %>%
  dplyr::select(-S.P,-S.N) %>%
  gather(Attribute,Count,Sig.N:E.N)

p2 <- ggplot(Spearman.HALLMARK3,aes(x=Count,y=reorder(HALLMARK,Count),fill=Attribute))+
  geom_col()+
  ggthemes::theme_few()+
  facet_wrap(~SCORE)+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(y="",x="No. of cancer types",fill="Type")+
  #scale_fill_discrete()+
  scale_fill_manual(breaks=c("E.N","E.P","Sig.N","Sig.P"),
                    labels=c("Insignificant negative correlation",
                             "Insignificant positive correlation",
                             "Significant negative correlation",
                             "Significant positive correlation"),
                    values =  c("#4DAF4A","#FF7F00","#377EB8","#E41A1C"))

Figure.F <- p1 %>% aplot::insert_right(p2,width = 1)

ggsave(Figure.F,filename = "BCAA.ssGSEA.HALLMARK.Spearman.pdf",height = 5,width = 8)
ggsave(Figure.F,filename = "BCAA.ssGSEA.HALLMARK.Spearman.Legend.pdf",height = 10,width = 8)

######## BCAA score PATHWAY Spearman ############
setwd("K:/TCGA/Anlysis/BCAA")
library(tidyverse)
PhenoData <- read.csv("../../mRNA.Phenotype.csv") %>%
  filter(Sample.Type %in% c("Primary Tumor",
                            "Primary Blood Derived Cancer - Peripheral Blood",
                            "Solid Tissue Normal")) %>%
  mutate(SampleType = ifelse(Sample.Type == "Solid Tissue Normal","Normal","Tumor")) %>%
  mutate(ProjectID=str_remove_all(Project.ID,"TCGA-")) %>%
  dplyr::select(ProjectID,SampleType,Sample.ID) %>%
  filter(SampleType=="Tumor")

BCAA.ssGSEA <- read.csv("BCAA.ssGSEA/BCAA.ssGSEA.csv") %>%
  mutate(barcode = substring(Sample.ID,1,12)) %>%
  filter(Sample.ID %in% PhenoData$Sample.ID)
BCAA.ID <- read.csv("BCAA.ID.csv")

Cancer.Score.PAS.Spearman <- data.frame()
files <- list.files(path = "../../GSCA/",pattern = "tsv.gz")
files <- files[str_detect(files,"PAS")]
for (file in files) {
  project = file %>% str_remove_all("\\..*")
  library(data.table)
  data <- fread(file.path("../../GSCA/",file),sep = "\t",check.names = F)
  Exp.PAS <- merge(data,BCAA.ssGSEA,by="barcode")
  for (index in unique(Exp.PAS$pathway)) {
    mid.data <- Exp.PAS %>% filter(pathway==index)
    for (gene in c("BCAA")) {
      test <- cor.test(mid.data[[gene]],mid.data[["score"]],
                       method = "spearman",
                       alternative = "two.sided",
                       exact = F)
      Cancer.Score.PAS.Spearman <- data.frame(Cancer=project,
                                              SYMBOL=gene,
                                              PATHWAY=index,
                                              Rho=test$estimate,
                                              Pvalue=test$p.value,
                                              Count=nrow(mid.data)) %>%
        rbind.data.frame(Cancer.Score.PAS.Spearman)
    }
  }
}
write.csv(Cancer.Score.PAS.Spearman,file = "BCAA.ssGSEA.PAS.Spearman.csv",row.names = F)

#### Figure ####
setwd("K:/TCGA/Anlysis/BCAA")
library(tidyverse)
Cancer.Score.PAS.Spearman <- read.csv(file = "BCAA.ssGSEA.PAS.Spearman.csv")
Cancer.Score.PAS.Spearman2 <- Cancer.Score.PAS.Spearman %>%
  mutate(PATHWAY=factor(PATHWAY,levels = rev(sort(unique(.$PATHWAY))))) %>%
  arrange(Cancer,SYMBOL) %>%
  rstatix::add_significance("Pvalue") %>%
  mutate(label=if_else(Pvalue.signif == "ns","",Pvalue.signif)) %>%
  mutate(Xaxis=paste(Cancer,"(",Count,")",sep = ""))

p21 <- ggplot()+
  geom_tile(data=Cancer.Score.PAS.Spearman2,
            aes(x=Xaxis,y=PATHWAY,fill=Rho),
            color="white")+
  geom_point(data=Cancer.Score.PAS.Spearman2 %>% 
               filter(Pvalue <= 0.05),
             aes(x=Xaxis,y=PATHWAY,size=-log10(Pvalue)),
             shape=1,color="black")+
  geom_point(data=Cancer.Score.PAS.Spearman2 %>% 
               filter(Pvalue > 0.05),
             aes(x=Xaxis,y=PATHWAY,size=-log10(Pvalue)),
             shape=4,color="white")+
  ggthemes::theme_few()+
  theme(legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9,angle = 90,vjust = 0.5,hjust = 1),
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="",y="")+
  scale_fill_gradient2("Spearman'CC",
                       low = "#053061",high="#67001F",mid = "white")
ggsave(p21,filename = "BCAA.ssGSEA.PAS.Spearman.pdf",height = 4,width = 8)

######## BCAA score Subtype ####
setwd("K:/TCGA/Anlysis/BCAA")
library(tidyverse)
BCAA.ssGSEA <- read.csv("BCAA.ssGSEA/BCAA.ssGSEA.csv")
BCAA.ID <- read.csv("BCAA.ID.csv")

PhenoData <- read.csv("../../mRNA.Phenotype.csv") %>%
  filter(Sample.Type %in% c("Primary Tumor",
                            "Primary Blood Derived Cancer - Peripheral Blood",
                            "Solid Tissue Normal")) %>%
  mutate(SampleType = ifelse(Sample.Type == "Solid Tissue Normal","Normal","Tumor")) %>%
  mutate(ProjectID=str_remove_all(Project.ID,"TCGA-")) %>%
  dplyr::select(ProjectID,SampleType,Sample.ID) %>%
  filter(SampleType=="Tumor")

BCAA.SCORE <- merge(BCAA.ssGSEA,PhenoData,by="Sample.ID") %>%
  data.frame(check.names = F) %>%
  mutate(barcode=substr(Sample.ID,1,12))

files <- list.files(path = "../../GSCA/",pattern = "gz$")
files <- files[str_detect(files,"subtype")]

Subtype.Data <- data.frame()
library(data.table)
for (file in files) {
  project <- file %>% str_remove_all("\\..*")
  Subtype.Data <- fread(file.path("../../GSCA/",file)) %>%
    merge(BCAA.SCORE,by="barcode") %>%
    mutate(Cancer=project) %>%
    rbind.data.frame(Subtype.Data)
}
write.csv(Subtype.Data,file = "BCAA.ssGSEA.SubType.Data.csv",row.names = F)
data2 <- compare_means(formula = BCAA~subtype,data = Subtype.Data,method = "kruskal.test",
                       group.by="Cancer") %>%
  data.frame(check.names = F)
write.csv(data2,file = "BCAA.ssGSEA.SubType.Data.Kruskal.csv",row.names = F)
#### Figure ####
Subtype.Data <- read.csv("BCAA.ssGSEA.SubType.Data.csv")
data2 <- read.csv("BCAA.ssGSEA.SubType.Data.Kruskal.csv")
Subtype.Data2 <- Subtype.Data %>%
  arrange(desc(Cancer),subtype) %>%
  mutate(Yaxis=paste(Cancer,subtype,sep = ":")) %>%
  mutate(Yaxis=factor(Yaxis,levels = unique(Yaxis))) %>%
  group_by(Yaxis) %>%
  summarise(Mean=mean(BCAA)) %>%
  mutate(Cancer=str_remove_all(Yaxis,":.*"))

Subtype.Data3 <- Subtype.Data %>%
  arrange(desc(Cancer),subtype) %>%
  mutate(Yaxis=paste(Cancer,subtype,sep = ":")) %>%
  mutate(Yaxis=factor(Yaxis,levels = unique(Yaxis)))

p10 <- ggplot(Subtype.Data3,aes(x=BCAA,y=Yaxis,fill=Cancer))+
  stat_summary(fun.data = "mean_cl_boot", # # 误差线
               geom = "errorbar",
               width = .4) +
  stat_summary(fun = "mean", geom = "point",shape=21)+
  labs(x="",y="")+
  #facet_wrap(~Phenotype,scale="free_y",ncol = 6)+
  ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9),
        axis.text.y = element_text(size=12),
        #axis.ticks.x = element_blank(),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  #scale_color_manual("",values = brewer.pal(9,"Set1")[c(2:1)])+
  #scale_fill_manual("",values = brewer.pal(9,"Set1")[c(2:1)])+
  #scale_x_continuous(trans = "log10")+
  #ggsci::scale_fill_npg()
  ggthemes::scale_fill_stata()+
  scale_x_continuous(limits = c(0,10),breaks = c(0,3,6,9))
p10.F <- p10 %>% aplot::insert_right(p10.R,width = 0.15)
ggsave(p10.F,filename = "BCAA.ssGSEA.Subtype.Tumor2.pdf",height = 5,width = 4)

p10<-
  ggplot(data=Subtype.Data2,
         aes(x=Mean,y=Yaxis,fill=Cancer))+
  geom_segment(aes(x=0,xend=Mean,y=Yaxis,yend=Yaxis),color="black")+
  geom_point(size=3,shape=21,color="black")+
  #geom_boxplot()+
  #geom_text(data=data2,aes(x=10.5,y=Cancer,label=p.signif))+
  #stat_ellipse(aes(fill=Group),type = "norm", geom ="polygon",alpha=0.2,color=NA)+
  #geom_beeswarm(aes(color=SampleType),size=1,dodge.width=0.75,priority = "ascending")+
  labs(x="",y="")+
  #facet_wrap(~Phenotype,scale="free_y",ncol = 6)+
  ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9),
        axis.text.y = element_text(size=12),
        #axis.ticks.x = element_blank(),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  #scale_color_manual("",values = brewer.pal(9,"Set1")[c(2:1)])+
  #scale_fill_manual("",values = brewer.pal(9,"Set1")[c(2:1)])+
  #scale_x_continuous(trans = "log10")+
  #ggsci::scale_fill_npg()
  ggthemes::scale_fill_stata()+
  scale_x_continuous(limits = c(0,10),breaks = c(0,3,6,9))

p10.R <- Subtype.Data2 %>%
  dplyr::select(Yaxis,Cancer) %>%
  unique() %>%
  mutate(Yaxis=factor(Yaxis,levels = unique(Yaxis))) %>%
  ggplot(aes(x=1,y=Yaxis,fill=Cancer))+
  geom_col(width=1)+
  theme_void()+
  #ggsci::scale_fill_nejm()
  ggthemes::scale_fill_stata()

p10.F <- p10 %>% aplot::insert_right(p10.R,width = 0.15)
ggsave(p10.F,filename = "BCAA.ssGSEA.Subtype.Tumor.pdf",height = 5,width = 4)

######## BCAA score Grade ####
setwd("K:/TCGA/Anlysis/BCAA")
library(tidyverse)
BCAA.ssGSEA <- read.csv("BCAA.ssGSEA/BCAA.ssGSEA.csv") %>%
  mutate(SampleName=substring(Sample.ID,1,15))
BCAA.ID <- read.csv("BCAA.ID.csv")

library(readxl)
data <- read_xlsx("../../PanCancer.DiverseIndex.xlsx",sheet = 18)
data2<-data[,-1]
BCAA.GRADE <- merge(BCAA.ssGSEA,data2,by="SampleName")

data2 <- compare_means(formula = BCAA~Grade,data = BCAA.GRADE,method = "kruskal.test",
                       group.by="Cancer") %>%
  data.frame(check.names = F)

write.csv(BCAA.GRADE,file = "BCAA.ssGSEA.Grade.Data.csv",row.names = F)
write.csv(data2,file = "BCAA.ssGSEA.Grade.Data.Kruskal.csv",row.names = F)

BCAA.GRADE <- read.csv("BCAA.ssGSEA.Grade.Data.csv")
data2 <- read.csv("BCAA.ssGSEA.Grade.Data.Kruskal.csv")
BCAA.GRADE2 <- BCAA.GRADE %>%
  mutate(Xaxis=paste(Cancer,Grade,sep = "_"))

p11<-
  ggplot(data=BCAA.GRADE2,
         aes(x=Xaxis,y=BCAA, #Cancer,#
             fill=Cancer))+
  #geom_boxplot()+
  #geom_text(data=data2,aes(x=10.5,y=Cancer,label=p.signif))+
  #stat_ellipse(aes(fill=Group),type = "norm", geom ="polygon",alpha=0.2,color=NA)+
  #geom_beeswarm(aes(color=SampleType),size=1,dodge.width=0.75,priority = "ascending")+
  stat_summary(fun.data = "mean_cl_boot", # # 误差线
             geom = "errorbar",
             width = .4) +
  stat_summary(fun = "mean", geom = "point",shape=21)+
  labs(x="",y="BCAA score")+
  #facet_wrap(~Phenotype,scale="free_y",ncol = 6)+
  ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9,angle = 90,vjust=0.5),
        axis.text.y = element_text(size=9),
        #axis.ticks.x = element_blank(),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  #scale_color_manual("",values = brewer.pal(9,"Set1")[c(2:1)])+
  #scale_fill_manual("",values = brewer.pal(9,"Set1")[c(2:1)])+
  #scale_x_continuous(trans = "log10")+
  #ggsci::scale_fill_npg()
  ggthemes::scale_fill_stata()
p11.T <- BCAA.GRADE2 %>%
  dplyr::select(Xaxis,Cancer) %>%
  unique() %>%
  #mutate(Xaxis=factor(Xaxis,levels = unique(Xaxis))) %>%
  ggplot(aes(x=Xaxis,y=1,fill=Cancer))+
  geom_col(width=1)+
  theme_void()+
  #ggsci::scale_fill_nejm()
  ggthemes::scale_fill_stata()+
  theme(legend.position = "bottom")

p11.F <- p11 %>% aplot::insert_top(p11.T,height = 0.15)
ggsave(p11.F,filename = "BCAA.ssGSEA.Grade.Tumor2.pdf",height = 2.5,width = 8)
ggsave(p11.F,filename = "BCAA.ssGSEA.Grade.Tumor2.Ledend.pdf",height = 12,width = 8)
######## BCAA score Stage ####
setwd("K:/TCGA/Anlysis/BCAA")
library(tidyverse)
BCAA.ssGSEA <- read.csv("BCAA.ssGSEA/BCAA.ssGSEA.csv") %>%
  mutate(SampleName=substring(Sample.ID,1,15))
BCAA.ID <- read.csv("BCAA.ID.csv")

library(readxl)
data <- read_xlsx("../../PanCancer.DiverseIndex.xlsx",sheet = 19)
data2<-data[,-1]
BCAA.STAGE <- merge(BCAA.ssGSEA,data2,by="SampleName")

data2 <- compare_means(formula = BCAA~Stage,
                       data = BCAA.STAGE,
                       method = "kruskal.test",
                       group.by="Cancer") %>%
  data.frame(check.names = F)

write.csv(BCAA.STAGE,file = "BCAA.ssGSEA.BCAA.STAGE.Data.csv",row.names = F)
write.csv(data2,file = "BCAA.ssGSEA.BCAA.STAGE.Data.Kruskal.csv",row.names = F)

BCAA.STAGE <- read.csv("BCAA.ssGSEA.BCAA.STAGE.Data.csv")
data2 <- read.csv("BCAA.ssGSEA.BCAA.STAGE.Data.Kruskal.csv")
BCAA.STAGE2 <- BCAA.STAGE %>%
  mutate(STAGE=str_remove_all(Stage,"Stage ")) %>%
  mutate(Xaxis=paste(Cancer,STAGE,sep = "_"))

p11<-
  ggplot(data=BCAA.STAGE2,
         aes(x=Xaxis,y=BCAA, #Cancer,#
             fill=Cancer))+
  #geom_boxplot()+
  #geom_text(data=data2,aes(x=10.5,y=Cancer,label=p.signif))+
  #stat_ellipse(aes(fill=Group),type = "norm", geom ="polygon",alpha=0.2,color=NA)+
  #geom_beeswarm(aes(color=SampleType),size=1,dodge.width=0.75,priority = "ascending")+
  stat_summary(fun.data = "mean_cl_boot", # # 误差线
               geom = "errorbar",
               width = .4) +
  stat_summary(fun = "mean", geom = "point",shape=21)+
  labs(x="",y="BCAA score")+
  #facet_wrap(~Phenotype,scale="free_y",ncol = 6)+
  ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_blank(),#element_text(size=9,angle = 90,vjust=0.5),
        axis.text.y = element_text(size=9),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  #scale_color_manual("",values = brewer.pal(9,"Set1")[c(2:1)])+
  #scale_fill_manual("",values = brewer.pal(9,"Set1")[c(2:1)])+
  #scale_x_continuous(trans = "log10")+
  #ggsci::scale_fill_npg()
  #ggthemes::scale_fill_stata()
  scale_fill_manual(values = Cancer.Col)

p11.T1 <- BCAA.STAGE2 %>%
  dplyr::select(Xaxis,Cancer) %>%
  unique() %>%
  #mutate(Xaxis=factor(Xaxis,levels = unique(Xaxis))) %>%
  ggplot(aes(x=Xaxis,y=1,fill=Cancer))+
  geom_col(width=1)+
  theme_void()+
  #ggsci::scale_fill_nejm()
  #ggthemes::scale_fill_stata()+
  scale_fill_manual(values = Cancer.Col)+
  theme(legend.position = "bottom")

p11.T2 <- BCAA.STAGE2 %>%
  dplyr::select(Xaxis,STAGE) %>%
  unique() %>%
  #mutate(Xaxis=factor(Xaxis,levels = unique(Xaxis))) %>%
  ggplot(aes(x=Xaxis,y=1,fill=STAGE))+
  geom_col(width=1)+
  theme_void()+
  ggsci::scale_fill_npg()+
  #ggthemes::scale_fill_stata()+
  #scale_fill_manual(values = Cancer.Col)+
  theme(legend.position = "bottom")

p11.F <- p11 %>% 
  aplot::insert_top(p11.T2,height = 0.1) %>%
  aplot::insert_top(p11.T1,height = 0.1)
ggsave(p11.F,filename = "BCAA.ssGSEA.STAGE.Tumor2.pdf",height = 2,width = 12)
ggsave(p11.F,filename = "BCAA.ssGSEA.STAGE.Tumor2.Ledend.pdf",height = 12,width = 12)

######## DEGS => GSEA ########
#### DEGs ####
setwd("K:/TCGA/Anlysis/BCAA")
library(tidyverse)
BCAA.ssGSEA <- read.csv("BCAA.ssGSEA/BCAA.ssGSEA.csv") %>%
  mutate(SampleName=substring(Sample.ID,1,15))
#BCAA.ssGSEA <- read.csv("BCAA.ssGSEA.csv") %>%
#  mutate(SampleName=substring(Sample.ID,1,15))
BCAA.ID <- read.csv("BCAA.ID.csv")
library(GSVA)
library(GSEABase)
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(limma)
library(tidyverse)
library(TCGAbiolinks)
library(DESeq2)
PhenoType <- read.csv("../../mRNA.Phenotype.csv")
PhenoType <- PhenoType %>% filter(Sample.Type == "Primary Tumor" | 
                                    Sample.Type == "Solid Tissue Normal" | 
                                    Sample.Type == "Primary Blood Derived Cancer - Peripheral Blood") %>% 
  dplyr::select(Project.ID,Sample.ID,Sample.Type) %>%
  mutate(Type=ifelse(Sample.Type=="Solid Tissue Normal","Normal","Tumor")) %>%
  mutate(Cancer=str_remove(Project.ID,"TCGA-"))


Hallmark <- read.gmt("../../msigdb/h.all.v7.5.1.symbols.gmt")
dir.create("DEGs.HALLMARK")
files <- list.files("../../mRNA.Pancancer.Exp.Count/",pattern = "csv$")
#files <- list.files("../mRNA.Pancancer.Exp.Count/",pattern = "csv$")
for (file in files) {
  project <- file %>%str_remove("PanCancer\\.") %>%
    str_remove("\\..*") %>% str_remove("TCGA-")
  print(project)
  #Count.data <- read.csv(file = file.path("../../mRNA.Pancancer.Exp.Count/",file),check.names = F,row.names = 1)
  Count.data <- read.csv(file = file.path("../mRNA.Pancancer.Exp.Count/",file),check.names = F,row.names = 1)
  filename <- paste("TCGA-",project,".mRNA.TPM.csv",sep = "")
  TPM.data <- read.csv(file = file.path("../mRNA.PanCancer.Exp/",filename),check.names = F,row.names = 1)
  #TPM.data <- read.csv(file = file.path("../../mRNA.PanCancer.Exp/",filename),check.names = F,row.names = 1)
  mid.Pheno <- PhenoType %>% filter(Type=="Tumor") %>%
    filter(Cancer==project)
  for (gene in BCAA.ID$ENSEMBL) {
    print(gene)
    SYMBOL=(BCAA.ID%>%filter(ENSEMBL==gene))$SYMBOL
    outfilename <- paste(project,SYMBOL,"DEGs.csv",sep = ".")
    if (!file.exists(file.path("DEGs.HALLMARK",outfilename))) {
      data2 <- TPM.data %>% dplyr::select(intersect(colnames(.),mid.Pheno$Sample.ID))
      GENE.EXP = data2[gene,] %>% t() %>% data.frame(check.names = F) %>%
        rownames_to_column("Sample.ID")
      GENE.EXP$Group=if_else(GENE.EXP[[gene]]>median(GENE.EXP[[gene]]),"high","low")
      DEG.data <- Count.data %>% rownames_to_column("ENSEMBL") %>%
        filter(ENSEMBL != gene) %>%
        column_to_rownames("ENSEMBL") %>%
        dplyr::select(GENE.EXP$Sample.ID)
      #High.Count <- DEG.data[,(GENE.EXP %>% filter(Group=="high"))$Sample.ID]
      #Low.Count <- DEG.data[,(GENE.EXP %>% filter(Group=="low"))$Sample.ID]
      countData <- DEG.data[rowMeans(DEG.data)>1,] 
      condition <- factor(GENE.EXP$Group,levels = c("high","low"))
      colData <- data.frame(row.names=colnames(countData), condition)
      
      dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)
      dds1 <- DESeq(dds)
      res <- results(dds1,contrast=c("condition","high","low"))
      res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
      
      write.csv(res1,file = file.path("DEGs.HALLMARK/",outfilename))
    }
  }
  
  mid.BCAA <- BCAA.ssGSEA %>% filter(Cancer==project) %>%
    filter(Sample.ID %in% mid.Pheno$Sample.ID)
  for (gene in c("BCAA")) {
    print(gene)
    outfilename <- paste(project,gene,"DEGs.csv",sep = ".")
    if (!file.exists(file.path("DEGs.HALLMARK",outfilename))) {
      mid.BCAA$Group=if_else(mid.BCAA[[gene]]>median(mid.BCAA[[gene]]),"high","low")
      DEG.data <- Count.data[,mid.BCAA$Sample.ID] %>% data.frame(check.names = F)
      countData <- DEG.data[rowMeans(DEG.data)>1,] 
      condition <- factor(mid.BCAA$Group,levels = c("high","low"))
      colData <- data.frame(row.names=colnames(countData), condition)
      
      dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)
      dds1 <- DESeq(dds)
      res <- results(dds1,contrast=c("condition","high","low"))
      res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
      #SYMBOL=(BCAA.ID%>%filter(ENSEMBL==gene))$SYMBOL
      write.csv(res1,file = file.path("DEGs.HALLMARK/",outfilename))
    }
  }
}

#### GSEA => HALLMARK ####
setwd("K:/TCGA/Anlysis/BCAA")
library(tidyverse)
library(GSVA)
library(GSEABase)
library(msigdbr)
library(fgsea)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(limma)
library(tidyverse)
library(TCGAbiolinks)
library(DESeq2)

Hallmark <- read.gmt("../../msigdb/h.all.v7.5.1.symbols.gmt")

GeneType <- read.table("../../../../Metastasis/TissueMicrobiome/Human.GRC38.GeneType.txt",header = T,sep = ",")
GeneType.Unique <- GeneType %>% 
  dplyr::select(Gene.stable.ID,Gene.name) %>% unique() %>%
  na.omit() %>%
  filter(Gene.name != "")
HALLMARK.GSEA <- data.frame()
files <- list.files("DEGs.HALLMARK/",pattern = "csv$")
for (file in files) {
  STRING <- file %>% str_split("\\.") %>% unlist()
  project = STRING[1]
  SYMBOL = STRING[2]
  print(project)
  data <- read.csv(file.path("DEGs.HALLMARK/",file),check.names = F,row.names = 1)
  data2 <- data %>% rownames_to_column("ENSEMBL") %>%
    merge(GeneType.Unique,by.x="ENSEMBL",by.y="Gene.stable.ID") %>%
    dplyr::select(-ENSEMBL)
  
  ge = data2$log2FoldChange
  names(ge) = data2$Gene.name
  ge = sort(ge,decreasing = T)
  
  GSEA.Res <- GSEA(ge, TERM2GENE = Hallmark,eps = 1e-100,pvalueCutoff =1)
  filename <- paste(project,SYMBOL,"GSEA.HALLMARK.Rdata",sep = ".")
  save(GSEA.Res,file = file.path("DEGs.HALLMARK/",filename))
  HALLMARK.GSEA <- GSEA.Res@result %>%
    mutate(Cancer=project,SYMBOL=SYMBOL) %>%
    rbind.data.frame(HALLMARK.GSEA)
}
#View(HALLMARK.GSEA)
write.csv(HALLMARK.GSEA,file = "DEGs.HALLMARK/DEGs.HALLMARK.GSEA.csv",row.names = F)
#### GSEA => HALLMARK => Figure ####
setwd("K:/TCGA/Anlysis/BCAA")
library(ggDoubleHeat)
library(ggplot2)
library(tidyverse)
library(readxl)
BCAA.ID <- read.csv("BCAA.ID.csv")
TERM <- c("HALLMARK_E2F_TARGETS",
          "HALLMARK_G2M_CHECKPOINT",
          "HALLMARK_MYC_TARGETS_V2",
          "HALLMARK_MYC_TARGETS_V1",
          "HALLMARK_MITOTIC_SPINDLE",
          "HALLMARK_NOTCH_SIGNALING",
          "HALLMARK_WNT_BETA_CATENIN_SIGNALING",
          "HALLMARK_MTORC1_SIGNALING",
          "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
          "HALLMARK_TGF_BETA_SIGNALING",
          "HALLMARK_KRAS_SIGNALING_UP",
          "HALLMARK_KRAS_SIGNALING_DN",
          "HALLMARK_IL6_JAK_STAT3_SIGNALING",
          "HALLMARK_INTERFERON_GAMMA_RESPONSE",
          "HALLMARK_INTERFERON_ALPHA_RESPONSE",
          "HALLMARK_INFLAMMATORY_RESPONSE",
          "HALLMARK_HEDGEHOG_SIGNALING",
          "HALLMARK_COMPLEMENT",
          "HALLMARK_DNA_REPAIR",
          "HALLMARK_HYPOXIA",
          "HALLMARK_APOPTOSIS",
          "HALLMARK_P53_PATHWAY",
          "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
          "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")
HALLMARK.GSEA <- read.csv("DEGs.HALLMARK/DEGs.HALLMARK.GSEA.csv")
HALLMARK.GSEA.Data <- HALLMARK.GSEA %>%
  filter(p.adjust <= 0.05) %>%
  mutate(Type=if_else(NES > 0,"High","Low")) %>%
  group_by(SYMBOL,ID,Type) %>%
  summarise(Count=n()) %>%
  ungroup() %>%
  spread(Type,Count,fill=0) %>%
  filter(ID %in% TERM) %>%
  mutate(ID=str_remove_all(ID,"HALLMARK_")) #%>%


BCAA.ID <- read.csv("BCAA.ID.csv")

HALLMARK.GSEA <- read.csv("DEGs.HALLMARK/DEGs.HALLMARK.GSEA.csv")
HALLMARK.GSEA.Data <- HALLMARK.GSEA %>%
  filter(p.adjust <= 0.05) %>%
  mutate(Type=if_else(NES > 0,"High","Low")) %>%
  group_by(SYMBOL,ID,Type) %>%
  summarise(Count=n()) %>%
  ungroup() %>%
  spread(Type,Count,fill=0) %>%
  #filter(ID %in% TERM) %>%
  #filter(str_detect(ID,"KEGG_") | str_detect(ID,"REACTOME_")) %>%
  #mutate(ID=str_remove_all(ID,"KEGG_"))
  filter(High >= 10 | Low >= 10)

High.Cancer <- HALLMARK.GSEA %>%
  filter(p.adjust <= 0.05) %>%
  filter(NES > 0) %>%
  dplyr::select(Cancer,SYMBOL,ID) %>%
  group_by(SYMBOL,ID) %>%
  mutate(High.Cancer=paste0(Cancer,collapse = ",")) %>% 
  dplyr::select(-Cancer) %>%
  unique()

Low.Cancer <- HALLMARK.GSEA %>%
  filter(p.adjust <= 0.05) %>%
  filter(NES < 0) %>%
  dplyr::select(Cancer,SYMBOL,ID) %>%
  group_by(SYMBOL,ID) %>%
  mutate(Low.Cancer=paste0(Cancer,collapse = ",")) %>% 
  dplyr::select(-Cancer) %>%
  unique()

HALLMARK.GSEA.Data2 <- HALLMARK.GSEA.Data %>% 
  merge(High.Cancer,by=c("ID","SYMBOL")) %>%
  merge(Low.Cancer,by=c("ID","SYMBOL")) %>%
  arrange(SYMBOL)
write.csv(HALLMARK.GSEA.Data2,"TableS6.csv",row.names = F)

#mutate(ID=factor(ID,levels = rev(TERM %>% str_remove_all("HALLMARK_"))))
#mutate(SYMBOL = factor(SYMBOL,levels = c("BCAA","GATOR1","GATOR2",BCAA.ID$SYMBOL)))
##是这个R包的小hug，如果不改顺序后面颜色会乱 , 将用于作图的数据集，按照用于绘图的X轴数据列的首字母从小到大进行序列重排
p12 <- ggplot(data = HALLMARK.GSEA.Data, 
              aes(x = SYMBOL, y = ID)) +
  geom_heat_tri(lower = High, upper = Low,
                lower_colors = RColorBrewer::brewer.pal(11,"RdBu")[6:1],
                upper_colors = RColorBrewer::brewer.pal(11,"RdBu")[6:11])+
  ggthemes::theme_few()+
  theme(legend.position = "top",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9,angle = 90,vjust = 0.5),
        axis.text.y = element_text(size=9),
        axis.title = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  geom_text(aes(label = Low), nudge_x = 0.22, nudge_y = 0.2,size= 3,color = c('#FF7F00'))+
  geom_text(aes(label = High), nudge_x = -0.22, nudge_y = -0.2,size= 3,color = c('#1B9E77'))
#scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11,"RdBu")[6:11])

ggsave(p12,filename = "DEGs.HALLMARK.GSEA.pdf",height = 6,width = 8)
ggsave(p12,filename = "DEGs.HALLMARK/DEGs.HALLMARK.GSEA.pdf",height = 6,width = 8)


#### GSEA => KEGG ####
setwd("K:/TCGA/Anlysis/BCAA")
library(tidyverse)
library(GSVA)
library(GSEABase)
library(msigdbr)
library(fgsea)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(limma)
library(tidyverse)
library(TCGAbiolinks)
library(DESeq2)

Hallmark <- read.gmt("../../msigdb/c2.all.v7.5.1.symbols.gmt")

GeneType <- read.table("../../../../Metastasis/TissueMicrobiome/Human.GRC38.GeneType.txt",header = T,sep = ",")
GeneType.Unique <- GeneType %>% 
  dplyr::select(Gene.stable.ID,Gene.name) %>% unique() %>%
  na.omit() %>%
  filter(Gene.name != "")
HALLMARK.GSEA <- data.frame()
files <- list.files("DEGs.HALLMARK/",pattern = "csv$")
files <- files[str_detect(files,".DEGs.csv")]
for (file in files) {
  STRING <- file %>% str_split("\\.") %>% unlist()
  project = STRING[1]
  SYMBOL = STRING[2]
  print(project)
  data <- read.csv(file.path("DEGs.HALLMARK/",file),check.names = F,row.names = 1)
  data2 <- data %>% rownames_to_column("ENSEMBL") %>%
    merge(GeneType.Unique,by.x="ENSEMBL",by.y="Gene.stable.ID") %>%
    dplyr::select(-ENSEMBL)
  
  ge = data2$log2FoldChange
  names(ge) = data2$Gene.name
  ge = sort(ge,decreasing = T)
  
  GSEA.Res <- GSEA(ge, TERM2GENE = Hallmark,eps = 1e-100,pvalueCutoff =1)
  filename <- paste(project,SYMBOL,"GSEA.KEGG.Rdata",sep = ".")
  save(GSEA.Res,file = file.path("DEGs.HALLMARK/",filename))
  HALLMARK.GSEA <- GSEA.Res@result %>%
    mutate(Cancer=project,SYMBOL=SYMBOL) %>%
    rbind.data.frame(HALLMARK.GSEA)
}
#View(HALLMARK.GSEA)
write.csv(HALLMARK.GSEA,file = "DEGs.HALLMARK/DEGs.KEGG.GSEA.csv",row.names = F)

#### GSEA => KEGG => Figure ####
setwd("K:/TCGA/Anlysis/BCAA")
library(ggDoubleHeat)
library(ggplot2)
library(tidyverse)
library(readxl)
BCAA.ID <- read.csv("BCAA.ID.csv")

HALLMARK.GSEA <- read.csv("DEGs.HALLMARK/DEGs.KEGG.GSEA.csv")
HALLMARK.GSEA.Data <- HALLMARK.GSEA %>%
  filter(p.adjust <= 0.05) %>%
  mutate(Type=if_else(NES > 0,"High","Low")) %>%
  group_by(SYMBOL,ID,Type) %>%
  summarise(Count=n()) %>%
  ungroup() %>%
  spread(Type,Count,fill=0) %>%
  #filter(ID %in% TERM) %>%
  filter(str_detect(ID,"KEGG_") | str_detect(ID,"REACTOME_")) %>%
  #mutate(ID=str_remove_all(ID,"KEGG_"))
  filter(High >= 10 | Low >= 10)

High.Cancer <- HALLMARK.GSEA %>%
  filter(p.adjust <= 0.05) %>%
  filter(NES > 0) %>%
  dplyr::select(Cancer,SYMBOL,ID) %>%
  group_by(SYMBOL,ID) %>%
  mutate(High.Cancer=paste0(Cancer,collapse = ",")) %>% 
  dplyr::select(-Cancer) %>%
  unique()
  
Low.Cancer <- HALLMARK.GSEA %>%
  filter(p.adjust <= 0.05) %>%
  filter(NES < 0) %>%
  dplyr::select(Cancer,SYMBOL,ID) %>%
  group_by(SYMBOL,ID) %>%
  mutate(Low.Cancer=paste0(Cancer,collapse = ",")) %>% 
  dplyr::select(-Cancer) %>%
  unique()
  
HALLMARK.GSEA.Data2 <- HALLMARK.GSEA.Data %>% 
  merge(High.Cancer,by=c("ID","SYMBOL")) %>%
  merge(Low.Cancer,by=c("ID","SYMBOL")) %>%
  arrange(SYMBOL)
write.csv(HALLMARK.GSEA.Data2,"TableS7.csv",row.names = F)

######## Risk score LASSO COX => glmnet => LASSO model ####
library(glmnet)
setwd("K:/TCGA/Anlysis/BCAA")
library(tidyverse)
library(survival)
library(survminer)
#BCAA.ssGSEA <- read.csv("BCAA.ssGSEA.csv") %>%
#  mutate(SampleName=substring(Sample.ID,1,15))
BCAA.ID <- read.csv("BCAA.ID.csv")

files <- list.files(path = "../../mRNA.PanCancer.Exp/",pattern = "csv$")
LASSO.Feature <- data.frame()
for (file in files) {
  project <- file %>% str_remove_all("\\..*") %>% str_remove_all("TCGA-")
  print(project)
  phenotype <- read.csv(paste("../../Clinical.XENA/",project,"/","TCGA-",project,".GDC_phenotype.tsv",sep = ""),sep = "\t",header = T) %>% 
    dplyr::filter(sample_type.samples == "Primary Tumor" | 
                    sample_type.samples == "Primary Blood Derived Cancer - Peripheral Blood")
  TPM <- read.csv(paste("../../mRNA.PanCancer.Exp/TCGA-",project,".mRNA.TPM.csv",sep = ""),check.names = F,row.names = 1)
  TPM.log <- log2(TPM+1)
  
  #### survival datasets
  survival.data <- read.csv(paste("../../Clinical.XENA/",project,"/TCGA-",project,".survival.tsv",sep = ""),sep = "\t",header = T,check.names = F)
  Inter.Samples <- intersect(intersect(phenotype$submitter_id.samples,colnames(TPM)),survival.data$sample)
  #### symbol expression
  TPM.ENSEMBL.Target <- TPM.log[BCAA.ID$ENSEMBL,Inter.Samples] %>% 
    t() %>% data.frame(check.names = F) %>% 
    mutate(SampleID = rownames(.)) %>% 
    merge(.,survival.data,by.x="SampleID",by.y="sample")
  
  x <- as.matrix(TPM.ENSEMBL.Target[,BCAA.ID$ENSEMBL]) #gene => column
  y <- Surv(as.double(TPM.ENSEMBL.Target$OS.time),as.double(TPM.ENSEMBL.Target$OS))
  
  Select.Variables <- c()
  for (i in 1:500) {
    print(i)
    #fit <- glmnet(x, y, family="cox")
    #plot(fit,xvar = "lambda")
    #using 10-fold CV to select lambda:
    cv.fit =cv.glmnet(x, y, family="cox", nfolds=5)
    #cv.fit$lambda.min
    fit_regL=cv.fit
    fit_regL_CVdev=cv.fit
    
    # final Lasso model:
    model_lasso_min <- glmnet(x=x, y=y,family="cox",lambda=cv.fit$lambda.min)
    fit_regL_coef <- coef(fit_regL, s=fit_regL_CVdev$lambda.min)
    ### Plot Cross-Validation LASSO model
    
    #par(mfrow=c(1,2))
    #plot(fit_regL_CVdev,las=1, main="Lasso fit CV")
    #abline(v=log(fit_regL_CVdev$lambda.min), col="orange", lwd=2)
    #text(log(fit_regL_CVdev$lambda.min), 1.4, paste("lambda.min=",round(fit_regL_CVdev$lambda.min,4),"\n",length(fit_regL_coef@i), " genes" ,sep=""), col="orange", cex=0.75, pos=4)
    
    ### Plot lambda fit
    #plot(fit_regL_CVdev, xvar="lambda", main="Lasso coefficient")
    #abline(v=log(fit_regL_CVdev$lambda.min), col="orange", lwd=2)
    ###############
    #summary(cv.fit)
    #print(cv.fit)
    #cv.fit$lambda.min
    #cv.fit$lambda.1se
    #coef(cv.fit, s = "lambda.min") #----------candidates
    #coef.min= coef(cv.fit, s = "lambda.min")
    select.varialbes = rownames(as.data.frame(which(coef(cv.fit, s = "lambda.min")[,1]!=0)))
    #select.varialbes = select.varialbes[-1] #remove intercept
    #select.varialbes
    Select.Variables<-c(Select.Variables,select.varialbes)
  }
  if (length(Select.Variables) != 0) {
    LASSO.Feature <- data.frame(Feature=Select.Variables,Cancer=project) %>%
      group_by(Cancer,Feature)%>%
      summarise(Count=n()) %>%
      rbind.data.frame(LASSO.Feature)
  }
}
#write.csv(LASSO.Feature,file = "Risk.LASSO.COX.csv",row.names = F)
######## Risk score Confirm Risk group Cancer ####
setwd("K:/TCGA/Anlysis/BCAA")
library(tidyverse)
library(survival)
library(survminer)
LASSO.Feature <- read.csv("Risk.LASSO.COX.csv")

Feature.Count <- LASSO.Feature %>%
  filter(Count >= 50) %>%
  group_by(Cancer) %>%
  summarise(Count=n()) %>%
  filter(Count >= 2)
Feature.Filter <- LASSO.Feature %>%
  filter(Count >= 50) %>%
  filter(Cancer %in% Feature.Count$Cancer)

CIndex <- data.frame()
RiskScore.Cancer <- data.frame()
files <- list.files("RiskGroup",pattern = "Risk.Model.rds")
Projects <- files %>% str_remove_all("\\..*")
dir.create("RiskGroup")
for (project in Projects) {
  print(project)
  phenotype <- read.csv(paste("../../Clinical.XENA/",project,"/","TCGA-",project,".GDC_phenotype.tsv",sep = ""),sep = "\t",header = T) %>% 
    dplyr::filter(sample_type.samples == "Primary Tumor" | 
                    sample_type.samples == "Primary Blood Derived Cancer - Peripheral Blood")
  TPM <- read.csv(paste("../../mRNA.PanCancer.Exp/TCGA-",project,".mRNA.TPM.csv",sep = ""),check.names = F,row.names = 1)
  TPM.log <- log2(TPM+1)
  
  #Confirm.Feature <- Feature.Filter %>%
  #  filter(Cancer==project)
  #### survival datasets
  survival.data <- read.csv(paste("../../Clinical.XENA/",project,"/TCGA-",project,".survival.tsv",sep = ""),sep = "\t",header = T,check.names = F)
  Inter.Samples <- intersect(intersect(phenotype$submitter_id.samples,colnames(TPM)),survival.data$sample)
  #### symbol expression
  #TPM.ENSEMBL.Target <- TPM.log[Confirm.Feature$Feature,Inter.Samples] %>% 
  #  t() %>% data.frame(check.names = F) %>% 
  #  mutate(SampleID = rownames(.)) %>% 
  #  merge(.,survival.data,by.x="SampleID",by.y="sample")
  
  ## formaula
  #formula <- as.formula(paste0("Surv(OS.time,OS)~",paste(Confirm.Feature$Feature,collapse = "+")))
  #multi_variate_cox <- coxph(formula,data=TPM.ENSEMBL.Target)
  #saveRDS(multi_variate_cox,file = file.path("RiskGroup",paste(project,".Risk.Model.rds",sep = "")))
  ## check PH hypothesis
  #ph_hypo_multi <- cox.zph(multi_variate_cox)
  #ph_hypo_table <- ph_hypo_multi$table[-nrow(ph_hypo_multi$table),]
  multi_variate_cox <- readRDS(file = file.path("RiskGroup",paste(project,".Risk.Model.rds",sep = "")))
  
  TPM.ENSEMBL.Target <- TPM.log[names(multi_variate_cox$coefficients),Inter.Samples] %>% 
    t() %>% data.frame(check.names = F) %>% 
    mutate(SampleID = rownames(.)) %>% 
    merge(.,survival.data,by.x="SampleID",by.y="sample")
  
  #riskscore <- function(survival_cancer_df,candidate_genes_for_cox,cox_model){
  #  library(dplyr)
  #  risk_score_table <- survival_cancer_df %>% 
  #    dplyr::select(candidate_genes_for_cox)
  #  for (each_sig_gene in colnames(risk_score_table)) {
  #    risk_score_table$each_sig_gene <- risk_score_table[,each_sig_gene]*(summary(cox_model)$coefficients[each_sig_gene,1])
  #  }
  #  risk_score_table <- cbind(risk_score_table, 'total_risk_score'=exp(rowSums(risk_score_table))) %>%
  #    cbind(survival_cancer_df[,c("SampleID",'OS.time','OS')])
  #  risk_score_table <- risk_score_table[,c("SampleID",'OS.time','OS',candidate_genes_for_cox,'total_risk_score')]
  #  return(risk_score_table)
  #}
  RiskScore.Cancer <- predict(multi_variate_cox,
                              type="risk",
                              newdata=TPM.ENSEMBL.Target[,names(multi_variate_cox$coefficients)]) %>%
    as.data.frame() %>% 
    magrittr::set_colnames("total_risk_score") %>%
    #mutate(TPM.ENSEMBL.Target$SampleID) %>%
    cbind.data.frame(TPM.ENSEMBL.Target) %>%
    data.frame(check.names = F) %>%
    dplyr::select(c("SampleID",'OS.time','OS',"total_risk_score")) %>%
    mutate(Cancer=project)%>%
    rbind.data.frame(RiskScore.Cancer)
  
  COX.Summary <- summary(multi_variate_cox)
  CIndex <- data.frame(C.Index = COX.Summary$concordance[1]) %>%
    mutate(Cancer=project) %>%
    rbind.data.frame(CIndex)
}
#write.csv(CIndex,file = "RiskGroup.COX.TCGA.CIndex.csv",row.names = F)
#write.csv(RiskScore.Cancer,file = "RiskGroup.COX.TCGA.RiskScore.csv",row.names = F)

#### Figure => C index + Gene count => 23 Cancers ####
BCAA.ID <- read.csv("BCAA.ID.csv")
LASSO.Feature <- read.csv("Risk.LASSO.COX.csv")

CIndex <- read.csv(file = "RiskGroup.COX.TCGA.CIndex.csv") %>%
  #arrange(C.Index) %>%
  mutate(Cancer=factor(Cancer,levels = rev(sort(unique(.$Cancer)))))
Immune.Col <- c("#293462","#a64942","#fe5f55","#fff1c1",
                "#5bd1d7","#348498","#004d61","#ff502f",
                "#e41749","#f5587b","#ff8a5c","#fff591",
                "#001871","#ff585d","#ffb549","#41b6e6",
                "#515bd4","#8134af","#dd2a7b","#feda77",
                "#96ceb4","#ffeead","#d9534f","#ffad60",
                "#05445c","#f2317f","#5c4f74","#040000",
                "#de4307","#f29c2b","#f6d04d","#8bc24c",
                "#fef4a9","#3b9a9c","#4bc2c5","#78fee0",
                "#fad3cf","#a696c8","#2470a0","#c50d66",
                "#f07810","#eec60a","#dd7777","#1687a7",
                "#014955","#dd0a35","#ED5485","#FFE869",
                "#bc8420","#F0B775","#D25565","#2E94B9")

p1.B <- ggplot(CIndex,aes(x=C.Index,y=Cancer,fill=Cancer))+
  geom_segment(aes(x=0.5,xend=C.Index,y=Cancer,yend=Cancer))+
  geom_point(shape=21,color="black",size=3)+
  ggthemes::theme_few()+
  theme(legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9,angle = 90,vjust = 0.3,hjust = 1),
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="C index",y="")+
  scale_fill_manual(values = Immune.Col)+
  scale_x_continuous(limits = c(0.5,1))+
  geom_vline(xintercept = seq(0.5,0.9,0.2),linetype=2)
  #ggthemes::scale_color_stata()
Feature.Count <- LASSO.Feature %>%
  filter(Count >= 50) %>%
  group_by(Cancer) %>%
  summarise(Count=n()) %>%
  filter(Count >= 2) %>%
  mutate(Cancer=factor(Cancer,levels = levels(CIndex$Cancer)))
library(RColorBrewer)
Cancer.Col <- c(brewer.pal(9,"Set1"),
                brewer.pal(8,"Set2")[1:7],
                brewer.pal(12,"Set3")[1:3],
                brewer.pal(8,"Dark2")[6],
                brewer.pal(12,"Set3")[4:12],
                brewer.pal(8,"Dark2")[1:4])
p1.C <- ggplot(Feature.Count,aes(Count,Cancer,fill=Cancer))+
  geom_segment(aes(x=0,xend=Count,y=Cancer,yend=Cancer))+
  geom_point(shape=21,color="black",size=3)+
  ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9,angle = 90,vjust = 0.3,hjust = 1),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="Gene count",y="")+
  #ggthemes::scale_color_stata()
  scale_fill_manual(values = Immune.Col)+
  scale_x_continuous(limits = c(0,16))

p1 <- p1.B %>% aplot::insert_right(p1.C)
ggsave(p1,filename = "RiskGroup.COX.TCGA.CIndex.Figure.pdf",height = 5,width = 3.2)

p1.B %>% aplot::insert_right(p1.C) %>%
  aplot::insert_right(p2)

p<-LASSO.Feature %>%
  filter(Count >= 50) %>%
  filter(Cancer %in% Feature.Count$Cancer) %>%
  merge(BCAA.ID,by.x="Feature",by.y="ENSEMBL") %>%
  mutate(Count=1) %>%
  #mutate(Cancer=factor(Cancer,levels = rev(sort(unique(.$Cancer))))) %>%
  ggplot(aes(Cancer,Count,fill=SYMBOL))+
  geom_col(color="black")+
  #geom_text(aes(label=Count),hjust=1)+
  ggthemes::theme_few()+
  theme(legend.position = "bottom",
        axis.text = element_text(color = "black"),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size=9,angle = 90,vjust = 0.5),
        axis.ticks.y = element_blank(),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="",y="",fill="")+
  #ggthemes::scale_color_stata()
  scale_fill_manual(values = Cancer.Col)#+
  #scale_x_continuous(limits = c(0,20))
ggsave(p,filename = "RiskGroup.COX.TCGA.CIndex.CancerWithGene.pdf",width = 3.5,height = 4.8)
ggsave(p,filename = "RiskGroup.COX.TCGA.CIndex.CancerWithGene2.pdf",width = 7,height = 4)


dataset.Sankey <- LASSO.Feature %>%
  filter(Count >= 50) %>%
  filter(Cancer %in% Feature.Count$Cancer) %>%
  merge(BCAA.ID,by.x="Feature",by.y="ENSEMBL") %>%
  mutate(Count=1)

write.csv(dataset.Sankey,file = "dataset.Sankey.csv",row.names = F)
dataset.Sankey <- read.csv("dataset.Sankey.csv")
library(sankeywheel)
sankeywheel(from = dataset.Sankey$Cancer,
            to = dataset.Sankey$SYMBOL,
            weight = dataset.Sankey$Count,
            title = "",
            type = "dependencywheel",
            theme = "sunset",
            width = "200%")

Data.Circle <- dataset.Sankey %>% 
  dplyr::select(Cancer,SYMBOL,Count) %>%
  spread(Cancer,Count,fill=0) %>%
  remove_rownames() %>%
  column_to_rownames("SYMBOL") %>%
  as.matrix()
library(circlize)
circos.par(gap.after = c(rep(5, nrow(mat)-1), 15, rep(5, ncol(mat)-1), 15),
           start.degree = 85)
chordDiagram(Data.Circle, order = c(rev(colnames(Data.Circle)), rev(rownames(Data.Circle))))
circos.clear()



Feature.Filter <- LASSO.Feature %>%
  filter(Count >= 50) %>%
  filter(Cancer %in% Feature.Count$Cancer) %>%
  merge(BCAA.ID,by.x="Feature",by.y="ENSEMBL") %>%
  group_by(SYMBOL) %>%
  summarise(Count=n()) %>%
  arrange(Count) %>%
  mutate(SYMBOL=factor(SYMBOL,levels = c(.$SYMBOL)))

p2 <- ggplot(Feature.Filter,aes(Count,SYMBOL,fill=SYMBOL))+
  geom_col()+
  geom_text(aes(label=Count),hjust=1)+
  ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9,angle = 90,vjust = 0.3,hjust = 1),
        axis.text.y = element_text(size=9),
        #axis.ticks.y = element_blank(),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="No. of cancer",y="")+
  #ggthemes::scale_color_stata()
  scale_fill_manual(values = Immune.Col)+
  scale_x_continuous(limits = c(0,20))

ggsave(p2,filename = "RiskGroup.COX.TCGA.CIndex.GeneForCancer.pdf",width = 5,height = 4)

#### Figure => TimeROC => 23 Cancers ######
setwd("K:/TCGA/Anlysis/BCAA")
library(tidyverse)
library(survminer)
library(survivalROC)
library(survival)
dir.create("riskScore.timeROC")
AUC.Data <- data.frame()
RiskScore.Cancer <- read.csv(file = "RiskGroup.COX.TCGA.RiskScore.csv")
for (project in unique(RiskScore.Cancer$Cancer)) {
  middata <- RiskScore.Cancer %>% filter(Cancer==project)
  library(timeROC)
  with(middata,
       ROC_riskscore <<- timeROC(T = OS.time,
                                 delta = OS,
                                 marker = total_risk_score,
                                 cause = 1,
                                 weighting = "marginal",
                                 times = c(365,1080,1800),
                                 ROC = TRUE,
                                 iid = TRUE)
  )
  #plot(ROC_riskscore, time = 365, col = "red", add = F,title = "")
  #plot(ROC_riskscore, time = 1080, col = "blue", add = T)
  #plot(ROC_riskscore, time = 1800, col = "purple", add = T)
  #legend("bottomright",c("1-Year","3-Year","5-Year"),col=c("red","blue","purple"),lty=1,lwd=2)
  #text(0.5,0.2,paste("1-Year AUC = ",round(ROC_riskscore$AUC[1],3)))
  #text(0.5,0.15,paste("3-Year AUC = ",round(ROC_riskscore$AUC[2],3)))
  #text(0.5,0.1,paste("5-Year AUC = ",round(ROC_riskscore$AUC[3],3)))
  dat = data.frame(fpr = as.numeric(ROC_riskscore$FP),
                   tpr = as.numeric(ROC_riskscore$TP),
                   time = rep(as.factor(c(365,1095,1825)),each = nrow(ROC_riskscore$TP)))
  write.csv(dat,file = file.path("riskScore.timeROC",paste("ROC.Data.",project,".csv",sep = "")),row.names = F)
  
  AUC.Data <- data.frame(year1=round(ROC_riskscore$AUC,2)[1],
             year3=round(ROC_riskscore$AUC,2)[2],
             year5=round(ROC_riskscore$AUC,2)[3]) %>%
    mutate(Cancer=project) %>%
    rbind.data.frame(AUC.Data)
  
  library(ggplot2)
  p <- ggplot() + 
    geom_line(data = dat,aes(x = fpr, y = tpr,color = time),size = 1) + 
    scale_color_manual(name = NULL,values = c("#E64B35FF", "#3C5488FF", "#00A087FF"),
                       labels = paste0("AUC of ",c(1,3,5),"-y survival: ",
                                       format(round(ROC_riskscore$AUC,2),nsmall = 2)))+
    #ggsci::scale_color_npg()+
    geom_line(aes(x=c(0,1),y=c(0,1)),color = "grey",linetype=2)+
    ggthemes::theme_few()+
    theme(#legend.position = "top",
      axis.text = element_text(color = "black"),
      axis.text.x = element_text(size=9,angle = 90,vjust = 0.3,hjust = 1),
      axis.text.y = element_text(size=9),
      #axis.ticks.x = element_blank(),
      axis.title = element_text(size=12),
      panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
    labs(x = "1 - Specificity",
         y = "Sensitivity")+
    coord_fixed()
  ggsave(p,filename = file.path("riskScore.timeROC",paste("ROC.",project,".pdf",sep = "")),
         height = 3,width = 5)
  #result <-with(new_dat, 
  #survivalROC(Stime=time,
  #            status=event,
  #            marker=fp,
  #            predict.time=365,
  #            method="KM"))
}
write.csv(AUC.Data,file = file.path("riskScore.timeROC","AUC.Pancancer.csv"),row.names = F)
#### Figure ####
AUC.Data <- read.csv(file.path("riskScore.timeROC","AUC.Pancancer.csv")) %>%
  gather(Survival,AUC,-Cancer) %>%
  mutate(Cancer=factor(Cancer,levels = rev(sort(unique(.$Cancer)))))


Cancer.Col <- c(brewer.pal(9,"Set1")[1:5],
                brewer.pal(8,"Set2")[1:7],
                brewer.pal(12,"Set3")[3:3],
                brewer.pal(8,"Dark2")[6],
                brewer.pal(12,"Set3")[4:6],
                brewer.pal(12,"Set3")[8:12],
                brewer.pal(8,"Dark2")[1:4])

p2 <- ggplot(AUC.Data,aes(x=AUC,y=reorder(Cancer,AUC),fill=Cancer))+
  geom_segment(aes(x=0,y=Cancer,xend=AUC,yend=Cancer))+
  facet_wrap(~Survival,ncol = 3)+
  geom_point(shape=21,color="black",size=3)+
  ggthemes::theme_few()+
  theme(legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9,angle = 90,vjust = 0.3,hjust = 1),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x = "AUROC",
       y = "")+
  scale_fill_manual(values = Immune.Col)+
  geom_vline(xintercept = seq(0.5,0.9,0.2),linetype=2)

p2.1 <- p1.B %>% aplot::insert_right(p1.C) %>%
  aplot::insert_right(p2,width = 3)
ggsave(p2.1,filename = "RiskGroup.COX.TCGA.timeROC.Figure.pdf",height = 5,width = 7)


######## Validation ######
#### PCAWG => linux #####
setwd("K:/TCGA/Anlysis/BCAA")
library(tidyverse)
library(survival)
library(survminer)
LASSO.Feature <- read.csv("Risk.LASSO.COX.csv")
BCAA.ID <- read.csv("BCAA.ID.csv")
Feature.Count <- LASSO.Feature %>%
  filter(Count >= 50) %>%
  group_by(Cancer) %>%
  summarise(Count=n()) %>%
  filter(Count >= 2)
Feature.Filter <- LASSO.Feature %>%
  filter(Count >= 50) %>%
  filter(Cancer %in% Feature.Count$Cancer)

library(readxl)
PCAWG <- read_xlsx("PCAWG/pcawg_donor_clinical_August2016_v9.xlsx")
Express <- read.csv("PCAWG/pcawg.RNA.TPM.csv")

Specimen <- read.csv("PCAWG/pcawg.sp_specimen_type.txt",sep = "\t") %>%
  filter(str_detect(dcc_specimen_type,"Primary"))
Survival.D <- read.csv("PCAWG/pcawg.survival_sp.txt",sep = "\t")
Sample.Project <- read.csv("PCAWG/pcawg.project_code_sp.txt",sep = "\t") %>%
  filter(!str_detect(dcc_project_code,"-US")) %>%
  filter(dcc_project_code %in% c("OV-AU","PRAD-CA","PRAD-UK",
                                 "LIRI-JP","PACA-AU","PACA-CA",
                                 "LAML-KR","ESAD-UK","LINC-JP",
                                 "LICA-FR","EOPC-DE","BRCA-UK",
                                 "BRCA-EU","BTCA-SG","GACA-CN")) %>%
  merge(data.frame(dcc_project_code = c("OV-AU","PRAD-CA","PRAD-UK",
                                        "LIRI-JP","PACA-AU","PACA-CA",
                                        "LAML-KR","ESAD-UK","LINC-JP",
                                        "LICA-FR","EOPC-DE","BRCA-UK",
                                        "BRCA-EU","BTCA-SG","GACA-CN"),
                   Cancer = c("OV","PRAD","PRAD",
                              "LIHC","PAAD","PAAD",
                              "LAML","ESCA","LIHC",
                              "LIHC","PRAD","BRCA","BRCA",
                              "CHOL","STAD")))

Phenotype <- merge(Specimen,Survival.D,by.x="icgc_specimen_id",by.y="xena_sample") %>%
  merge(Sample.Project,by = "icgc_specimen_id")

for (project1 in unique(Phenotype$dcc_project_code)) {
  print(project1)
  Target.Data <- Phenotype %>% filter(dcc_project_code == project1)
  project2 <- Target.Data$Cancer %>% unique()
  if (project2 %in% Feature.Filter$Cancer) {
    middata <- Feature.Filter %>% filter(Cancer == project2)
    Exp <- Express %>% filter(X %in% middata$Feature) %>%
      remove_rownames() %>%
      column_to_rownames("X") %>%
      t() %>% data.frame(check.names = F) %>%
      rownames_to_column("icgc_specimen_id") %>%
      merge(Target.Data,by="icgc_specimen_id")
    print(nrow(Exp))
    if (nrow(Exp) >= 0) {
      saveRDS(Exp,file = paste("BCAA.riskScore.",project1,".",project2,".RDS",sep = ""))
    }
  }
}

###### Validation => Function ####
setwd("K:/TCGA/Anlysis/BCAA")
library(tidyverse)
library(survival)
library(survminer)
LASSO.Feature <- read.csv("Risk.LASSO.COX.csv")
BCAA.ID <- read.csv("BCAA.ID.csv")
Feature.Count <- LASSO.Feature %>%
  filter(Count >= 50) %>%
  group_by(Cancer) %>%
  summarise(Count=n()) %>%
  filter(Count >= 2)

Feature.Filter <- LASSO.Feature %>%
  filter(Count >= 50) %>%
  filter(Cancer %in% Feature.Count$Cancer) %>%
  merge(BCAA.ID,by.x = "Feature",by.y = "ENSEMBL")

Cal_function <- function(dataset,feature,project,timeVector){
  return.list <- list()
  formula <- as.formula(paste0("Surv(OS.time,OS)~",
                               paste(feature,collapse = "+")))
  dataset$OS <- as.numeric(as.double(dataset$OS))
  dataset$OS.time <- as.numeric(as.double(dataset$OS.time))
  multi_variate_cox <- coxph(formula,data=dataset)
  dataset$total_risk_score <- predict(multi_variate_cox,
                                      type="risk",
                                      newdata=dataset[,feature])
  return.list[[1]] = dataset
  COX.Summary <- summary(multi_variate_cox)
  return.list[[2]] = COX.Summary$concordance[1]
  library(timeROC)
  with(dataset,
       ROC_riskscore <<- timeROC(T = OS.time,
                                 delta = OS,
                                 marker = total_risk_score,
                                 cause = 1,
                                 weighting = "marginal",
                                 times = timeVector,# c(12,36,60),
                                 ROC = TRUE,
                                 iid = TRUE)
  )
  return.list[[3]] = ROC_riskscore
  
  res.cut <- surv_cutpoint(dataset,
                           time = paste("OS",".time",sep = ""),
                           event = "OS",
                           variables = "total_risk_score")
  dat <- surv_categorize(res.cut) #%>% data.frame() %>% mutate(Sample=SNHGs.ICI.data$Sample.ID)
  fit <- survfit(Surv(OS.time, OS) ~ total_risk_score,
                 data = dat)
  index="OS"
  diff=survdiff(as.formula(paste("Surv(",index,".time,",index,") ~ ","total_risk_score",sep = "")),data = dat)
  pValue=1-pchisq(diff$chisq,df=1)
  return.list[[4]] = pValue
  library(ggthemes)
  theme.sur <- ggthemes::theme_few()+
    theme(legend.position = "none",
          axis.text = element_text(color = "black"),
          axis.text.x = element_text(size=15),
          axis.text.y = element_text(size=15),
          axis.title = element_text(size=18),
          panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))
  
  p<-ggsurvplot(fit,pval =TRUE, data = dat, 
                surv.median.line = "hv",
                legend.title = "Risk",
                conf.int.style = "step",
                xlab = "Time in days",
                #break.time.by = 500,
                risk.table = "abs_pct",
                #risk.table.y.text.col = T,
                #risk.table.y.text = FALSE,
                legend.labs = c("High", "Low"),
                #pval = TRUE,
                conf.int = TRUE,
                palette = "Set1",
                #ggtheme = theme.sur,
                risk.table.y.text.col = T,
                risk.table.y.text = F,
                ggtheme = theme.sur)
  print(p)
  #dev.off()
  library(export)
  graph2pdf(file=paste("RiskGroup/Validation.OS.",project,".RiskScore.BestCutoff.KM.pdf",sep = ''),height = 5,width = 4)
  names(return.list) <- c("Dataset","C.index","timeROC","KM.Pvalue")
  return(return.list)
}

#### LIHC => LIRI-JP ####
LIHC.Feature <- Feature.Filter %>% filter(Cancer == "LIHC")
Exp.Samples <- read.csv("../../../../HepG2-CBX2/Validate/ICGC-LIRI-JP/HCCDB18.sample.txt",sep = "\t",check.names = F,row.names = 1) %>%
  t() %>% data.frame(check.names = F) %>% 
  filter(TYPE=="HCC") %>%
  rownames_to_column("SampleID") %>%
  dplyr::select(SampleID,TYPE,PATIENT_ID)
Exp.Clinical <- read.csv("../../../../HepG2-CBX2/Validate/ICGC-LIRI-JP/HCCDB18.patient.txt",sep = "\t",check.names = F,row.names = 1) %>%
  t() %>% data.frame(check.names = F) %>%
  rownames_to_column("PATIENT_ID") %>%
  mutate(OS=if_else(STATUS=="Dead",1,0)) %>%
  mutate(OS.time=SUR) %>%
  dplyr::select(PATIENT_ID,PATIENT,OS,OS.time) %>%
  merge(Exp.Samples,by="PATIENT_ID")

Exp <- read.csv("../../../../HepG2-CBX2/Validate/ICGC-LIRI-JP/HCCDB18_mRNA_level3.txt",sep = "\t",check.names = F)

Exp.Data <- Exp[,-1] %>% data.frame(check.names = F) %>%
  remove_rownames() %>%
  filter(Symbol %in% LIHC.Feature$SYMBOL) %>%
  column_to_rownames("Symbol") %>%
  t() %>% data.frame(check.names = F) %>%
  rownames_to_column("SampleID") %>%
  merge(Exp.Clinical,by="SampleID")

DATA <- Cal_function(Exp.Data,LIHC.Feature$SYMBOL,"LIHC.LIRI-JP",c(12,36,60))
saveRDS(DATA,file = "riskScore.Validation.LIHC.LIRI-JP.RDS")
#### COAD => GSE28722 ####
library(data.table)
COAD.Feature <- Feature.Filter %>% filter(Cancer == "COAD")
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE28722.gene_expression_array.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  filter(Gene_Symbol %in% COAD.Feature$SYMBOL) %>%
  column_to_rownames("Gene_Symbol") %>%
  t() %>% data.frame(check.names = F)
#Exp.Data <- log2(Exp.Data+1)
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE28722.clinical.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  mutate(SampleID=str_remove_all(Patient,".*, ")) %>%
  filter(Sample_Type == "Tumor") %>%
  dplyr::select(Accession,OS_Time,OS_Status) %>%
  magrittr::set_colnames(c("Accession","OS.time","OS")) %>%
  merge(Exp.Data %>% rownames_to_column("Accession"),by="Accession")

DATA <- Cal_function(Exp.Data,COAD.Feature$SYMBOL,"COAD.GSE28722",c(365,1080,1800))
saveRDS(DATA,file = "riskScore.Validation.COAD.GSE28722.RDS")
#### COAD => GSE103479 ####
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE103479.gene_expression_array.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  filter(Gene_Symbol %in% COAD.Feature$SYMBOL) %>%
  column_to_rownames("Gene_Symbol") %>%
  t() %>% data.frame(check.names = F)
#Exp.Data <- log2(Exp.Data+1)
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE103479.clinical.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  mutate(SampleID=str_remove_all(Patient,".*, ")) %>%
  filter(Sample_Type == "Tumor") %>%
  dplyr::select(Accession,OS_Time,OS_Status) %>%
  magrittr::set_colnames(c("Accession","OS.time","OS")) %>%
  merge(Exp.Data %>% rownames_to_column("Accession"),by="Accession")

DATA <- Cal_function(Exp.Data,COAD.Feature$SYMBOL,"COAD.GSE103479",c(365,1080,1800))
saveRDS(DATA,file = "riskScore.Validation.COAD.GSE103479.RDS")
#### LUSC => GSE157009 ####
LUSC.Feature <- Feature.Filter %>% filter(Cancer == "LUSC")
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE157009.gene_expression_array.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  filter(Gene_Symbol %in% LUSC.Feature$SYMBOL) %>%
  column_to_rownames("Gene_Symbol") %>%
  t() %>% data.frame(check.names = F)
#Exp.Data <- log2(Exp.Data+1)
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE157009.clinical.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  mutate(SampleID=str_remove_all(Patient,".*, ")) %>%
  filter(Sample_Type == "Tumor") %>%
  dplyr::select(Accession,OS_Time,OS_Status) %>%
  magrittr::set_colnames(c("Accession","OS.time","OS")) %>%
  merge(Exp.Data %>% rownames_to_column("Accession"),by="Accession")

DATA <- Cal_function(Exp.Data,LUSC.Feature$SYMBOL,"LUSC.GSE157009",c(365,1080,1800))
saveRDS(DATA,file = "riskScore.Validation.LUSC.GSE157009.RDS")
#### LUSC => GSE37745 ####
LUSC.Feature <- Feature.Filter %>% filter(Cancer == "LUSC")
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE37745.gene_expression_array.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  filter(Gene_Symbol %in% LUSC.Feature$SYMBOL) %>%
  column_to_rownames("Gene_Symbol") %>%
  t() %>% data.frame(check.names = F)
#Exp.Data <- log2(Exp.Data+1)
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE37745.clinical.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  mutate(SampleID=str_remove_all(Patient,".*, ")) %>%
  filter(Sample_Type == "Tumor") %>%
  dplyr::select(Accession,OS_Time,OS_Status) %>%
  magrittr::set_colnames(c("Accession","OS.time","OS")) %>%
  merge(Exp.Data %>% rownames_to_column("Accession"),by="Accession")

DATA <- Cal_function(Exp.Data,LUSC.Feature$SYMBOL,"LUSC.GSE37745",c(365,1080,1800))
saveRDS(DATA,file = "riskScore.Validation.LUSC.GSE37745.RDS")
#### LUSC => GSE50081 ####
LUSC.Feature <- Feature.Filter %>% filter(Cancer == "LUSC")
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE50081.gene_expression_array.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  filter(Gene_Symbol %in% LUSC.Feature$SYMBOL) %>%
  column_to_rownames("Gene_Symbol") %>%
  t() %>% data.frame(check.names = F)
#Exp.Data <- log2(Exp.Data+1)
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE50081.clinical.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  mutate(SampleID=str_remove_all(Patient,".*, ")) %>%
  filter(Sample_Type == "Tumor") %>%
  dplyr::select(Accession,OS_Time,OS_Status) %>%
  magrittr::set_colnames(c("Accession","OS.time","OS")) %>%
  merge(Exp.Data %>% rownames_to_column("Accession"),by="Accession")

DATA <- Cal_function(Exp.Data,LUSC.Feature$SYMBOL,"LUSC.GSE50081",c(365,1080,1800))
saveRDS(DATA,file = "riskScore.Validation.LUSC.GSE50081.RDS")
#### BLCA => GSE31684 ####
BLCA.Feature <- Feature.Filter %>% filter(Cancer == "BLCA")
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE31684.gene_expression_array.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  filter(Gene_Symbol %in% BLCA.Feature$SYMBOL) %>%
  column_to_rownames("Gene_Symbol") %>%
  t() %>% data.frame(check.names = F)
#Exp.Data <- log2(Exp.Data+1)
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE31684.clinical.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  mutate(SampleID=str_remove_all(Patient,".*, ")) %>%
  filter(Sample_Type == "Tumor") %>%
  dplyr::select(Accession,OS_Time,OS_Status) %>%
  magrittr::set_colnames(c("Accession","OS.time","OS")) %>%
  merge(Exp.Data %>% rownames_to_column("Accession"),by="Accession")

DATA <- Cal_function(Exp.Data,BLCA.Feature$SYMBOL,"BLCA.GSE31684",c(365,1080,1800))
saveRDS(DATA,file = "riskScore.Validation.BLCA.GSE31684.RDS")
#### BLCA => GSE69795 ####
BLCA.Feature <- Feature.Filter %>% filter(Cancer == "BLCA")
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE69795.gene_expression_array.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  filter(Gene_Symbol %in% BLCA.Feature$SYMBOL) %>%
  column_to_rownames("Gene_Symbol") %>%
  t() %>% data.frame(check.names = F)
#Exp.Data <- log2(Exp.Data+1)
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE69795.clinical.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  mutate(SampleID=str_remove_all(Patient,".*, ")) %>%
  filter(Sample_Type == "Tumor") %>%
  dplyr::select(Accession,OS_Time,OS_Status) %>%
  magrittr::set_colnames(c("Accession","OS.time","OS")) %>%
  merge(Exp.Data %>% rownames_to_column("Accession"),by="Accession")

DATA <- Cal_function(Exp.Data,BLCA.Feature$SYMBOL,"BLCA.GSE69795",c(365,1080,1800))
saveRDS(DATA,file = "riskScore.Validation.BLCA.GSE69795.RDS")
#### ESCA => GSE53625 ####
ESCA.Feature <- Feature.Filter %>% filter(Cancer == "ESCA")
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE53625.gene_expression_array.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  filter(Gene_Symbol %in% ESCA.Feature$SYMBOL) %>%
  column_to_rownames("Gene_Symbol") %>%
  t() %>% data.frame(check.names = F)
#Exp.Data <- log2(Exp.Data+1)
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE53625.clinical.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  mutate(SampleID=str_remove_all(Patient,".*, ")) %>%
  filter(Sample_Type == "Tumor") %>%
  dplyr::select(Accession,OS_Time,OS_Status) %>%
  magrittr::set_colnames(c("Accession","OS.time","OS")) %>%
  merge(Exp.Data %>% rownames_to_column("Accession"),by="Accession")

DATA <- Cal_function(Exp.Data,ESCA.Feature$SYMBOL,"ESCA.GSE53625",c(365,1080,1800))
saveRDS(DATA,file = "riskScore.Validation.ESCA.GSE53625.RDS")

#### ESCA => GSE53624 ####
ESCA.Feature <- Feature.Filter %>% filter(Cancer == "ESCA")
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE53624.gene_expression_array.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  filter(Gene_Symbol %in% ESCA.Feature$SYMBOL) %>%
  column_to_rownames("Gene_Symbol") %>%
  t() %>% data.frame(check.names = F)
#Exp.Data <- log2(Exp.Data+1)
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE53624.clinical.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  mutate(SampleID=str_remove_all(Patient,".*, ")) %>%
  filter(Sample_Type == "Tumor") %>%
  dplyr::select(Accession,OS_Time,OS_Status) %>%
  magrittr::set_colnames(c("Accession","OS.time","OS")) %>%
  merge(Exp.Data %>% rownames_to_column("Accession"),by="Accession")

DATA <- Cal_function(Exp.Data,ESCA.Feature$SYMBOL,"ESCA.GSE53624",c(365,1080,1800))
saveRDS(DATA,file = "riskScore.Validation.ESCA.GSE53624.RDS")
#### HNSC => GSE65858 ####
HNSC.Feature <- Feature.Filter %>% filter(Cancer == "HNSC")
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE65858.gene_expression_array.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  filter(Gene_Symbol %in% HNSC.Feature$SYMBOL) %>%
  column_to_rownames("Gene_Symbol") %>%
  t() %>% data.frame(check.names = F)
#Exp.Data <- log2(Exp.Data+1)
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE65858.clinical.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  mutate(SampleID=str_remove_all(Patient,".*, ")) %>%
  filter(Sample_Type == "Tumor") %>%
  dplyr::select(Accession,OS_Time,OS_Status) %>%
  magrittr::set_colnames(c("Accession","OS.time","OS")) %>%
  merge(Exp.Data %>% rownames_to_column("Accession"),by="Accession")

DATA <- Cal_function(Exp.Data,HNSC.Feature$SYMBOL,"HNSC.GSE65858",c(365,1080,1800))
saveRDS(DATA,file = "riskScore.Validation.HNSC.GSE65858.RDS")
#### KIRC => GSE167573 ####
KIRC.Feature <- Feature.Filter %>% filter(Cancer == "KIRC")
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE167573.gene_expression_RNAseq.Normalized_Counts.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  filter(Gene_Symbol %in% KIRC.Feature$SYMBOL) %>%
  column_to_rownames("Gene_Symbol") %>%
  t() %>% data.frame(check.names = F)
Exp.Data <- log2(Exp.Data+1)
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE167573.clinical.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  mutate(SampleID=str_remove_all(Patient,".*, ")) %>%
  filter(Sample_Type == "Tumor") %>%
  dplyr::select(SampleID,OS_Time,OS_Status) %>%
  magrittr::set_colnames(c("SampleID","OS.time","OS")) %>%
  merge(Exp.Data %>% rownames_to_column("SampleID"),by="SampleID")
Exp.Data <- na.omit(Exp.Data)

DATA <- Cal_function(Exp.Data,KIRC.Feature$SYMBOL,"KIRC.GSE167573",c(365,1080,1800))
saveRDS(DATA,file = "riskScore.Validation.KIRC.GSE167573.RDS")

#### LAML => GSE71014 ####
LAML.Feature <- Feature.Filter %>% filter(Cancer == "LAML")
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE71014.gene_expression_array.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  filter(Gene_Symbol %in% LAML.Feature$SYMBOL) %>%
  column_to_rownames("Gene_Symbol") %>%
  t() %>% data.frame(check.names = F)
#Exp.Data <- log2(Exp.Data+1)
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE71014.clinical.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  mutate(SampleID=str_remove_all(Patient,".*, ")) %>%
  filter(Sample_Type == "Tumor") %>%
  dplyr::select(Accession,OS_Time,OS_Status) %>%
  magrittr::set_colnames(c("Accession","OS.time","OS")) %>%
  merge(Exp.Data %>% rownames_to_column("Accession"),by="Accession")

DATA <- Cal_function(Exp.Data,LAML.Feature$SYMBOL,"LAML.GSE71014",c(365,1080,1800))
saveRDS(DATA,file = "riskScore.Validation.LAML.GSE71014.RDS")
#### LAML => GSE106291 ####
LAML.Feature <- Feature.Filter %>% filter(Cancer == "LAML")
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE106291.gene_expression_RNAseq.Counts.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  filter(Gene_Symbol %in% LAML.Feature$SYMBOL) %>%
  column_to_rownames("Gene_Symbol") %>%
  t() %>% data.frame(check.names = F)
Exp.Data <- log2(Exp.Data+1)
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE106291.clinical.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  mutate(SampleID=str_remove_all(Patient,".*, ")) %>%
  filter(Sample_Type == "Tumor") %>%
  dplyr::select(SampleID,OS_Time,OS_Status) %>%
  magrittr::set_colnames(c("SampleID","OS.time","OS")) %>%
  mutate(SampleID=str_remove_all(SampleID,"-")) %>%
  merge(Exp.Data %>% rownames_to_column("SampleID") %>% 
          mutate(SampleID=str_remove_all(SampleID,"-")),
        by="SampleID")

DATA <- Cal_function(Exp.Data,LAML.Feature$SYMBOL,"LAML.GSE106291",c(365,1080,1800))
saveRDS(DATA,file = "riskScore.Validation.LAML.GSE106291.RDS")
#### LAML => GSE146173 ####
LAML.Feature <- Feature.Filter %>% filter(Cancer == "LAML")
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE146173.gene_expression_RNAseq.Salmon_Estimated_Counts.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  filter(Gene_Symbol %in% LAML.Feature$SYMBOL) %>%
  column_to_rownames("Gene_Symbol") %>%
  t() %>% data.frame(check.names = F)
Exp.Data <- log2(Exp.Data+1)
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE146173.clinical.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  mutate(SampleID=str_remove_all(Patient,".*, ")) %>%
  filter(Sample_Type == "Tumor") %>%
  dplyr::select(SampleID,OS_Time,OS_Status) %>%
  magrittr::set_colnames(c("SampleID","OS.time","OS")) %>%
  mutate(SampleID=str_remove_all(SampleID,"-") %>% str_remove_all(":.*")) %>%
  merge(Exp.Data %>% rownames_to_column("SampleID") %>% 
          mutate(SampleID=str_remove_all(SampleID,"-")),
        by="SampleID")

DATA <- Cal_function(Exp.Data,LAML.Feature$SYMBOL,"LAML.GSE146173",c(365,1080,1800))
saveRDS(DATA,file = "riskScore.Validation.LAML.GSE146173.RDS")

#### LAML => GSE12417 ####
LAML.Feature <- Feature.Filter %>% filter(Cancer == "LAML")
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE12417.gene_expression_array.GPL570.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  filter(Gene_Symbol %in% LAML.Feature$SYMBOL) %>%
  column_to_rownames("Gene_Symbol") %>%
  t() %>% data.frame(check.names = F)
#Exp.Data <- log2(Exp.Data+1)
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE12417.clinical.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  mutate(SampleID=str_remove_all(Patient,".*, ")) %>%
  filter(Sample_Type == "Tumor") %>%
  dplyr::select(Accession,OS_Time,OS_Status) %>%
  magrittr::set_colnames(c("Accession","OS.time","OS")) %>%
  merge(Exp.Data %>% rownames_to_column("Accession"),by="Accession")

DATA <- Cal_function(Exp.Data,LAML.Feature$SYMBOL,"LAML.GSE12417",c(365,1080,1800))
saveRDS(DATA,file = "riskScore.Validation.LAML.GSE12417.RDS")
#### LGG  => GSE121720 ####
LGG.Feature <- Feature.Filter %>% filter(Cancer == "LGG")
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE121720.gene_expression_RNAseq.TPM.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  filter(Gene_Symbol %in% LGG.Feature$SYMBOL) %>%
  column_to_rownames("Gene_Symbol") %>%
  t() %>% data.frame(check.names = F)
Exp.Data <- log2(Exp.Data+1)
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE121720.clinical.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  mutate(SampleID=str_remove_all(Patient," .*")) %>%
  filter(Sample_Type == "Tumor") %>%
  dplyr::select(SampleID,OS_Time,OS_Status) %>%
  magrittr::set_colnames(c("SampleID","OS.time","OS")) %>%
  merge(Exp.Data %>% rownames_to_column("SampleID"),by="SampleID")

DATA <- Cal_function(Exp.Data,LGG.Feature$SYMBOL,"LGG.GSE121720",c(365,1080,1800))
saveRDS(DATA,file = "riskScore.Validation.LGG.GSE121720.RDS")

#### LGG  => GSE4412 ####
LGG.Feature <- Feature.Filter %>% filter(Cancer == "LGG")
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE4412.gene_expression_array.GPL96.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  filter(Gene_Symbol %in% LGG.Feature$SYMBOL) %>%
  column_to_rownames("Gene_Symbol") %>%
  t() %>% data.frame(check.names = F)
#Exp.Data <- log2(Exp.Data+1)
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE4412.clinical.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  mutate(SampleID=str_remove_all(Patient,".*, ")) %>%
  filter(Sample_Type == "Tumor") %>%
  dplyr::select(Accession,OS_Time,OS_Status) %>%
  magrittr::set_colnames(c("Accession","OS.time","OS")) %>%
  merge(Exp.Data %>% rownames_to_column("Accession"),by="Accession")

DATA <- Cal_function(Exp.Data,LGG.Feature$SYMBOL,"LGG.GSE4412",c(365,1080,1800))
saveRDS(DATA,file = "riskScore.Validation.LGG.GSE4412.RDS")

#### LGG  => GSE7696 ####
LGG.Feature <- Feature.Filter %>% filter(Cancer == "LGG")
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE7696.gene_expression_array.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  filter(Gene_Symbol %in% LGG.Feature$SYMBOL) %>%
  column_to_rownames("Gene_Symbol") %>%
  t() %>% data.frame(check.names = F)
#Exp.Data <- log2(Exp.Data+1)
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE7696.clinical.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  mutate(SampleID=str_remove_all(Patient,".*, ")) %>%
  filter(Sample_Type == "Tumor") %>%
  dplyr::select(Accession,OS_Time,OS_Status) %>%
  magrittr::set_colnames(c("Accession","OS.time","OS")) %>%
  merge(Exp.Data %>% rownames_to_column("Accession"),by="Accession")

DATA <- Cal_function(Exp.Data,LGG.Feature$SYMBOL,"LGG.GSE7696",c(365,1080,1800))
saveRDS(DATA,file = "riskScore.Validation.LGG.GSE7696.RDS")
#### LIHC => GSE76427 ####
LIHC.Feature <- Feature.Filter %>% filter(Cancer == "LIHC")
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE76427.gene_expression_array.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  filter(Gene_Symbol %in% LIHC.Feature$SYMBOL) %>%
  column_to_rownames("Gene_Symbol") %>%
  t() %>% data.frame(check.names = F)
#Exp.Data <- log2(Exp.Data+1)
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE76427.clinical.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  mutate(SampleID=str_remove_all(Patient,".*, ")) %>%
  filter(Sample_Type == "Tumor") %>%
  dplyr::select(Accession,OS_Time,OS_Status) %>%
  magrittr::set_colnames(c("Accession","OS.time","OS")) %>%
  merge(Exp.Data %>% rownames_to_column("Accession"),by="Accession")

DATA <- Cal_function(Exp.Data,LIHC.Feature$SYMBOL,"LIHC.GSE76427",c(365,1080,1800))
saveRDS(DATA,file = "riskScore.Validation.LIHC.GSE76427.RDS")
#### LIHC => GSE116174 ####
LIHC.Feature <- Feature.Filter %>% filter(Cancer == "LIHC")
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE116174.gene_expression_array.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  filter(Gene_Symbol %in% LIHC.Feature$SYMBOL) %>%
  column_to_rownames("Gene_Symbol") %>%
  t() %>% data.frame(check.names = F)
#Exp.Data <- log2(Exp.Data+1)
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE116174.clinical.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  mutate(SampleID=str_remove_all(Patient,".*, ")) %>%
  dplyr::select(Accession,OS_Time,OS_Status) %>%
  magrittr::set_colnames(c("Accession","OS.time","OS")) %>%
  merge(Exp.Data %>% rownames_to_column("Accession"),by="Accession")

DATA <- Cal_function(Exp.Data,LIHC.Feature$SYMBOL,"LIHC.GSE116174",c(365,1080,1800))
saveRDS(DATA,file = "riskScore.Validation.LIHC.GSE116174.RDS")
#### LIHC => GSE144269 ####
LIHC.Feature <- Feature.Filter %>% filter(Cancer == "LIHC")
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE144269.gene_expression_RNAseq.Counts.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  filter(Gene_Symbol %in% LIHC.Feature$SYMBOL) %>%
  column_to_rownames("Gene_Symbol") %>%
  t() %>% data.frame(check.names = F)
Exp.Data <- log2(Exp.Data+1)
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE144269.clinical.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  mutate(SampleID=str_remove_all(Patient,".*\\[") %>% str_remove_all("\\]")) %>%
  filter(Sample_Type == "Tumor") %>%
  dplyr::select(SampleID,OS_Time,OS_Status) %>%
  magrittr::set_colnames(c("SampleID","OS.time","OS")) %>%
  merge(Exp.Data %>% rownames_to_column("SampleID"),by="SampleID")

DATA <- Cal_function(Exp.Data,LIHC.Feature$SYMBOL,"LIHC.GSE144269",c(365,1080,1800))
saveRDS(DATA,file = "riskScore.Validation.LIHC.GSE144269.RDS")
#### LUAD => GSE42127 ####
LUAD.Feature <- Feature.Filter %>% filter(Cancer == "LUAD")
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE42127.gene_expression_array.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  filter(Gene_Symbol %in% LUAD.Feature$SYMBOL) %>%
  column_to_rownames("Gene_Symbol") %>%
  t() %>% data.frame(check.names = F)
#Exp.Data <- log2(Exp.Data+1)
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE42127.clinical.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  mutate(SampleID=str_remove_all(Patient,".*, ")) %>%
  filter(Sample_Type == "Tumor") %>%
  dplyr::select(Accession,OS_Time,OS_Status) %>%
  magrittr::set_colnames(c("Accession","OS.time","OS")) %>%
  merge(Exp.Data %>% rownames_to_column("Accession"),by="Accession")

DATA <- Cal_function(Exp.Data,LUAD.Feature$SYMBOL,"LUAD.GSE42127",c(365,1080,1800))
saveRDS(DATA,file = "riskScore.Validation.LUAD.GSE42127.RDS")
#### LUAD => GSE30219 ####
LUAD.Feature <- Feature.Filter %>% filter(Cancer == "LUAD")
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE30219_expression.txt",header = T) %>%
  data.frame(check.names = F) %>%
  filter(`Gene Symbol` %in% LUAD.Feature$SYMBOL) %>%
  column_to_rownames("Gene Symbol") %>%
  t() %>% data.frame(check.names = F)
#Exp.Data <- log2(Exp.Data+1)
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE30219_clinical.txt",header = T) %>%
  data.frame(check.names = F) %>%
  mutate(OS=if_else(status=="DEAD",1,0),
         OS.time=`disease free survival in months`) %>%
  filter(source=="Lung Tumour") %>%
  filter(OS.time!="na") %>%
  filter(OS!="na") %>%
  dplyr::select(accession,OS.time,OS) %>%
  merge(Exp.Data %>% rownames_to_column("accession"),by="accession")

DATA <- Cal_function(Exp.Data,LUAD.Feature$SYMBOL,"LUAD.GSE30219",c(12,36,60))
saveRDS(DATA,file = "riskScore.Validation.LUAD.GSE30219.RDS")
#### LUAD => GSE31210 ####
LUAD.Feature <- Feature.Filter %>% filter(Cancer == "LUAD")
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE31210_expression.txt",header = T) %>%
  data.frame(check.names = F) %>%
  filter(`Gene Symbol` %in% LUAD.Feature$SYMBOL) %>%
  column_to_rownames("Gene Symbol") %>%
  t() %>% data.frame(check.names = F)
#Exp.Data <- log2(Exp.Data+1)
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE31210_clinical.txt",header = T) %>%
  data.frame(check.names = F) %>%
  mutate(OS=if_else(death=="dead",1,0),
         OS.time=`days before death/censor`) %>%
  #filter(source=="Lung Tumour") %>%
  filter(OS.time!="na") %>%
  filter(OS!="na") %>%
  dplyr::select(accession,OS.time,OS) %>%
  filter(OS.time!="na") %>%
  filter(OS!="na") %>%
  mutate(OS=as.numeric(as.character(OS))) %>%
  mutate(OS.time=as.numeric(as.character(OS.time))) %>%
  merge(Exp.Data %>% rownames_to_column("accession"),by="accession")

DATA <- Cal_function(Exp.Data,LUAD.Feature$SYMBOL,"LUAD.GSE31210",c(365,1080,1800))
saveRDS(DATA,file = "riskScore.Validation.LUAD.GSE31210.RDS")
#### MESO => GSE17118 ####
MESO.Feature <- Feature.Filter %>% filter(Cancer == "MESO")
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE17118.gene_expression_array.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  filter(Gene_Symbol %in% MESO.Feature$SYMBOL) %>%
  column_to_rownames("Gene_Symbol") %>%
  t() %>% data.frame(check.names = F)
#Exp.Data <- log2(Exp.Data+1)
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE17118.clinical.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  mutate(SampleID=str_remove_all(Patient,".*, ")) %>%
  filter(Sample_Type == "Tumor") %>%
  dplyr::select(Accession,OS_Time,OS_Status) %>%
  magrittr::set_colnames(c("Accession","OS.time","OS")) %>%
  merge(Exp.Data %>% rownames_to_column("Accession"),by="Accession")

DATA <- Cal_function(Exp.Data,MESO.Feature$SYMBOL,"MESO.GSE17118",c(365,1080,1800))
saveRDS(DATA,file = "riskScore.Validation.MESO.GSE17118.RDS")
#### OV   => GSE138866 ####
OV.Feature <- Feature.Filter %>% filter(Cancer == "OV")
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE138866.gene_expression_RNAseq.DESeq2_Normalized_Counts.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  filter(Gene_Symbol %in% OV.Feature$SYMBOL) %>%
  column_to_rownames("Gene_Symbol") %>%
  t() %>% data.frame(check.names = F)
Exp.Data <- log2(Exp.Data+1)
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE138866.clinical.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  mutate(SampleID=str_remove_all(Patient,".*, ")) %>%
  filter(Sample_Type == "Tumor") %>%
  dplyr::select(SampleID,OS_Time,OS_Status) %>%
  magrittr::set_colnames(c("SampleID","OS.time","OS")) %>%
  merge(Exp.Data %>% rownames_to_column("SampleID"),by="SampleID")

DATA <- Cal_function(Exp.Data,OV.Feature$SYMBOL,"OV.GSE138866",c(365,1080,1800))
saveRDS(DATA,file = "riskScore.Validation.OV.GSE138866.RDS")
#### OV   => GSE63885 ####
OV.Feature <- Feature.Filter %>% filter(Cancer == "OV")
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE63885.gene_expression_array.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  filter(Gene_Symbol %in% OV.Feature$SYMBOL) %>%
  column_to_rownames("Gene_Symbol") %>%
  t() %>% data.frame(check.names = F)
#Exp.Data <- log2(Exp.Data+1)
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE63885.clinical.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  mutate(SampleID=str_remove_all(Patient,".*, ")) %>%
  filter(Sample_Type == "Tumor") %>%
  dplyr::select(Accession,OS_Time,OS_Status) %>%
  magrittr::set_colnames(c("Accession","OS.time","OS")) %>%
  merge(Exp.Data %>% rownames_to_column("Accession"),by="Accession")

DATA <- Cal_function(Exp.Data,OV.Feature$SYMBOL,"OV.GSE63885",c(365,1080,1800))
saveRDS(DATA,file = "riskScore.Validation.OV.GSE63885.RDS")
#### OV   => GSE17260 ####
OV.Feature <- Feature.Filter %>% filter(Cancer == "OV")
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE17260.gene_expression_array.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  filter(Gene_Symbol %in% OV.Feature$SYMBOL) %>%
  column_to_rownames("Gene_Symbol") %>%
  t() %>% data.frame(check.names = F)
#Exp.Data <- log2(Exp.Data+1)
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE17260.clinical.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  mutate(SampleID=str_remove_all(Patient,".*, ")) %>%
  filter(Sample_Type == "Tumor") %>%
  dplyr::select(Accession,OS_Time,OS_Status) %>%
  magrittr::set_colnames(c("Accession","OS.time","OS")) %>%
  merge(Exp.Data %>% rownames_to_column("Accession"),by="Accession")

DATA <- Cal_function(Exp.Data,OV.Feature$SYMBOL,"OV.GSE17260",c(365,1080,1800))
saveRDS(DATA,file = "riskScore.Validation.OV.GSE17260.RDS")
#### OV   => GSE73614 ####
OV.Feature <- Feature.Filter %>% filter(Cancer == "OV")
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE73614.gene_expression_array.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  filter(Gene_Symbol %in% OV.Feature$SYMBOL) %>%
  column_to_rownames("Gene_Symbol") %>%
  t() %>% data.frame(check.names = F)
#Exp.Data <- log2(Exp.Data+1)
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE73614.clinical.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  mutate(SampleID=str_remove_all(Patient,".*, ")) %>%
  filter(Sample_Type == "Tumor") %>%
  dplyr::select(Accession,OS_Time,OS_Status) %>%
  magrittr::set_colnames(c("Accession","OS.time","OS")) %>%
  merge(Exp.Data %>% rownames_to_column("Accession"),by="Accession")

DATA <- Cal_function(Exp.Data,OV.Feature$SYMBOL,"OV.GSE73614",c(365,1080,1800))
saveRDS(DATA,file = "riskScore.Validation.OV.GSE73614.RDS")
#### OV   => GSE26712 ####
OV.Feature <- Feature.Filter %>% filter(Cancer == "OV")
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE26712.gene_expression_array.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  filter(Gene_Symbol %in% OV.Feature$SYMBOL) %>%
  column_to_rownames("Gene_Symbol") %>%
  t() %>% data.frame(check.names = F)
#Exp.Data <- log2(Exp.Data+1)
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE26712.clinical.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  mutate(SampleID=str_remove_all(Patient,".*, ")) %>%
  filter(Sample_Type == "Tumor") %>%
  dplyr::select(Accession,OS_Time,OS_Status) %>%
  magrittr::set_colnames(c("Accession","OS.time","OS")) %>%
  merge(Exp.Data %>% rownames_to_column("Accession"),by="Accession")

DATA <- Cal_function(Exp.Data,OV.Feature$SYMBOL,"OV.GSE26712",c(365,1080,1800))
saveRDS(DATA,file = "riskScore.Validation.OV.GSE26712.RDS")
#### OV   => GSE18520 ####
OV.Feature <- Feature.Filter %>% filter(Cancer == "OV")
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE18520.gene_expression_array.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  filter(Gene_Symbol %in% OV.Feature$SYMBOL) %>%
  column_to_rownames("Gene_Symbol") %>%
  t() %>% data.frame(check.names = F)
#Exp.Data <- log2(Exp.Data+1)
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE18520.clinical.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  mutate(SampleID=str_remove_all(Patient,".*, ")) %>%
  filter(Sample_Type == "Tumor") %>%
  dplyr::select(Accession,OS_Time,OS_Status) %>%
  magrittr::set_colnames(c("Accession","OS.time","OS")) %>%
  merge(Exp.Data %>% rownames_to_column("Accession"),by="Accession")

DATA <- Cal_function(Exp.Data,OV.Feature$SYMBOL,"OV.GSE18520",c(365,1080,1800))
saveRDS(DATA,file = "riskScore.Validation.OV.GSE18520.RDS")
#### OV   => GSE30161 ####
OV.Feature <- Feature.Filter %>% filter(Cancer == "OV")
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE30161.gene_expression_array.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  filter(Gene_Symbol %in% OV.Feature$SYMBOL) %>%
  column_to_rownames("Gene_Symbol") %>%
  t() %>% data.frame(check.names = F)
#Exp.Data <- log2(Exp.Data+1)
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE30161.clinical.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  mutate(SampleID=str_remove_all(Patient,".*, ")) %>%
  filter(Sample_Type == "Tumor") %>%
  dplyr::select(Accession,OS_Time,OS_Status) %>%
  magrittr::set_colnames(c("Accession","OS.time","OS")) %>%
  merge(Exp.Data %>% rownames_to_column("Accession"),by="Accession")

DATA <- Cal_function(Exp.Data,OV.Feature$SYMBOL,"OV.GSE30161",c(365,1080,1800))
saveRDS(DATA,file = "riskScore.Validation.OV.GSE30161.RDS")
#### PAAD => GSE62452 ####
PAAD.Feature <- Feature.Filter %>% filter(Cancer == "PAAD")
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE62452.gene_expression_array.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  filter(Gene_Symbol %in% PAAD.Feature$SYMBOL) %>%
  column_to_rownames("Gene_Symbol") %>%
  t() %>% data.frame(check.names = F)
#Exp.Data <- log2(Exp.Data+1)
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE62452.clinical.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  mutate(SampleID=str_remove_all(Patient,".*, ")) %>%
  filter(Sample_Type == "Tumor") %>%
  dplyr::select(Accession,OS_Time,OS_Status) %>%
  magrittr::set_colnames(c("Accession","OS.time","OS")) %>%
  merge(Exp.Data %>% rownames_to_column("Accession"),by="Accession")

DATA <- Cal_function(Exp.Data,PAAD.Feature$SYMBOL,"PAAD.GSE62452",c(365,1080,1800))
saveRDS(DATA,file = "riskScore.Validation.PAAD.GSE62452.RDS")
#### PAAD => GSE28735 ####
PAAD.Feature <- Feature.Filter %>% filter(Cancer == "PAAD")
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE28735.gene_expression_array.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  filter(Gene_Symbol %in% PAAD.Feature$SYMBOL) %>%
  column_to_rownames("Gene_Symbol") %>%
  t() %>% data.frame(check.names = F)
#Exp.Data <- log2(Exp.Data+1)
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE28735.clinical.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  mutate(SampleID=str_remove_all(Patient,".*, ")) %>%
  filter(Sample_Type == "Tumor") %>%
  dplyr::select(Accession,OS_Time,OS_Status) %>%
  magrittr::set_colnames(c("Accession","OS.time","OS")) %>%
  merge(Exp.Data %>% rownames_to_column("Accession"),by="Accession")

DATA <- Cal_function(Exp.Data,PAAD.Feature$SYMBOL,"PAAD.GSE28735",c(365,1080,1800))
saveRDS(DATA,file = "riskScore.Validation.PAAD.GSE28735.RDS")
#### PAAD => GSE85916 ####
PAAD.Feature <- Feature.Filter %>% filter(Cancer == "PAAD")
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE85916.gene_expression_array.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  filter(Gene_Symbol %in% PAAD.Feature$SYMBOL) %>%
  column_to_rownames("Gene_Symbol") %>%
  t() %>% data.frame(check.names = F)
#Exp.Data <- log2(Exp.Data+1)
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE85916.clinical.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  mutate(SampleID=str_remove_all(Patient,".*, ")) %>%
  filter(Sample_Type == "Tumor") %>%
  dplyr::select(Accession,OS_Time,OS_Status) %>%
  magrittr::set_colnames(c("Accession","OS.time","OS")) %>%
  merge(Exp.Data %>% rownames_to_column("Accession"),by="Accession")

DATA <- Cal_function(Exp.Data,PAAD.Feature$SYMBOL,"PAAD.GSE85916",c(365,1080,1800))
saveRDS(DATA,file = "riskScore.Validation.PAAD.GSE85916.RDS")
#### PAAD => GSE79668 ####
PAAD.Feature <- Feature.Filter %>% filter(Cancer == "PAAD")
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE79668.gene_expression_RNAseq.Counts.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  filter(Gene_Symbol %in% PAAD.Feature$SYMBOL) %>%
  column_to_rownames("Gene_Symbol") %>%
  t() %>% data.frame(check.names = F)
Exp.Data <- log2(Exp.Data+1)
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE79668.clinical.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  mutate(SampleID=str_remove_all(Patient,"_SL.*")%>%str_remove_all("_T$")) %>%
  filter(Sample_Type == "Tumor") %>%
  dplyr::select(SampleID,OS_Time,OS_Status) %>%
  magrittr::set_colnames(c("SampleID","OS.time","OS")) %>%
  merge(Exp.Data %>% rownames_to_column("SampleID"),by="SampleID")

DATA <- Cal_function(Exp.Data,PAAD.Feature$SYMBOL,"PAAD.GSE79668",c(365,1080,1800))
saveRDS(DATA,file = "riskScore.Validation.PAAD.GSE79668.RDS")
#### PRAD => GSE107299 ####
PRAD.Feature <- Feature.Filter %>% filter(Cancer == "PRAD")
Exp.Data <- readRDS("../AnimoAcidSensing/Validation.Survival/CPC-Gene_eSet.RDS")
Sur.Data <- Exp.Data@phenoData@data
Exp.Data <- log2(Exp.Data@assayData$exprs+1) %>% data.frame(check.names = F) %>%
  rownames_to_column("Gene_Symbol") %>%
  filter(Gene_Symbol %in% PRAD.Feature$Feature) %>%
  data.frame(check.names = F) %>% remove_rownames() %>%
  column_to_rownames("Gene_Symbol")
Exp.Data <- Exp.Data %>% t() %>%
  data.frame(check.names = F) %>%
  rownames_to_column("sample_id") %>%
  merge(Sur.Data,by="sample_id") %>%
  mutate(OS=bcr_status,
         OS.time=time_to_bcr) %>%
  filter(sample_type == "Primary")
DATA <- Cal_function(Exp.Data,PRAD.Feature$Feature,"PRAD.GSE107299",c(12,36,60))
saveRDS(DATA,file = "riskScore.Validation.PRAD.GSE107299.RDS")

#### PRAD => GSE21034####
PRAD.Feature <- Feature.Filter %>% filter(Cancer == "PRAD")
Exp.Data <- readRDS("../AnimoAcidSensing/Validation.Survival/Taylor_eSet.RDS")
Sur.Data <- Exp.Data@phenoData@data
Exp.Data <- log2(Exp.Data@assayData$exprs+1) %>% data.frame(check.names = F) %>%
  rownames_to_column("Gene_Symbol") %>%
  filter(Gene_Symbol %in% PRAD.Feature$Feature) %>%
  data.frame(check.names = F) %>% remove_rownames() %>%
  column_to_rownames("Gene_Symbol")
Exp.Data <- Exp.Data %>% t() %>%
  data.frame(check.names = F) %>%
  rownames_to_column("Accession") %>%
  merge(Sur.Data %>% rownames_to_column("Accession"),by="Accession") %>%
  mutate(OS=bcr_status,
         OS.time=time_to_bcr) %>%
  filter(sample_type == "Primary")
DATA <- Cal_function(Exp.Data,PRAD.Feature$Feature,"PRAD.GSE21034",c(12,36,60))
saveRDS(DATA,file = "riskScore.Validation.PRAD.GSE21034.RDS")

#### PRAD => DKFZ ####
PRAD.Feature <- Feature.Filter %>% filter(Cancer == "PRAD")
Exp.Data <- readRDS("../AnimoAcidSensing/Validation.Survival/DKFZ_eSet.RDS")
Sur.Data <- Exp.Data@phenoData@data
Exp.Data <- log2(Exp.Data@assayData$exprs+1) %>% data.frame(check.names = F) %>%
  rownames_to_column("Gene_Symbol") %>%
  filter(Gene_Symbol %in% PRAD.Feature$Feature) %>%
  data.frame(check.names = F) %>% remove_rownames() %>%
  column_to_rownames("Gene_Symbol") %>%
  t() %>%
  data.frame(check.names = F) %>%
  rownames_to_column("Accession") %>%
  merge(Sur.Data %>% rownames_to_column("Accession"),by="Accession") %>%
  mutate(OS=bcr_status,
         OS.time=time_to_bcr) %>%
  filter(sample_type == "Primary")
DATA <- Cal_function(Exp.Data,PRAD.Feature$Feature,"PRAD.DKFZ",c(12,36,60))
saveRDS(DATA,file = "riskScore.Validation.PRAD.DKFZ.RDS")

#### PRAD => GSE54460 ####
PRAD.Feature <- Feature.Filter %>% filter(Cancer == "PRAD")
Exp.Data <- readRDS("../AnimoAcidSensing/Validation.Survival/GSE54460_eSet.RDS")
Sur.Data <- Exp.Data@phenoData@data
Exp.Data <- Exp.Data@assayData$exprs %>% data.frame(check.names = F) %>%
  rownames_to_column("Gene_Symbol") %>%
  filter(Gene_Symbol %in% PRAD.Feature$Feature) %>%
  data.frame(check.names = F) %>% remove_rownames() %>%
  column_to_rownames("Gene_Symbol")
Exp.Data <- Exp.Data %>%
  t() %>%
  data.frame(check.names = F) %>%
  rownames_to_column("Accession") %>%
  merge(Sur.Data %>% rownames_to_column("Accession"),by="Accession") %>%
  mutate(OS=bcr_status,
         OS.time=time_to_bcr) %>%
  filter(sample_type == "Primary")
DATA <- Cal_function(Exp.Data,PRAD.Feature$Feature,"PRAD.GSE54460",c(12,36,60))
saveRDS(DATA,file = "riskScore.Validation.PRAD.GSE54460.RDS")
#### PRAD => GSE70768 ####
PRAD.Feature <- Feature.Filter %>% filter(Cancer == "PRAD")
Exp.Data <- readRDS("../AnimoAcidSensing/Validation.Survival/Cambridge_eSet.RDS")
Sur.Data <- Exp.Data@phenoData@data
Exp.Data <- Exp.Data@assayData$exprs %>% data.frame(check.names = F) %>%
  rownames_to_column("Gene_Symbol") %>%
  filter(Gene_Symbol %in% PRAD.Feature$Feature) %>%
  data.frame(check.names = F) %>% remove_rownames() %>%
  column_to_rownames("Gene_Symbol") %>%
  t() %>%
  data.frame(check.names = F) %>%
  rownames_to_column("Accession") %>%
  merge(Sur.Data %>% rownames_to_column("Accession"),by="Accession") %>%
  mutate(OS=bcr_status,
         OS.time=time_to_bcr) %>%
  filter(sample_type == "Primary")
DATA <- Cal_function(Exp.Data,PRAD.Feature$Feature,"PRAD.GSE70768",c(12,36,60))
saveRDS(DATA,file = "riskScore.Validation.PRAD.GSE70768.RDS")
#### PRAD => GSE70769 ####
PRAD.Feature <- Feature.Filter %>% filter(Cancer == "PRAD")
Exp.Data <- readRDS("../AnimoAcidSensing/Validation.Survival/Stockholm_eSet.RDS")
Sur.Data <- Exp.Data@phenoData@data
Exp.Data <- Exp.Data@assayData$exprs %>% data.frame(check.names = F) %>%
  rownames_to_column("Gene_Symbol") %>%
  filter(Gene_Symbol %in% PRAD.Feature$Feature) %>%
  data.frame(check.names = F) %>% remove_rownames() %>%
  column_to_rownames("Gene_Symbol") %>%
  t() %>%
  data.frame(check.names = F) %>%
  rownames_to_column("Accession") %>%
  merge(Sur.Data %>% rownames_to_column("Accession"),by="Accession") %>%
  mutate(OS=bcr_status,
         OS.time=time_to_bcr) %>%
  filter(sample_type == "Primary")
DATA <- Cal_function(Exp.Data,PRAD.Feature$Feature,"PRAD.GSE70769",c(12,36,60))
saveRDS(DATA,file = "riskScore.Validation.PRAD.GSE70769.RDS")
#### PRAD => GSE70767 ####
PRAD.Feature <- Feature.Filter %>% filter(Cancer == "PRAD")
Exp.Data <- readRDS("../AnimoAcidSensing/Validation.Survival/CancerMap_eSet.RDS")
Sur.Data <- Exp.Data@phenoData@data
Exp.Data <- Exp.Data@assayData$exprs %>% data.frame(check.names = F) %>%
  rownames_to_column("Gene_Symbol") %>%
  filter(Gene_Symbol %in% PRAD.Feature$Feature) %>%
  data.frame(check.names = F) %>% remove_rownames() %>%
  column_to_rownames("Gene_Symbol") %>%
  t() %>%
  data.frame(check.names = F) %>%
  rownames_to_column("Accession") %>%
  merge(Sur.Data %>% rownames_to_column("Accession"),by="Accession") %>%
  mutate(OS=bcr_status,
         OS.time=time_to_bcr) %>%
  filter(str_detect(sample_type,"Tumor"))
DATA <- Cal_function(Exp.Data,PRAD.Feature$Feature,"PRAD.GSE70767",c(12,36,60))
saveRDS(DATA,file = "riskScore.Validation.PRAD.GSE70767.RDS")
#### PRAD => E-MTAB-6128 ####
PRAD.Feature <- Feature.Filter %>% filter(Cancer == "PRAD")
Exp.Data <- readRDS("../AnimoAcidSensing/Validation.Survival/CIT_eSet.RDS")
Sur.Data <- Exp.Data@phenoData@data
Exp.Data <- Exp.Data@assayData$exprs %>% data.frame(check.names = F) %>%
  rownames_to_column("Gene_Symbol") %>%
  filter(Gene_Symbol %in% PRAD.Feature$Feature) %>%
  data.frame(check.names = F) %>% remove_rownames() %>%
  column_to_rownames("Gene_Symbol") %>%
  t() %>%
  data.frame(check.names = F) %>%
  rownames_to_column("Accession") %>%
  merge(Sur.Data %>% rownames_to_column("Accession"),by="Accession") %>%
  mutate(OS=bcr_status,
         OS.time=time_to_bcr) %>%
  filter(str_detect(sample_type,"Primary"))
DATA <- Cal_function(Exp.Data,PRAD.Feature$Feature,"PRAD.E-MTAB-6128",c(12,36,60))
saveRDS(DATA,file = "riskScore.Validation.PRAD.E-MTAB-6128.RDS")
#### READ => GSE87211 ####
READ.Feature <- Feature.Filter %>% filter(Cancer == "READ")
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE87211.gene_expression_array.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  filter(Gene_Symbol %in% READ.Feature$SYMBOL) %>%
  column_to_rownames("Gene_Symbol") %>%
  t() %>% data.frame(check.names = F)
#Exp.Data <- log2(Exp.Data+1)
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE87211.clinical.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  mutate(SampleID=str_remove_all(Patient,".*, ")) %>%
  filter(Sample_Type == "Tumor") %>%
  dplyr::select(Accession,OS_Time,OS_Status) %>%
  magrittr::set_colnames(c("Accession","OS.time","OS")) %>%
  merge(Exp.Data %>% rownames_to_column("Accession"),by="Accession")

DATA <- Cal_function(Exp.Data,READ.Feature$SYMBOL,"READ.GSE87211",c(365,1080,1800))
saveRDS(DATA,file = "riskScore.Validation.READ.GSE87211.RDS")

#### READ => GSE103479 ####
READ.Feature <- Feature.Filter %>% filter(Cancer == "READ")
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE103479.gene_expression_array.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  filter(Gene_Symbol %in% READ.Feature$SYMBOL) %>%
  column_to_rownames("Gene_Symbol") %>%
  t() %>% data.frame(check.names = F)
#Exp.Data <- log2(Exp.Data+1)
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE103479.clinical.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  mutate(SampleID=str_remove_all(Patient,".*, ")) %>%
  filter(Sample_Type == "Tumor") %>%
  dplyr::select(Accession,OS_Time,OS_Status) %>%
  magrittr::set_colnames(c("Accession","OS.time","OS")) %>%
  merge(Exp.Data %>% rownames_to_column("Accession"),by="Accession")

DATA <- Cal_function(Exp.Data,READ.Feature$SYMBOL,"READ.GSE103479",c(365,1080,1800))
saveRDS(DATA,file = "riskScore.Validation.READ.GSE103479.RDS")
#### READ => GSE12945 ####
READ.Feature <- Feature.Filter %>% filter(Cancer == "READ")
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE12945.gene_expression_array.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  filter(Gene_Symbol %in% READ.Feature$SYMBOL) %>%
  column_to_rownames("Gene_Symbol") %>%
  t() %>% data.frame(check.names = F)
#Exp.Data <- log2(Exp.Data+1)
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE12945.clinical.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  mutate(SampleID=str_remove_all(Patient,".*, ")) %>%
  filter(Sample_Type == "Tumor") %>%
  dplyr::select(Accession,OS_Time,OS_Status) %>%
  magrittr::set_colnames(c("Accession","OS.time","OS")) %>%
  merge(Exp.Data %>% rownames_to_column("Accession"),by="Accession")

DATA <- Cal_function(Exp.Data,READ.Feature$SYMBOL,"READ.GSE12945",c(365,1080,1800))
saveRDS(DATA,file = "riskScore.Validation.READ.GSE12945.RDS")
#### SKCM => GSE22153 ####
SKCM.Feature <- Feature.Filter %>% filter(Cancer == "SKCM")
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE22153.gene_expression_array.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  filter(Gene_Symbol %in% SKCM.Feature$SYMBOL) %>%
  column_to_rownames("Gene_Symbol") %>%
  t() %>% data.frame(check.names = F)
#Exp.Data <- log2(Exp.Data+1)
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE22153.clinical.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  mutate(SampleID=str_remove_all(Patient,".*, ")) %>%
  filter(Sample_Type == "Tumor") %>%
  dplyr::select(Accession,OS_Time,OS_Status) %>%
  magrittr::set_colnames(c("Accession","OS.time","OS")) %>%
  merge(Exp.Data %>% rownames_to_column("Accession"),by="Accession")

DATA <- Cal_function(Exp.Data,SKCM.Feature$SYMBOL,"SKCM.GSE22153",c(365,1080,1800))
saveRDS(DATA,file = "riskScore.Validation.SKCM.GSE22153.RDS")
#### SKCM => GSE65904 ####
SKCM.Feature <- Feature.Filter %>% filter(Cancer == "SKCM")
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE65904_expression.txt",header = T) %>%
  data.frame(check.names = F) %>%
  filter(Symbol %in% SKCM.Feature$SYMBOL) %>%
  column_to_rownames("Symbol") %>%
  t() %>% data.frame(check.names = F)
#Exp.Data <- log2(Exp.Data+1)
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE65904_clinical.txt",header = T) %>%
  data.frame(check.names = F) %>%
  mutate(OS=`disease specific survival (1=death, 0=alive)`,
         OS.time=`disease specific survival in days`) %>%
  #mutate(SampleID=str_remove_all(Patient,".*, ")) %>%
  #filter(Sample_Type == "Tumor") %>%
  dplyr::select(accession,OS.time,OS) %>%
  #magrittr::set_colnames(c("accession","OS.time","OS")) %>%
  merge(Exp.Data %>% rownames_to_column("accession"),by="accession")

DATA <- Cal_function(Exp.Data,SKCM.Feature$SYMBOL,"SKCM.GSE65904",c(365,1080,1800))
saveRDS(DATA,file = "riskScore.Validation.SKCM.GSE65904.RDS")
#### SKCM => GSE19234 ####
SKCM.Feature <- Feature.Filter %>% filter(Cancer == "SKCM")
library(data.table)
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE19234.gene_expression_array.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  filter(Gene_Symbol %in% SKCM.Feature$SYMBOL) %>%
  column_to_rownames("Gene_Symbol") %>%
  t() %>% data.frame(check.names = F)
#Exp.Data <- log2(Exp.Data+1)
Exp.Data <- fread("../AnimoAcidSensing/Validation.Survival/GSE19234.clinical.tsv.gz",header = T) %>%
  data.frame(check.names = F) %>%
  mutate(SampleID=str_remove_all(Patient,".*, ")) %>%
  filter(Sample_Type == "Tumor") %>%
  dplyr::select(Accession,OS_Time,OS_Status) %>%
  magrittr::set_colnames(c("Accession","OS.time","OS")) %>%
  merge(Exp.Data %>% rownames_to_column("Accession"),by="Accession")

DATA <- Cal_function(Exp.Data,SKCM.Feature$SYMBOL,"SKCM.GSE19234",c(365,1080,1800))
saveRDS(DATA,file = "riskScore.Validation.SKCM.GSE19234.RDS")
##### Figure #####
files <- list.files(".",pattern = "^riskScore")
files <- files[str_detect(files,"RDS$")]
Cindex <- data.frame()
AUC.Data <- data.frame()
KM.P <- data.frame()
for (file in files) {
  project <- (file %>% str_split("\\.") %>% unlist())[3]
  code <- (file %>% str_split("\\.") %>% unlist())[4]
  print(code)
  data <- readRDS(file)
  AUC.Data <- data.frame(year1=round(data$timeROC$AUC,2)[1],
                         year3=round(data$timeROC$AUC,2)[2],
                         year5=round(data$timeROC$AUC,2)[3]) %>%
    mutate(Cancer=project,
           Code=code) %>%
    rbind.data.frame(AUC.Data)
  
  Cindex <- data.frame(C.Index = data$C.index) %>%
    mutate(Cancer=project,
           Code=code) %>%
    rbind.data.frame(Cindex)
  
  KM.P <- data.frame(KM.Pvalue = data$KM.Pvalue) %>%
    mutate(Cancer=project,
           Code=code) %>%
    rbind.data.frame(KM.P)
  
  dat = data.frame(fpr = as.numeric(data$timeROC$FP),
                   tpr = as.numeric(data$timeROC$TP),
                   time = rep(as.factor(data$timeROC$times),
                              each = nrow(data$timeROC$TP)))
  library(ggplot2)
  p <- ggplot() + 
    geom_line(data = dat,aes(x = fpr, y = tpr,color = time),size = 1) + 
    scale_color_manual(name = NULL,values = c("#E64B35FF", "#3C5488FF", "#00A087FF"),
                       labels = paste0("AUC of ",c(1,3,5),"-y survival: ",
                                       format(round(data$timeROC$AUC,2),nsmall = 2)))+
    #ggsci::scale_color_npg()+
    geom_line(aes(x=c(0,1),y=c(0,1)),color = "grey",linetype=2)+
    ggthemes::theme_few()+
    theme(#legend.position = "top",
      axis.text = element_text(color = "black"),
      axis.text.x = element_text(size=9,angle = 90,vjust = 0.3,hjust = 1),
      axis.text.y = element_text(size=9),
      #axis.ticks.x = element_blank(),
      axis.title = element_text(size=12),
      panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
    labs(x = "1 - Specificity",
         y = "Sensitivity")+
    coord_fixed()
  ggsave(p,filename = file.path("riskScore.timeROC",paste("ROC.",project,"_",code,".pdf",sep = "")),
         height = 3,width = 5)
}
library(RColorBrewer)
Cancer.Col <- c(brewer.pal(9,"Set1"),
                brewer.pal(8,"Set2")[1:7],
                brewer.pal(12,"Set3")[1:3],
                brewer.pal(8,"Dark2")[6],
                brewer.pal(12,"Set3")[4:12],
                brewer.pal(8,"Dark2")[1:4])
Cindex <- Cindex %>% arrange(desc(Cancer),Code) %>%
  mutate(Code=factor(Code,levels = .$Code))
middata <- Cindex %>% group_by(Cancer) %>% 
  summarise(Count=n()) %>%
  arrange(desc(Cancer))
RECT<-data.frame(xmax=accumulate(middata$Count,sum)+0.5,
           xmin=accumulate(middata$Count,sum)-middata$Count+0.5,
           Cancer=middata$Cancer)
write.csv(Cindex,file = "riskScore.Validation.Cindex.csv",row.names = F)
P1 <- ggplot(data = Cindex,aes(C.Index,Code,fill=Cancer))+
  geom_segment(aes(x=0.5,xend=C.Index,
                   y=Code,yend=Code),
               color="black")+
  geom_point(size=3,shape=21,color="black")+
  ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9,angle = 90,vjust = 0.3,hjust = 1),
        axis.text.y = element_text(size=9),
        #axis.ticks.x = element_blank(),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x = "C index",
       y = "")+
  annotate("rect",
           xmin=-Inf,xmax=Inf,
           ymin=RECT$xmin,
           ymax=RECT$xmax,
           fill=rev(Cancer.Col[1:16]),
           alpha = 0.4
  )+
  scale_fill_manual(values = Cancer.Col)+
  scale_color_manual(values = Cancer.Col)+
  geom_vline(xintercept = seq(0.6,1,0.2),linetype=2)
  

P1.R <-  ggplot(Cindex,aes(1,Code,fill=Cancer))+
  geom_col(width = 1)+
  theme_void()+
  #theme(legend.position = "bottom")+
  scale_fill_manual(values = Cancer.Col)
library(tidyverse)
AUC.Data <- AUC.Data %>% gather(Year,AUROC,year1:year5) #%>%
  #arrange(desc(Cancer),Code) #%>%
  #mutate(Code=factor(Code,levels = levels(Cindex$Code)))
write.csv(AUC.Data,file = "riskScore.Validation.AUROC.csv",row.names = F)

P2 <- ggplot(AUC.Data,aes(AUROC,Code,fill=Cancer))+
  geom_segment(aes(x=0.2,
                   xend=AUROC,
                   y=Code,
                   yend=Code),color="black")+
  geom_point(size=3,shape=21,color="black")+
  ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9,angle = 90,vjust = 0.3,hjust = 1),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x = "AUROC",
       y = "")+
  facet_wrap(~Year,ncol = 3)+
  scale_fill_manual(values = Cancer.Col)+
  scale_color_manual(values = Cancer.Col)+
  geom_vline(xintercept = seq(0.6,1,0.2),linetype=2)

KM.P <- KM.P %>% #gather(Year,AUROC,year1:year5) %>%
  arrange(desc(Cancer),Code) %>%
  mutate(Code=factor(Code,levels = levels(Cindex$Code)))

P3 <- ggplot()+
  geom_segment(data=KM.P,aes(x=0,xend=-log10(KM.Pvalue),y=Code,yend=Code),
               size=1,color="black")+
  geom_point(data=KM.P,
             aes(x=-log10(KM.Pvalue),y=Code,
                 fill=-log10(KM.Pvalue),size=-log10(KM.Pvalue)),
             color="white",stroke=1,shape=21)+
  geom_point(data=KM.P %>% filter(KM.Pvalue <= 0.05),
             aes(x=-log10(KM.Pvalue),y=Code,
                 fill=-log10(KM.Pvalue),size=-log10(KM.Pvalue)),
             color="black",stroke=1,shape=21)+
  scale_size_continuous(range = c(1,5))+
  ggthemes::theme_few()+
  theme(#legend.position = "bottom",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9,angle = 90,vjust = 0.3,hjust = 1),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x = "-log10(P)",
       y = "")+
  annotate("rect",
           xmin=-Inf,xmax=Inf,
           ymin=RECT$xmin,
           ymax=RECT$xmax,
           fill=rev(Cancer.Col[1:16]),
           alpha = 0.4
  )+
  scale_x_continuous(limits = c(0,14))+
  geom_vline(xintercept = c(-log10(0.05)),linetype=2)+
  #scale_color_manual(values = Cancer.Col)+
  scale_fill_gradientn(colours = brewer.pal(11,"RdBu")[5:1])

P.F <- P1 %>% aplot::insert_left(P1.R,width = 0.1) %>%
  aplot::insert_right(P2,width = 3) %>%
  aplot::insert_right(P3,width = 1)
ggsave(P.F,filename = "riskScore.Validation.pdf",height = 8,width = 12)
ggsave(P.F,filename = "riskScore.Validation.Legend.pdf",height = 8,width = 12)


######## Risk score Risk Group => Survival ####
setwd("K:/TCGA/Anlysis/BCAA")
library(tidyverse)
RiskScore.Cancer <- read.csv(file = "RiskGroup.COX.TCGA.RiskScore.csv")
Risk.Cutoff <- data.frame()
Survival.data <- data.frame()
for (project in unique(RiskScore.Cancer$Cancer)) {
  mid.Risk <- RiskScore.Cancer %>% filter(Cancer == project)
  
  res.cut <- surv_cutpoint(mid.Risk,
                           time = paste("OS",".time",sep = ""),
                           event = "OS",
                           variables = "total_risk_score")
  #res.cat <- surv_categorize(res.cut)
  Risk.Cutoff <- data.frame(CutPoint=res.cut$cutpoint[1]) %>%
    mutate(Cancer=project) %>%
    rbind.data.frame(Risk.Cutoff)
  theme.sur <- ggthemes::theme_few()+
    theme(legend.position = "none",
          axis.text = element_text(color = "black"),
          axis.text.x = element_text(size=15),
          axis.text.y = element_text(size=15),
          axis.title = element_text(size=18),
          panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))
  
  dat <- surv_categorize(res.cut) #%>% data.frame() %>% mutate(Sample=SNHGs.ICI.data$Sample.ID)
  fit <- survfit(Surv(OS.time, OS) ~ total_risk_score,
                 data = dat)
  library(ggthemes)
  p<-ggsurvplot(fit,pval =TRUE, data = dat, 
                surv.median.line = "hv",
                legend.title = "Risk",
                conf.int.style = "step",
                xlab = "Time in days",
                #break.time.by = 500,
                risk.table = "abs_pct",
                #risk.table.y.text.col = T,
                #risk.table.y.text = FALSE,
                legend.labs = c("High", "Low"),
                #pval = TRUE,
                conf.int = TRUE,
                palette = "Set1",
                #ggtheme = theme.sur,
                risk.table.y.text.col = T,
                risk.table.y.text = F,
                ggtheme = theme.sur)
  print(p)
  #dev.off()
  library(export)
  graph2pdf(file=paste("RiskGroup/OS.",project,".RiskScore.BestCutoff.KM.pdf",sep = ''),height = 5,width = 4)
  
  index="OS"
  diff=survdiff(as.formula(paste("Surv(",index,".time,",index,") ~ ","total_risk_score",sep = "")),data = dat)
  pValue=1-pchisq(diff$chisq,df=1)
  
  cox <- coxph(as.formula(paste("Surv(OS.time, OS) ~ ","total_risk_score",sep = "")), data = mid.Risk)
  coxSummary = summary(cox)
  
  Survival.data=rbind(Survival.data,
                      data.frame(ENSEMBL="RiskScore",
                                 KM.pvalue=pValue,
                                 HR=coxSummary$conf.int[,"exp(coef)"],
                                 HR.95L=coxSummary$conf.int[,"lower .95"],
                                 HR.95H=coxSummary$conf.int[,"upper .95"],
                                 HR.pvalue=coxSummary$coefficients[,"Pr(>|z|)"],
                                 Cancer=project))
}
write.csv(Survival.data,file = "RiskGroup.COX.TCGA.RiskScore.Survival.csv",row.names = F)
write.csv(Risk.Cutoff,file = "RiskGroup.COX.TCGA.RiskScore.CutOff.csv",row.names = F)

######## Risk score => TIDE ####
setwd("K:/TCGA/Anlysis/BCAA")
library(tidyverse)
RiskScore.Cancer <- read.csv(file = "RiskGroup.COX.TCGA.RiskScore.csv")

PhenoType <- read.csv("../../mRNA.Phenotype.csv")
PhenoType <- PhenoType %>% filter(Sample.Type == "Primary Tumor" | 
                                    Sample.Type == "Solid Tissue Normal" | 
                                    Sample.Type == "Primary Blood Derived Cancer - Peripheral Blood") %>% 
  dplyr::select(Project.ID,Sample.ID,Sample.Type) %>%
  mutate(Type=ifelse(Sample.Type=="Solid Tissue Normal","Normal","Tumor")) %>%
  mutate(Cancer=str_remove(Project.ID,"TCGA-")) %>%
  filter(Type=="Tumor")

TIDE.data <- read.csv("../../TIDE/TIDE.csv")
Risk.Cutoff <- read.csv(file = "RiskGroup.COX.TCGA.RiskScore.CutOff.csv")

RiskScore.TIDE <- read.csv(file = "RiskGroup.COX.TCGA.RiskScore.csv") %>%
  filter(SampleID %in% PhenoType$Sample.ID) %>%
  mutate(Sample.ID = substring(SampleID,1,12)) %>%
  merge(TIDE.data,by.x="Sample.ID",by.y="SampleID") %>%
  #merge(Risk.Cutoff,by.x = "Cancer.x",by.y="Cancer") %>%
  group_by(Cancer.x) %>%
  mutate(cutpoint=median(total_risk_score)) %>%
  mutate(Group = if_else(total_risk_score >= cutpoint,"High","Low"))

write.csv(RiskScore.TIDE,file = "RiskGroup.COX.TCGA.RiskScore.TIDE.csv",row.names = F)

library(ggpubr)
data2 <- compare_means(formula = Exclusion~Group,data = RiskScore.TIDE,method = "wilcox.test",
              group.by="Cancer.x")
write.csv(data2,file = "RiskGroup.COX.TCGA.RiskScore.TIDE.Exclusion.Wilcox.csv")

data3 <- compare_means(formula = Dysfunction~Group,data = RiskScore.TIDE,method = "wilcox.test",
                       group.by="Cancer.x")
write.csv(data3,file = "RiskGroup.COX.TCGA.RiskScore.TIDE.Dysfunction.Wilcox.csv")

library(RColorBrewer)
p1<-ggplot()+
  geom_boxplot(data=RiskScore.TIDE,
               aes(x=Exclusion,y=reorder(Cancer.x,Exclusion), #Cancer,#
                   fill=Group))+
  geom_text(data=data2,aes(x=0.3,y=Cancer.x,label=p.signif))+
  #stat_ellipse(aes(fill=Group),type = "norm", geom ="polygon",alpha=0.2,color=NA)+
  #geom_beeswarm(aes(color=SampleType),size=1,dodge.width=0.75,priority = "ascending")+
  labs(x="",y="")+
  #facet_wrap(~Phenotype,scale="free_y",ncol = 6)+
  ggthemes::theme_few()+
  theme(#legend.position = "top",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9,angle = 90,vjust = 0.3,hjust = 1),
        axis.text.y = element_text(size=9),
        #axis.ticks.x = element_blank(),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  scale_color_manual("",values = brewer.pal(9,"Set1")[c(1:2)])+
  scale_fill_manual("",values = brewer.pal(9,"Set1")[c(1:2)])+
  scale_x_continuous(limits = c(-0.3,0.4))+
  ggtitle("Exclusion")+
  theme(plot.title = element_text(hjust = 0.5))

p2<-ggplot()+
  geom_boxplot(data=RiskScore.TIDE,
               aes(x=Dysfunction,y=reorder(Cancer.x,Dysfunction), #Cancer,#
                   fill=Group))+
  geom_text(data=data3,aes(x=0.5,y=Cancer.x,label=p.signif))+
  #stat_ellipse(aes(fill=Group),type = "norm", geom ="polygon",alpha=0.2,color=NA)+
  #geom_beeswarm(aes(color=SampleType),size=1,dodge.width=0.75,priority = "ascending")+
  labs(x="",y="")+
  #facet_wrap(~Phenotype,scale="free_y",ncol = 6)+
  ggthemes::theme_few()+
  theme(legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9,angle = 90,vjust = 0.3,hjust = 1),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  scale_color_manual("",values = brewer.pal(9,"Set1")[c(1:2)])+
  scale_fill_manual("",values = brewer.pal(9,"Set1")[c(1:2)])+
  scale_x_continuous(limits = c(-0.4,0.6))+
  ggtitle("Dysfunction")+
  theme(plot.title = element_text(hjust = 0.5))

p12 <- p1 %>% aplot::insert_right(p2,width = 1)
ggsave(p12,filename = "RiskGroup.COX.TCGA.RiskScore.TIDE.pdf",height = 4,width = 5)

######## Risk score => IMMUNECELL ####
setwd("K:/TCGA/Anlysis/BCAA")
library(tidyverse)
RiskScore.Cancer <- read.csv(file = "RiskGroup.COX.TCGA.RiskScore.csv")
PhenoType <- read.csv("../../mRNA.Phenotype.csv")
PhenoType <- PhenoType %>% filter(Sample.Type == "Primary Tumor" | 
                                    Sample.Type == "Solid Tissue Normal" | 
                                    Sample.Type == "Primary Blood Derived Cancer - Peripheral Blood") %>% 
  dplyr::select(Project.ID,Sample.ID,Sample.Type) %>%
  mutate(Type=ifelse(Sample.Type=="Solid Tissue Normal","Normal","Tumor")) %>%
  mutate(Cancer=str_remove(Project.ID,"TCGA-")) %>%
  filter(Type=="Tumor")

TIMER2 <- read.csv("../../TIMER2/infiltration_estimation_for_tcga.csv",
                   check.names = F) %>%
  remove_rownames() %>% column_to_rownames("cell_type") %>%
  dplyr::select(ends_with("CIBERSORT")) %>%
  rownames_to_column("Sample.ID") %>%
  gather(CellType,Abundance,-Sample.ID) %>%
  mutate(CellType=str_remove_all(CellType,"_CIBERSORT"))

Risk.Cutoff <- read.csv(file = "RiskGroup.COX.TCGA.RiskScore.CutOff.csv")

RiskScore.TIMER2 <- read.csv(file = "RiskGroup.COX.TCGA.RiskScore.csv") %>%
  filter(SampleID %in% PhenoType$Sample.ID) %>%
  mutate(Sample.ID = substring(SampleID,1,15)) %>%
  merge(TIMER2,by="Sample.ID") %>%
  #merge(Risk.Cutoff,by = "Cancer") %>%
  group_by(Cancer,CellType) %>%
  mutate(cutpoint=median(total_risk_score)) %>%
  mutate(Group = if_else(total_risk_score >= cutpoint,"High","Low"))

write.csv(RiskScore.TIMER2,file = "RiskGroup.COX.TCGA.RiskScore.CIBERSORT.csv",row.names = F)

data2 <- compare_means(formula = Abundance~Group,data = RiskScore.TIMER2,method = "wilcox.test",
              group.by=c("Cancer","CellType"))

write.csv(data2,file = "RiskGroup.COX.TCGA.RiskScore.CIBERSORT.Wilcox.csv",row.names = F)

RiskScore.TIMER2 <- read.csv("RiskGroup.COX.TCGA.RiskScore.CIBERSORT.csv")
data2 <- read.csv("RiskGroup.COX.TCGA.RiskScore.CIBERSORT.Wilcox.csv")

RiskScore.TIMER2.FigureD <- RiskScore.TIMER2 %>%
  group_by(Cancer,CellType) %>%
  mutate(cutpoint=median(total_risk_score)) %>%
  mutate(Group = if_else(total_risk_score >= cutpoint,"High","Low")) %>%
  group_by(Cancer,CellType,Group) %>%
  summarise(MeanValue=mean(Abundance)) %>%
  ungroup() %>% data.frame(check.names = F) %>%
  merge(data2 %>% dplyr::select(Cancer,CellType,p.signif),
        by=c("Cancer","CellType")) %>%
  mutate(p.signif=if_else(p.signif=="ns","",p.signif)) %>%
  #filter(Attribute %in% c("IPS","TMB","Neoantigen","MSI","HRD","MATH","Ploidy","Purity","RNAss")) %>%
  spread(Group,MeanValue) %>%
  mutate(Enrich=if_else(p.signif == "","n.s.",if_else(High>Low,"High","Low"))) %>%
  gather(Group,Score,High:Low) %>%
  arrange(Cancer,CellType) %>%
  group_by(Cancer,CellType) %>%
  mutate(Enrich2=if_else(Enrich=="n.s.","n.s.",if_else(Enrich==Group,"Higher","Lower"))) %>%
  mutate(CellType=factor(CellType,levels = unique(sort(CellType)))) %>%
  mutate(Cancer=factor(Cancer,levels = rev(sort(unique(.$Cancer))))) %>%
  mutate(p.signif=factor(p.signif,levels = c("*","**","***","****",""))) %>%
  mutate(Xaxis = paste(CellType,Group,sep = "_"))
write.csv(RiskScore.TIMER2.FigureD,file = "RiskGroup.COX.TCGA.RiskScore.CIBERSORT.FigureData.csv",row.names = F)

p1 <- ggplot(RiskScore.TIMER2.FigureD,aes(x=Xaxis,y=Cancer))+
  geom_point(aes(size=Score,color=p.signif,shape=Enrich2),
             stroke=1.2)+
  scale_size_continuous("Mean",range = c(2,4))+
  #facet_wrap(~Attribute,ncol = 11)+
  labs(x="",y="")+
  #facet_wrap(~Attribute,scale="free_x",ncol = 15)+
  ggthemes::theme_few()+
  theme(legend.position = "top",
    axis.text = element_text(color = "black"),
    axis.text.x = element_blank(), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    axis.ticks.x = element_blank(),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  scale_color_manual("P value",values = c(brewer.pal(9,"Set1")[1:4],"#F7F7F7"))+
  scale_shape_manual("Score",values = c(1,3,4))+
  geom_vline(xintercept = seq(0.5,44,2),linetype=2)

p2 <- ggplot(RiskScore.TIMER2.FigureD %>% data.frame() %>%
         dplyr::select(Xaxis,Group) %>%
         unique(),
       aes(x=Xaxis,y=1,fill=Group))+
  geom_col(width = 1)+
  theme_void()+
  scale_fill_manual("Risk score",values = c("#fe5f55","#004d61"))


p3 <- ggplot(RiskScore.TIMER2.FigureD %>% data.frame() %>%
               dplyr::select(Xaxis,CellType) %>%
               unique(),
             aes(x=Xaxis,y=1,fill=CellType))+
  geom_col(width = 1)+
  theme_void()+
  scale_fill_manual("Cell type",values = Cancer.Col)

p123 <- p1 %>% aplot::insert_top(p2,height = 0.05) %>%
  aplot::insert_top(p3,height = 0.05)
ggsave(p123,filename = "RiskGroup.COX.TCGA.RiskScore.CIBERSORT.FigureData.pdf",height = 4,width = 12)
ggsave(p123,filename = "RiskGroup.COX.TCGA.RiskScore.CIBERSORT.FigureData.Legend.pdf",height = 16,width = 12)


for (project in unique(RiskScore.TIMER2$Cancer)) {
  middata <- RiskScore.TIMER2 %>% filter(Cancer==project)
  Label.Data <- data2 %>% filter(Cancer==project)
  p <- ggplot()+
    geom_boxplot(data=middata,
                 aes(x=Abundance,y=CellType, #Cancer,#
                     fill=Group))+
    geom_text(data=Label.Data,
              aes(x=max(middata$Abundance)+0.03,y=CellType,
                  label=p.signif))+
    #stat_ellipse(aes(fill=Group),type = "norm", geom ="polygon",alpha=0.2,color=NA)+
    #geom_beeswarm(aes(color=SampleType),size=1,dodge.width=0.75,priority = "ascending")+
    labs(x="",y="")+
    #facet_wrap(~Phenotype,scale="free_y",ncol = 6)+
    ggthemes::theme_few()+
    theme(#legend.position = "top",
      axis.text = element_text(color = "black"),
      axis.text.x = element_text(size=9,angle = 90,vjust = 0.3,hjust = 1),
      axis.text.y = element_text(size=9),
      #axis.ticks.x = element_blank(),
      axis.title = element_text(size=12),
      panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
    scale_color_manual("Risk score",values = brewer.pal(9,"Set1")[c(1:2)])+
    scale_fill_manual("Risk score",values = brewer.pal(9,"Set1")[c(1:2)])+
    scale_x_continuous(limits = c(0,max(middata$Abundance)+0.1))
  ggsave(p,filename = paste("RiskGroup/RiskGroup.COX.",project,".RiskScore.CIBERSORT.pdf",sep = ""),height = 5,width = 5)
}


  #scale_x_continuous(trans = "log2")
  #scale_x_continuous(limits = c(-0.3,0.4))+
  #ggtitle("Exclusion")+
  #theme(plot.title = element_text(hjust = 0.5))
######## Risk score => INDEX ####
setwd("K:/TCGA/Anlysis/BCAA")
library(tidyverse)
RiskScore.Cancer <- read.csv(file = "RiskGroup.COX.TCGA.RiskScore.csv") %>%
  mutate(ID=substring(SampleID,1,15))

PhenoType <- read.csv("../../mRNA.Phenotype.csv")
PhenoType <- PhenoType %>% filter(Sample.Type == "Primary Tumor" | 
                                    Sample.Type == "Solid Tissue Normal" | 
                                    Sample.Type == "Primary Blood Derived Cancer - Peripheral Blood") %>% 
  dplyr::select(Project.ID,Sample.ID,Sample.Type) %>%
  mutate(Type=ifelse(Sample.Type=="Solid Tissue Normal","Normal","Tumor")) %>%
  mutate(Cancer=str_remove(Project.ID,"TCGA-")) %>%
  filter(Type=="Tumor")

Risk.Index <- data.frame()

data <- read_xlsx("../../PanCancer.DiverseIndex.xlsx",sheet = 3)
data2<-data[,-1]
Risk.Index <- merge(RiskScore.Cancer,data2,by="ID") %>% unique() %>%
  data.frame() %>%
  filter(SampleID %in% PhenoType$Sample.ID) %>%
  gather(Attribute,Score,MHC:IPS) %>%
  dplyr::select(ID,SampleID,Attribute,Score,Cancer,total_risk_score) %>%
  rbind.data.frame(Risk.Index)

for (index in c(4,6,7,10,5,8,9,13)) {
  data <- read_xlsx("../../PanCancer.DiverseIndex.xlsx",sheet = index)
  data2<-data[,-1]
  Risk.Index <- merge(RiskScore.Cancer,data2,by.x="ID",by.y="SampleName") %>% 
    unique() %>%
    data.frame() %>%
    filter(SampleID %in% PhenoType$Sample.ID) %>%
    mutate(Attribute=colnames(data2)[1]) %>%
    mutate(Score=.[[colnames(data2)[1]]]) %>%
    dplyr::select(ID,SampleID,Attribute,Score,Cancer,total_risk_score) %>%
  rbind.data.frame(Risk.Index)
}

Risk.Index2 <- Risk.Index %>%
  group_by(Cancer,Attribute) %>%
  mutate(cutpoint=median(total_risk_score)) %>%
  mutate(Group = if_else(total_risk_score >= cutpoint,"High","Low"))

write.csv(Risk.Index,file = "RiskGroup.COX.TCGA.RiskScore.Index.csv",row.names = F)

data2 <- compare_means(formula = Score~Group,data = Risk.Index2,method = "wilcox.test",
                       group.by=c("Cancer","Attribute"))

write.csv(data2,file = "RiskGroup.COX.TCGA.RiskScore.Index.Wilcox.csv",row.names = F) 

#### Figure ####
Risk.Index <- read.csv("RiskGroup.COX.TCGA.RiskScore.Index.csv")
data2 <- read.csv("RiskGroup.COX.TCGA.RiskScore.Index.Wilcox.csv")
Risk.Index2 <- Risk.Index %>%
  group_by(Cancer,Attribute) %>%
  mutate(cutpoint=median(total_risk_score)) %>%
  mutate(Group = if_else(total_risk_score >= cutpoint,"High","Low")) %>%
  group_by(Cancer,Attribute,Group) %>%
  summarise(MeanValue=mean(Score)) %>%
  ungroup() %>% data.frame(check.names = F) %>%
  merge(data2 %>% dplyr::select(Cancer,Attribute,p.signif),
        by=c("Cancer","Attribute")) %>%
  mutate(p.signif=if_else(p.signif=="ns","",p.signif)) %>%
  filter(Attribute %in% c("IPS","TMB","Neoantigen","MSI","HRD","MATH","Ploidy","Purity","RNAss")) %>%
  spread(Group,MeanValue) %>%
  mutate(Enrich=if_else(p.signif == "","n.s.",if_else(High>Low,"High","Low"))) %>%
  gather(Group,Score,High:Low) %>%
  arrange(Cancer,Attribute) %>%
  group_by(Cancer,Attribute) %>%
  mutate(Enrich2=if_else(Enrich=="n.s.","n.s.",if_else(Enrich==Group,"Higher","Lower")))
  #mutate(Attribute=factor(Attribute,levels = c("IPS","TMB","Neoantigen","MSI","HRD","MATH","Ploidy","Purity","RNAss"))) %>%
  #mutate(Cancer=factor(Cancer,levels = rev(sort(unique(.$Cancer))))) %>%
  #mutate(p.signif=factor(p.signif,levels = c("*","**","***","****","")))


library(ggpubr)
RiskScore.TIDE <- read.csv(file = "RiskGroup.COX.TCGA.RiskScore.TIDE.csv") %>%
  gather(Attribute,Score,Exclusion:Dysfunction) %>%
  dplyr::select(Cancer.x,Attribute,Score,Group) %>%
  magrittr::set_colnames(c("Cancer","Attribute","Score","Group"))
data2 <- compare_means(formula = Score~Group,data = RiskScore.TIDE,method = "wilcox.test",
                       group.by=c("Cancer","Attribute"))

Risk.Index3 <- RiskScore.TIDE %>%
  #group_by(Cancer,Attribute) %>%
  #mutate(cutpoint=median(total_risk_score)) %>%
  #mutate(Group = if_else(total_risk_score >= cutpoint,"High","Low")) %>%
  group_by(Cancer,Attribute,Group) %>%
  summarise(MeanValue=mean(Score)) %>%
  ungroup() %>% data.frame(check.names = F) %>%
  merge(data2 %>% dplyr::select(Cancer,Attribute,p.signif),by=c("Cancer","Attribute")) %>%
  mutate(p.signif=if_else(p.signif=="ns","",p.signif)) %>%
  #filter(Attribute %in% c("IPS","TMB","Neoantigen","MSI","HRD","MATH","Ploidy","Purity","RNAss")) %>%
  spread(Group,MeanValue) %>%
  mutate(Enrich=if_else(p.signif == "","n.s.",if_else(High>Low,"High","Low"))) %>%
  gather(Group,Score,High:Low) %>%
  arrange(Cancer,Attribute) %>%
  group_by(Cancer,Attribute) %>%
  mutate(Enrich2=if_else(Enrich=="n.s.","n.s.",if_else(Enrich==Group,"Higher","Lower")))

Risk.Index4 <- rbind(Risk.Index2,Risk.Index3) %>%
  data.frame(check.names = F) %>%
  mutate(Attribute=factor(Attribute,levels = c("IPS","TMB","Exclusion","Dysfunction","Neoantigen","MSI","HRD","MATH","Ploidy","Purity","RNAss"))) %>%
  mutate(Cancer=factor(Cancer,levels = rev(sort(unique(.$Cancer))))) %>%
  mutate(p.signif=factor(p.signif,levels = c("*","**","***","****","")))
  

p <- ggplot(Risk.Index4,aes(x=Group,y=Cancer))+
  geom_point(aes(size=Score,color=p.signif,shape=Enrich2))+
  scale_size_continuous("Mean",range = c(2,5))+
  facet_wrap(~Attribute,ncol = 11)+
  labs(x="",y="")+
  #facet_wrap(~Attribute,scale="free_x",ncol = 15)+
  ggthemes::theme_few()+
  theme(#legend.position = "top",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  scale_color_manual("P value",values = c(brewer.pal(9,"Set1")[1:4],"grey"))+
  scale_shape_manual("Score",values = c(15:17))
  #scale_fill_manual("Risk score",values = brewer.pal(9,"Set1")[c(1:2)])
ggsave(p,filename = "RiskGroup.COX.TCGA.RiskScore.Index.All.pdf",height = 4,width = 10)


for (attr in unique(Risk.Index2$Attribute)) {
  middata <- Risk.Index2 %>% filter(Attribute == attr) %>%
    mutate(Cancer=factor(Cancer,levels = rev(sort(unique(.$Cancer)))))
  data3 <- data2 %>% filter(Attribute == attr)
  p1 <- ggplot()+
    geom_boxplot(data=middata,
                 aes(x=Score,y=Cancer, #Cancer,#
                     fill=Group))+
    #stat_ellipse(aes(fill=Group),type = "norm", geom ="polygon",alpha=0.2,color=NA)+
    #geom_beeswarm(aes(color=SampleType),size=1,dodge.width=0.75,priority = "ascending")+
    labs(x="",y="")+
    #facet_wrap(~Attribute,scale="free_x",ncol = 15)+
    ggthemes::theme_few()+
    theme(#legend.position = "top",
      axis.text = element_text(color = "black"),
      axis.text.x = element_text(size=9,angle = 90,vjust = 0.3,hjust = 1),
      axis.text.y = element_text(size=9),
      #axis.ticks.x = element_blank(),
      axis.title = element_text(size=12),
      panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
    scale_color_manual("Risk score",values = brewer.pal(9,"Set1")[c(1:2)])+
    scale_fill_manual("Risk score",values = brewer.pal(9,"Set1")[c(1:2)])#+
    #scale_x_continuous(trans = "log2")
  p2 <- ggplot(data3,aes(x=0.5,y=Cancer))+
    geom_text(aes(x=0.5,y=Cancer,label=p.signif),vjust=0.5)+
    theme_void()
  p12 <- p1 %>% aplot::insert_right(p2,width = 0.1)
  ggsave(p12,filename = paste("RiskGroup.COX.TCGA.RiskScore.Index.",attr,".pdf",sep = ""),
         height = 4.5,width = 4)
}

######## Risk score => ESTIMATE ####
setwd("K:/TCGA/Anlysis/BCAA")
library(tidyverse)
RiskScore.Cancer <- read.csv(file = "RiskGroup.COX.TCGA.RiskScore.csv") %>%
  mutate(ID=substring(SampleID,1,15))

PhenoType <- read.csv("../../mRNA.Phenotype.csv")
PhenoType <- PhenoType %>% filter(Sample.Type == "Primary Tumor" | 
                                    Sample.Type == "Solid Tissue Normal" | 
                                    Sample.Type == "Primary Blood Derived Cancer - Peripheral Blood") %>% 
  dplyr::select(Project.ID,Sample.ID,Sample.Type) %>%
  mutate(Type=ifelse(Sample.Type=="Solid Tissue Normal","Normal","Tumor")) %>%
  mutate(Cancer=str_remove(Project.ID,"TCGA-")) %>%
  filter(Type=="Tumor")

library(readxl)
data <- read_xlsx("../../PanCancer.DiverseIndex.xlsx",sheet = 2)
data2<-data[,-1]
Index.data2 <- merge(RiskScore.Cancer,data2,by="ID") %>% unique() %>%
  data.frame() %>%
  filter(SampleID %in% PhenoType$Sample.ID) %>%
  gather(Attribute,Score,StromalScore:ESTIMATEScore) %>%
  dplyr::select(ID,SampleID,Attribute,Score,Cancer,total_risk_score) %>%
  group_by(Cancer,Attribute) %>%
  mutate(cutpoint=median(total_risk_score)) %>%
  mutate(Group = if_else(total_risk_score >= cutpoint,"High","Low"))

write.csv(Index.data2,file = "RiskGroup.COX.TCGA.RiskScore.ESTIMATE.csv",row.names = F)

data2 <- compare_means(formula = Score~Group,data = Index.data2,method = "wilcox.test",
                       group.by=c("Cancer","Attribute"))

write.csv(data2,file = "RiskGroup.COX.TCGA.RiskScore.ESTIMATE.Wilcox.csv",row.names = F) 
#### Figure ####
for (attr in unique(Index.data2$Attribute)) {
  middata <- Index.data2 %>% filter(Attribute == attr) %>%
    mutate(Cancer=factor(Cancer,levels = rev(sort(unique(.$Cancer)))))
  data3 <- data2 %>% filter(Attribute == attr)
  p1 <- ggplot()+
    geom_boxplot(data=middata,
                 aes(x=Score,y=Cancer, #Cancer,#
                     fill=Group))+
    #stat_ellipse(aes(fill=Group),type = "norm", geom ="polygon",alpha=0.2,color=NA)+
    #geom_beeswarm(aes(color=SampleType),size=1,dodge.width=0.75,priority = "ascending")+
    labs(x="",y="")+
    #facet_wrap(~Attribute,scale="free_x",ncol = 15)+
    ggthemes::theme_few()+
    theme(#legend.position = "top",
      axis.text = element_text(color = "black"),
      axis.text.x = element_text(size=9,angle = 90,vjust = 0.3,hjust = 1),
      axis.text.y = element_text(size=9),
      #axis.ticks.x = element_blank(),
      axis.title = element_text(size=12),
      panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
    scale_color_manual("Risk score",values = brewer.pal(9,"Set1")[c(1:2)])+
    scale_fill_manual("Risk score",values = brewer.pal(9,"Set1")[c(1:2)])#+
  #scale_x_continuous(trans = "log2")
  p2 <- ggplot(data3,aes(x=0.5,y=Cancer))+
    geom_text(aes(x=0.5,y=Cancer,label=p.signif),vjust=0.5)+
    theme_void()
  p12 <- p1 %>% aplot::insert_right(p2,width = 0.1)
  ggsave(p12,filename = paste("RiskGroup.COX.TCGA.RiskScore.ESTIMATE.",attr,".pdf",sep = ""),
         height = 4.5,width = 4.5)
}
######## Risk score => Immune score #####
setwd("K:/TCGA/Anlysis/BCAA")
library(tidyverse)
RiskScore.Cancer <- read.csv(file = "RiskGroup.COX.TCGA.RiskScore.csv")
PhenoType <- read.csv("../../mRNA.Phenotype.csv")
PhenoType <- PhenoType %>% filter(Sample.Type == "Primary Tumor" | 
                                    Sample.Type == "Solid Tissue Normal" | 
                                    Sample.Type == "Primary Blood Derived Cancer - Peripheral Blood") %>% 
  dplyr::select(Project.ID,Sample.ID,Sample.Type) %>%
  mutate(Type=ifelse(Sample.Type=="Solid Tissue Normal","Normal","Tumor")) %>%
  mutate(Cancer=str_remove(Project.ID,"TCGA-")) %>%
  filter(Type=="Tumor")

TIMER2 <- read.csv("../../TIMER2/infiltration_estimation_for_tcga.csv",
                   check.names = F) %>%
  remove_rownames() %>% column_to_rownames("cell_type") %>%
  dplyr::select(contains("score")) %>%
  rownames_to_column("Sample.ID") %>%
  gather(CellType,Abundance,-Sample.ID)

RiskScore.TIMER2 <- read.csv(file = "RiskGroup.COX.TCGA.RiskScore.csv") %>%
  filter(SampleID %in% PhenoType$Sample.ID) %>%
  mutate(Sample.ID = substring(SampleID,1,15)) %>%
  merge(TIMER2,by="Sample.ID") %>%
  #merge(Risk.Cutoff,by = "Cancer") %>%
  group_by(Cancer,CellType) %>%
  mutate(cutpoint=median(total_risk_score)) %>%
  mutate(Group = if_else(total_risk_score >= cutpoint,"High","Low")) %>%
  ungroup() %>%
  filter(!str_detect(CellType,"cytotoxicity")) %>%
  mutate(CellType=str_remove_all(CellType,"_.*"))

write.csv(RiskScore.TIMER2,file = "RiskGroup.COX.TCGA.RiskScore.ImmuneScore.XCELL.csv",row.names = F)

data2 <- compare_means(formula = Abundance~Group,data = RiskScore.TIMER2,method = "wilcox.test",
                       group.by=c("Cancer","CellType"))

write.csv(data2,file = "RiskGroup.COX.TCGA.RiskScore.ImmuneScore.XCELL.Wilcox.csv",row.names = F) 
data2 <- read.csv("RiskGroup.COX.TCGA.RiskScore.ImmuneScore.XCELL.Wilcox.csv")
data2 %>% filter(Cancer=="BLCA")
#### Figure ####
for (attr in unique(RiskScore.TIMER2$CellType)) {
  middata <- RiskScore.TIMER2 %>% filter(CellType == attr) %>%
    mutate(Cancer=factor(Cancer,levels = rev(sort(unique(.$Cancer)))))
  data3 <- data2 %>% filter(CellType == attr)
  p1 <- ggplot()+
    geom_boxplot(data=middata,
                 aes(x=Abundance,y=Cancer, #Cancer,#
                     fill=Group))+
    #stat_ellipse(aes(fill=Group),type = "norm", geom ="polygon",alpha=0.2,color=NA)+
    #geom_beeswarm(aes(color=SampleType),size=1,dodge.width=0.75,priority = "ascending")+
    labs(x="",y="")+
    #facet_wrap(~Attribute,scale="free_x",ncol = 15)+
    ggthemes::theme_few()+
    theme(#legend.position = "top",
      axis.text = element_text(color = "black"),
      axis.text.x = element_text(size=9,angle = 90,vjust = 0.3,hjust = 1),
      axis.text.y = element_text(size=9),
      #axis.ticks.x = element_blank(),
      axis.title = element_text(size=12),
      panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
    scale_color_manual("Risk score",values = brewer.pal(9,"Set1")[c(1:2)])+
    scale_fill_manual("Risk score",values = brewer.pal(9,"Set1")[c(1:2)])+
    ggtitle(Hmisc::capitalize(attr))+
    theme(plot.title = element_text(hjust = 0.5))
  #scale_x_continuous(trans = "log2")
  p2 <- ggplot(data3,aes(x=0.5,y=Cancer))+
    geom_text(aes(x=0.5,y=Cancer,label=p.signif),vjust=0.5)+
    theme_void()
  p12 <- p1 %>% aplot::insert_right(p2,width = 0.1)
  ggsave(p12,filename = paste("RiskGroup.COX.TCGA.RiskScore.ImmuneScore.XCELL.",Hmisc::capitalize(attr),".pdf",sep = ""),
         height = 4.5,width = 4.5)
}

#### BLCA ####
setwd("K:/TCGA/Anlysis/BCAA")
library(tidyverse)
RiskScore.Cancer <- read.csv(file = "RiskGroup.COX.TCGA.RiskScore.csv")
PhenoType <- read.csv("../../mRNA.Phenotype.csv")
PhenoType <- PhenoType %>% filter(Sample.Type == "Primary Tumor" | 
                                    Sample.Type == "Solid Tissue Normal" | 
                                    Sample.Type == "Primary Blood Derived Cancer - Peripheral Blood") %>% 
  dplyr::select(Project.ID,Sample.ID,Sample.Type) %>%
  mutate(Type=ifelse(Sample.Type=="Solid Tissue Normal","Normal","Tumor")) %>%
  mutate(Cancer=str_remove(Project.ID,"TCGA-")) %>%
  filter(Type=="Tumor")
TIMER2 <- read.csv("../../TIMER2/infiltration_estimation_for_tcga.csv",
                   check.names = F) %>%
  remove_rownames() %>% column_to_rownames("cell_type") %>%
  dplyr::select(contains("score")) %>%
  rownames_to_column("Sample.ID") %>%
  gather(CellType,Abundance,-Sample.ID)

RiskScore.TIMER2 <- read.csv(file = "RiskGroup.COX.TCGA.RiskScore.csv") %>%
  filter(SampleID %in% PhenoType$Sample.ID) %>%
  mutate(Sample.ID = substring(SampleID,1,15)) %>%
  merge(TIMER2,by="Sample.ID") %>%
  #merge(Risk.Cutoff,by = "Cancer") %>%
  group_by(Cancer,CellType) %>%
  mutate(cutpoint=median(total_risk_score)) %>%
  mutate(Group = if_else(total_risk_score >= cutpoint,"High","Low")) %>%
  ungroup() %>%
  filter(!str_detect(CellType,"cytotoxicity")) %>%
  mutate(CellType=str_remove_all(CellType,"_.*")) %>%
  filter(Cancer=="BLCA") %>%
  mutate(Group=factor(Group,levels = c("Low","High")))
library(ggbeeswarm)
p12 <- ggplot(RiskScore.TIMER2,aes(x=Abundance,y=Group,fill=Group))+
  geom_violin(width=1.4) +
  geom_boxplot(width=0.2, color="grey", alpha=1)+
  #geom_jitter()+
  #stat_ellipse(aes(fill=Group),type = "norm", geom ="polygon",alpha=0.2,color=NA)+
  #geom_beeswarm(aes(color=Group),size=1,dodge.width=0.75,priority = "ascending")+
  facet_grid(CellType~.)+
  #stat_compare_means(aes(group=CellType),comparisons = list(c("Low","High")))+
  labs(x="Score",y="")+
  #facet_wrap(~Attribute,scale="free_x",ncol = 15)+
  ggthemes::theme_few()+
  theme(legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9,angle = 90,vjust = 0.3,hjust = 1),
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  #scale_color_manual("Risk score",values = brewer.pal(9,"Set1")[c(1:2)])+
  scale_fill_manual("Risk score",values = brewer.pal(9,"Set1")[c(2:1)])#+
  #scale_x_continuous(trans = "log2")
  #ggtitle(Hmisc::capitalize(attr))+
  #theme(plot.title = element_text(hjust = 0.5))
ggsave(p12,filename = "RiskGroup.COX.TCGA.RiskScore.ImmuneScore.XCELL.BLCA.pdf",height = 4,width = 2)


######## Risk score => Immune score+ESTIMATE ######
RiskScore.TIMER2 <- read.csv(file = "RiskGroup.COX.TCGA.RiskScore.ImmuneScore.XCELL.csv")
data2 <- read.csv(file = "RiskGroup.COX.TCGA.RiskScore.ImmuneScore.XCELL.Wilcox.csv") 
Wilcox.Data <- data2 %>% 
  dplyr::select(Cancer,CellType,group1,group2,p.signif) %>%
  gather(Name,Group,group1:group2) %>%
  dplyr::select(-Name)

RiskScore.TIMER2.FigureD1 <- RiskScore.TIMER2 %>%
  group_by(Cancer,CellType) %>%
  mutate(cutpoint=median(total_risk_score)) %>%
  mutate(Group = if_else(total_risk_score >= cutpoint,"High","Low")) %>%
  group_by(Cancer,CellType,Group) %>%
  summarise(MeanValue=mean(Abundance)) %>%
  ungroup() %>% data.frame(check.names = F) %>%
  merge(data2 %>% dplyr::select(Cancer,CellType,p.signif),
        by=c("Cancer","CellType")) %>%
  mutate(p.signif=if_else(p.signif=="ns","",p.signif)) %>%
  #filter(Attribute %in% c("IPS","TMB","Neoantigen","MSI","HRD","MATH","Ploidy","Purity","RNAss")) %>%
  spread(Group,MeanValue) %>%
  mutate(Enrich=if_else(p.signif == "","n.s.",if_else(High>Low,"High","Low"))) %>%
  gather(Group,Score,High:Low) %>%
  arrange(Cancer,CellType) %>%
  group_by(Cancer,CellType) %>%
  mutate(Enrich2=if_else(Enrich=="n.s.","n.s.",if_else(Enrich==Group,"Higher","Lower")))

colnames(RiskScore.TIMER2.FigureD1)[2]="Attribute"

Index.data2 <- read.csv(file = "RiskGroup.COX.TCGA.RiskScore.ESTIMATE.csv")
data2 <- read.csv(file = "RiskGroup.COX.TCGA.RiskScore.ESTIMATE.Wilcox.csv") 

RiskScore.TIMER2.FigureD2 <- Index.data2 %>%
  group_by(Cancer,Attribute) %>%
  mutate(cutpoint=median(total_risk_score)) %>%
  mutate(Group = if_else(total_risk_score >= cutpoint,"High","Low")) %>%
  group_by(Cancer,Attribute,Group) %>%
  summarise(MeanValue=mean(Score)) %>%
  ungroup() %>% data.frame(check.names = F) %>%
  merge(data2 %>% dplyr::select(Cancer,Attribute,p.signif),
        by=c("Cancer","Attribute")) %>%
  mutate(p.signif=if_else(p.signif=="ns","",p.signif)) %>%
  #filter(Attribute %in% c("IPS","TMB","Neoantigen","MSI","HRD","MATH","Ploidy","Purity","RNAss")) %>%
  spread(Group,MeanValue) %>%
  mutate(Enrich=if_else(p.signif == "","n.s.",if_else(High>Low,"High","Low"))) %>%
  gather(Group,Score,High:Low) %>%
  arrange(Cancer,Attribute) %>%
  group_by(Cancer,Attribute) %>%
  mutate(Enrich2=if_else(Enrich=="n.s.","n.s.",if_else(Enrich==Group,"Higher","Lower")))

Risk.Index4 <- rbind(RiskScore.TIMER2.FigureD1,RiskScore.TIMER2.FigureD2) %>%
  data.frame(check.names = F) %>%
  mutate(Attribute=factor(Attribute,
                          levels = c("ImmuneScore","StromalScore","ESTIMATEScore",
                                     "immune score","stroma score","microenvironment score"))) %>%
  mutate(Cancer=factor(Cancer,levels = rev(sort(unique(.$Cancer))))) %>%
  mutate(p.signif=factor(p.signif,levels = c("*","**","***","****","")))

p <- ggplot(Risk.Index4,aes(x=Group,y=Cancer))+
  geom_point(aes(size=Score,color=p.signif,shape=Enrich2))+
  scale_size_continuous("Mean",range = c(2,5))+
  facet_wrap(~Attribute,ncol = 11)+
  labs(x="",y="")+
  #facet_wrap(~Attribute,scale="free_x",ncol = 15)+
  ggthemes::theme_few()+
  theme(#legend.position = "top",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  scale_color_manual("P value",values = c(brewer.pal(9,"Set1")[1:4],"grey"))+
  scale_shape_manual("Score",values = c(15:17))
#scale_fill_manual("Risk score",values = brewer.pal(9,"Set1")[c(1:2)])
ggsave(p,filename = "RiskGroup.COX.TCGA.RiskScore.Estimate.All.pdf",height = 4.5,width = 7)

######## Risk score => HALLMARK ####
setwd("K:/TCGA/Anlysis/BCAA")
library(tidyverse)
RiskScore.Cancer <- read.csv(file = "RiskGroup.COX.TCGA.RiskScore.csv")
PhenoType <- read.csv("../../mRNA.Phenotype.csv")
PhenoType <- PhenoType %>% filter(Sample.Type == "Primary Tumor" | 
                                    Sample.Type == "Solid Tissue Normal" | 
                                    Sample.Type == "Primary Blood Derived Cancer - Peripheral Blood") %>% 
  dplyr::select(Project.ID,Sample.ID,Sample.Type) %>%
  mutate(Type=ifelse(Sample.Type=="Solid Tissue Normal","Normal","Tumor")) %>%
  mutate(Cancer=str_remove(Project.ID,"TCGA-")) %>%
  filter(Type=="Tumor")

TERM <- c("HALLMARK_E2F_TARGETS",
          "HALLMARK_G2M_CHECKPOINT",
          "HALLMARK_MYC_TARGETS_V2",
          "HALLMARK_MYC_TARGETS_V1",
          "HALLMARK_MITOTIC_SPINDLE",
          "HALLMARK_NOTCH_SIGNALING",
          "HALLMARK_WNT_BETA_CATENIN_SIGNALING",
          "HALLMARK_MTORC1_SIGNALING",
          "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
          "HALLMARK_TGF_BETA_SIGNALING",
          "HALLMARK_KRAS_SIGNALING_UP",
          "HALLMARK_KRAS_SIGNALING_DN",
          "HALLMARK_IL6_JAK_STAT3_SIGNALING",
          "HALLMARK_INTERFERON_GAMMA_RESPONSE",
          "HALLMARK_INTERFERON_ALPHA_RESPONSE",
          "HALLMARK_INFLAMMATORY_RESPONSE",
          "HALLMARK_HEDGEHOG_SIGNALING",
          "HALLMARK_COMPLEMENT",
          "HALLMARK_DNA_REPAIR",
          "HALLMARK_HYPOXIA",
          "HALLMARK_APOPTOSIS",
          "HALLMARK_P53_PATHWAY",
          "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
          "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")

files <- list.files("../../HALLMARK.ssGSEA.Pancancer/")
RiskScore.Hallmark <- data.frame()
for (file in files) {
  project <- file %>% str_remove_all("\\..*") %>% str_remove_all("TCGA-")
  print(project)
  data <- read.csv(file.path("../../HALLMARK.ssGSEA.Pancancer/",file),
                   check.names = F,row.names = 1) %>%
    rownames_to_column("HALLMARK") %>%
    gather(SampleID,Score,-HALLMARK)
  RiskScore.Hallmark <- RiskScore.Cancer %>% 
    data.frame(check.names = F) %>%
    filter(Cancer==project) %>%
    mutate(Group=ifelse(total_risk_score>=median(.$total_risk_score),"High","Low")) %>%
    merge(data,by="SampleID") %>%
    group_by(HALLMARK,Group) %>%
    #mutate(Scale.Score=scale(Score)) %>%
    mutate(MeanValue=mean(Score),
           MedianValue=median(Score),
           CVValue=sd(Score)/MeanValue*100) %>%
    ungroup() %>%
    filter(HALLMARK %in% TERM) %>%
    rbind.data.frame(RiskScore.Hallmark)
    
}
write.csv(RiskScore.Hallmark,file = "RiskGroup.COX.TCGA.RiskScore.HALLMARKscore.csv",row.names = F)

data2 <- compare_means(formula = Score~Group,data = RiskScore.Hallmark,method = "wilcox.test",
                       group.by=c("Cancer","HALLMARK"))

write.csv(data2,file = "RiskGroup.COX.TCGA.RiskScore.HALLMARKscore.Wilcox.csv",row.names = F) 

#### Figure ####
Wilcox.Data <- data2 %>% dplyr::select(Cancer,HALLMARK,group1,group2,p.signif) %>%
  gather(Name,Group,group1:group2) %>%
  dplyr::select(-Name) %>%
  mutate(HALLMARK=str_remove(HALLMARK,"HALLMARK_"))
RiskScore.Hallmark <- read.csv(file = "RiskGroup.COX.TCGA.RiskScore.HALLMARKscore.csv")
RiskScore.Hallmark2 <- RiskScore.Hallmark  %>%
  dplyr::select(Cancer,Group,HALLMARK,MeanValue,MedianValue,CVValue) %>%
  arrange(Cancer,Group,HALLMARK) %>%
  unique() %>%
  mutate(Xaxis=paste(Cancer,Group,sep = "_")) %>%
  mutate(HALLMARK=str_remove(HALLMARK,"HALLMARK_")) %>%
  merge(Wilcox.Data,by=c("Cancer","HALLMARK","Group")) %>%
  dplyr::select(Cancer,HALLMARK,Group,MedianValue,p.signif) %>%
  spread(Group,MedianValue) %>%
  mutate(Attributes=if_else(p.signif == "ns","ns",
                            if_else(High>Low,"High","Low"))) %>%
  gather(Group,MedianValue,High:Low) %>%
  mutate(Xaxis=paste(Cancer,Group,sep = "_")) %>%
  mutate(HALLMARK=factor(HALLMARK,levels = rev(sort(unique(.$HALLMARK)))))
  
write.csv(RiskScore.Hallmark2,file = "RiskGroup.COX.TCGA.RiskScore.HALLMARKscore.FigureData.csv",row.names = F)


p1 <- ggplot()+
  geom_point(data = RiskScore.Hallmark2,
             aes(x=Xaxis,y=HALLMARK,fill=Attributes,
                 size=MedianValue),
             shape=21,color="white",stroke=1.2)+
  geom_point(data = RiskScore.Hallmark2 %>% filter(p.signif != "ns"),
             aes(x=Xaxis,y=HALLMARK,fill=Attributes,
                 size=MedianValue),
             shape=21,color="black",stroke=1.2)+
  ggthemes::theme_few()+
  scale_size_continuous("Median of\nssGSEA score",range = c(1,5))+
  theme(#legend.position = "top",
    axis.text = element_text(color = "black"),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=12),
    axis.ticks.x = element_blank(),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  #scale_color_manual("Risk score",values = brewer.pal(9,"Set1")[c(1:2)])+
  #scale_fill_manual("Risk score",values = brewer.pal(9,"Set1")[c(1:3)])
  geom_vline(aes(xintercept=c(seq(2.5,44.5,2))),
             color="black",linetype="dashed")+
  ggsci::scale_fill_nejm()+
  labs(x="",y="",fill="Higher in group")

Cancer.Col <- c(brewer.pal(9,"Set1"),
                brewer.pal(8,"Set2")[1:7],
                brewer.pal(8,"Dark2")[1:6],
                brewer.pal(12,"Set3")[11:1])
p2 <- ggplot(data = RiskScore.Hallmark2 %>% 
               dplyr::select(Xaxis,Cancer) %>%
               arrange(Cancer,Xaxis) %>%
               unique(),
             aes(x=Xaxis,y=1,fill=Cancer))+
  geom_col(width = 1)+
  theme_void()+
  scale_fill_manual(values = Cancer.Col)
  #ggthemes::scale_fill_colorblind()
p3 <- ggplot(data = RiskScore.Hallmark2 %>% 
               dplyr::select(Xaxis,Group) %>%
               arrange(Group,Xaxis) %>%
               unique(),
             aes(x=Xaxis,y=1,fill=Group))+
  geom_col(width = 1)+
  theme_void()+
  scale_fill_manual("Risk score",values = c("#fe5f55","#004d61"))

p123 <- p1 %>% aplot::insert_top(p3,height = 0.05) %>%
  aplot::insert_top(p2,height = 0.05)
            
ggsave(p123,filename = "RiskGroup.COX.TCGA.RiskScore.HALLMARKscore.Figure.pdf",height = 5,width = 12)
ggsave(p123,filename = "RiskGroup.COX.TCGA.RiskScore.HALLMARKscore.Figure.Legend.pdf",height = 12,width = 12)
######## Risk score => KEGG pathway #######
setwd("K:/TCGA/Anlysis/BCAA")
library(tidyverse)
library(GSVA)
library(GSVA)
library(GSEABase)
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
#gene_set <- read.gmt("../../msigdb/c2.cp.kegg.v7.5.1.symbols.gmt")
gene_set <- getGmt("../../msigdb/c2.cp.kegg.v7.5.1.symbols.gmt")
#geneset <- gene_set %>% split(x = .$gene, f = .$term)
dirs <- list.dirs("../../Clinical.XENA/",full.names = T)
dirs <- dirs[str_detect(dirs,"//")]
dir.create("KEGG.ssGSEA.Pancancer")
for (dir in dirs) {
  project <- str_remove(dir,".*/")
  print(project)
  phenotype <- read.csv(paste(dir,"/TCGA-",project,".GDC_phenotype.tsv",sep = ""),sep = "\t",header = T) %>% 
    dplyr::filter(sample_type.samples == "Primary Tumor" | sample_type.samples == "Primary Blood Derived Cancer - Peripheral Blood")
  TPM <- read.csv(paste("../../mRNA.PanCancer.Exp/TCGA-",project,".mRNA.TPM.csv",sep = ""),check.names = F,row.names = 1)
  TPM.log <- log2(TPM+1)
  
  mRNA.Exp.Target <- TPM #[,intersect(phenotype$submitter_id.samples,colnames(TPM.log))] %>% data.frame(check.names = F)
  
  G_list <- read.csv("../LIHC/biomaRt.GeneID.tansfer.csv")
  
  dat.2 <- mRNA.Exp.Target %>% data.frame(check.names = F) %>%
    rownames_to_column("ensembl_gene_id") %>% data.frame(check.names = F) %>%
    merge(.,G_list[,1:2],by="ensembl_gene_id") %>% 
    data.frame(check.names = F) %>%
    na.omit() %>% dplyr::select(-ensembl_gene_id)
  expr_mean=aggregate(.~hgnc_symbol,mean,data=dat.2) #calculate mean for same symbol
  dat.3 <- expr_mean %>% remove_rownames() %>% column_to_rownames("hgnc_symbol")
  
  res.ssgsea <- gsva(as.matrix(dat.3), gene_set, method = "ssgsea", kcdf = "Gaussian", min.sz = 5,abs.ranking=TRUE)
  write.csv(res.ssgsea,file = paste("KEGG.ssGSEA.Pancancer/TCGA-",project,".ssGSEA.csv",sep = ""))
}

TERM <- gene_set$term[str_detect(gene_set$term,"METABOLISM$") | 
  str_detect(gene_set$term,"BIOSYNTHESIS$") |
  str_detect(gene_set$term,"DEGRADATION$") |
  str_detect(gene_set$term,"KEGG_CELL_CYCLE")|
  str_detect(gene_set$term,"KEGG_DNA_REPLICATION")|
  str_detect(gene_set$term,"KEGG_RNA_POLYMERASE")|
  str_detect(gene_set$term,"KEGG_NOTCH_SIGNALING_PATHWAY")|
  str_detect(gene_set$term,"KEGG_HEDGEHOG_REPLICATION")|
  str_detect(gene_set$term,"KEGG_BASE_EXCISION_REPAIR")|
  str_detect(gene_set$term,"KEGG_NUCLEOTIDE_EXCISION_REPAIR")|
  str_detect(gene_set$term,"KEGG_MISMATCH_REPAIR")|
  str_detect(gene_set$term,"KEGG_HOMOLOGOUS_RECOMBINATION")|
  str_detect(gene_set$term,"KEGG_MAPK_SIGNALING_PATHWAY")] %>%
  sort() %>% unique() %>% as.character() 
TERM <- TERM[TERM != "KEGG_GLYCOSYLPHOSPHATIDYLINOSITOL_GPI_ANCHOR_BIOSYNTHESIS"]

setwd("K:/TCGA/Anlysis/BCAA")
library(tidyverse)
RiskScore.Cancer <- read.csv(file = "RiskGroup.COX.TCGA.RiskScore.csv")
PhenoType <- read.csv("../../mRNA.Phenotype.csv")
PhenoType <- PhenoType %>% filter(Sample.Type == "Primary Tumor" | 
                                    Sample.Type == "Solid Tissue Normal" | 
                                    Sample.Type == "Primary Blood Derived Cancer - Peripheral Blood") %>% 
  dplyr::select(Project.ID,Sample.ID,Sample.Type) %>%
  mutate(Type=ifelse(Sample.Type=="Solid Tissue Normal","Normal","Tumor")) %>%
  mutate(Cancer=str_remove(Project.ID,"TCGA-")) %>%
  filter(Type=="Tumor")

files <- list.files(path = "KEGG.ssGSEA.Pancancer/",pattern = "csv$")
RiskScore.KEGG <- data.frame()
for (file in files) {
  project <- file %>% str_remove_all("\\..*") %>% str_remove_all("TCGA-")
  print(project)
  mid.pheno <- PhenoType %>% filter(Cancer==project)
  data <- read.csv(file.path("KEGG.ssGSEA.Pancancer/",file),
                   check.names = F,row.names = 1) %>%
    rownames_to_column("HALLMARK") %>%
    gather(SampleID,Score,-HALLMARK) %>%
    filter(SampleID %in% mid.pheno$Sample.ID)
  
  RiskScore.KEGG <- RiskScore.Cancer %>% 
    data.frame(check.names = F) %>%
    filter(Cancer==project) %>%
    mutate(Group=ifelse(total_risk_score>=median(.$total_risk_score),"High","Low")) %>%
    merge(data,by="SampleID") %>%
    group_by(HALLMARK,Group) %>%
    #mutate(Scale.Score=scale(Score)) %>%
    mutate(MeanValue=mean(Score),
           MedianValue=median(Score),
           CVValue=sd(Score)/MeanValue*100) %>%
    ungroup() %>%
    filter(HALLMARK %in% TERM) %>%
    rbind.data.frame(RiskScore.KEGG)
  
}
write.csv(RiskScore.KEGG,file = "RiskGroup.COX.TCGA.RiskScore.KEGGscore.csv",row.names = F)
library(ggpubr)
data2 <- compare_means(formula = Score~Group,data = RiskScore.KEGG,method = "wilcox.test",
                       group.by=c("Cancer","HALLMARK"))

write.csv(data2,file = "RiskGroup.COX.TCGA.RiskScore.KEGGscore.Wilcox.csv",row.names = F) 
#### Figure ####
Wilcox.Data <- data2 %>% dplyr::select(Cancer,HALLMARK,group1,group2,p.signif) %>%
  gather(Name,Group,group1:group2) %>%
  dplyr::select(-Name) %>%
  mutate(HALLMARK=str_remove(HALLMARK,"KEGG_"))
RiskScore.KEGG <- read.csv(file = "RiskGroup.COX.TCGA.RiskScore.KEGGscore.csv")
RiskScore.KEGG2 <- RiskScore.KEGG  %>%
  dplyr::select(Cancer,Group,HALLMARK,MeanValue,MedianValue,CVValue) %>%
  arrange(Cancer,Group,HALLMARK) %>%
  unique() %>%
  mutate(Xaxis=paste(Cancer,Group,sep = "_")) %>%
  mutate(HALLMARK=str_remove(HALLMARK,"KEGG_")) %>%
  merge(Wilcox.Data,by=c("Cancer","HALLMARK","Group")) %>%
  dplyr::select(Cancer,HALLMARK,Group,MedianValue,p.signif) %>%
  spread(Group,MedianValue) %>%
  mutate(Attributes=if_else(p.signif == "ns","ns",
                            if_else(High>Low,"High","Low"))) %>%
  gather(Group,MedianValue,High:Low) %>%
  mutate(Xaxis=paste(Cancer,Group,sep = "_")) %>%
  mutate(HALLMARK=factor(HALLMARK,levels = rev(sort(unique(.$HALLMARK)))))

write.csv(RiskScore.KEGG2,file = "RiskGroup.COX.TCGA.RiskScore.KEGGscore.FigureData.csv",row.names = F)

p1 <- ggplot()+
  geom_point(data = RiskScore.KEGG2,
             aes(x=Xaxis,y=HALLMARK,fill=Attributes,
                 size=MedianValue),
             shape=21,color="white",stroke=1.2)+
  geom_point(data = RiskScore.KEGG2 %>% filter(p.signif != "ns"),
             aes(x=Xaxis,y=HALLMARK,fill=Attributes,
                 size=MedianValue),
             shape=21,color="black",stroke=1.2)+
  ggthemes::theme_few()+
  scale_size_continuous("Median of\nssGSEA score",range = c(1,4))+
  theme(#legend.position = "top",
    axis.text = element_text(color = "black"),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=9),
    axis.ticks.x = element_blank(),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  #scale_color_manual("Risk score",values = brewer.pal(9,"Set1")[c(1:2)])+
  #scale_fill_manual("Risk score",values = brewer.pal(9,"Set1")[c(1:3)])
  geom_vline(aes(xintercept=c(seq(2.5,44.5,2))),
             color="black",linetype="dashed")+
  ggsci::scale_fill_nejm()+
  labs(x="",y="",fill="Higher in group")

Cancer.Col <- c(brewer.pal(9,"Set1"),
                brewer.pal(8,"Set2")[1:7],
                brewer.pal(8,"Dark2")[1:6],
                brewer.pal(12,"Set3")[11:1])
p2 <- ggplot(data = RiskScore.KEGG2 %>% 
               dplyr::select(Xaxis,Cancer) %>%
               arrange(Cancer,Xaxis) %>%
               unique(),
             aes(x=Xaxis,y=1,fill=Cancer))+
  geom_col(width = 1)+
  theme_void()+
  scale_fill_manual(values = Cancer.Col)
#ggthemes::scale_fill_colorblind()
p3 <- ggplot(data = RiskScore.KEGG2 %>% 
               dplyr::select(Xaxis,Group) %>%
               arrange(Group,Xaxis) %>%
               unique(),
             aes(x=Xaxis,y=1,fill=Group))+
  geom_col(width = 1)+
  theme_void()+
  scale_fill_manual("Risk score",values = c("#fe5f55","#004d61"))

p123 <- p1 %>% aplot::insert_top(p3,height = 0.05) %>%
  aplot::insert_top(p2,height = 0.05)

ggsave(p123,filename = "RiskGroup.COX.TCGA.RiskScore.KEGGscore.Figure.pdf",height = 10,width = 12)


######## Risk score => MUTATION + CNV ####
#### Mutation => Cancer #####
setwd("K:/TCGA/Anlysis/BCAA")
#Driver <- read.csv("mutation_download_tab.txt",sep = "\t",check.names = F)
library(tidyverse)
RiskScore.Cancer <- read.csv(file = "RiskGroup.COX.TCGA.RiskScore.csv")
PhenoType <- read.csv("../../mRNA.Phenotype.csv")
PhenoType <- PhenoType %>% filter(Sample.Type == "Primary Tumor" | 
                                    Sample.Type == "Solid Tissue Normal" | 
                                    Sample.Type == "Primary Blood Derived Cancer - Peripheral Blood") %>% 
  dplyr::select(Project.ID,Sample.ID,Sample.Type) %>%
  mutate(Type=ifelse(Sample.Type=="Solid Tissue Normal","Normal","Tumor")) %>%
  mutate(Cancer=str_remove(Project.ID,"TCGA-")) %>%
  filter(Type=="Tumor")

dirs <- list.dirs("../../Clinical.XENA/")
dirs <- dirs[str_detect(dirs,"//")]
dir.create("RiskGenome")
for (dir in dirs) {
  project <- str_remove_all(dir,".*/")
  print(project)
  filename <- paste(project,"/TCGA-",project,".mutect2_snv.tsv",sep = "")
  tmp = read.csv(file.path("../../Clinical.XENA/",filename),sep = "\t",header=TRUE)
  colnames(tmp) =c( "Tumor_Sample_Barcode", "Hugo_Symbol", 
                    "Chromosome", "Start_Position", 
                    "End_Position", "Reference_Allele", "Tumor_Seq_Allele2", 
                    "HGVSp_Short" , 'effect' ,"Consequence",
                    "vaf" )
  tmp$Entrez_Gene_Id =1
  tmp$Center ='ucsc'
  tmp$NCBI_Build ='GRCh38'
  tmp$NCBI_Build ='GRCh38'
  tmp$Strand ='+'
  tmp$Variant_Classification = tmp$effect
  tail(sort(table(tmp$Variant_Classification )))
  tmp$Tumor_Seq_Allele1 = tmp$Reference_Allele
  tmp$Variant_Type = ifelse(
    tmp$Reference_Allele %in% c('A','C','T','G') & tmp$Tumor_Seq_Allele2 %in% c('A','C','T','G'),
    'SNP','INDEL'
  )
  tmp <- tmp %>% mutate(Variant_Classification = 
                      if_else(Variant_Classification == "inframe_insertion","In_Frame_Ins",
                              if_else(Variant_Classification == "inframe_deletion","In_Frame_Del",
                                      if_else(Variant_Classification == "frameshift_variant","Frame_Shift_Var",
                                              if_else(Variant_Classification == "missense_variant","Missense_Mutation",
                                                      if_else(Variant_Classification == "synonymous_variant","Nonsense_Mutation",
                                                              if_else(Variant_Classification == "3_prime_UTR_variant" | 
                                                                        Variant_Classification == "5_prime_UTR_variant","Prime_UTR_Var",
                                                                      if_else(Variant_Classification == "intron_variant","Intron_Var",
                                                                              if_else(str_detect(Variant_Classification,";"),"Multi_Hit",
                                                                                      if_else(str_detect(Variant_Classification,"splice"),"Splice_Site","Multi_Hit"))))))))))
  #DriverGene <- Driver %>% filter(cancer_type_abbr == project)
  #DriverGeneList <- DriverGene$driver_gene %>% str_split(",") %>% unlist() %>% unique()
  
  if (c(project) %in% RiskScore.Cancer$Cancer) {
    RiskGroup <- RiskScore.Cancer %>% filter(Cancer==project) %>%
      mutate(Group=if_else(total_risk_score >= median(total_risk_score),"High","Low"))
    
    High.Risk <- RiskGroup %>% filter(Group=="High")
    Low.Risk <- RiskGroup %>% filter(Group=="Low")
    
    HighRisk.Sample <- tmp %>% 
      filter(Tumor_Sample_Barcode %in% High.Risk$SampleID)
    
    LowRisk.Sample <- tmp %>% 
      filter(Tumor_Sample_Barcode %in% Low.Risk$SampleID)
    write.csv(HighRisk.Sample,file = file.path("RiskGenome",paste(project,".HighRisk.csv",sep = "")),row.names = F)
    write.csv(LowRisk.Sample,file = file.path("RiskGenome",paste(project,".LowRisk.csv",sep = "")),row.names = F)
    
    library(RColorBrewer)
    #vc_cols = c("#fff1c1",
    #            "#5bd1d7","#348498","#ff502f","#004d61",
    #            "#e41749","#f5587b","#ff8a5c","#8bc24c")
    vc_cols <- rev(pal_npg("nrc")(9)) #ggthemes::ggthemes_data$stata$colors$schemes$economist$value[14:6]
    names(vc_cols) <- c("In_Frame_Ins",
                        "Missense_Mutation",
                        "In_Frame_Del",
                        "Frame_Shift_Var",
                        "Nonsense_Mutation",
                        "Prime_UTR_Var",
                        "Intron_Var",
                        "Splice_Site",
                        "Multi_Hit"
    )
    #DriverGeneList2 <- DriverGeneList[DriverGeneList %in% HighRisk.Sample$Hugo_Symbol]
    library(maftools)
    HighRisk = read.maf(maf = HighRisk.Sample,
                        vc_nonSyn=c("In_Frame_Ins",
                                    "In_Frame_Del",
                                    "Frame_Shift_Var",
                                    "Missense_Mutation",
                                    "Nonsense_Mutation",
                                    "Prime_UTR_Var",
                                    "Intron_Var",
                                    "Splice_Site"
                        ))
    LowRisk = read.maf(maf = LowRisk.Sample,
                       vc_nonSyn=c("In_Frame_Ins",
                                   "In_Frame_Del",
                                   "Frame_Shift_Var",
                                   "Missense_Mutation",
                                   "Nonsense_Mutation",
                                   "Prime_UTR_Var",
                                   "Intron_Var",
                                   "Splice_Site"
                       ))
    
    library(magrittr)
    laml <- read.maf(maf = tmp %>% 
                       filter(Tumor_Sample_Barcode %in% RiskGroup$SampleID),
                     clinicalData = RiskGroup %>% 
                       dplyr::select(SampleID,Group) %>%
                       set_colnames(c("Tumor_Sample_Barcode","Riskscore")),
                     vc_nonSyn=c("In_Frame_Ins",
                                 "In_Frame_Del",
                                 "Frame_Shift_Var",
                                 "Missense_Mutation",
                                 "Nonsense_Mutation",
                                 "Prime_UTR_Var",
                                 "Intron_Var",
                                 "Splice_Site"
                     ))
    fab.ce = clinicalEnrichment(maf = laml, clinicalFeature = 'Riskscore')
    #Sig.Table <- fab.ce$groupwise_comparision[p_value < 0.05]
    write.csv(fab.ce$groupwise_comparision,file = file.path("RiskGenome/",paste("High.vs.Low.",project,".csv",sep = "")),row.names = F)
    GeneList <- laml@gene.summary$Hugo_Symbol[1:15]
    coOncoplot(m1 = HighRisk, 
               m2 = LowRisk, 
               m1Name = 'High risk', 
               m2Name = 'Low risk', 
               genes = GeneList,
               #genes = 
               removeNonMutated = F,
               colors = vc_cols)
    library(export)
    graph2pdf(file=file.path("RiskGenome/",
                             paste("High.vs.Low.",project,".Oncoplot.pdf",sep = "")),
              height = 3.5,width = 8)
  }
}

#### Mutation => Pancancer ####
setwd("K:/TCGA/Anlysis/BCAA")
Driver <- read.csv("mutation_download_tab.txt",sep = "\t",check.names = F)
library(tidyverse)
RiskScore.Cancer <- read.csv(file = "RiskGroup.COX.TCGA.RiskScore.csv")
  
PhenoType <- read.csv("../../mRNA.Phenotype.csv")
PhenoType <- PhenoType %>% filter(Sample.Type == "Primary Tumor" | 
                                    Sample.Type == "Solid Tissue Normal" | 
                                    Sample.Type == "Primary Blood Derived Cancer - Peripheral Blood") %>% 
  dplyr::select(Project.ID,Sample.ID,Sample.Type) %>%
  mutate(Type=ifelse(Sample.Type=="Solid Tissue Normal","Normal","Tumor")) %>%
  mutate(Cancer=str_remove(Project.ID,"TCGA-")) %>%
  filter(Type=="Tumor")

dirs <- list.dirs("../../Clinical.XENA/")
dirs <- dirs[str_detect(dirs,"//")]
#dir.create("RiskGenome")
SNP.Data <- data.frame()
SNP.Clinical <- data.frame()
for (dir in dirs) {
  project <- str_remove_all(dir,".*/")
  print(project)
  filename <- paste(project,"/TCGA-",project,".mutect2_snv.tsv",sep = "")
  tmp = read.csv(file.path("../../Clinical.XENA/",filename),sep = "\t",header=TRUE)
  colnames(tmp) =c( "Tumor_Sample_Barcode", "Hugo_Symbol", 
                    "Chromosome", "Start_Position", 
                    "End_Position", "Reference_Allele", "Tumor_Seq_Allele2", 
                    "HGVSp_Short" , 'effect' ,"Consequence",
                    "vaf" )
  tmp$Entrez_Gene_Id =1
  tmp$Center ='ucsc'
  tmp$NCBI_Build ='GRCh38'
  tmp$NCBI_Build ='GRCh38'
  tmp$Strand ='+'
  tmp$Variant_Classification = tmp$effect
  tail(sort(table(tmp$Variant_Classification )))
  tmp$Tumor_Seq_Allele1 = tmp$Reference_Allele
  tmp$Variant_Type = ifelse(
    tmp$Reference_Allele %in% c('A','C','T','G') & tmp$Tumor_Seq_Allele2 %in% c('A','C','T','G'),
    'SNP','INDEL'
  )
  tmp <- tmp %>% mutate(Variant_Classification = 
                          if_else(Variant_Classification == "inframe_insertion","In_Frame_Ins",
                                  if_else(Variant_Classification == "inframe_deletion","In_Frame_Del",
                                          if_else(Variant_Classification == "frameshift_variant","Frame_Shift_Var",
                                                  if_else(Variant_Classification == "missense_variant","Missense_Mutation",
                                                          if_else(Variant_Classification == "synonymous_variant","Nonsense_Mutation",
                                                                  if_else(Variant_Classification == "3_prime_UTR_variant" | 
                                                                            Variant_Classification == "5_prime_UTR_variant","Prime_UTR_Var",
                                                                          if_else(Variant_Classification == "intron_variant","Intron_Var",
                                                                                  if_else(str_detect(Variant_Classification,";"),"Multi_Hit",
                                                                                          if_else(str_detect(Variant_Classification,"splice"),"Splice_Site","Multi_Hit"))))))))))
  if (c(project) %in% RiskScore.Cancer$Cancer) {
    RiskGroup <- RiskScore.Cancer %>% filter(Cancer==project) %>%
      mutate(Group=if_else(total_risk_score >= median(total_risk_score),"High","Low"))
    
    SNP.Data <- tmp %>% filter(Tumor_Sample_Barcode %in% RiskGroup$SampleID) %>%
      rbind.data.frame(SNP.Data)
    SNP.Clinical <- data.frame(Tumor_Sample_Barcode=tmp$Tumor_Sample_Barcode %>% unique(),
                               Cancer=project) %>%
      filter(Tumor_Sample_Barcode %in% RiskGroup$SampleID) %>%
      merge(RiskGroup[,c("SampleID","Group")],
            by.x="Tumor_Sample_Barcode",
            by.y="SampleID") %>%
      rbind.data.frame(SNP.Clinical)
    }
}
write.csv(SNP.Data,file = "RiskGenome/SNP.Data.csv",row.names = F)
write.csv(SNP.Clinical,file = "RiskGenome/SNP.Clinical.csv",row.names = F)

setwd("K:/TCGA/Anlysis/BCAA")
library(tidyverse)
SNP.Data <- read.csv(file = "RiskGenome/SNP.Data.csv")
SNP.Clinical <- read.csv(file = "RiskGenome/SNP.Clinical.csv")
library(ggsci)
vc_cols <- rev(pal_npg("nrc")(9)) #ggthemes::ggthemes_data$stata$colors$schemes$economist$value[14:6]
names(vc_cols) <- c("In_Frame_Ins",
                    "Missense_Mutation",
                    "In_Frame_Del",
                    "Frame_Shift_Var",
                    "Nonsense_Mutation",
                    "Prime_UTR_Var",
                    "Intron_Var",
                    "Splice_Site",
                    "Multi_Hit"
)

Immune.Col <- c("#293462","#a64942","#fe5f55","#fff1c1",
                "#5bd1d7","#348498","#004d61","#ff502f",
                "#e41749","#f5587b","#ff8a5c","#fff591",
                "#001871","#ff585d","#ffb549","#41b6e6",
                "#515bd4","#8134af","#dd2a7b","#feda77",
                "#96ceb4","#ffeead","#d9534f","#ffad60",
                "#05445c","#f2317f","#5c4f74","#040000",
                "#de4307","#f29c2b","#f6d04d","#8bc24c",
                "#fef4a9","#3b9a9c","#4bc2c5","#78fee0",
                "#fad3cf","#a696c8","#2470a0","#c50d66",
                "#f07810","#eec60a","#dd7777","#1687a7",
                "#014955","#dd0a35","#ED5485","#FFE869",
                "#bc8420","#F0B775","#D25565","#2E94B9")

Cancer.Col = Immune.Col[c(11:43)]
names(Cancer.Col) = unique(SNP.Clinical$Cancer)
Cancer.Col <- list(Cancer=Cancer.Col)

Group.Col <- Immune.Col[2:1]
names(Group.Col) = unique(SNP.Clinical$Group)
Group.Col <- list(Group=Group.Col)

SNP.Clinical2 <- SNP.Clinical %>%
  arrange(Cancer,Group)
High.D <- SNP.Clinical %>% filter(Group=="High")
Low.D <- SNP.Clinical %>% filter(Group=="Low")
library(maftools)
tcga <- read.maf(maf = SNP.Data,clinicalData = SNP.Clinical,
                 vc_nonSyn=c("In_Frame_Ins",
                             "In_Frame_Del",
                             "Frame_Shift_Var",
                             "Missense_Mutation",
                             "Nonsense_Mutation",
                             "Prime_UTR_Var",
                             "Intron_Var",
                             "Splice_Site"
                 ))


GeneList <- tcga@gene.summary$Hugo_Symbol[1:100]

tcga.high <- read.maf(maf = SNP.Data %>% 
                        filter(Tumor_Sample_Barcode %in% High.D$Tumor_Sample_Barcode),
                      clinicalData = High.D,
                      vc_nonSyn=c("In_Frame_Ins",
                                  "In_Frame_Del",
                                  "Frame_Shift_Var",
                                  "Missense_Mutation",
                                  "Nonsense_Mutation",
                                  "Prime_UTR_Var",
                                  "Intron_Var",
                                  "Splice_Site"
                      ))

maftools::oncoplot(maf = tcga.high,bgCol = "white",
                   colors = vc_cols,
                   additionalFeatureCol = "white",
                   showTitle =F,
                   annoBorderCol = "white",
                   sortByAnnotation = T,
                   drawRowBar = F,
                   genes = GeneList,
                   #sampleOrder = SNP.Clinical2$Tumor_Sample_Barcode,
                   clinicalFeatures=c("Cancer"))

tcga.Low <- read.maf(maf = SNP.Data %>% 
                        filter(Tumor_Sample_Barcode %in% Low.D$Tumor_Sample_Barcode),
                      clinicalData = Low.D,
                      vc_nonSyn=c("In_Frame_Ins",
                                  "In_Frame_Del",
                                  "Frame_Shift_Var",
                                  "Missense_Mutation",
                                  "Nonsense_Mutation",
                                  "Prime_UTR_Var",
                                  "Intron_Var",
                                  "Splice_Site"
                      ))

maftools::oncoplot(maf = tcga.Low,bgCol = "white",
                   colors = vc_cols,
                   additionalFeatureCol = "white",
                   showTitle =F,
                   annoBorderCol = "white",
                   sortByAnnotation = T,
                   drawRowBar = F,
                   genes = GeneList,
                   #sampleOrder = SNP.Clinical2$Tumor_Sample_Barcode,
                   clinicalFeatures=c("Cancer"))
#### Mutation => Fisher test => Pancancer => Figure ####
High.ED <- tcga.high@gene.summary
Low.ED <- tcga.Low@gene.summary
Gene.Fisher <- data.frame()
for (gene in GeneList) {
  Count1 <- (High.ED[High.ED$Hugo_Symbol == gene,])$AlteredSamples
  Count2 <- (Low.ED[Low.ED$Hugo_Symbol == gene,])$AlteredSamples
  test <- fisher.test(matrix(c(Count1,3203-Count1,Count2,3183-Count2),nrow = 2))
  Gene.Fisher <- data.frame(Pvalue=test$p.value) %>%
    mutate(Cancer="Pooled",SYMBOL=gene) %>%
    rbind.data.frame(Gene.Fisher)
}
write.csv(Gene.Fisher,file = "RiskScore.Top10.SNP.csv")
Gene.Fisher$SYMBOL <- factor(Gene.Fisher$SYMBOL,levels = Gene.Fisher$SYMBOL)
p <- ggplot(Gene.Fisher,aes(x=-log10(Pvalue),y=SYMBOL,color=SYMBOL))+
  geom_segment(aes(x=0,xend=-log10(Pvalue),y=SYMBOL,yend=SYMBOL))+
  geom_point(size=3)+
  ggthemes::theme_few()+
  scale_size_continuous("Median of\nssGSEA score",range = c(1,5))+
  theme(legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9,angle = 90,vjust = 0.5,hjust = 1),
    axis.text.y = element_text(size=12),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  #scale_color_manual("Risk score",values = brewer.pal(9,"Set1")[c(1:2)])+
  #scale_fill_manual("Risk score",values = brewer.pal(9,"Set1")[c(1:3)])
  geom_vline(aes(xintercept=-log10(0.05)),color="black",linetype="dashed")+
  #ggsci::scale_color_nejm()
  ggthemes::scale_color_stata(scheme = "s1rcolor")+
  labs(y="")
ggsave(p,filename = "RiskGenome.High.vs.Low.PanCancer.pdf",height = 4,width = 2.5)

#### Mutation => Fisher test => Cancer => Figure ####
files <- list.files("RiskGenome/",pattern = "csv$")
files <- files[str_detect(files,"vs")]
Fisher.Cancer <- data.frame()
for (file in files) {
  project=(file %>% str_split("\\.") %>% unlist())[4]
  Fisher.Cancer <- read.csv(file.path("RiskGenome/",file)) %>%
    filter(p_value <= 0.05) %>%
    #dplyr::select(p_value,Hugo_Symbol) %>%
    mutate(MutatedCount = str_remove_all(n_mutated_group1," of .*")) %>%
    mutate(Cancer=project) %>%
    rbind.data.frame(Fisher.Cancer)
}
Fisher.Cancer <- Fisher.Cancer %>% unique()
write.csv(Fisher.Cancer,file = "RiskGenome.High.vs.Low.Cancer.csv",row.names = F)
Fisher.Cancer2 <- Fisher.Cancer %>%
  data.frame() %>%
  group_by(Cancer) %>%
  dplyr::top_n(-3,p_value) %>%
  arrange(Cancer) %>%
  mutate(Hugo_Symbol=factor(Hugo_Symbol,levels = rev(unique(.$Hugo_Symbol))))

p2 <- ggplot(Fisher.Cancer2,aes(x=Cancer,y=Hugo_Symbol,
                          color=-log10(p_value),
                          size=-log10(p_value)))+
  geom_point()+
  ggthemes::theme_few()+
  scale_size_continuous("-log10(Fisher's P)",range = c(1,5))+
  theme(#legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9,angle = 90,vjust = 0.3,hjust = 1),
        axis.text.y = element_text(size=9),
        #axis.ticks.x = element_blank(),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  #scale_color_manual("Risk score",values = brewer.pal(9,"Set1")[c(1:2)])+
  #scale_fill_manual("Risk score",values = brewer.pal(9,"Set1")[c(1:3)])
  #geom_vline(aes(xintercept=-log10(0.05)),color="black",linetype="dashed")+
  #ggsci::scale_color_nejm()
  #ggthemes::scale_color_stata(scheme = "s1rcolor")+
  scale_color_gradientn("-log10(Fisher's P)",colours = brewer.pal(11,"RdBu")[5:1])+
  labs(y="",x="")
ggsave(p2,filename = "RiskGenome.High.vs.Low.Cancer.pdf",height = 8,width = 5)

###### CNV #####
setwd("K:/TCGA/Anlysis/BCAA")
#Driver <- read.csv("mutation_download_tab.txt",sep = "\t",check.names = F)
library(tidyverse)
RiskScore.Cancer <- read.csv(file = "RiskGroup.COX.TCGA.RiskScore.csv")
PhenoType <- read.csv("../../mRNA.Phenotype.csv")
PhenoType <- PhenoType %>% filter(Sample.Type == "Primary Tumor" | 
                                    Sample.Type == "Solid Tissue Normal" | 
                                    Sample.Type == "Primary Blood Derived Cancer - Peripheral Blood") %>% 
  dplyr::select(Project.ID,Sample.ID,Sample.Type) %>%
  mutate(Type=ifelse(Sample.Type=="Solid Tissue Normal","Normal","Tumor")) %>%
  mutate(Cancer=str_remove(Project.ID,"TCGA-")) %>%
  filter(Type=="Tumor")

for (project in unique(RiskScore.Cancer$Cancer)) {
  print(project)
  RiskScore <- RiskScore.Cancer %>% filter(Cancer==project) %>%
    mutate(Group=if_else(total_risk_score > median(total_risk_score),"High","Low")) %>%
    dplyr::select(SampleID,Group) %>%
    filter(SampleID %in% PhenoType$Sample.ID)
  library(data.table)
  data <- fread(file.path("../../GSCA/",paste(project,".cnv.tsv.gz",sep = ""))) %>%
    dplyr::select(-entrez,-cancer_type) %>% remove_rownames() %>%
    column_to_rownames("symbol") %>%
    t() %>% data.frame(check.names = F) %>%
    rownames_to_column("Sample.ID") %>%
    mutate(SampleID=substring(Sample.ID,1,16)) %>%
    dplyr::select(-Sample.ID)
  
  HighRisk.Data <- merge(RiskScore,data,by="SampleID") %>%
    data.frame(check.names = F) %>%
    filter(Group=="High") %>%
    remove_rownames() %>%
    column_to_rownames("SampleID") %>%
    dplyr::select(-Group)
  HighRisk.Count <- nrow(HighRisk.Data)
  
  LowRisk.Data <- merge(RiskScore,data,by="SampleID") %>%
    data.frame(check.names = F) %>%
    filter(Group=="Low") %>%
    remove_rownames() %>%
    column_to_rownames("SampleID") %>%
    dplyr::select(-Group)
  LowRisk.Count <- nrow(LowRisk.Data)
  
  Amplification.Fisher <- data.frame()
  Deletion.Fisher <- data.frame()
  
  for (name in colnames(HighRisk.Data)) {
    print(name)
    middata1 <- HighRisk.Data[[name]]
    middata1[middata1<0]=0
    middata1[middata1==2]=1
    
    middata2 <- HighRisk.Data[[name]]
    middata2[middata2>0]=0
    middata2[middata2==-2]=-1
    
    middata3 <- LowRisk.Data[[name]]
    middata3[middata3<0]=0
    middata3[middata3==2]=1
    
    middata4 <- LowRisk.Data[[name]]
    middata4[middata4>0]=0
    middata4[middata4==-2]=-1
    #LowRisk.Deletion <- data.frame(SYMBOL=name,Count1=abs(sum(middata4))) %>%
    #  mutate(Cancer=project) %>%
    #  mutate(Count2=LowRisk.Count-Count1) %>%
    #  rbind.data.frame(LowRisk.Deletion)
    
    #test <- fisher.test(matrix(c(sum(middata1),
    #                             HighRisk.Count-sum(middata1),
    #                             sum(middata3),
    #                             LowRisk.Count-sum(middata3)),nrow = 2))
    #Amplification.Fisher <- data.frame(SYMBOL=name,
    #                                   Fisher.P=test$p.value,
    #                                   OR=test$estimate) %>%
    #  mutate(High.CNV=sum(middata1),Low.CNV=sum(middata3)) %>%
    #  mutate(High.Sample=HighRisk.Count,Low.Sample=LowRisk.Count) %>%
    #  mutate(Cancer=project,CNVType="Amplification") %>%
    #  rbind.data.frame(Amplification.Fisher)
    
    test <- fisher.test(matrix(c(abs(sum(middata2)),
                                 HighRisk.Count-abs(sum(middata2)),
                                 abs(sum(middata4)),
                                 LowRisk.Count-abs(sum(middata4))),nrow = 2))
    
    Deletion.Fisher <- data.frame(SYMBOL=name,
                                  Fisher.P=test$p.value,
                                  OR=test$estimate) %>%
      mutate(High.CNV=abs(sum(middata2)),Low.CNV=abs(sum(middata4))) %>%
      mutate(High.Sample=HighRisk.Count,Low.Sample=LowRisk.Count) %>%
      mutate(Cancer=project,CNVType="Deletion") %>%
      rbind.data.frame(Deletion.Fisher)
    
  }
  write.csv(Deletion.Fisher,file = file.path("RiskGenome/",paste("RiskGenome.CNV.Deletion.Fisher.",project,".csv",sep = "")),row.names = F)
  #write.csv(Amplification.Fisher,file = file.path("RiskGenome/",paste("RiskGenome.CNV.Amplification.Fisher.",project,".csv",sep = "")),row.names = F)
  
}
#write.csv(Deletion.Fisher,file = "RiskGenome.CNV.Deletion.Fisher.Cancer.csv",row.names = F)
#write.csv(Amplification.Fisher,file = "RiskGenome.CNV.Amplification.Fisher.Cancer.csv",row.names = F)
Pancancer.A.CNV <- data.frame()
Pancancer.D.CNV <- data.frame()
for (project in unique(RiskScore.Cancer$Cancer)) {
  Pancancer.D.CNV <- read.csv(file.path("RiskGenome/",paste("RiskGenome.CNV.Deletion.Fisher.",project,".csv",sep = ""))) %>%
    rbind.data.frame(Pancancer.D.CNV)
  Pancancer.A.CNV <- read.csv(file.path("RiskGenome/",paste("RiskGenome.CNV.Amplification.Fisher.",project,".csv",sep = ""))) %>%
    rbind.data.frame(Pancancer.A.CNV)
}
### amplification ###
Pancancer.A.CNV2 <- Pancancer.A.CNV %>%
  group_by(SYMBOL) %>%
  summarise(High.CNV=sum(High.CNV),
            Low.CNV=sum(Low.CNV),
            High.Sample=sum(High.Sample),
            Low.Sample=sum(Low.Sample)) %>%
  mutate(Total.CNV=High.CNV+Low.CNV) %>%
  arrange(desc(Total.CNV))

Pancancer.A.Fisher <- data.frame()
for (gene in Pancancer.A.CNV2$SYMBOL) {
  print(gene)
  middata<-Pancancer.A.CNV2 %>% filter(SYMBOL==gene)
  test <- fisher.test(matrix(c(middata$High.CNV,
                               middata$Low.CNV,
                               middata$High.Sample-middata$High.CNV,
                               middata$Low.Sample-middata$Low.CNV),nrow = 2))
  if (test$p.value <= 0.05) {
    if (middata$High.CNV > middata$Low.CNV) {
      Pancancer.A.Fisher <- data.frame(SYMBOL=gene,Pvalue=test$p.value,OR=test$estimate,
                                       Enrich="RiskHigh") %>%
        rbind.data.frame(Pancancer.A.Fisher)
    }else{
      Pancancer.A.Fisher <- data.frame(SYMBOL=gene,Pvalue=test$p.value,OR=test$estimate,
                                       Enrich="RiskLow") %>%
        rbind.data.frame(Pancancer.A.Fisher)
    }
  }else{
    Pancancer.A.Fisher <- data.frame(SYMBOL=gene,Pvalue=test$p.value,OR=test$estimate,
                                     Enrich="NS") %>%
      rbind.data.frame(Pancancer.A.Fisher)
  }
}
write.csv(Pancancer.A.Fisher,file = "RiskGenome/RiskGenome.CNV.Amplification.Fisher.Pancancer.csv",row.names = F)
### Deletion ###
Pancancer.D.CNV2 <- Pancancer.D.CNV %>%
  group_by(SYMBOL) %>%
  summarise(High.CNV=sum(High.CNV),
            Low.CNV=sum(Low.CNV),
            High.Sample=sum(High.Sample),
            Low.Sample=sum(Low.Sample)) %>%
  mutate(Total.CNV=High.CNV+Low.CNV) %>%
  arrange(desc(Total.CNV))

Pancancer.D.CNV2 <- Pancancer.D.CNV %>%
  group_by(SYMBOL) %>%
  summarise(High.CNV=sum(High.CNV),
            Low.CNV=sum(Low.CNV),
            High.Sample=sum(High.Sample),
            Low.Sample=sum(Low.Sample)) %>%
  mutate(Total.CNV=High.CNV+Low.CNV) %>%
  arrange(desc(Total.CNV))

Pancancer.D.Fisher <- data.frame()
for (gene in Pancancer.D.CNV2$SYMBOL[]) {
  print(gene)
  middata<-Pancancer.D.CNV2 %>% filter(SYMBOL==gene)
  test <- fisher.test(matrix(c(middata$High.CNV,
                               middata$Low.CNV,
                               middata$High.Sample-middata$High.CNV,
                               middata$Low.Sample-middata$Low.CNV),nrow = 2))
  if (test$p.value <= 0.05) {
    if (middata$High.CNV > middata$Low.CNV) {
      Pancancer.D.Fisher <- data.frame(SYMBOL=gene,Pvalue=test$p.value,OR=test$estimate,
                                       Enrich="RiskHigh") %>%
        rbind.data.frame(Pancancer.D.Fisher)
    }else{
      Pancancer.D.Fisher <- data.frame(SYMBOL=gene,Pvalue=test$p.value,OR=test$estimate,
                                       Enrich="RiskLow") %>%
        rbind.data.frame(Pancancer.D.Fisher)
    }
  }else{
    Pancancer.D.Fisher <- data.frame(SYMBOL=gene,Pvalue=test$p.value,OR=test$estimate,
                                     Enrich="NS") %>%
      rbind.data.frame(Pancancer.D.Fisher)
  }
}
write.csv(Pancancer.D.Fisher,file = "RiskGenome/RiskGenome.CNV.Deletion.Fisher.Pancancer.csv",row.names = F)
#### Figure => Pancancer ####
Pancancer.A.Fisher <- read.csv(file = "RiskGenome/RiskGenome.CNV.Amplification.Fisher.Pancancer.csv")
Pancancer.D.Fisher <- read.csv(file = "RiskGenome/RiskGenome.CNV.Deletion.Fisher.Pancancer.csv")
Pancancer.A.D.Fisher <- Pancancer.A.Fisher %>%
  filter(Pvalue <= 0.05) %>%
  mutate(CNVType="Amplification") %>%
  rbind.data.frame(Pancancer.D.Fisher %>% filter(Pvalue <= 0.05) %>% mutate(CNVType="Deletion")) %>%
  remove_rownames() %>%
  arrange(Pvalue)
write.csv(Pancancer.A.D.Fisher,file = "RiskGenome/RiskGenome.CNV.Deletion.Amplification.Fisher.Sig.Pancancer.csv",row.names = F)

Pancancer.A.D.Fisher.FigureData <- Pancancer.A.D.Fisher %>%
  filter(SYMBOL %in% c("TP53","MYC","PTEN","RB1","EGFR","BRAF","APC","TTN","PIK3CA")) %>%
  arrange(CNVType,Pvalue) %>%
  mutate(SYMBOL=factor(SYMBOL,levels = .$SYMBOL))

p <- ggplot()+
  geom_point(data=Pancancer.A.D.Fisher.FigureData,
             aes(x=CNVType,y=SYMBOL,size=-log10(Pvalue),fill=Enrich),
             shape=21,color="black",stroke=1)+
  ggthemes::theme_few()+
  scale_size_continuous("-log10(Fisher's P)",range = c(1,5))+
  theme(#legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12),#,angle = 90,vjust = 0.5,hjust = 1
        axis.text.y = element_text(size=12),
        #axis.ticks.x = element_blank(),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  #scale_color_manual("Risk score",values = brewer.pal(9,"Set1")[c(1:2)])+
  #scale_fill_manual("Risk score",values = brewer.pal(9,"Set1")[c(1:3)])
  geom_vline(aes(xintercept=-log10(0.05)),color="black",linetype="dashed")+
  #ggsci::scale_color_nejm()
  ggthemes::scale_color_stata(scheme = "s1rcolor")+
  #labs(y="")
  ggsci::scale_fill_nejm()+
  labs(x="",y="",fill="Higher in group")+
  scale_x_discrete(label=c("Amp.","Del."))
ggsave(p,filename = "RiskGenome.CNV.Deletion.Amplification.Fisher.Sig.Pancancer.pdf",height = 2.5,width = 3.5)

Pancancer.A.D.Fisher.Count <- Pancancer.A.D.Fisher %>%
  group_by(Enrich,CNVType) %>%
  summarise(Count=n())
write.csv(Pancancer.A.D.Fisher.Count,file = "RiskGenome/RiskGenome.CNV.Deletion.Amplification.Fisher.Sig.Pancancer.Count.csv",row.names = F)

p2 <- ggplot(Pancancer.A.D.Fisher.Count,aes(x=Enrich,y=Count,fill=CNVType))+
  geom_col()+
  scale_y_log10()+
  ggthemes::theme_few()+
  #scale_size_continuous("-log10(Fisher's P)",range = c(1,5))+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=12),#,angle = 90,vjust = 0.5,hjust = 1
    axis.text.y = element_text(size=12),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  #scale_color_manual("Risk score",values = brewer.pal(9,"Set1")[c(1:2)])+
  #scale_fill_manual("Risk score",values = brewer.pal(9,"Set1")[c(1:3)])
  #geom_vline(aes(xintercept=-log10(0.05)),color="black",linetype="dashed")+
  #ggsci::scale_color_nejm()
  ggthemes::scale_color_stata(scheme = "s1rcolor")+
  #labs(y="")
  ggsci::scale_fill_nejm()+
  labs(x="",y="Frequency",fill="CNV type")+
  scale_x_discrete(label=c("High risk","Low risk"))
ggsave(p2,filename = "RiskGenome.CNV.Deletion.Amplification.Fisher.Sig.Pancancer.Count.pdf",height = 2.5,width = 3.5)
#### Figure => Cancer => Amplification+Deletion => IntOGen Driver gene ####
setwd("K:/TCGA/Anlysis/BCAA")
library(tidyverse)
Drivergene <- read.csv("IntOGen-DriverGenes.tsv",sep = "\t")
RiskScore.Cancer <- read.csv(file = "RiskGroup.COX.TCGA.RiskScore.csv")
Pancancer.A.CNV <- data.frame()
Pancancer.D.CNV <- data.frame()
for (project in unique(RiskScore.Cancer$Cancer)) {
  Pancancer.D.CNV <- read.csv(file.path("RiskGenome/",paste("RiskGenome.CNV.Deletion.Fisher.",project,".csv",sep = ""))) %>%
    rbind.data.frame(Pancancer.D.CNV)
  Pancancer.A.CNV <- read.csv(file.path("RiskGenome/",paste("RiskGenome.CNV.Amplification.Fisher.",project,".csv",sep = ""))) %>%
    rbind.data.frame(Pancancer.A.CNV)
}

Cancer.A.CNV <- Pancancer.A.CNV %>%
  filter(Fisher.P <= 0.05) %>%
  mutate(Enrich=if_else(High.CNV > Low.CNV,"High","Low"))
write.csv(Cancer.A.CNV,file = "RiskGenome/RiskGenome.CNV.Amplification.Fisher.Cancer.csv",row.names = F)

Cancer.D.CNV <- Pancancer.D.CNV %>%
  filter(Fisher.P <= 0.05) %>%
  mutate(Enrich=if_else(High.CNV > Low.CNV,"High","Low"))
write.csv(Cancer.D.CNV,file = "RiskGenome/RiskGenome.CNV.Deletion.Fisher.Cancer.csv",row.names = F)

Cancer.CNV <- Cancer.A.CNV %>%
  rbind.data.frame(Cancer.D.CNV) %>%
  filter(SYMBOL %in% Drivergene$Symbol[1:60]) %>%
  #group_by(Cancer,CNVType) %>%
  #top_n(-2,Fisher.P) %>%
  mutate(Cancer=factor(Cancer,levels = rev(sort(unique(.$Cancer))))) %>%
  #arrange(Cancer,CNVType,Fisher.P) %>%
  mutate(SYMBOL=factor(SYMBOL,levels = Drivergene$Symbol[1:60]))

library(RColorBrewer)
p2 <- ggplot(Cancer.CNV,aes(x=SYMBOL,y=Cancer,
                            color=Enrich,
                            size=-log10(Fisher.P),
                            shape=CNVType))+
  geom_point(stroke=1.2)+
  scale_shape_manual(values=c(1,4))+
  ggthemes::theme_few()+
  scale_size_continuous("-log10(Fisher's P)",range = c(1,4))+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9,angle = 90,vjust = 0.3,hjust = 1),
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  #scale_color_manual("Risk score",values = brewer.pal(9,"Set1")[c(1:2)])+
  #scale_fill_manual("Risk score",values = brewer.pal(9,"Set1")[c(1:3)])
  #geom_vline(aes(xintercept=-log10(0.05)),color="black",linetype="dashed")+
  ggsci::scale_color_nejm()+
  #ggthemes::scale_color_stata(scheme = "s1rcolor")+
  #scale_color_gradientn("-log10(Fisher's P)",colours = brewer.pal(11,"RdBu")[4:1])+
  labs(y="",x="")
#p2.F <- p2 %>% aplot::insert_top(p2.Top,height = 0.2)
ggsave(p2,filename = "RiskGenome.CNV.Amplification.Fisher.Cancer.Examples.pdf",height = 4.5,width = 12)

######## Risk score => ImmuneGenes ###########
setwd("K:/TCGA/Anlysis/BCAA")
library(tidyverse)
RiskScore.Cancer <- read.csv(file = "RiskGroup.COX.TCGA.RiskScore.csv")
PhenoType <- read.csv("../../mRNA.Phenotype.csv")
PhenoType <- PhenoType %>% filter(Sample.Type == "Primary Tumor" | 
                                    Sample.Type == "Solid Tissue Normal" | 
                                    Sample.Type == "Primary Blood Derived Cancer - Peripheral Blood") %>% 
  dplyr::select(Project.ID,Sample.ID,Sample.Type) %>%
  mutate(Type=ifelse(Sample.Type=="Solid Tissue Normal","Normal","Tumor")) %>%
  mutate(Cancer=str_remove(Project.ID,"TCGA-")) %>%
  filter(Type=="Tumor")
ICI.P.G.E <- data.frame(Protein = c("PD-1","PD-L1","PD-L2","CTLA-4","TIGIT","LAG-3","TIM-3","ARG1","CD163","CD68","FOXP3"),
                        SYMBOL = c("PDCD1","CD274","CD273","CTLA4","TIGIT","LAG3","HAVCR2","ARG1","CD163","CD68","FOXP3"),
                        ENSEMBL = c("ENSG00000188389","ENSG00000120217","ENSG00000197646",
                                      "ENSG00000163599","ENSG00000181847","ENSG00000089692",
                                      "ENSG00000135077","ENSG00000118520","ENSG00000177575",
                                      "ENSG00000129226","ENSG00000049768"))

dirs <- list.dirs("../../Clinical.XENA/",full.names = T)
dirs <- dirs[str_detect(dirs,"//")]

Exp.Risk <- data.frame()
for (dir in dirs) {
  project <- str_remove(dir,".*/")
  print(project)
  phenotype <- read.csv(paste(dir,"/TCGA-",project,".GDC_phenotype.tsv",sep = ""),sep = "\t",header = T) %>% 
    dplyr::filter(sample_type.samples == "Primary Tumor" | sample_type.samples == "Primary Blood Derived Cancer - Peripheral Blood")
  TPM <- read.csv(paste("../../mRNA.PanCancer.Exp/TCGA-",project,".mRNA.TPM.csv",sep = ""),check.names = F,row.names = 1)
  TPM.log <- log2(TPM+1)

  Exp.Risk <- TPM.log[ICI.P.G.E$ENSEMBL,] %>%
    data.frame(check.names = F) %>%
    rownames_to_column("ENSEMBL") %>%
    merge(ICI.P.G.E,by="ENSEMBL") %>%
    dplyr::select(-ENSEMBL,-Protein) %>%
    column_to_rownames("SYMBOL") %>%
    t() %>%
    data.frame(check.names = F) %>%
    rownames_to_column("Sample.ID") %>%
    filter(Sample.ID %in% PhenoType$Sample.ID) %>%
    mutate(Cancer=project) %>%
    rbind.data.frame(Exp.Risk)
}
Exp.Risk2 <- merge(Exp.Risk,RiskScore.Cancer,by.x="Sample.ID",by.y="SampleID")
write.csv(Exp.Risk2,file = "RiskGroup.COX.TCGA.RiskScore.ImmuneGenes.csv",row.names = F)
#### Figure => Correlation ####
setwd("K:/TCGA/Anlysis/BCAA")
library(tidyverse)
Exp.Risk2 <- read.csv("RiskGroup.COX.TCGA.RiskScore.ImmuneGenes.csv")
### PDCD1
library(RColorBrewer)
Cancer.Col <- c(brewer.pal(9,"Set1"),
                brewer.pal(8,"Set2")[1:7],
                brewer.pal(12,"Set3")[1:3],
                brewer.pal(8,"Dark2")[6],
                brewer.pal(12,"Set3")[4:12],
                brewer.pal(8,"Dark2")[1:4])
colors <- Cancer.Col[1:23]
names(colors)=sort(unique(Exp.Risk2$Cancer.x))
Spearman.Data <- data.frame()
for (project in unique(Exp.Risk2$Cancer.x)) {
  middata <- Exp.Risk2 %>% 
    filter(Cancer.x == project) %>%
    filter(SYMBOL=="PDCD1")
  test <- cor.test(middata[["TPM"]],
                   middata$total_risk_score,
                   method = "spearman",
                   exact = F)
  Spearman.Data <- data.frame(Rho=test$estimate,
                              Pvalue=test$p.value,
                              Cancer=project) %>%
    rbind.data.frame(Spearman.Data)
}
write.csv(Spearman.Data,file = "RiskGroup.COX.TCGA.RiskScore.PDCD1.RiskScore.Spearman.csv",row.names = F)
### Figure
suppressMessages(library(rstatix))
Spearman.Data <- Spearman.Data %>% add_significance("Pvalue")

p1 <- ggplot()+
  geom_segment(data = Spearman.Data,
               aes(x=Cancer,xend=Cancer,y=0,yend=Rho),
               color="black")+
  geom_point(data = Spearman.Data,
             aes(x=reorder(Cancer,Rho),y=Rho,fill=Cancer,size=-log10(Pvalue)),
             shape=21,color="white",stroke=1.1)+
  geom_point(data = Spearman.Data %>% filter(Pvalue <= 0.05),
             aes(x=reorder(Cancer,Rho),y=Rho,fill=Cancer,size=-log10(Pvalue)),
             shape=21,color="black",stroke=1.1)+
  scale_size_continuous(range = c(2,5))+
  #geom_text(data=Spearman.Data,aes(x=Cancer,y=Rho,label=Pvalue.signif))+
  ggthemes::theme_few()+
  #scale_size_continuous("-log10(Fisher's P)",range = c(1,4))+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9,angle = 90,vjust = 0.5,hjust = 1),
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  #scale_color_manual("Risk score",values = brewer.pal(9,"Set1")[c(1:2)])+
  #scale_fill_manual("Risk score",values = brewer.pal(9,"Set1")[c(1:3)])
  #geom_vline(aes(xintercept=-log10(0.05)),color="black",linetype="dashed")+
  #ggsci::scale_color_nejm()+
  #ggthemes::scale_color_stata(scheme = "s1rcolor")+
  #scale_color_gradientn("-log10(Fisher's P)",colours = brewer.pal(11,"RdBu")[4:1])+
  labs(y="Spearman'CC",x="")+
  scale_fill_manual(values = colors)+
  geom_hline(yintercept = c(0),linetype=2)+
  scale_y_continuous(limits = c(-0.3,0.65))
ggsave(p1,filename = "RiskGroup.COX.TCGA.RiskScore.PDCD1.Spearman.pdf",height = 2,width = 7)
ggsave(p1,filename = "RiskGroup.COX.TCGA.RiskScore.PDCD1.Spearman.Legend.pdf",height = 10,width = 7)
### FOXP3
Spearman.Data <- data.frame()
for (project in unique(Exp.Risk2$Cancer.x)) {
  middata <- Exp.Risk2 %>% 
    filter(Cancer.x == project) %>%
    filter(SYMBOL=="FOXP3")
  test <- cor.test(middata[["TPM"]],
                   middata$total_risk_score,
                   method = "spearman",
                   exact = F)
  Spearman.Data <- data.frame(Rho=test$estimate,
                              Pvalue=test$p.value,
                              Cancer=project) %>%
    rbind.data.frame(Spearman.Data)
}
write.csv(Spearman.Data,file = "RiskGroup.COX.TCGA.RiskScore.FOXP3.RiskScore.Spearman.csv",row.names = F)
### Figure
suppressMessages(library(rstatix))
Spearman.Data <- Spearman.Data %>% add_significance("Pvalue")

p1 <- ggplot()+
  geom_segment(data = Spearman.Data,
               aes(x=Cancer,xend=Cancer,y=0,yend=Rho),
               color="black")+
  geom_point(data = Spearman.Data,
             aes(x=reorder(Cancer,Rho),y=Rho,fill=Cancer,size=-log10(Pvalue)),
             shape=21,color="white",stroke=1.1)+
  geom_point(data = Spearman.Data %>% filter(Pvalue <= 0.05),
             aes(x=reorder(Cancer,Rho),y=Rho,fill=Cancer,size=-log10(Pvalue)),
             shape=21,color="black",stroke=1.1)+
  scale_size_continuous(range = c(2,5))+
  #geom_text(data=Spearman.Data,aes(x=Cancer,y=Rho,label=Pvalue.signif))+
  ggthemes::theme_few()+
  #scale_size_continuous("-log10(Fisher's P)",range = c(1,4))+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9,angle = 90,vjust = 0.5,hjust = 1),
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  #scale_color_manual("Risk score",values = brewer.pal(9,"Set1")[c(1:2)])+
  #scale_fill_manual("Risk score",values = brewer.pal(9,"Set1")[c(1:3)])
  #geom_vline(aes(xintercept=-log10(0.05)),color="black",linetype="dashed")+
  #ggsci::scale_color_nejm()+
  #ggthemes::scale_color_stata(scheme = "s1rcolor")+
  #scale_color_gradientn("-log10(Fisher's P)",colours = brewer.pal(11,"RdBu")[4:1])+
  labs(y="Spearman'CC",x="")+
  scale_fill_manual(values = colors)+
  geom_hline(yintercept = c(0),linetype=2)+
  scale_y_continuous(limits = c(-0.3,0.65))
ggsave(p1,filename = "RiskGroup.COX.TCGA.RiskScore.FOXP3.Spearman.pdf",height = 2,width = 7)
ggsave(p1,filename = "RiskGroup.COX.TCGA.RiskScore.FOXP3.Spearman.Legend.pdf",height = 10,width = 7)
#### Figure ####
setwd("K:/TCGA/Anlysis/BCAA")
library(tidyverse)
Exp.Risk2 <- read.csv("RiskGroup.COX.TCGA.RiskScore.ImmuneGenes.csv")
library(ggpubr)
Exp.Risk2 <- Exp.Risk2 %>%
  gather(SYMBOL,TPM,FOXP3:CD273) %>%
  group_by(Cancer.x) %>%
  mutate(cutpoint=median(total_risk_score)) %>%
  mutate(Group = if_else(total_risk_score >= cutpoint,"High","Low")) %>%
  ungroup()

data2 <- compare_means(formula = TPM~Group,data = Exp.Risk2,method = "wilcox.test",
                       group.by=c("Cancer.x","SYMBOL"))
write.csv(data2,file="RiskGroup.COX.TCGA.RiskScore.ImmuneGenes.Wilcoxon.csv",row.names = F)

Exp.Risk3 <- Exp.Risk2 %>%
  group_by(Cancer.x,SYMBOL,Group) %>%
  summarise(MeanValue=mean(TPM)) %>%
  ungroup() %>% data.frame(check.names = F) %>%
  merge(data2 %>% dplyr::select(Cancer.x,SYMBOL,p.signif),
        by=c("Cancer.x","SYMBOL")) %>%
  mutate(p.signif=if_else(p.signif=="ns","",p.signif)) %>%
  spread(Group,MeanValue) %>%
  mutate(Enrich=if_else(p.signif == "","n.s.",if_else(High>Low,"High","Low"))) %>%
  gather(Group,Score,High:Low) %>%
  arrange(Cancer.x,SYMBOL) %>%
  group_by(Cancer.x,SYMBOL) %>%
  mutate(Enrich2=if_else(Enrich=="n.s.","n.s.",if_else(Enrich==Group,"Higher","Lower"))) %>%
  mutate(Cancer.x=factor(Cancer.x,levels = sort(unique(.$Cancer.x)))) %>%
  mutate(p.signif=factor(p.signif,levels = c("*","**","***","****",""))) %>%
  mutate(SYMBOL = factor(SYMBOL,levels = rev(unique(sort(.$SYMBOL))))) %>%
  mutate(Xaxis=paste(Cancer.x,Group,sep = "."))

p1 <- ggplot(Exp.Risk3,aes(x=Xaxis,y=SYMBOL))+
  geom_point(aes(size=Score,color=p.signif,shape=Enrich2),
             stroke=1)+
  scale_size_continuous("Mean of log2(TPM+1)",range = c(1,4))+
  #facet_wrap(~Attribute,ncol = 11)+
  labs(x="",y="")+
  #facet_wrap(~Attribute,scale="free_x",ncol = 15)+
  ggthemes::theme_few()+
  theme(#legend.position = "top",
    axis.text = element_text(color = "black"),
    axis.text.x = element_blank(), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    axis.ticks.x = element_blank(),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  scale_color_manual("P value",values = c(brewer.pal(9,"Set1")[1:4],"#F7F7F7"))+
  scale_shape_manual("Score",values = c(1,3,4))+
  geom_vline(xintercept = seq(0.5,44.5,2),linetype=2)
#scale_fill_manual("Risk score",values = brewer.pal(9,"Set1")[c(1:2)])

p2 <- ggplot(Exp.Risk3 %>% data.frame() %>%
               dplyr::select(Xaxis,Group) %>%
               unique(),
             aes(x=Xaxis,y=1,fill=Group))+
  geom_col(width = 1)+
  theme_void()+
  scale_fill_manual("Risk score",values = c("#fe5f55","#004d61"))


p3 <- ggplot(Exp.Risk3 %>% data.frame() %>%
               dplyr::select(Xaxis,Cancer.x) %>%
               unique(),
             aes(x=Xaxis,y=1,fill=Cancer.x))+
  geom_col(width = 1)+
  theme_void()+
  scale_fill_manual("Cell type",values = Cancer.Col)

p123 <- p1 %>% aplot::insert_top(p2,height = 0.05) %>%
  aplot::insert_top(p3,height = 0.05)
ggsave(p123,filename = "RiskGroup.COX.TCGA.RiskScore.ImmuneGene.FigureData.pdf",height = 3,width = 12)
ggsave(p123,filename = "RiskGroup.COX.TCGA.RiskScore.ImmuneGene.FigureData.Legend.pdf",height = 16,width = 10)




ggsave(p,filename = "RiskGroup.COX.TCGA.RiskScore.Estimate.All.pdf",height = 4.5,width = 7)




######## Risk score => Pre+Post => NR+R ###########
setwd("K:/TCGA/Anlysis/BCAA")
library(tidyverse)
library(survival)
library(survminer)
LASSO.Feature <- read.csv("Risk.LASSO.COX.csv")
BCAA.ID <- read.csv("BCAA.ID.csv")
Feature.Count <- LASSO.Feature %>%
  filter(Count >= 50) %>%
  group_by(Cancer) %>%
  summarise(Count=n()) %>%
  filter(Count >= 2)

Feature.Filter <- LASSO.Feature %>%
  filter(Count >= 50) %>%
  filter(Cancer %in% Feature.Count$Cancer) %>%
  merge(BCAA.ID,by.x = "Feature",by.y = "ENSEMBL")

Cal_function <- function(dataset,feature,project,timeVector){
  return.list <- list()
  formula <- as.formula(paste0("Surv(OS.time,OS)~",
                               paste(feature,collapse = "+")))
  dataset$OS <- as.numeric(as.double(dataset$OS))
  dataset$OS.time <- as.numeric(as.double(dataset$OS.time))
  multi_variate_cox <- coxph(formula,data=dataset)
  dataset$total_risk_score <- predict(multi_variate_cox,
                                      type="risk",
                                      newdata=dataset[,feature])
  return.list[[1]] = dataset
  COX.Summary <- summary(multi_variate_cox)
  return.list[[2]] = COX.Summary$concordance[1]
  library(timeROC)
  with(dataset,
       ROC_riskscore <<- timeROC(T = OS.time,
                                 delta = OS,
                                 marker = total_risk_score,
                                 cause = 1,
                                 weighting = "marginal",
                                 times = timeVector,# c(12,36,60),
                                 ROC = TRUE,
                                 iid = TRUE)
  )
  return.list[[3]] = ROC_riskscore
  
  res.cut <- surv_cutpoint(dataset,
                           time = paste("OS",".time",sep = ""),
                           event = "OS",
                           variables = "total_risk_score")
  print(res.cut)
  dat <- surv_categorize(res.cut) #%>% data.frame() %>% mutate(Sample=SNHGs.ICI.data$Sample.ID)
  fit <- survfit(Surv(OS.time, OS) ~ total_risk_score,
                 data = dat)
  index="OS"
  diff=survdiff(as.formula(paste("Surv(",index,".time,",index,") ~ ","total_risk_score",sep = "")),data = dat)
  pValue=1-pchisq(diff$chisq,df=1)
  return.list[[4]] = pValue
  library(ggthemes)
  theme.sur <- ggthemes::theme_few()+
    theme(legend.position = "none",
          axis.text = element_text(color = "black"),
          axis.text.x = element_text(size=15),
          axis.text.y = element_text(size=15),
          axis.title = element_text(size=18),
          panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))
  
  p<-ggsurvplot(fit,pval =TRUE, data = dat, 
                surv.median.line = "hv",
                legend.title = "Risk",
                conf.int.style = "step",
                xlab = "Time in days",
                #break.time.by = 500,
                risk.table = "abs_pct",
                #risk.table.y.text.col = T,
                #risk.table.y.text = FALSE,
                legend.labs = c("High", "Low"),
                #pval = TRUE,
                conf.int = TRUE,
                palette = "Set1",
                #ggtheme = theme.sur,
                risk.table.y.text.col = T,
                risk.table.y.text = F,
                ggtheme = theme.sur)
  print(p)
  #dev.off()
  library(export)
  graph2pdf(file=paste("RiskGroup/Validation.OS.",project,".RiskScore.BestCutoff.KM.pdf",sep = ''),height = 5,width = 4)
  names(return.list) <- c("Dataset","C.index","timeROC","KM.Pvalue")
  return(return.list)
}

###### TIGER => GBM-PRJNA482620 #### 
BCAA.ID <- read.csv("BCAA.ID.csv")
GBM.Feature <- Feature.Filter %>% filter(Cancer == "LGG")
Exp <- read.csv("../../TIGER.Database/GBM-PRJNA482620.Response.tsv",sep = "\t")
Clinical <- read.csv("../../TIGER.Database/GBM-PRJNA482620.Response (1).tsv",sep = "\t")
#BCAA.ID$SYMBOL %in% Exp$GENE_SYMBOL
Exp.Clinical <- Exp %>% 
  filter(GENE_SYMBOL %in% GBM.Feature$SYMBOL) %>%
  remove_rownames() %>% column_to_rownames("GENE_SYMBOL") %>%
  t() %>% data.frame(check.names = F) %>%
  rownames_to_column("sample_id") %>%
  merge(Clinical,by="sample_id") %>%
  mutate(OS=if_else(vital.status=="Dead",1,0),
         OS.time=as.numeric(overall.survival..days.))
DATA <- Cal_function(Exp.Clinical,GBM.Feature$SYMBOL,"GBM-PRJNA482620",c(365,1080,1800))
saveRDS(DATA,file = "ICI.riskScore.Validation.GBM-PRJNA482620.Rds")
#### Figure ####
p1 <- ggplot(DATA$Dataset %>% filter(total_risk_score<=10),
       aes(x=response,y=total_risk_score,fill=response))+
  geom_violin(cex=1.3)+
  geom_boxplot(width=0.4,color="white")+
  labs(x="",y="Risk score")+
  #facet_wrap(~Attribute,scale="free_x",ncol = 15)+
  ggthemes::theme_few()+
  theme(legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=12), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  stat_compare_means(method = "wilcox.test",label.x = 1.2,label = "p.format",size=5)+
  #scale_x_discrete(labels=c("NR","R"))+
  ggsci::scale_fill_npg()
ggsave(p1,filename = "ICI.riskScore.Validation.GBM-PRJNA482620.Boxplot.pdf",height = 2.5,width = 2.5)

DATA$Dataset %>% mutate(Group=if_else(total_risk_score > median(total_risk_score),"High","Low")) -> dataset2
fisher.test(table(dataset2$Group,dataset2$response))
###### TIGER => Melanoma-GSE91061 ####
SKCM.Feature <- Feature.Filter %>% filter(Cancer == "SKCM")
Exp <- read.csv("../../TIGER.Database/Melanoma-GSE91061.Response.tsv",sep = "\t")
Clinical <- read.csv("../../TIGER.Database/Melanoma-GSE91061.Response (1).tsv",sep = "\t")

Exp.Clinical <- Exp %>% 
  filter(GENE_SYMBOL %in% GBM.Feature$SYMBOL) %>%
  remove_rownames() %>% column_to_rownames("GENE_SYMBOL") %>%
  t() %>% data.frame(check.names = F) %>%
  rownames_to_column("sample_id") %>%
  merge(Clinical,by="sample_id") %>%
  mutate(OS=if_else(vital.status=="Dead",1,0),
         OS.time=as.numeric(overall.survival..days.))
DATA <- Cal_function(Exp.Clinical,SKCM.Feature$SYMBOL,"Melanoma-GSE91061",c(365,1080,1800))
saveRDS(DATA,file = "ICI.riskScore.Validation.Melanoma-GSE91061.Rds")
#### Figure ####
DATA$Dataset %>% mutate(RESPONSE = if_else(response %in% c("CR","PR"),"R",
                                           if_else(response %in% c("PD","SD"),"NR","UNK"))) %>%
  mutate(Group=if_else(total_risk_score > median(total_risk_score),"High","Low")) -> dataset2
  
p1 <- ggplot(dataset2 , #%>% filter(Treatment=="ON") %>% filter(total_risk_score<=10)
             aes(x=RESPONSE,y=total_risk_score,fill=RESPONSE))+
  geom_violin(cex=1.3)+
  geom_boxplot(width=0.2,color="white")+
  labs(x="",y="Risk score")+
  #facet_wrap(~Attribute,scale="free_x",ncol = 15)+
  ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12), #,angle = 90,vjust = 0.3,hjust = 1
        axis.text.y = element_text(size=9),
        #axis.ticks.x = element_blank(),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  stat_compare_means(comparisons = list(c("NR","R")),label.x = 1,label = "p.format",size=5,
                      method="wilcox.test")+
  #scale_x_discrete(labels=c("NR","R"))+
  ggsci::scale_fill_npg()
ggsave(p1,filename = "ICI.riskScore.Validation.Melanoma-GSE91061.Boxplot.pdf",height = 2.5,width = 2.5)

fisher.test(table(dataset2$Group,dataset2$RESPONSE))

p2 <- ggplot(dataset2 , #%>% filter(Treatment=="ON") %>% filter(total_risk_score<=10)
             aes(x=Treatment,y=total_risk_score,fill=Treatment))+
  geom_violin(cex=1.3)+
  geom_boxplot(width=0.2,color="white")+
  labs(x="",y="Risk score")+
  #facet_wrap(~Attribute,scale="free_x",ncol = 15)+
  ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12), #,angle = 90,vjust = 0.3,hjust = 1
        axis.text.y = element_text(size=9),
        #axis.ticks.x = element_blank(),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  stat_compare_means(label.x = 1,label = "p.format",size=5,
                     method="wilcox.test")+
  #scale_x_discrete(labels=c("NR","R"))+
  ggsci::scale_fill_npg()
ggsave(p2,filename = "ICI.riskScore.Validation.Melanoma-GSE91061.Boxplot.PREON.pdf",height = 2.5,width = 2.5)

###### TIGER => Melanoma-phs000452 #### 
SKCM.Feature <- Feature.Filter %>% filter(Cancer == "SKCM")
Exp <- read.csv("../../TIGER.Database/Melanoma-phs000452.Response.tsv",sep = "\t")
Clinical <- read.csv("../../TIGER.Database/Melanoma-phs000452.Response (1).tsv",sep = "\t")

Exp.Clinical <- Exp %>% 
  filter(GENE_SYMBOL %in% SKCM.Feature$SYMBOL) %>%
  remove_rownames() %>% column_to_rownames("GENE_SYMBOL") %>%
  t() %>% data.frame(check.names = F) %>%
  rownames_to_column("sample_id") %>%
  merge(Clinical,by="sample_id") %>%
  mutate(OS=if_else(vital.status=="Dead",1,0),
         OS.time=as.numeric(overall.survival..days.))
DATA <- Cal_function(Exp.Clinical,SKCM.Feature$SYMBOL,"Melanoma-phs000452",c(365,1080,1800))
saveRDS(DATA,file = "ICI.riskScore.Validation.Melanoma-phs000452.Rds")
#### Figure ####
DATA$Dataset %>% #mutate(RESPONSE = if_else(response %in% c("CR","PR"),"R",if_else(response %in% c("PD","SD"),"NR","UNK"))) %>%
  mutate(Group=if_else(total_risk_score > median(total_risk_score),"High","Low")) -> dataset2

p1 <- ggplot(dataset2 %>% filter(Treatment=="EDT"), # %>% filter(total_risk_score<=10)
             aes(x=response_NR,y=total_risk_score,fill=response_NR))+
  geom_violin(cex=1.3)+
  geom_boxplot(width=0.2,color="white")+
  labs(x="",y="Risk score")+
  #facet_wrap(~Attribute,scale="free_x",ncol = 15)+
  ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12), #,angle = 90,vjust = 0.3,hjust = 1
        axis.text.y = element_text(size=9),
        #axis.ticks.x = element_blank(),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  stat_compare_means(comparisons = list(c("N","R")),label.x = 1,label = "p.format",size=5,
                     method="wilcox.test")+
  #scale_x_discrete(labels=c("NR","R"))+
  ggsci::scale_fill_npg()
ggsave(p1,filename = "ICI.riskScore.Validation.Melanoma-phs000452.Boxplot.pdf",height = 2.5,width = 2.5)

fisher.test(table(dataset2$Group,dataset2$response_NR))

###### TIGER => Melanoma-PRJEB23709 #### 
SKCM.Feature <- Feature.Filter %>% filter(Cancer == "SKCM")
Exp <- read.csv("../../TIGER.Database/Melanoma-PRJEB23709.Response.tsv",sep = "\t")
Clinical <- read.csv("../../TIGER.Database/Melanoma-PRJEB23709.Response (1).tsv",sep = "\t")

Exp.Clinical <- Exp %>% 
  filter(GENE_SYMBOL %in% SKCM.Feature$SYMBOL) %>%
  remove_rownames() %>% column_to_rownames("GENE_SYMBOL") %>%
  t() %>% data.frame(check.names = F) %>%
  rownames_to_column("sample_id") %>%
  merge(Clinical,by="sample_id") %>%
  mutate(OS=if_else(vital.status=="Dead",1,0),
         OS.time=as.numeric(overall.survival..days.))
DATA <- Cal_function(Exp.Clinical,SKCM.Feature$SYMBOL,"Melanoma-PRJEB23709",c(365,1080,1800))
saveRDS(DATA,file = "ICI.riskScore.Validation.Melanoma-PRJEB23709.Rds")
#### Figure ####
DATA$Dataset %>% #mutate(RESPONSE = if_else(response %in% c("CR","PR"),"R",if_else(response %in% c("PD","SD"),"NR","UNK"))) %>%
  mutate(Group=if_else(total_risk_score > median(total_risk_score),"High","Low")) -> dataset2

p1 <- ggplot(dataset2 , #%>% filter(Treatment=="ON") %>% filter(total_risk_score<=10)
             aes(x=response_NR,y=total_risk_score,fill=response_NR))+
  geom_violin(cex=1.3)+
  geom_boxplot(width=0.2,color="white")+
  labs(x="",y="Risk score")+
  #facet_wrap(~Attribute,scale="free_x",ncol = 15)+
  ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12), #,angle = 90,vjust = 0.3,hjust = 1
        axis.text.y = element_text(size=9),
        #axis.ticks.x = element_blank(),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  stat_compare_means(comparisons = list(c("N","R")),label.x = 1,label = "p.format",size=5,
                     method="wilcox.test")+
  #scale_x_discrete(labels=c("NR","R"))+
  ggsci::scale_fill_npg()
ggsave(p1,filename = "ICI.riskScore.Validation.Melanoma-phs000452.Boxplot.pdf",height = 2.5,width = 2.5)

fisher.test(table(dataset2$Group,dataset2$response_NR))

###### TIGER => RCC-Braun_2020 #####
KIRC.Feature <- Feature.Filter %>% filter(Cancer == "KIRC")
Exp <- read.csv("../../TIGER.Database/RCC-Braun_2020.Response.tsv",sep = "\t")
Clinical <- read.csv("../../TIGER.Database/RCC-Braun_2020.Response (1).tsv",sep = "\t")
Exp.Clinical <- Exp %>% 
  filter(GENE_SYMBOL %in% KIRC.Feature$SYMBOL) %>%
  remove_rownames() %>% column_to_rownames("GENE_SYMBOL") %>%
  t() %>% data.frame(check.names = F) %>%
  rownames_to_column("sample_id") %>%
  merge(Clinical,by="sample_id") %>%
  mutate(OS=if_else(vital.status=="Dead",1,0),
         OS.time=as.numeric(overall.survival..days.))
DATA <- Cal_function(Exp.Clinical,KIRC.Feature$SYMBOL,"RCC-Braun_2020",c(365,1080,1800))
saveRDS(DATA,file = "ICI.riskScore.Validation.RCC-Braun_2020.Rds")

DATA$Dataset %>% #mutate(RESPONSE = if_else(response %in% c("CR","PR"),"R",if_else(response %in% c("PD","SD"),"NR","UNK"))) %>%
  mutate(Group=if_else(total_risk_score > median(total_risk_score),"High","Low")) -> dataset2

p1 <- ggplot(dataset2 , #%>% filter(Treatment=="ON") %>% filter(total_risk_score<=10)
             aes(x=response_NR,y=total_risk_score,fill=response_NR))+
  geom_violin(cex=1.3)+
  geom_boxplot(width=0.2,color="white")+
  labs(x="",y="Risk score")+
  #facet_wrap(~Attribute,scale="free_x",ncol = 15)+
  ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12), #,angle = 90,vjust = 0.3,hjust = 1
        axis.text.y = element_text(size=9),
        #axis.ticks.x = element_blank(),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  stat_compare_means(comparisons = list(c("N","R")),label.x = 1,label = "p.format",size=5,
                     method="wilcox.test")+
  #scale_x_discrete(labels=c("NR","R"))+
  ggsci::scale_fill_npg()
ggsave(p1,filename = "ICI.riskScore.Validation.RCC-Braun_2020.Boxplot.pdf",height = 2.5,width = 2.5)

###### IMvigor210 #####
BLCA.Feature <- Feature.Filter %>% filter(Cancer == "BLCA")

load("../../IMvigor210CoreBiologies/IMvigor210CoreBiologies.Rdata")
# annoData是基因的注释信息
# phenoData是样本信息
# expreSet 是表达矩阵，行是基因31286，列是样本348
SYMBOLS <- annoData %>% 
  filter(symbol %in% KIRC.Feature$SYMBOL) %>%
  dplyr::select(entrez_id,symbol)
Exp.Clinical <- expreSet[SYMBOLS$entrez_id,] %>% data.frame(check.names = F) %>%
  rownames_to_column("entrez_id") %>%
  merge(SYMBOLS,by="entrez_id") %>%
  dplyr::select(-entrez_id) %>%
  remove_rownames() %>% column_to_rownames("symbol") %>%
  t() %>% data.frame(check.names = F) %>%
  rownames_to_column("sample_id") %>%
  merge(phenoData %>% rownames_to_column("sample_id"),by="sample_id") %>%
  mutate(OS=censOS,
         OS.time=as.numeric(os))

DATA <- Cal_function(Exp.Clinical,BLCA.Feature$SYMBOL,"IMvigor210",c(1,3,5))
saveRDS(DATA,file = "ICI.riskScore.Validation.IMvigor210.Rds")

DATA$Dataset %>% #mutate(RESPONSE = if_else(response %in% c("CR","PR"),"R",if_else(response %in% c("PD","SD"),"NR","UNK"))) %>%
  mutate(Group=if_else(total_risk_score > 0.9881938,"High","Low")) %>%
  mutate(RESPONSE=if_else(`Best Confirmed Overall Response` %in% c("CR","PR"),"R",
                          if_else(`Best Confirmed Overall Response` %in% c("PD","SD"),"NR","UNK"))) -> dataset2


p1 <- ggplot(dataset2 %>% filter(RESPONSE != "UNK") , #%>% filter(Treatment=="ON") %>% filter(total_risk_score<=10)
             aes(x=RESPONSE,y=total_risk_score,fill=RESPONSE))+
  geom_violin(cex=1.3)+
  geom_boxplot(width=0.2,color="white")+
  labs(x="",y="Risk score")+
  #facet_wrap(~Attribute,scale="free_x",ncol = 15)+
  ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12), #,angle = 90,vjust = 0.3,hjust = 1
        axis.text.y = element_text(size=9),
        #axis.ticks.x = element_blank(),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  stat_compare_means(comparisons = list(c("NR","R")),label.x = 1,label = "p.format",size=5,
                     method="wilcox.test")+
  #scale_x_discrete(labels=c("NR","R"))+
  ggsci::scale_fill_npg()
ggsave(p1,filename = "ICI.riskScore.Validation.IMvigor210.Boxplot.pdf",height = 2.5,width = 2.5)

dataset3 <- dataset2 %>% filter(RESPONSE != "UNK")
fisher.test(table(dataset3$Group,dataset3$RESPONSE))
fisher.test(matrix(c(80,20,76,24),ncol = 2))

######## Drugs => druggable gene mRNA => 46 genes >= 10 Cancers ###########
setwd("K:/TCGA/Anlysis/BCAA")
library(tidyverse)
Drug.Gene <- data.table::fread("../../DrugsDGIDB/genes.tsv")
GeneType <- read.table("../../../../Metastasis/TissueMicrobiome/Human.GRC38.GeneType.txt",header = T,sep = ",")
GeneType.Unique <- GeneType %>% dplyr::select(Gene.stable.ID,HGNC.symbol,Gene.type,Gene.name,NCBI.gene..formerly.Entrezgene..ID) %>% unique()
GeneType.Unique %>% filter(Gene.type == "protein_coding") %>%
  dplyr::select(Gene.stable.ID,HGNC.symbol) %>%
  filter(HGNC.symbol %in% Drug.Gene$gene_name) %>%
  na.omit() %>%
  mutate(HGNC.symbol=as.character(HGNC.symbol)) %>%
  filter(HGNC.symbol != "") %>%
  magrittr::set_colnames(c("ENSEMBL","SYMBOL")) -> druggene

RiskScore.Cancer <- read.csv(file = "RiskGroup.COX.TCGA.RiskScore.csv")
PhenoType <- read.csv("../../mRNA.Phenotype.csv")
PhenoType <- PhenoType %>% filter(Sample.Type == "Primary Tumor" | 
                                    Sample.Type == "Solid Tissue Normal" | 
                                    Sample.Type == "Primary Blood Derived Cancer - Peripheral Blood") %>% 
  dplyr::select(Project.ID,Sample.ID,Sample.Type) %>%
  mutate(Type=ifelse(Sample.Type=="Solid Tissue Normal","Normal","Tumor")) %>%
  mutate(Cancer=str_remove(Project.ID,"TCGA-")) %>%
  filter(Type=="Tumor")

dirs <- list.dirs("../../Clinical.XENA/",full.names = T)
dirs <- dirs[str_detect(dirs,"//")]

for (dir in dirs) {
  Exp.Spearman <- data.frame()
  project <- str_remove(dir,".*/")
  if (project %in% RiskScore.Cancer$Cancer) {
    print(project)
    phenotype <- read.csv(paste(dir,"/TCGA-",project,".GDC_phenotype.tsv",sep = ""),sep = "\t",header = T) %>% 
      dplyr::filter(sample_type.samples == "Primary Tumor" | sample_type.samples == "Primary Blood Derived Cancer - Peripheral Blood")
    TPM <- read.csv(paste("../../mRNA.PanCancer.Exp/TCGA-",project,".mRNA.TPM.csv",sep = ""),check.names = F,row.names = 1)
    TPM.log <- log2(TPM+1)
    
    Exp.Risk <- TPM.log[intersect(rownames(TPM.log),druggene$ENSEMBL) %>% sort()%>%unique(),] %>%
      data.frame(check.names = F) %>%
      #rownames_to_column("ENSEMBL") %>%
      #merge(druggene,by="ENSEMBL") %>%
      #dplyr::select(-ENSEMBL) %>%
      #column_to_rownames("SYMBOL") %>%
      t() %>%
      data.frame(check.names = F) %>%
      rownames_to_column("Sample.ID") %>%
      filter(Sample.ID %in% PhenoType$Sample.ID) %>%
      #mutate(Cancer=project) %>%
      merge(RiskScore.Cancer,by.x="Sample.ID",by.y="SampleID")
    for (gene in intersect(rownames(TPM.log),druggene$ENSEMBL) %>% sort()%>%unique()) {
      print(gene)
      test <- cor.test(Exp.Risk[[gene]],Exp.Risk$total_risk_score,
                       method = "spearman",exact = F)
      Exp.Spearman <- data.frame(ENSEMBL=gene,
                                 Rho=test$estimate,
                                 Pvalue=test$p.value,
                                 Cancer=project) %>%
        rbind.data.frame(Exp.Spearman)
    }
  }
  write.csv(Exp.Spearman,file = paste("Drugs.DGI.",project,".csv",sep = ""),row.names = F)
}
#Exp.Spearman %>% filter(Cancer!="LIHC") %>% write.csv("Drugs.DGI.Part1.csv",row.names = F)
library(tidyverse)
Exp.Spearman <- data.frame()
files <- list.files(".",pattern = "^Drugs")
files <- files[str_detect(files,"DGI")]
for (file in files) {
  data <- read.csv(file)
  Exp.Spearman <- rbind.data.frame(Exp.Spearman,data)
}
write.csv(Exp.Spearman,"DRUGs.DGI.All.csv",row.names = F)
#### Figure ####
Exp.Spearman <- read.csv("DRUGs.DGI.All.csv")
n=length(unique(Exp.Spearman$Cancer))
dim(Exp.Spearman)
Selec.Gene <- Exp.Spearman %>% 
  filter(Pvalue <= 0.05 & abs(Rho) >= 0.3) %>%
  mutate(Type=if_else(Rho>0,"P","N")) %>%
  filter(! ENSEMBL %in% BCAA.ID$ENSEMBL) %>%
  #filter(Type=="P") %>%
  group_by(ENSEMBL,Type) %>%
  summarise(Count=n()) %>%
  arrange(desc(Count)) %>%
  filter(Count >= 10) %>%
  #top_n(100,Count) %>%
  arrange(Count) %>%
  merge(druggene,by="ENSEMBL")
write.csv(Selec.Gene,"CERES.SelectGene.csv",row.names = F)

p1 <- Selec.Gene %>% 
  #top_n(7,Count) %>%
  filter(SYMBOL %in% CERES.Score1$SYMBOL) %>%
  ggplot(aes(x=Count,y=SYMBOL,fill=SYMBOL))+
  geom_col()+
  labs(x="No. of cancers",y="")+
  #facet_wrap(~Attribute,scale="free_x",ncol = 15)+
  ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9,angle = 90,vjust = 0.5,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
        axis.text.y = element_text(size=9),
        #axis.ticks.x = element_blank(),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  ggthemes::scale_fill_stata()+
  geom_text(aes(label=Count))+
  coord_flip()

ggsave(p1,filename = "DRUGs.DGI.Positive.pdf",width = 3,height = 2)

######## CRISPR dependency score (CERES) ############
CERES.data <- read.csv("../../DepMap/CRISPR_(Project_Score,_CERES).csv")
CERES.Target <- CERES.data %>% dplyr::select(intersect(Selec.Gene$SYMBOL,colnames(CERES.data)))
CERES.Score1 <- CERES.Target %>%
  gather(SYMBOL,Score) %>%
  filter(SYMBOL %in% Selec.Gene$SYMBOL) %>%
  group_by(SYMBOL) %>%
  summarise(Score.M = median(Score)) %>%
  filter(Score.M <= -1)

CERES.Score2 <- CERES.Target %>%
  gather(SYMBOL,Score) %>%
  filter(SYMBOL %in% Selec.Gene$SYMBOL) %>%
  filter(SYMBOL %in% CERES.Score1$SYMBOL)

p2 <- ggplot(CERES.Score2,aes(x=reorder(SYMBOL,Score),y=Score,fill=SYMBOL))+
  geom_boxplot()+
  #geom_col()+
  labs(x="",y="CERES")+
  #facet_wrap(~Attribute,scale="free_x",ncol = 15)+
  ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9,angle = 90,vjust = 0.5,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
        axis.text.y = element_text(size=9),
        #axis.ticks.x = element_blank(),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  ggthemes::scale_fill_stata()+
  geom_hline(yintercept = c(-1),linetype=2)
ggsave(p2,filename = "DRUGs.DGI.Positive.13.SYMBOL.CERES.pdf",width = 4,height = 2)

######## RISK SCORE => DEGs #########
setwd("K:/TCGA/Anlysis/BCAA")
library(tidyverse)
RiskScore.Cancer <- read.csv(file = "RiskGroup.COX.TCGA.RiskScore.csv") %>%
  group_by(Cancer) %>%
  mutate(cutpoint=median(total_risk_score)) %>%
  mutate(Group=if_else(total_risk_score >= cutpoint,"High","Low")) %>%
  ungroup()
PhenoType <- read.csv("../../mRNA.Phenotype.csv")
PhenoType <- PhenoType %>% filter(Sample.Type == "Primary Tumor" | 
                                    Sample.Type == "Solid Tissue Normal" | 
                                    Sample.Type == "Primary Blood Derived Cancer - Peripheral Blood") %>% 
  dplyr::select(Project.ID,Sample.ID,Sample.Type) %>%
  mutate(Type=ifelse(Sample.Type=="Solid Tissue Normal","Normal","Tumor")) %>%
  mutate(Cancer=str_remove(Project.ID,"TCGA-")) %>%
  filter(Type=="Tumor")

dir.create("DEGs.RiskScore")
files <- list.files("../../mRNA.Pancancer.Exp.Count/",pattern = "csv$")
#files <- list.files("../mRNA.Pancancer.Exp.Count/",pattern = "csv$")
for (file in files) {
  project <- file %>%str_remove("PanCancer\\.") %>%
    str_remove("\\..*") %>% str_remove("TCGA-")
  print(project)
  if (project %in% RiskScore.Cancer$Cancer) {
    #Count.data <- read.csv(file = file.path("../../mRNA.Pancancer.Exp.Count/",file),check.names = F,row.names = 1)
    Count.data <- read.csv(file = file.path("../../mRNA.Pancancer.Exp.Count/",file),check.names = F,row.names = 1)
    #filename <- paste("TCGA-",project,".mRNA.TPM.csv",sep = "")
    #TPM.data <- read.csv(file = file.path("../../mRNA.PanCancer.Exp/",filename),check.names = F,row.names = 1)
    #TPM.data <- read.csv(file = file.path("../../mRNA.PanCancer.Exp/",filename),check.names = F,row.names = 1)
    Count.data <- Count.data %>% t() %>% data.frame(check.names = F) %>%
      rownames_to_column("SampleID") %>%
      filter(SampleID %in% RiskScore.Cancer$SampleID) %>%
      column_to_rownames("SampleID") %>%
      t() %>% data.frame(check.names = F)
    countData <- Count.data[rowMeans(Count.data)>1,] 
    colData <- data.frame(SampleID=colnames(Count.data)) %>%
      merge(RiskScore.Cancer %>% dplyr::select(SampleID,Group),by="SampleID") %>%
      column_to_rownames("SampleID")
    colData$Group <- factor(colData$Group,levels = c("High","Low"))
    countData <- countData %>% dplyr::select(rownames(colData))
    
    #condition <- factor(GENE.EXP$Group,levels = c("high","low"))
    #colData <- data.frame(row.names=colnames(countData), condition)
    library(DESeq2)
    dds <- DESeqDataSetFromMatrix(countData = countData, 
                                  colData = colData, 
                                  design = ~ Group)
    dds1 <- DESeq(dds)
    res <- results(dds1,contrast=c("Group","High","Low"))
    res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
    outfilename <- paste(project,".DEG.RiskScore.csv",sep = "")
    write.csv(res1,file = file.path("DEGs.RiskScore/",outfilename))
  }
}

#### Dataset2 ####
Exp.Data <- data.frame()
files <- list.files("DEGs.RiskScore/",pattern = "csv$")
for (file in files) {
  project <- file %>% str_remove_all("\\..*")
  print(project)
  Exp.Data <- read.csv(file.path("DEGs.RiskScore/",file)) %>%
    mutate(Cancer=project) %>%
    rbind.data.frame(Exp.Data)
}

GeneType <- read.table("../../../../Metastasis/TissueMicrobiome/Human.GRC38.GeneType.txt",header = T,sep = ",")
GeneType.Unique <- GeneType %>% dplyr::select(Gene.stable.ID,HGNC.symbol,Gene.type,Gene.name,NCBI.gene..formerly.Entrezgene..ID) %>% unique()
GeneType.Unique %>% #filter(Gene.type == "protein_coding") %>%
  dplyr::select(Gene.stable.ID,HGNC.symbol) %>%
  filter(HGNC.symbol %in% Drug.Gene$gene_name) %>%
  na.omit() %>%
  mutate(HGNC.symbol=as.character(HGNC.symbol)) %>%
  filter(HGNC.symbol != "") %>%
  magrittr::set_colnames(c("ENSEMBL","SYMBOL")) -> druggene

Up.Gene <- Exp.Data %>% filter(padj <= 0.05) %>%
  filter(log2FoldChange >= 1) %>%
  #mutate(Type=if_else(log2FoldChange >= 1,"Up","Down")) %>%
  group_by(X) %>%
  summarise(Count=n()) %>%
  arrange(desc(Count)) %>%
  filter(Count >= 5) %>%
  merge(druggene,by.x="X",by.y="ENSEMBL") %>%
  top_n(100,Count)


Down.Gene <- Exp.Data %>% filter(padj <= 0.05) %>%
  filter(log2FoldChange <= -1) %>%
  #mutate(Type=if_else(log2FoldChange >= 1,"Up","Down")) %>%
  group_by(X) %>%
  summarise(Count=n()) %>%
  arrange(desc(Count)) %>%
  filter(Count >= 5) %>%
  merge(druggene,by.x="X",by.y="ENSEMBL")

write.csv(Up.Gene %>% dplyr::select(SYMBOL),"CMAP.UP.SYMBOL2.csv",row.names = F)

write.csv(Down.Gene %>% dplyr::select(SYMBOL),"CMAP.DOWN.SYMBOL2.csv",row.names = F)

GPL96=read.table("GPL96-57554.txt",  header = TRUE,fill = T,sep = "\t",
                 comment.char = "#",
                 stringsAsFactors = FALSE,
                 quote = "")

Up.Gene.Trans <- GPL96[GPL96$Gene.Symbol %in% Up.Gene$SYMBOL,c(1,11)]
write.csv(Up.Gene.Trans,"CMAP.UP.SYMBOL.csv",row.names = F)
Down.Gene.Trans <- GPL96[GPL96$Gene.Symbol %in% Down.Gene$SYMBOL,c(1,11)]
write.csv(Down.Gene.Trans,"CMAP.DOWN.SYMBOL.csv",row.names = F)

#### Confirm SYMBOLS => SPIED3 + CMAP ####
GeneType <- read.table("../../../../Metastasis/TissueMicrobiome/Human.GRC38.GeneType.txt",header = T,sep = ",")
GeneType.Unique <- GeneType %>% dplyr::select(Gene.stable.ID,HGNC.symbol,Gene.type,Gene.name,NCBI.gene..formerly.Entrezgene..ID) %>% unique()
GeneType.Unique %>% filter(Gene.type == "protein_coding") %>%
  dplyr::select(Gene.stable.ID,HGNC.symbol) %>%
  #filter(HGNC.symbol %in% Drug.Gene$gene_name) %>%
  na.omit() %>%
  mutate(HGNC.symbol=as.character(HGNC.symbol)) %>%
  filter(HGNC.symbol != "") %>%
  magrittr::set_colnames(c("ENSEMBL","SYMBOL")) -> druggene

files <- list.files("DEGs.RiskScore/",pattern = "csv$")
files <- files[str_detect(files,"DEG")]
for (file in files) {
  project <- file %>% str_remove_all("\\..*")
  print(project)
  data1 <- read.csv(file.path("DEGs.RiskScore/",file)) %>%
    mutate(Cancer=project) %>%
    filter(padj <= 0.05) %>%
    filter(abs(log2FoldChange) >= 1) %>%
    merge(druggene,by.x="X",by.y="ENSEMBL") %>% 
    dplyr::select(SYMBOL,log2FoldChange) %>%
    arrange(desc(log2FoldChange)) %>% 
    filter()
  data2 <- read.csv(file.path("DEGs.RiskScore/",file)) %>%
    mutate(Cancer=project) %>%
    filter(padj <= 0.05) %>%
    filter(abs(log2FoldChange) >= 1) %>%
    merge(druggene,by.x="X",by.y="ENSEMBL") %>% 
    dplyr::select(SYMBOL,log2FoldChange) %>%
    arrange(desc(log2FoldChange)) %>% 
    filter()
    write.table(file = file.path("DEGs.RiskScore/",
                               paste(project,".Sig.CMAP.txt",sep = "")),
              row.names = F,col.names = F,sep = "\t",quote = F)
}

#### CMAP => DRUGS ####
setwd("K:/TCGA/Anlysis/BCAA/CMAP")
library(data.table)
library(tidyverse)
files <- list.files("./",recursive = T)
files <- files[str_detect(files,"query_result")]
CMAP.GRUGS <- data.frame()
for (file in files) {
  project = file %>% str_remove("\\..*")
  print(project)
  data <- read.table(file,
                     skip = 2,sep = "\t",header=T)
  data <- data %>% filter(pert_type == "trt_cp") %>%
    mutate(fdr_q_nlog10=as.numeric(fdr_q_nlog10)) %>%
    arrange(desc(fdr_q_nlog10)) %>%
    filter(fdr_q_nlog10 >= 10)
  CMAP.GRUGS <- rbind.data.frame(CMAP.GRUGS,data%>%mutate(Cancer=project))
  write.csv(data,file = paste("CMAP.",project,".csv",sep = ""),row.names = F)
}
Anno <- fread("CMAP.Drugs.Repurposing_Hub_export.txt")
Class <- CMAP.GRUGS %>% #merge(Anno,by.x="pert_iname",by.y="Name",all.x = T) %>%
  group_by(moa) %>%
  summarise(Count=n()) %>%
  arrange(desc(Count)) %>%
  filter(moa != -666) %>%
  top_n(10,Count)
write.csv(CMAP.GRUGS,file = "CMAP.Drugs.csv",row.names = F)

p1 <- ggplot(Class,aes(x=Count,y=reorder(moa,Count),fill=moa))+
  geom_col()+
  ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9,angle = 90,vjust = 0.5,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
        axis.text.y = element_text(size=9),
        #axis.ticks.x = element_blank(),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  ggthemes::scale_fill_stata()+
  geom_text(aes(label=Count),hjust=1)+
  labs(y="",x="Frequence")
ggsave(p1,filename = "CMAP.pdf",height = 3,width = 3)

# c("Fluorouracil", "Cisplatin", "Oxaliplatin", "Sorafenib","Docetaxel", "Irinotecan", "Capecitabine")

###### oncoPredict ##############
setwd("K:/TCGA/Anlysis/BCAA")
library(tidyverse)
library(oncoPredict)
library(data.table)
library(gtools)
library(reshape2)
library(ggpubr)
library(ggplot2)
dir.create("oncoPredict")
options("expressions"=20000)
#memory.limit(size=8000000)
dir1 = "../../oncoPredict/DataFiles/DataFiles/Training Data/"
GDSC2_Expr = readRDS(file=file.path(dir1,'GDSC2_Expr (RMA Normalized and Log Transformed).rds'))
GDSC2_Res = readRDS(file = file.path(dir1,"GDSC2_Res.rds"))
GDSC2_Res <- exp(GDSC2_Res) 

CTRP2_Expr = readRDS(file=file.path(dir1,'CTRP2_Expr (TPM, not log transformed).rds'))
CTRP2_Res = readRDS(file = file.path(dir1,"CTRP2_Res.rds"))
CTRP2_Res <- exp(CTRP2_Res)

#testExpr<- GDSC2_Expr[,sample(1:ncol(GDSC2_Expr),10)]
#colnames(testExpr)=paste0('test',colnames(testExpr))
#testExpr[1:4,1:4] 
RiskScore.Cancer <- read.csv(file = "RiskGroup.COX.TCGA.RiskScore.csv")
PhenoType <- read.csv("../../mRNA.Phenotype.csv")
PhenoType <- PhenoType %>% filter(Sample.Type == "Primary Tumor" | 
                                    Sample.Type == "Solid Tissue Normal" | 
                                    Sample.Type == "Primary Blood Derived Cancer - Peripheral Blood") %>% 
  dplyr::select(Project.ID,Sample.ID,Sample.Type) %>%
  mutate(Type=ifelse(Sample.Type=="Solid Tissue Normal","Normal","Tumor")) %>%
  mutate(Cancer=str_remove(Project.ID,"TCGA-")) %>%
  filter(Type=="Tumor")

dirs <- list.dirs("../../Clinical.XENA/",full.names = T)
dirs <- dirs[str_detect(dirs,"//")]

for (dir in dirs) {
  project <- str_remove(dir,".*/")
  print(project)
  if (project %in% RiskScore.Cancer$Cancer) {
    #phenotype <- read.csv(paste(dir,"/TCGA-",project,".GDC_phenotype.tsv",sep = ""),sep = "\t",header = T) %>% 
    #  dplyr::filter(sample_type.samples == "Primary Tumor" | sample_type.samples == "Primary Blood Derived Cancer - Peripheral Blood")
    TPM <- read.csv(paste("../../mRNA.PanCancer.Exp/TCGA-",project,".mRNA.TPM.csv",sep = ""),check.names = F,row.names = 1)
    #TPM.log <- log2(TPM+1)
    G_list <- read.csv("../LIHC/biomaRt.GeneID.tansfer.csv")
    mRNA.Exp.Target <- TPM %>% #[,intersect(phenotype$submitter_id.samples,colnames(TPM.log))] %>% data.frame(check.names = F)
      t() %>% data.frame(check.names = F) %>%
      rownames_to_column("SampleID") %>%
      filter(SampleID %in% PhenoType$Sample.ID) %>%
      filter(SampleID %in% RiskScore.Cancer$SampleID) %>%
      remove_rownames() %>%
      column_to_rownames("SampleID") %>%
      t() %>% data.frame(check.names = F)
    
    dat.2 <- mRNA.Exp.Target %>% data.frame(check.names = F) %>%
      rownames_to_column("ensembl_gene_id") %>% data.frame(check.names = F) %>%
      merge(.,G_list[,1:2],by="ensembl_gene_id") %>% 
      data.frame(check.names = F) %>%
      na.omit() %>% dplyr::select(-ensembl_gene_id)
    expr_mean=aggregate(.~hgnc_symbol,mean,data=dat.2) #calculate mean for same symbol
    dat.3 <- expr_mean %>% remove_rownames() %>% column_to_rownames("hgnc_symbol")
    if(! file.exists(paste("oncoPredict/GDSC2.",project,".drugSensitivity.csv",sep = ""))){
      calcPhenotype(trainingExprData = GDSC2_Expr,
                    trainingPtype = GDSC2_Res,
                    testExprData = as.matrix(dat.3),
                    batchCorrect = 'eb',  #   "eb" for ComBat  
                    powerTransformPhenotype = TRUE,
                    removeLowVaryingGenes = 0.2,
                    minNumSamples = 10, 
                    printOutput = TRUE, 
                    removeLowVaringGenesFrom = 'homogenizeData' )
      
      testPtype <- fread('./calcPhenotype_Output/DrugPredictions.csv', data.table = F)
      write.csv(testPtype,file = paste("oncoPredict/GDSC2.",project,".drugSensitivity.csv",sep = ""),quote = F)
      
    }
    if(! file.exists(paste("oncoPredict/CTRP2.",project,".drugSensitivity.csv",sep = ""))){
      calcPhenotype(trainingExprData = CTRP2_Expr,
                    trainingPtype = CTRP2_Res,
                    testExprData = as.matrix(dat.3),
                    batchCorrect = 'eb',  #   "eb" for ComBat  
                    powerTransformPhenotype = TRUE,
                    removeLowVaryingGenes = 0.2,
                    minNumSamples = 10, 
                    printOutput = TRUE, 
                    removeLowVaringGenesFrom = 'rawData' )
      
      testPtype <- fread('./calcPhenotype_Output/DrugPredictions.csv', data.table = F)
      write.csv(testPtype,file = paste("oncoPredict/CTRP2.",project,".drugSensitivity.csv",sep = ""),quote = F)
    }
  print("DONE")
  }
}  

calcPhenotype(trainingExprData = GDSC2_Expr,
              trainingPtype = GDSC2_Res,
              testExprData = as.matrix(dat.3),
              batchCorrect = 'eb',  #   "eb" for ComBat  
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 10, 
              printOutput = TRUE, 
              removeLowVaringGenesFrom = 'rawData' )
library(data.table)
testPtype <- fread('./calcPhenotype_Output/DrugPredictions.csv', data.table = F)


###### GDSC2 ######
setwd("K:/TCGA/Anlysis/BCAA")
library(tidyverse)

RiskScore.Cancer <- read.csv(file = "RiskGroup.COX.TCGA.RiskScore.csv") %>%
  group_by(Cancer) %>%
  mutate(cutpoint=median(total_risk_score)) %>%
  mutate(Group=if_else(total_risk_score >= cutpoint,"High","Low")) %>%
  ungroup()

files <- list.files("oncoPredict/GDSC/",pattern = "csv$")
Drugs.Wilcox <- data.frame()
for (file in files) {
  project = (file %>% str_split("\\.") %>% unlist())[2]
  if (project %in% RiskScore.Cancer$Cancer) {
    print(project)
    library(rstatix)
    data <- read.csv(file.path("oncoPredict/GDSC/",file),row.names = 1) %>%
      merge(RiskScore.Cancer,by.x = "V1",by.y="SampleID") %>%
      gather(Drugs,Sensitivity,Camptothecin_1003:JQ1_2172) 
    data2 <- data %>% na.omit() %>%
      group_by(Group,Drugs) %>%
      summarise(S.M = median(Sensitivity)) %>%
      ungroup() %>%
      spread(Group,S.M)
    library(ggpubr)
    Drugs.Wilcox <- compare_means(Sensitivity~Group,data = data,group.by = "Drugs",method = "wilcox") %>%
      mutate(Cancer=project) %>% 
      merge(data2,by="Drugs") %>%
      rbind.data.frame(Drugs.Wilcox)
  }
}
write.csv(Drugs.Wilcox,file = "Pancancer.GDSC2.Drugs.Sensitivity.csv",row.names = F)
#### Figure ####
Drugs.Wilcox2 <- Drugs.Wilcox %>% 
  filter(p.signif != "ns") %>%
  mutate(Label = if_else(High >= Low,"High","Low")) %>%
  mutate(Drugs2 = str_remove_all(Drugs,"_.*")) %>%  #filter(Drugs2 == "Docetaxel")
  group_by(Drugs,Label) %>%
  summarise(Count=n()) %>%
  ungroup() %>%
  arrange(desc(Count)) %>% #filter(Drugs == "Docetaxel")
  filter(Count >= 10)
#Phenotype <- data.frame(Drugs.Wilcox2$Drugs,Rename=c("BRD Inhibitor","PDGFR tyrosine kinase receptor inhibitor",
#                                                     "K-Ras(G12C)-inhibitor-12",
#                                                     "Protein Kinase inhibitor",
#                                                     "Mitochondrial inhibitor",
#                                                     "ATM Kinase Inhibitor",
#                                                     "Protein tyrosine kinases inhibitor",
#                                                     "Unknown",
#                                                     "RAF inhibitor",
#                                                     ))
p1 <- ggplot(data=Drugs.Wilcox %>% mutate(Label = if_else(p.signif=="ns","N.S.",
                                                    if_else(High >= Low,"High","Low"))) %>%
         filter(Drugs %in% Drugs.Wilcox2$Drugs),
       aes(x=Cancer,y=Drugs,color=Label,shape=Label)
       )+
  geom_point(size=2)+
  ggthemes::theme_few()+
  theme(#legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9,angle = 90,vjust = 0.5,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
        axis.text.y = element_text(size=9),
        #axis.ticks.x = element_blank(),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  #ggsci::scale_color_npg()+
  scale_color_manual(values = c("#E64B35FF","#3C5488FF","#B09C85FF"))+
  scale_shape_manual(values = c(1,3,4))+
  labs(x="",y="",color="Drug sensitivity\nfor risk group",
       shape="Drug sensitivity\nfor risk group")
p2 <- ggplot(data=Drugs.Wilcox %>% mutate(Label = if_else(p.signif=="ns","N.S.",
                                                          if_else(High >= Low,"High","Low"))) %>%
               filter(Drugs %in% Drugs.Wilcox2$Drugs) %>%
               group_by(Drugs,Label) %>%
               summarise(Count=n()) %>%
               mutate(Label=factor(Label,levels = rev(c("High","Low","N.S.")))),
             aes(x=Count,y=Drugs,fill=Label))+
  geom_col()+
  scale_fill_manual(values = rev(c("#E64B35FF","#3C5488FF","#B09C85FF")))+
  labs(x="No. of cancer",y="")+
  ggthemes::theme_few()+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9,angle = 90,vjust = 0.5,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))
p<-p1 %>% aplot::insert_right(p2,width = 0.3)
ggsave(p,filename = "Drugs.GDSC2.Grugs.Sensitivity.pdf",height = 5,width = 7)

###### PRISM + CTRP => OncoPredict ##########
setwd("K:/TCGA/Anlysis/BCAA")
library(tidyverse)
library(oncoPredict)
library(data.table)
library(gtools)
library(reshape2)
library(ggpubr)
library(ggplot2)
#library(impute)
#library(SimDesign) # Information used to suppress the output of the susceptibility prediction process
prism.auc <- read.csv("DepMap/Drug_sensitivity_AUC_(PRISM_Repurposing_Secondary_Screen).csv",
                      check.names = F,stringsAsFactors = F,row.names = 1)
ctrp.auc <- read.csv("DepMap/Drug_sensitivity_AUC_(CTD^2).csv",
                     check.names = F,stringsAsFactors = F,row.names = 1)
#celline.anno <- read.csv("DepMap/CBX2 Expression Public 22Q4.csv")
## Remove drugs with more than 20% missing values
prism.auc <- prism.auc[,apply(prism.auc,2,function(x) sum(is.na(x))) < 0.2*nrow(prism.auc)]
ctrp.auc <- ctrp.auc[,apply(ctrp.auc,2,function(x) sum(is.na(x))) < 0.2*nrow(ctrp.auc)]
### Load expression profiles of CCLE cell lines as training set
ccle.exp <- read.csv("DepMap/DepMap.22Q2.TPM.Processed.csv",
                     check.names = F,row.names = 1)

PhenoType <- read.csv("../mRNA.Phenotype.csv")
PhenoType <- PhenoType %>% filter(Sample.Type == "Primary Tumor" |
                                    Sample.Type == "Solid Tissue Normal" |
                                    Sample.Type == "Primary Blood Derived Cancer - Peripheral Blood") %>%
  dplyr::select(Project.ID,Sample.ID,Sample.Type) %>%
  mutate(Type=ifelse(Sample.Type=="Solid Tissue Normal","Normal","Tumor")) %>%
  mutate(Cancer=str_remove(Project.ID,"TCGA-")) %>%
  filter(Type=="Tumor")

dirs <- list.dirs("../Clinical.XENA/",full.names = T)
dirs <- dirs[str_detect(dirs,"//")]
for (dir in dirs) {
  project <- str_remove(dir,".*/")
  print(project)
  TPM <- read.csv(paste("../mRNA.PanCancer.Exp/TCGA-",project,".mRNA.TPM.csv",sep = ""),check.names = F,row.names = 1)
  #TPM.log <- log2(TPM+1)
  G_list <- read.csv("biomaRt.GeneID.tansfer.csv")
  mRNA.Exp.Target <- TPM %>% #[,intersect(phenotype$submitter_id.samples,colnames(TPM.log))] %>% data.frame(check.names = F)
    t() %>% data.frame(check.names = F) %>%
    rownames_to_column("SampleID") %>%
    filter(SampleID %in% PhenoType$Sample.ID) %>%
    #filter(SampleID %in% RiskScore.Cancer$SampleID) %>%
    remove_rownames() %>%
    column_to_rownames("SampleID") %>%
    t() %>% data.frame(check.names = F)
  
  dat.2 <- mRNA.Exp.Target %>% data.frame(check.names = F) %>%
    rownames_to_column("ensembl_gene_id") %>% data.frame(check.names = F) %>%
    merge(.,G_list[,1:2],by="ensembl_gene_id") %>%
    data.frame(check.names = F) %>%
    na.omit() %>% dplyr::select(-ensembl_gene_id)
  expr_mean=aggregate(.~hgnc_symbol,mean,data=dat.2) #calculate mean for same symbol
  dat.3 <- expr_mean %>% remove_rownames() %>% column_to_rownames("hgnc_symbol")
  
  keepgene <- apply(dat.3, 1, mad) > 0.5
  testExpr <- log2(dat.3[keepgene,] + 1)
  
  comgene <- intersect(rownames(ccle.exp),rownames(testExpr)) 
  trainExpr <- as.matrix(ccle.exp[comgene,])
  testExpr <- testExpr[comgene,]
  
  if (! file.exists(paste("OncoPredict/CTRP.AUC.",project,".csv",sep = ""))) {
    trainPtype <- ctrp.auc
    ### The loop is slow, please be patient
    calcPhenotype(trainingExprData = as.matrix(trainExpr),
                  trainingPtype = as.matrix(trainPtype),
                  testExprData = as.matrix(testExpr),  ##"testExprData" must be a matrix.
                  powerTransformPhenotype = F,
                  selection = 2,batchCorrect = "eb",
                  printOutput = T,removeLowVaringGenesFrom = "rawData",
                  minNumSamples = 10)
    ctrp.pred.auc <- fread('./calcPhenotype_Output/DrugPredictions.csv', data.table = F)
    write.csv(ctrp.pred.auc,file=paste("OncoPredict/CTRP.AUC.",project,".csv",sep = ""),quote=F, col.names=T, row.names=T)
  }
  if (! file.exists(paste("OncoPredict/PRISM.AUC.",project,".csv",sep = ""))) {
    trainPtype <- prism.auc
    ### The loop is slow, please be patient
    calcPhenotype(trainingExprData = as.matrix(trainExpr),
                  trainingPtype = as.matrix(trainPtype),
                  testExprData = as.matrix(testExpr),  ##"testExprData" must be a matrix.
                  powerTransformPhenotype = F,
                  removeLowVaringGenesFrom = "rawData",
                  printOutput = T,
                  selection = 2,batchCorrect = "eb",
                  minNumSamples = 10)
    prism.pred.auc <- fread('./calcPhenotype_Output/DrugPredictions.csv', data.table = F)
    write.csv(prism.pred.auc,file=paste("OncoPredict/PRISM.AUC.",project,".csv",sep = ""),quote=F, col.names=T, row.names=T)
  }
}

calcPhenotype(trainingExprData = trainExpr,
              trainingPtype = CTRP2_Res,
              testExprData = as.matrix(testExpr),
              batchCorrect = 'eb',  #   "eb" for ComBat  
              powerTransformPhenotype = F,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 10, 
              printOutput = TRUE)


###### PRISM #######
setwd("K:/TCGA/Anlysis/BCAA/Depmap")
files <- list.files(".",pattern = "csv$")
files <- files[str_detect(files,"PRISM")]
RiskScore.Cancer <- read.csv(file = "../RiskGroup.COX.TCGA.RiskScore.csv") %>%
  group_by(Cancer) %>%
  mutate(cutpoint=median(total_risk_score)) %>%
  mutate(Group=if_else(total_risk_score >= cutpoint,"High","Low")) %>%
  ungroup()
PRISM.RiskScore <- data.frame()
for (file in files) {
  project = (file %>% str_split("\\.") %>% unlist())[3]
  print(project)
  if (project %in% RiskScore.Cancer$Cancer) {
    data <- read.csv(file,row.names = 1,check.names = F) %>%
      magrittr::set_colnames(str_remove_all(colnames(.),"\\(.*")) %>%
      merge(RiskScore.Cancer,by.x="V1",by.y="SampleID")
    High.G <- data %>% filter(Group == "High")
    Low.G <- data %>% filter(Group == "Low")
    for (drug in colnames(data)[2:432]) {
      if (length(na.omit(data[[drug]])) > 0) {
        test <- wilcox.test(High.G[[drug]],Low.G[[drug]],alternative = "two.sided",
                            exact = F)
        PRISM.RiskScore <- data.frame(Pvalue=test$p.value,
                                      High.Mean=mean(High.G[[drug]]),
                                      High.Median=median(High.G[[drug]]),
                                      Low.Mean=mean(Low.G[[drug]]),
                                      Low.Median=median(Low.G[[drug]]),
                                      Drugs = drug,
                                      Cancer=project) %>%
          rbind.data.frame(PRISM.RiskScore)
      }
    }
  }
}
write.csv(PRISM.RiskScore,"Drugs.Depmap.PRISM.Wilcox.csv",row.names = F)

PRISM.RiskScore.S <- PRISM.RiskScore %>% filter(Pvalue <= 0.05) %>%
  mutate(Label = if_else(High.Median>Low.Median,"High","Low")) %>%
  group_by(Drugs,Label) %>%
  summarise(Count=n()) %>%
  arrange(desc(Count))%>%
  mutate(Drugs=str_remove_all(Drugs," "))

#### Figure ####
PRISM.RiskScore.S <- PRISM.RiskScore %>% #filter(Pvalue <= 0.05) %>%
  mutate(Label = if_else(Pvalue>0.05,"N.S.",
                         if_else(High.Median>Low.Median,"High","Low"))) %>%
  group_by(Drugs,Label) %>%
  summarise(Count=n()) %>%
  arrange(desc(Count))%>%
  filter(Label!="N.S.") %>%
  mutate(Drugs=str_remove_all(Drugs," ")) %>%
  filter(Count > 12)

PRISM.RiskScore2 <- PRISM.RiskScore %>%
  mutate(Drugs=str_remove_all(Drugs," ")) %>%
  mutate(Label = if_else(Pvalue>0.05,"N.S.",
                         if_else(High.Median>Low.Median,"High","Low"))) %>%
  filter(Drugs %in% PRISM.RiskScore.S$Drugs)

PRISM.RiskScore3 <- PRISM.RiskScore %>% #filter(Pvalue <= 0.05) %>%
  mutate(Label = if_else(Pvalue>0.05,"N.S.",
                         if_else(High.Median>Low.Median,"High","Low"))) %>%
  group_by(Drugs,Label) %>%
  summarise(Count=n()) %>%
  arrange(desc(Count))%>%
  #filter(Label!="N.S.") %>%
  mutate(Drugs=str_remove_all(Drugs," ")) %>%
  filter(Drugs %in% PRISM.RiskScore.S$Drugs)

p1 <- ggplot(data=PRISM.RiskScore2,
               aes(x=Cancer,y=Drugs,color=Label,shape=Label)
  )+
  geom_point(size=2)+
  ggthemes::theme_few()+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9,angle = 90,vjust = 0.5,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  #ggsci::scale_color_npg()+
  scale_color_manual(values = c("#E64B35FF","#3C5488FF","#B09C85FF"))+
  scale_shape_manual(values = c(1,3,4))+
  labs(x="",y="",color="Predicted AUC\nfor risk group",
       shape="Predicted AUC\nfor risk group")
p2 <- ggplot(data=PRISM.RiskScore3 %>%
               mutate(Label=factor(Label,levels = rev(c("High","Low","N.S.")))),
             aes(x=Count,y=Drugs,fill=Label))+
  geom_col()+
  scale_fill_manual(values = rev(c("#E64B35FF","#3C5488FF","#B09C85FF")))+
  labs(x="No. of cancer",y="")+
  ggthemes::theme_few()+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9,angle = 90,vjust = 0.5,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))
p<-p1 %>% aplot::insert_right(p2,width = 0.3)
ggsave(p,filename = "Drugs.PRISM.Drugs.Predicted.AUC.pdf",height = 6,width = 9)

  
  
  
  
###### CTRP  #######
setwd("K:/TCGA/Anlysis/BCAA/Depmap")
files <- list.files(".",pattern = "csv$")
files <- files[str_detect(files,"CTRP")]
RiskScore.Cancer <- read.csv(file = "../RiskGroup.COX.TCGA.RiskScore.csv") %>%
  group_by(Cancer) %>%
  mutate(cutpoint=median(total_risk_score)) %>%
  mutate(Group=if_else(total_risk_score >= cutpoint,"High","Low")) %>%
  ungroup()
CTRP.RiskScore <- data.frame()
for (file in files) {
  project = (file %>% str_split("\\.") %>% unlist())[3]
  print(project)
  if (project %in% RiskScore.Cancer$Cancer) {
    data <- read.csv(file,row.names = 1,check.names = F) %>%
      magrittr::set_colnames(str_remove_all(colnames(.),"\\(.*")) %>%
      merge(RiskScore.Cancer,by.x="V1",by.y="SampleID")
    High.G <- data %>% filter(Group == "High")
    Low.G <- data %>% filter(Group == "Low")
    for (drug in colnames(data)[2:407]) {
      if (length(na.omit(data[[drug]])) > 0) {
        test <- wilcox.test(High.G[[drug]],Low.G[[drug]],alternative = "two.sided",
                            exact = F)
        CTRP.RiskScore <- data.frame(Pvalue=test$p.value,
                                      High.Mean=mean(High.G[[drug]]),
                                      High.Median=median(High.G[[drug]]),
                                      Low.Mean=mean(Low.G[[drug]]),
                                      Low.Median=median(Low.G[[drug]]),
                                      Drugs = drug,
                                      Cancer=project) %>%
          rbind.data.frame(CTRP.RiskScore)
      }
    }
  }
}
write.csv(CTRP.RiskScore,"Drugs.Depmap.CTRP.Wilcox.csv",row.names = F)

CTRP.RiskScore.S <- CTRP.RiskScore %>% filter(Pvalue <= 0.05) %>%
  mutate(Label = if_else(High.Median>Low.Median,"High","Low")) %>%
  group_by(Drugs,Label) %>%
  summarise(Count=n()) %>%
  arrange(desc(Count)) %>%
  mutate(Drugs=str_remove_all(Drugs," "))
#### Figure ####
CTRP.RiskScore.S <- CTRP.RiskScore %>% #filter(Pvalue <= 0.05) %>%
  mutate(Label = if_else(Pvalue>0.05,"N.S.",
                         if_else(High.Median>Low.Median,"High","Low"))) %>%
  group_by(Drugs,Label) %>%
  summarise(Count=n()) %>%
  arrange(desc(Count))%>%
  filter(Label!="N.S.") %>%
  mutate(Drugs=str_remove_all(Drugs," ")) %>%
  filter(Count > 12)

CTRP.RiskScore2 <- CTRP.RiskScore %>%
  mutate(Drugs=str_remove_all(Drugs," ")) %>%
  mutate(Label = if_else(Pvalue>0.05,"N.S.",
                         if_else(High.Median>Low.Median,"High","Low"))) %>%
  filter(Drugs %in% c(CTRP.RiskScore.S$Drugs,"axitinib",
                      "docetaxel:tanespimycin","linsitinib","sorafenib"))

CTRP.RiskScore3 <- CTRP.RiskScore %>% #filter(Pvalue <= 0.05) %>%
  mutate(Label = if_else(Pvalue>0.05,"N.S.",
                         if_else(High.Median>Low.Median,"High","Low"))) %>%
  group_by(Drugs,Label) %>%
  summarise(Count=n()) %>%
  arrange(desc(Count))%>%
  #filter(Label!="N.S.") %>%
  mutate(Drugs=str_remove_all(Drugs," ")) %>%
  filter(Drugs %in% c(CTRP.RiskScore.S$Drugs,"axitinib",
                      "docetaxel:tanespimycin","linsitinib","sorafenib"))

p1 <- ggplot(data=CTRP.RiskScore2,
             aes(x=Cancer,y=Drugs,color=Label,shape=Label)
)+
  geom_point(size=2)+
  ggthemes::theme_few()+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9,angle = 90,vjust = 0.5,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  #ggsci::scale_color_npg()+
  scale_color_manual(values = c("#E64B35FF","#3C5488FF","#B09C85FF"))+
  scale_shape_manual(values = c(1,3,4))+
  labs(x="",y="",color="Predicted AUC\nfor risk group",
       shape="Predicted AUC\nfor risk group")
p2 <- ggplot(data=CTRP.RiskScore3 %>%
               mutate(Label=factor(Label,levels = rev(c("High","Low","N.S.")))),
             aes(x=Count,y=Drugs,fill=Label))+
  geom_col()+
  scale_fill_manual(values = rev(c("#E64B35FF","#3C5488FF","#B09C85FF")))+
  labs(x="No. of cancer",y="")+
  ggthemes::theme_few()+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9,angle = 90,vjust = 0.5,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))
p<-p1 %>% aplot::insert_right(p2,width = 0.3)
ggsave(p,filename = "Drugs.CTRP.Drugs.Predicted.AUC.pdf",height = 6,width = 9)

CTRP <- read.csv("Drugs.Depmap.CTRP.Wilcox.csv") %>%
  filter(Pvalue < 0.05) %>%
  mutate(Database="CTRP")

PRISM <- read.csv("Drugs.Depmap.PRISM.Wilcox.csv") %>%
  filter(Pvalue < 0.05) %>%
  mutate(Database="PRISM")

GDSC2 <- read.csv("../Pancancer.GDSC2.Drugs.Sensitivity.csv") %>%
  filter(p.adj < 0.05) %>%
  mutate(Database="GDSC2")
write.csv(GDSC2,file="GDSC2.csv")

rbind(CTRP,PRISM) %>% write.csv("TableS22.csv")

####################
###### Heterogeneity => T cell #######
setwd("K:/TCGA/ZhangZeMin.Pancancer.T")
library(data.table)
library(tidyverse)
library(Seurat)
library(SeuratObject)
library(magrittr)
# AUCell、UCell和singscore + ssGSEA
BCAA.ID <- read.csv("../Anlysis/BCAA/BCAA.ID.csv")
files <- list.files(".",pattern = ".gz$")
files <- files[str_detect(files,"CD8")|str_detect(files,"CD4")]
metadata <- fread("GSE156728_metadata.txt.gz")
CD4.CD8.Bulk.BCAA <- data.frame()
CD4.CD8.scRNA.ssGSEA.BCAA <- data.frame()
CD4.CD8.scRNA.AddModuleScore.BCAA <- data.frame()
library(ggsci)
color = c(pal_d3("category20")(20),pal_d3("category20b")(20),pal_d3("category20c")(20),pal_d3("category10")(10))

for (file in files) {
  print(file)
  project=(file %>%str_remove_all("\\..*") %>%str_split("_") %>% unlist())[2]
  Cell <- (file %>%str_split("\\.") %>% unlist())[2]
  data1 <- fread(file) %>% data.frame() %>% #filter(V1 %in% BCAA.ID$SYMBOL) %>%
    remove_rownames() %>% column_to_rownames("V1") %>%
    t() %>% data.frame(check.names = F) %>%
    rownames_to_column("cellID") %>%
    merge(metadata %>% dplyr::select(cellID,loc,meta.cluster,cancerType),
          by="cellID") %>%
    column_to_rownames("cellID") %>% dplyr::select(-cancerType)
  
  data2 <- aggregate(.~loc+meta.cluster,mean,data=data1) %>%
    mutate(axis=paste(loc,meta.cluster,sep = "_")) %>%
    dplyr::select(-loc,-meta.cluster) %>%
    remove_rownames() %>% column_to_rownames("axis") %>%
    t() %>% data.frame(check.names = F)
  
  library(GSVA)
  BCAA=list()
  BCAA[["BCAA"]]<-BCAA.ID$SYMBOL
  #data2 <- fread(file) %>% remove_rownames() %>% column_to_rownames("V1")
  res.ssgsea <- gsva(as.matrix(data2), 
                     gset.idx.list=BCAA, 
                     method = "ssgsea", 
                     kcdf = "Poisson")
  CD4.CD8.Bulk.BCAA <- data.frame(cellType=colnames(res.ssgsea),
                                  BCAA=res.ssgsea[1,],Cancer=project) %>%
    rbind.data.frame(CD4.CD8.Bulk.BCAA)
  
  data2 <- fread(file) %>% remove_rownames() %>% column_to_rownames("V1")
  res.ssgsea <- gsva(as.matrix(data2), 
                     gset.idx.list=BCAA, 
                     method = "ssgsea", 
                     kcdf = "Poisson")
  CD4.CD8.scRNA.ssGSEA.BCAA <- data.frame(cellType=colnames(res.ssgsea),
                                  BCAA=res.ssgsea[1,],Cancer=project) %>%
    rbind.data.frame(CD4.CD8.scRNA.ssGSEA.BCAA)
  
  sce <- CreateSeuratObject(data2)
  sce[["percent.mt"]] <- PercentageFeatureSet(
    sce,
    pattern = "^MT-"
  )
  
  meta <- metadata %>% filter(cancerType==project)
  #meta<-data.table::fread("GSE210347_meta.txt",data.table = F) %>%
  #  dplyr::select(cellname,SampleID,seurat_clusters,cluster,celltype,tissue,group)
  data <- sce@meta.data
  data <- rownames_to_column(data,var = "barcodes")
  meta1 <- merge(data,
                    meta,
                    by.x='barcodes',
                    by.y='cellID',
                    all=T)
  meta1 <- column_to_rownames(meta1,var = "barcodes")
  sce <- AddMetaData(sce,meta1)
  sce %<>% subset(subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20) %>% #nolint
    NormalizeData(normalization.method = "LogNormalize") %>%
    FindVariableFeatures(selection.method = "vst")
  sce %<>% ScaleData()
  sce %<>% RunPCA(object = ., features = VariableFeatures(object = .))
  sce <- RunUMAP(sce,  dims = 1:30, reduction = "pca")
  #sce <- AddMetaData(sce,metadata = cancer.metadata)
  p1 <- DimPlot(sce,group.by = "meta.cluster",
          split.by = "loc",cols = color)+
    theme(panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
          plot.title = element_blank())
  ggsave(p1,filename = paste(project,Cell,"loc.cluster.pdf",sep = "."),width = 12,height = 4)
  p1 <- DimPlot(sce,group.by = "meta.cluster",cols = color)+
    theme(panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
          plot.title = element_blank())
  ggsave(p1,filename = paste(project,Cell,"cluster.pdf",sep = "."),width = 8,height = 4)
  p1 <- FeaturePlot(sce,features = c("BCAT1"))+
    theme(panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
    scale_color_gradientn(colors = c("grey90",brewer.pal(11,"RdBu")[5:1]))
  ggsave(p1,filename = paste(project,Cell,"BCAT1.pdf",sep = "."),height = 3,width = 3.5)
  p1 <- FeaturePlot(sce,features = c("BCAT2"))+
    theme(panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
    scale_color_gradientn(colors = c("grey90",brewer.pal(11,"RdBu")[5:1]))
  ggsave(p1,filename = paste(project,Cell,"BCAT2.pdf",sep = "."),height = 3,width = 3.5)

  BCAA <- c("ACAD8","ACADSB","ACAT1","ALDH6A1","AUH","BCAT1","BCAT2","BCKDHA","BCKDHB","DBT","DLD","ECHS1","HIBADH","HIBCH","HSD17B10","IVD","MCCC1","MCCC2")
  gene <- as.list(BCAA)
  sce <- AddModuleScore(sce,features = gene,name = 'BCAA')
  saveRDS(sce, paste(project,Cell,"sce.Preprocess.Rds",sep = "."))
  p1 <- FeaturePlot(sce,features = c("BCAA1"))+
    theme(panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
    scale_color_gradientn(colors = c("grey90",brewer.pal(11,"RdBu")[5:1]))+
    labs(title="BCAA score")
  ggsave(p1,filename = paste(project,Cell,"BCAA.score.pdf",sep = "."),height = 3,width = 3.5)
  
  mydata<- FetchData(sce,vars = c("BCAA1")) %>% mutate(Cancer=project)
  #write.csv(mydata,"FetchData.UMAP.BCAA.csv")
  CD4.CD8.scRNA.AddModuleScore.BCAA <- rbind.data.frame(CD4.CD8.scRNA.AddModuleScore.BCAA,mydata)
}
#CD4.CD8.Bulk.BCAA <- data.frame()
#CD4.CD8.scRNA.ssGSEA.BCAA <- data.frame()
#CD4.CD8.scRNA.AddModuleScore.BCAA <- data.frame()
write.csv(CD4.CD8.Bulk.BCAA,file = "Pancancer.CD4.CD8.Bulk.BCAA.csv",row.names = F)
write.csv(CD4.CD8.scRNA.ssGSEA.BCAA,file = "Pancancer.CD4.CD8.scRNA.ssGSEA.BCAA",row.names = F)
write.csv(CD4.CD8.scRNA.AddModuleScore.BCAA,file = "Pancancer.CD4.CD8.scRNA.AddModuleScore.BCAA.csv",row.names = F)

#### Figure ####
setwd("K:/TCGA/Anlysis/BCAA/Hetero/T")
data <- read.csv("Pancancer.CD4.CD8.Bulk.BCAA.csv") %>%
  mutate(CellType=str_remove(cellType,".*_")) %>%
  mutate(Loc = str_remove(cellType,"_.*")) %>%
  mutate(Yaxis=paste(CellType,Loc,sep = "_"))

p1 <- ggplot(data,aes(y=Cancer,x=Yaxis))+
  geom_tile(aes(fill=BCAA),color="white")+
  ggthemes::theme_few()+
  theme(legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.y = element_text(size=9), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="",y="")+
  scale_fill_gradientn(colors = brewer.pal(9,"Reds")) # YlOrRd Blues BuGn

library(ggsci)
color = c(pal_d3("category20")(20),pal_d3("category20b")(20),pal_d3("category20c")(20),pal_d3("category10")(10))
mid <- data %>% dplyr::select(Yaxis,CellType) %>% unique() %>%
  group_by(CellType) %>% summarise(Count=n())

p2 <- ggplot(data %>% dplyr::select(Yaxis,CellType) %>% unique(),
             aes(y=1,x=Yaxis))+
  geom_col(aes(fill=CellType),width = 1)+
  #ggthemes::theme_few()+
  #theme_void()+
  ggthemes::theme_few()+
  theme(legend.position = "bottom",
    #axis.text = element_text(color = "black"),
    #axis.text.x = element_text(size=9,angle = 90,vjust = 0.5,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  geom_vline(xintercept = accumulate(mid$Count,sum)+0.5)+
  scale_fill_manual(values = color)+
  labs(x="",y="")
  #scale_fill_gradientn(colors = brewer.pal(9,"Reds")) # YlOrRd Blues BuGn
mid <- data %>% dplyr::select(Loc,CellType) %>% unique() %>%
  group_by(CellType) %>% summarise(Count=n())


p3 <- ggplot(data %>% dplyr::select(Yaxis,Loc) %>% unique(),
             aes(y=1,x=Yaxis))+
  geom_col(aes(fill=Loc),width = 1)+
  ggthemes::theme_few()+
  theme(legend.position = "bottom",
        #axis.text = element_text(color = "black"),
        #axis.text.x = element_text(size=9,angle = 90,vjust = 0.5,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  #ggsci::scale_color_nejm()+
  scale_fill_manual(values = rev(c("#DC0000FF","#F39B7FFF","#00A087FF")))+
  theme(legend.position = "bottom")+
  geom_vline(xintercept = accumulate(mid$Count,sum)+0.5)+
  labs(x="",y="")
  #scale_fill_manual(values = color)
p123 <- p1 %>% aplot::insert_top(p2,height = 0.1) %>%
  aplot::insert_top(p3,height = 0.1)
ggsave(p123,filename = "Hetero.T.BCAAscore.ssGSEA.pdf",height = 3,width = 18)
ggsave(p123,filename = "Hetero.T.BCAAscore.ssGSEA.Legend.pdf",height = 10,width = 12)


#### ESCA ####
setwd("K:/TCGA/Anlysis/BCAA/Hetero/T")
data <- read.csv("Pancancer.CD4.CD8.Bulk.BCAA.csv") %>%
  mutate(CellType=str_remove(cellType,".*_")) %>%
  mutate(Loc = str_remove(cellType,"_.*")) %>%
  mutate(Yaxis=paste(CellType,Loc,sep = "_")) %>%
  filter(Cancer=="ESCA")

p1 <- ggplot(data %>% mutate(Loc=factor(Loc,levels = c("T","N")))
             ,aes(x=BCAA,y=Yaxis))+
  geom_segment(aes(x=0,xend=BCAA,y=Yaxis,yend=Yaxis),color="black")+
  geom_point(aes(color=Loc),size=3)+
  ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=7,angle = 90,vjust = 0.5,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
        axis.text.y = element_text(size=8),
        #axis.ticks.x = element_blank(),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="BCAA score",y="")+
  geom_hline(yintercept = seq(2,100,2)+0.5,linetype=2)+
  scale_color_npg()

ggsave(p1,filename = "ESCA.BCAAscore.pdf",height = 8,width = 4)




#### BCAA score => AUCell ####
setwd("K:/TCGA/Anlysis/BCAA/Hetero/T")
files <- list.files(".",pattern = "AUCell.*.csv")
AUCell.Stat <- data.frame()
for (file in files) {
  data <- read.csv(file) %>% mutate(CellType=str_extract(meta.cluster,"\\..*\\.")) %>%
    mutate(CellType=str_remove_all(CellType,"\\.c.?.?\\.") %>% str_remove_all("\\."))
  string <- file %>% str_split("\\.") %>% unlist()
  project <- string[2]
  Cell <- string[3]
  print(project)
  if (project != "FTC") {
    if (project %in% c("BCL","MM")) {
      middata <- data %>% group_by(loc,CellType) %>% summarise(Count=n()) %>%
        arrange(CellType) %>% ungroup() %>%
        spread(loc,Count) %>% filter(P>=3) %>% filter(T>=3)
      
      AUCell.Stat <- data %>%
        filter(CellType %in% middata$CellType) %>%
        group_by(CellType) %>%
        rstatix::t_test(AUC~loc,detailed = T) %>%
        mutate(Cancer=project,Cell=Cell) %>%
        plyr::rbind.fill(AUCell.Stat)
    }else{
      middata <- data %>% group_by(loc,CellType) %>% summarise(Count=n()) %>%
        arrange(CellType) %>% ungroup() %>%
        spread(loc,Count) %>% filter(N>=3) %>% filter(T>=3)
      
      AUCell.Stat <- data %>%
        filter(CellType %in% middata$CellType) %>%
        group_by(CellType) %>%
        rstatix::t_test(AUC~loc,detailed = T) %>%
        mutate(Cancer=project,Cell=Cell) %>%
        plyr::rbind.fill(AUCell.Stat)
    }
  }
}
write.csv(AUCell.Stat,"AUCell.BCAA.T.Stat.csv",row.names = F)

#### BCAA SYMBOLs => T test ####
setwd("K:/TCGA/Anlysis/BCAA/Hetero/T")
data <- read.csv("BCAA.SYMBOLS.T.Expression.csv") %>% mutate(CellType=str_extract(meta.cluster,"\\..*\\.")) %>%
  mutate(CellType=str_remove_all(CellType,"\\.c.?.?\\.") %>% str_remove_all("\\.")) %>%
  mutate(CellType=paste(FCell,CellType,sep = "_"))
SYMBOLS.Stat <- data.frame()
for (project in unique(data$Cancer)) {
  #data <- read.csv(file) 
  #string <- file %>% str_split("\\.") %>% unlist()
  #project <- string[2]
  #Cell <- string[3]
  print(project)
  mid <- data %>% filter(Cancer==project) %>% 
    gather(SYMBOL,TPM,ACAD8:MCCC2)
  if (project != "FTC") {
    if (project %in% c("BCL","MM")) {
      middata <- mid %>% group_by(loc,CellType,SYMBOL) %>% summarise(Count=n()) %>%
        arrange(CellType) %>% ungroup() %>%
        spread(loc,Count) %>% filter(P>=3) %>% filter(T>=3)
      
      SYMBOLS.Stat <- mid %>%
        filter(CellType %in% middata$CellType) %>%
        group_by(CellType,SYMBOL) %>%
        rstatix::t_test(TPM~loc,detailed = T) %>%
        mutate(Cancer=project) %>%
        plyr::rbind.fill(SYMBOLS.Stat)
    }else{
      middata <- mid %>% group_by(loc,CellType,SYMBOL) %>% summarise(Count=n()) %>%
        arrange(CellType) %>% ungroup() %>%
        spread(loc,Count) %>% filter(N>=3) %>% filter(T>=3)
      
      SYMBOLS.Stat <- mid %>%
        filter(CellType %in% middata$CellType) %>%
        group_by(CellType,SYMBOL) %>%
        rstatix::t_test(TPM~loc,detailed = T) %>%
        mutate(Cancer=project) %>%
        plyr::rbind.fill(SYMBOLS.Stat)
    }
  }
}
write.csv(SYMBOLS.Stat,"SYMBOLS.BCAA.T.Stat.csv",row.names = F)

###### Heterogeneity => T cell => CMP ######
setwd("K:/TCGA/ZhangZeMin.Pancancer.T")
library(data.table)
library(tidyverse)
library(Seurat)
library(SeuratObject)
library(magrittr)
# AUCell、UCell和singscore + ssGSEA
countToFpkm <- function(counts, effLen) 
{ 
  N <- sum(counts) 
  exp( log(counts) + log(1e9) - log(effLen) - log(N) ) 
}

#Gene.Length <- read.csv("../gene.length.txt",sep = "\t",header = F)
library(vegan)
#BCAA.ID <- read.csv("../Anlysis/BCAA/BCAA.ID.csv")
files <- list.files(".",pattern = ".gz$")
files <- files[str_detect(files,"CD8")|str_detect(files,"CD4")]
metadata <- fread("GSE156728_metadata.txt.gz")
for (file in files) {
  print(file)
  project=(file %>%str_remove_all("\\..*") %>%str_split("_") %>% unlist())[2]
  Cell <- (file %>%str_split("\\.") %>% unlist())[2]
  if (!file.exists(paste(project,Cell,"CPM.csv",sep = "."))) {
    data1 <- fread(file) %>% data.frame() %>% #filter(V1 %in% BCAA.ID$SYMBOL) %>%
      remove_rownames() %>% column_to_rownames("V1") #%>%
    #t() %>% data.frame(check.names = F) %>%
    #rownames_to_column("cellID") %>%
    #merge(metadata %>% dplyr::select(cellID,loc,meta.cluster,cancerType),
    #      by="cellID") %>%
    #column_to_rownames("cellID") %>% dplyr::select(-cancerType)
    #data <- vegan::decostand(data1,method = "total",MARGIN = 2) * 1e6
    cpm <- apply(data1 ,2, function(x) { x/sum(x)*1000000 })
    write.csv(data,file = paste(project,Cell,"CPM.csv",sep = "."))
  }
}

#### Datasets ####
setwd("K:/TCGA/ZhangZeMin.Pancancer.T")
library(data.table)
library(tidyverse)
library(Seurat)
library(SeuratObject)
library(magrittr)
BCAA.ID <- read.csv("../Anlysis/BCAA/BCAA.ID.csv")
files <- list.files("CPM/",pattern = "csv$")
files <- files[str_detect(files,"CPM") & !str_detect(files,"Pancancer")]
metadata <- fread("GSE156728_metadata.txt.gz")
BCAA.T <- data.frame()
for (file in files) {
  print(file)
  #BCAA.T <- data.frame()
  BCAA.T <- read.csv(file.path("CPM",file),row.names = 1) %>%
    rownames_to_column("SYMBOL") %>%
    filter(SYMBOL %in% BCAA.ID$SYMBOL) %>%
    t() %>% data.frame(check.names = F) %>%
    rownames_to_column("cellID") %>%
    merge(metadata,by="cellID") %>%
    rbind.data.frame(BCAA.T)
}
write.csv(BCAA.T,"Pancancer.BCAA.SYMBOL.CPM.csv",row.names = F)
#### Figure => Data => AveExp + PerExp ####
setwd("K:/TCGA/Anlysis/BCAA/")
BCAA.ID <- read.csv("BCAA.ID.csv") %>%
  arrange(SYMBOL)
BCAA.T <- read.csv("Pancancer.BCAA.SYMBOL.CPM.csv")
library(data.table)
library(tidyverse)
library(Seurat)
library(SeuratObject)
library(magrittr)
colnames(BCAA.T)[2:19]=BCAA.ID$SYMBOL

BCAA.T.AveExp <- BCAA.T %>% filter(loc != "L") %>%
  gather(SYMBOL,Expression,ACAD8:MCCC2) %>%
  mutate(cell=str_remove_all(meta.cluster,"\\..*")) %>%
  mutate(type=str_extract(meta.cluster,"\\..*\\.") %>% 
           str_remove_all("\\.c..\\.") %>% str_remove_all("\\.")) %>%
  mutate(cluster=paste(cell,type,sep = "_")) %>%
  group_by(cancerType,cluster,loc,SYMBOL) %>%
  summarise(AveExp = mean(Expression)) %>%
  ungroup()

BCAA.T.PerExp <- BCAA.T %>% filter(loc != "L") %>%
  gather(SYMBOL,Expression,ACAD8:MCCC2) %>%
  mutate(cell=str_remove_all(meta.cluster,"\\..*")) %>%
  mutate(type=str_extract(meta.cluster,"\\..*\\.") %>% 
           str_remove_all("\\.c..\\.") %>% str_remove_all("\\.")) %>%
  mutate(cluster=paste(cell,type,sep = "_")) %>%
  group_by(cancerType,cluster,loc,SYMBOL) %>%
  mutate(TotalCount=sum(Expression>=0)) %>%
  mutate(ExpCount = sum(Expression >0)) %>%
  dplyr::select(cancerType,loc,cluster,TotalCount,ExpCount) %>%
  unique() %>%
  mutate(PerExp = ExpCount*100/TotalCount)
write.csv(BCAA.T.AveExp,file = "BCAA.SYMBOLS.T.CLUSTER.AVEEXP.csv",row.names = F)
write.csv(BCAA.T.PerExp,file = "BCAA.SYMBOLS.T.CLUSTER.PEREXP.csv",row.names = F)
#### Figure => CD4 ####
library(ggsci)
color = c(pal_d3("category20")(20),pal_d3("category20b")(20),pal_d3("category20c")(20),pal_d3("category10")(10))

BCAA.T.D <- merge(BCAA.T.AveExp,BCAA.T.PerExp,
                  by=c("cancerType","loc","cluster","SYMBOL")) %>%
  data.frame(check.names = F) %>%
  filter(str_detect(cluster,"CD4")) %>%
  #filter(SYMBOL=="BCAT1") %>%
  mutate(mid = paste(cancerType,loc,sep = "_")) %>%
  mutate(Xaxis=paste(cluster,loc,cancerType,sep = "_")) %>%
  arrange(cluster,mid) %>%
  mutate(Xaxis=factor(Xaxis,levels = unique(.$Xaxis)))
p1 <- ggplot(BCAA.T.D, 
             aes(x=Xaxis,y=SYMBOL,fill=log(AveExp)))+
  geom_tile(color="white")+
  ggthemes::theme_few()+
  theme(legend.position = "bottom",
        axis.text = element_text(color = "black"),
        axis.text.y = element_text(size=9), #,angle = 90,vjust = 0.3,hjust = 1
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="",y="")+
  scale_fill_gradient2(high = brewer.pal(11,"RdBu")[1],
                       low=brewer.pal(11,"RdBu")[11])

p2 <- ggplot(BCAA.T.D,aes(x=Xaxis,y=SYMBOL,fill=PerExp))+
  geom_tile(color="white")+
  ggthemes::theme_few()+
  theme(legend.position = "bottom",
        axis.text = element_text(color = "black"),
        axis.text.y = element_text(size=9), #,angle = 90,vjust = 0.3,hjust = 1
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="",y="")+
  scale_fill_gradientn(colours = brewer.pal(9,"Blues"))

p3 <- ggplot(BCAA.T.D %>% 
               dplyr::select(Xaxis,cancerType) %>% unique(),
             aes(y=1,x=Xaxis))+
  geom_col(aes(fill=cancerType),width = 1)+
  ggthemes::theme_few()+
  theme(legend.position = "bottom",
        #axis.text = element_text(color = "black"),
        #axis.text.x = element_text(size=9,angle = 90,vjust = 0.5,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  #ggsci::scale_color_nejm()+
  ggthemes::scale_fill_stata()+
  #scale_fill_manual(values = rev(c("#DC0000FF","#F39B7FFF","#00A087FF")))+
  theme(legend.position = "bottom")+
  #geom_vline(xintercept = accumulate(mid$Count,sum)+0.5)+
  labs(x="",y="") 

p4 <- ggplot(BCAA.T.D %>% 
               dplyr::select(Xaxis,cluster) %>% unique(),
             aes(y=1,x=Xaxis))+
  geom_col(aes(fill=cluster),width = 1)+
  ggthemes::theme_few()+
  theme(legend.position = "bottom",
        #axis.text = element_text(color = "black"),
        #axis.text.x = element_text(size=9,angle = 90,vjust = 0.5,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  #ggsci::scale_color_nejm()+
  #ggthemes::scale_fill_stata()+
  scale_fill_manual(values = color)+
  theme(legend.position = "bottom")+
  #geom_vline(xintercept = accumulate(mid$Count,sum)+0.5)+
  labs(x="",y="")

p5 <- ggplot(BCAA.T.D %>% 
               dplyr::select(Xaxis,loc) %>% unique(),
             aes(y=1,x=Xaxis))+
  geom_col(aes(fill=loc),width = 1)+
  ggthemes::theme_few()+
  theme(legend.position = "bottom",
        #axis.text = element_text(color = "black"),
        #axis.text.x = element_text(size=9,angle = 90,vjust = 0.5,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  #ggsci::scale_color_nejm()+
  #ggthemes::scale_fill_stata()+
  scale_fill_manual(values = rev(c("#DC0000FF","#F39B7FFF","#00A087FF")))+
  theme(legend.position = "bottom")+
  #geom_vline(xintercept = accumulate(mid$Count,sum)+0.5)+
  labs(x="",y="")

p12345 <- p1 %>% aplot::insert_bottom(p2,height = 1) %>%
  aplot::insert_top(p3,height = 0.1)%>%
  aplot::insert_top(p4,height = 0.1)%>%
  aplot::insert_top(p5,height = 0.1)
ggsave(p12345,filename = "Htero.SYMBOL.T.CD4.pdf",height = 8,width = 14)
#### Figure => CD8 #### 
BCAA.T.D <- merge(BCAA.T.AveExp,BCAA.T.PerExp,
                  by=c("cancerType","loc","cluster","SYMBOL")) %>%
  data.frame(check.names = F) %>%
  filter(str_detect(cluster,"CD8")) %>%
  #filter(SYMBOL=="BCAT1") %>%
  mutate(mid = paste(cancerType,loc,sep = "_")) %>%
  mutate(Xaxis=paste(cluster,loc,cancerType,sep = "_")) %>%
  arrange(cluster,mid) %>%
  mutate(Xaxis=factor(Xaxis,levels = unique(.$Xaxis)))
p1 <- ggplot(BCAA.T.D, 
             aes(x=Xaxis,y=SYMBOL,fill=log(AveExp)))+
  geom_tile(color="white")+
  ggthemes::theme_few()+
  theme(legend.position = "bottom",
        axis.text = element_text(color = "black"),
        axis.text.y = element_text(size=9), #,angle = 90,vjust = 0.3,hjust = 1
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="",y="")+
  scale_fill_gradient2(high = brewer.pal(11,"RdBu")[1],
                       low=brewer.pal(11,"RdBu")[11])

p2 <- ggplot(BCAA.T.D,aes(x=Xaxis,y=SYMBOL,fill=PerExp))+
  geom_tile(color="white")+
  ggthemes::theme_few()+
  theme(legend.position = "bottom",
        axis.text = element_text(color = "black"),
        axis.text.y = element_text(size=9), #,angle = 90,vjust = 0.3,hjust = 1
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="",y="")+
  scale_fill_gradientn(colours = brewer.pal(9,"Blues"))

p3 <- ggplot(BCAA.T.D %>% 
               dplyr::select(Xaxis,cancerType) %>% unique(),
             aes(y=1,x=Xaxis))+
  geom_col(aes(fill=cancerType),width = 1)+
  ggthemes::theme_few()+
  theme(legend.position = "bottom",
        #axis.text = element_text(color = "black"),
        #axis.text.x = element_text(size=9,angle = 90,vjust = 0.5,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  #ggsci::scale_color_nejm()+
  ggthemes::scale_fill_stata()+
  #scale_fill_manual(values = rev(c("#DC0000FF","#F39B7FFF","#00A087FF")))+
  theme(legend.position = "bottom")+
  #geom_vline(xintercept = accumulate(mid$Count,sum)+0.5)+
  labs(x="",y="") 

p4 <- ggplot(BCAA.T.D %>% 
               dplyr::select(Xaxis,cluster) %>% unique(),
             aes(y=1,x=Xaxis))+
  geom_col(aes(fill=cluster),width = 1)+
  ggthemes::theme_few()+
  theme(legend.position = "bottom",
        #axis.text = element_text(color = "black"),
        #axis.text.x = element_text(size=9,angle = 90,vjust = 0.5,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  #ggsci::scale_color_nejm()+
  #ggthemes::scale_fill_stata()+
  scale_fill_manual(values = color)+
  theme(legend.position = "bottom")+
  #geom_vline(xintercept = accumulate(mid$Count,sum)+0.5)+
  labs(x="",y="")

p5 <- ggplot(BCAA.T.D %>% 
               dplyr::select(Xaxis,loc) %>% unique(),
             aes(y=1,x=Xaxis))+
  geom_col(aes(fill=loc),width = 1)+
  ggthemes::theme_few()+
  theme(legend.position = "bottom",
        #axis.text = element_text(color = "black"),
        #axis.text.x = element_text(size=9,angle = 90,vjust = 0.5,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  #ggsci::scale_color_nejm()+
  #ggthemes::scale_fill_stata()+
  scale_fill_manual(values = rev(c("#DC0000FF","#F39B7FFF","#00A087FF")))+
  theme(legend.position = "bottom")+
  #geom_vline(xintercept = accumulate(mid$Count,sum)+0.5)+
  labs(x="",y="")

p12345 <- p1 %>% aplot::insert_bottom(p2,height = 1) %>%
  aplot::insert_top(p3,height = 0.1)%>%
  aplot::insert_top(p4,height = 0.1)%>%
  aplot::insert_top(p5,height = 0.1)
ggsave(p12345,filename = "Htero.SYMBOL.T.CD8.pdf",height = 8,width = 14)



#### Example => CD4 ####
BCAA.T.D <- merge(BCAA.T.AveExp,BCAA.T.PerExp,
                  by=c("cancerType","loc","cluster","SYMBOL")) %>%
  data.frame(check.names = F) %>%
  filter(str_detect(cluster,"CD4")) %>%
  #filter(SYMBOL=="BCAT1") %>%
  mutate(mid = paste(cancerType,loc,sep = "_")) %>%
  mutate(Xaxis=paste(cluster,loc,cancerType,sep = "_")) %>%
  arrange(cluster,mid) %>%
  mutate(Xaxis=factor(Xaxis,levels = unique(.$Xaxis)))
Exa.BCAT1 <- BCAA.T.D %>% filter(SYMBOL == "BCAT1") %>%
  mutate(Xaxis=paste(cancerType,loc,sep = "_"))
p1 <- ggplot(Exa.BCAT1,aes(x=Xaxis,y=cluster,
                           color=log10(AveExp),size=PerExp))+
  geom_point()+
  ggthemes::theme_few()+
  theme(legend.position = "bottom",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9,angle = 90,vjust = 0.5,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
        axis.text.y = element_text(size=9),
        #axis.ticks = element_blank(),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  #scale_color_gradientn(colors = brewer.pal(9,"Reds"))+
  scale_color_gradient2(high = brewer.pal(11,"RdBu")[1],
                       low=brewer.pal(11,"RdBu")[11],mid = "white")+
  scale_size_continuous(breaks=c(25,50,75,100))+
  labs(x="",y="")
ggsave(p1,filename = "Htero.SYMBOL.T.CD4.BCAT1.pdf",height = 4,width = 6)
#### BCAT1 => CD8 ####
BCAA.T.D <- merge(BCAA.T.AveExp,BCAA.T.PerExp,
                  by=c("cancerType","loc","cluster","SYMBOL")) %>%
  data.frame(check.names = F) %>%
  filter(str_detect(cluster,"CD8")) %>%
  #filter(SYMBOL=="BCAT1") %>%
  mutate(mid = paste(cancerType,loc,sep = "_")) %>%
  mutate(Xaxis=paste(cluster,loc,cancerType,sep = "_")) %>%
  arrange(cluster,mid) %>%
  mutate(Xaxis=factor(Xaxis,levels = unique(.$Xaxis)))
Exa.BCAT1 <- BCAA.T.D %>% filter(SYMBOL == "BCAT1") %>%
  mutate(Xaxis=paste(cancerType,loc,sep = "_"))
p1 <- ggplot(Exa.BCAT1,aes(x=Xaxis,y=cluster,
                           color=log10(AveExp),size=PerExp))+
  geom_point()+
  ggthemes::theme_few()+
  theme(legend.position = "bottom",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9,angle = 90,vjust = 0.5,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
        axis.text.y = element_text(size=9),
        #axis.ticks = element_blank(),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  #scale_color_gradientn(colors = brewer.pal(9,"Reds"))+
  scale_color_gradient2(high = brewer.pal(11,"RdBu")[1],
                        low=brewer.pal(11,"RdBu")[11],mid = "white")+
  scale_size_continuous(breaks=c(25,50,75,100))+
  labs(x="",y="")
ggsave(p1,filename = "Htero.SYMBOL.T.CD8.BCAT1.pdf",height = 4,width = 6)
###### Heterogeneity => myeloid cell #######
setwd("K:/TCGA/ZhangZeMin.Pancancer.Myeloid/")
library(data.table)
library(tidyverse)
library(Seurat)
library(SeuratObject)
library(magrittr)
# AUCell、UCell和singscore + ssGSEA
BCAA.ID <- read.csv("../Anlysis/BCAA/BCAA.ID.csv")
files <- list.files(".",pattern = ".gz$")
files <- files[str_detect(files,"expression")]
library(ggsci)
color = c(pal_d3("category20")(20),pal_d3("category20b")(20),pal_d3("category20c")(20),pal_d3("category10")(10))
#metadata <- fread("GSE156728_metadata.txt.gz")
Bulk.BCAA <- data.frame()

for (file in files) {
  print(file)
  project=(file %>% str_split("_") %>% unlist())[2]
  metadata <- fread(paste("GSE154763_",project,"_metadata.csv.gz",sep = ""))
  #Cell <- (file %>%str_split("\\.") %>% unlist())[2]
  data1 <- fread(file) %>% data.frame() %>% #filter(V1 %in% BCAA.ID$SYMBOL) %>%
    #remove_rownames() %>% column_to_rownames("index") %>%
    #t() %>% data.frame(check.names = F) %>%
    #rownames_to_column("cellID") %>%
    merge(metadata %>% dplyr::select(index,tissue,MajorCluster),
          by="index") %>%
    column_to_rownames("index")
  
  data2 <- aggregate(.~tissue+MajorCluster,mean,data=data1) %>%
    mutate(axis=paste(tissue,MajorCluster,sep = "_")) %>%
    dplyr::select(-tissue,-MajorCluster) %>%
    remove_rownames() %>% column_to_rownames("axis") %>%
    t() %>% data.frame(check.names = F)
  
  library(GSVA)
  BCAA=list()
  BCAA[["BCAA"]]<-BCAA.ID$SYMBOL
  #data2 <- fread(file) %>% remove_rownames() %>% column_to_rownames("V1")
  res.ssgsea <- gsva(as.matrix(data2), 
                     gset.idx.list=BCAA, 
                     method = "ssgsea", 
                     kcdf = "Gaussian")
  Bulk.BCAA <- data.frame(cellType=colnames(res.ssgsea),
                          BCAA=res.ssgsea[1,],Cancer=project) %>%
    rbind.data.frame(Bulk.BCAA)
  
  data2 <- fread(file) %>% remove_rownames() %>% column_to_rownames("index") %>% t()
  res.ssgsea <- gsva(as.matrix(data1), 
                     gset.idx.list=BCAA, 
                     method = "ssgsea", 
                     kcdf = "Gaussian")
  scRNA.ssGSEA.BCAA <- data.frame(cellType=colnames(res.ssgsea),
                                  BCAA=res.ssgsea[1,],Cancer=project) %>%
    rbind.data.frame(scRNA.ssGSEA.BCAA)
  
  sce <- CreateSeuratObject(data1)
  sce[["percent.mt"]] <- PercentageFeatureSet(
    sce,
    pattern = "^MT-"
  )
  
  #meta <- metadata %>% filter(cancerType==project)
  metadata <- fread(paste("GSE154763_",project,"_metadata.csv.gz",sep = ""))
  #meta<-data.table::fread("GSE210347_meta.txt",data.table = F) %>%
  #  dplyr::select(cellname,SampleID,seurat_clusters,cluster,celltype,tissue,group)
  data <- sce@meta.data
  data <- rownames_to_column(data,var = "barcodes")
  meta1 <- merge(data,
                 meta,
                 by.x='barcodes',
                 by.y='index',
                 all=T)
  meta1 <- column_to_rownames(meta1,var = "barcodes")
  sce <- AddMetaData(sce,meta1)
  sce %<>% #subset(subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20) %>% #nolint
    NormalizeData(normalization.method = "LogNormalize") %>%
    FindVariableFeatures(selection.method = "vst")
  sce %<>% ScaleData()
  #sce %<>% RunPCA(object = ., features = VariableFeatures(object = .))
  #sce <- RunUMAP(sce,  dims = 1:30, reduction = "pca")
  #sce <- AddMetaData(sce,metadata = cancer.metadata)
  #p1 <- DimPlot(sce,group.by = "meta.cluster",split.by = "loc",cols = color)+
  #  theme(panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
  #        plot.title = element_blank())
  #ggsave(p1,filename = paste(project,Cell,"loc.cluster.pdf",sep = "."),width = 12,height = 4)
  #p1 <- DimPlot(sce,group.by = "meta.cluster",cols = color)+
  #  theme(panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
  #        plot.title = element_blank())
  #ggsave(p1,filename = paste(project,Cell,"cluster.pdf",sep = "."),width = 8,height = 4)
  #p1 <- FeaturePlot(sce,features = c("BCAT1"))+
  #  theme(panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  #  scale_color_gradientn(colors = c("grey90",brewer.pal(11,"RdBu")[5:1]))
  #ggsave(p1,filename = paste(project,Cell,"BCAT1.pdf",sep = "."),height = 3,width = 3.5)
  #p1 <- FeaturePlot(sce,features = c("BCAT2"))+
  #  theme(panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  #  scale_color_gradientn(colors = c("grey90",brewer.pal(11,"RdBu")[5:1]))
  #ggsave(p1,filename = paste(project,Cell,"BCAT2.pdf",sep = "."),height = 3,width = 3.5)
  
  #BCAA <- c("ACAD8","ACADSB","ACAT1","ALDH6A1","AUH","BCAT1","BCAT2","BCKDHA","BCKDHB","DBT","DLD","ECHS1","HIBADH","HIBCH","HSD17B10","IVD","MCCC1","MCCC2")
  #gene <- as.list(BCAA)
  #gene=list
  sce <- AddModuleScore(sce,features = BCAA,name = 'BCAA')
  saveRDS(sce, paste(project,Cell,"sce.Preprocess.Rds",sep = "."))
  #p1 <- FeaturePlot(sce,features = c("BCAA1"))+
  #  theme(panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  #  scale_color_gradientn(colors = c("grey90",brewer.pal(11,"RdBu")[5:1]))+
  #  labs(title="BCAA score")
  #ggsave(p1,filename = paste(project,Cell,"BCAA.score.pdf",sep = "."),height = 3,width = 3.5)
  
  mydata<- FetchData(sce,vars = c("BCAA1")) %>% mutate(Cancer=project)
  #write.csv(mydata,paste("FetchData.UMAP.BCAA.csv"))
  scRNA.AddModuleScore.BCAA <- rbind.data.frame(scRNA.AddModuleScore.BCAA,mydata)
}
write.csv(Bulk.BCAA,file = "Pancancer.Myeloid.Bulk.BCAA.csv",row.names = F)
write.csv(scRNA.ssGSEA.BCAA,file = "Pancancer.Myeloid.scRNA.ssGSEA.BCAA",row.names = F)
write.csv(scRNA.AddModuleScore.BCAA,file = "Pancancer.Myeloid.scRNA.AddModuleScore.BCAA.csv",row.names = F)


#### Figure ####
setwd("K:/TCGA/Anlysis/BCAA/Hetero/Myeloid")
library(data.table)
library(tidyverse)
library(Seurat)
library(SeuratObject)
library(magrittr) 
files <- list.files(".",pattern = "csv$")
files <- files[!str_detect(files,"Hetero")]
BCAT1.BCAT2 <- data.frame()
for (file in files) {
  project = file %>% str_remove_all("_.*") 
  print(project)
  BCAT1.BCAT2 <- read.csv(file) %>%
    dplyr::select(Sub_Cluster,Global_Cluster,Tissue,Global_UMAP_1,Global_UMAP_2,BCAT1,BCAT2,CancerType) %>%
    mutate(Cancer=project) %>%
    rbind.data.frame(BCAT1.BCAT2)
}
write.csv(BCAT1.BCAT2,"Hetero.Myeloid.BCAT1.BCAT2.csv",row.names = F)
for (project in unique(BCAT1.BCAT2$Cancer)) {
  P1 <- ggplot(BCAT1.BCAT2 %>% arrange(BCAT1) %>%
                 filter(Cancer==project) %>%
                 filter(Tissue != "L"),
               aes(Global_UMAP_1,Global_UMAP_2,
                   color=BCAT1))+
    geom_point(size=1)+
    #facet_wrap(~Tissue,ncol=3)+
    ggthemes::theme_few()+
    theme(legend.position = "bottom",
          axis.text = element_text(color = "black"),
          axis.text.y = element_text(size=9), #,angle = 90,vjust = 0.3,hjust = 1
          axis.text.x = element_text(size=9),
          #axis.ticks.x = element_blank(),
          axis.title = element_text(size=12),
          panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
    #labs(x="",y="")+
    #scale_color_gradientn(colors = brewer.pal(9,"Reds"))
    scale_color_gradientn(colours = brewer.pal(9,"Reds"))
  
  P2 <- ggplot(BCAT1.BCAT2 %>% arrange(BCAT1) %>%
                 filter(Cancer==project) %>%
                 filter(Tissue != "L"),
               aes(Global_UMAP_1,Global_UMAP_2,
                   color=Global_Cluster))+
    geom_point(size=0.5)+
    facet_wrap(~Tissue,ncol=3)+
    ggthemes::theme_few()+
    theme(legend.position = "bottom",
          axis.text = element_text(color = "black"),
          axis.text.y = element_text(size=9), #,angle = 90,vjust = 0.3,hjust = 1
          axis.text.x = element_text(size=9),
          #axis.ticks.x = element_blank(),
          axis.title = element_text(size=12),
          panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
    #labs(x="",y="")+
    #scale_color_gradientn(colors = brewer.pal(9,"Reds"))
    scale_color_manual(values = color)
  library(cowplot)
  
  ggsave(plot_grid(P1, P2, ncol = 1, align = "v"),
         filename = paste("Hetero.",project,"Myeloid.BCAT1.Cluster.pdf",sep = ""),height =7,width = 6)
}
#### Percent Expressed + Average expression ####
files <- list.files(".",pattern = "csv$")
files <- files[!str_detect(files,"Hetero")]
BCAT1.BCAT2.All <- data.frame()
for (file in files) {
  project = file %>% str_remove_all("_.*") 
  print(project)
  BCAT1.BCAT2.All <- read.csv(file) %>%
    dplyr::select(Sub_Cluster,Global_Cluster,Tissue,Global_UMAP_1,Global_UMAP_2,ACAD8:MCCC2,CancerType) %>%
    gather(SYMBOL,Expression,ACAD8:MCCC2) %>%
    mutate(Cancer=project) %>%
    rbind.data.frame(BCAT1.BCAT2.All)
}
#BCAT1.BCAT2 <- read.csv("Hetero.Myeloid.BCAT1.BCAT2.csv")
write.csv(BCAT1.BCAT2.All,"Hetero.BCAT1.BCAT2.All.csv",row.names = F)

BCAT1.BCAT2.AveExp <- BCAT1.BCAT2.All %>%
  filter(Tissue != "L") %>%
  group_by(Cancer,Tissue,Global_Cluster,SYMBOL) %>%
  summarise(AveExp = mean(Expression)) %>% ungroup()

BCAT1.BCAT2.PerExp <- BCAT1.BCAT2.All %>%
  filter(Tissue != "L") %>%
  group_by(Cancer,Tissue,Global_Cluster,SYMBOL) %>%
  mutate(TotalCount=sum(Expression>=0)) %>%
  mutate(ExpCount = sum(Expression >0)) %>%
  dplyr::select(Cancer,Tissue,Global_Cluster,TotalCount,ExpCount) %>%
  unique() %>%
  mutate(PerExp = ExpCount*100/TotalCount)

BCAT1.BCAT2 <- merge(BCAT1.BCAT2.AveExp,BCAT1.BCAT2.PerExp,
                     by=c("Cancer","Tissue","Global_Cluster","SYMBOL")) %>%
  data.frame(check.names = F) %>%
  #filter(SYMBOL=="BCAT1") %>%
  mutate(mid = paste(Cancer,Tissue,sep = "_")) %>%
  mutate(Xaxis=paste(Global_Cluster,Tissue,Cancer,sep = "_")) %>%
  arrange(Global_Cluster,mid) %>%
  mutate(Xaxis=factor(Xaxis,levels = unique(.$Xaxis)))
  
  
p1 <- ggplot(BCAT1.BCAT2,aes(x=Xaxis,y=SYMBOL,fill=AveExp))+
  geom_tile(color="white")+
  ggthemes::theme_few()+
  theme(legend.position = "bottom",
        axis.text = element_text(color = "black"),
        axis.text.y = element_text(size=9), #,angle = 90,vjust = 0.3,hjust = 1
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="",y="")+
  scale_fill_gradientn(colours = brewer.pal(9,"Reds"))

p2 <- ggplot(BCAT1.BCAT2,aes(x=Xaxis,y=SYMBOL,fill=PerExp))+
  geom_tile(color="white")+
  ggthemes::theme_few()+
  theme(legend.position = "bottom",
        axis.text = element_text(color = "black"),
        axis.text.y = element_text(size=9), #,angle = 90,vjust = 0.3,hjust = 1
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="",y="")+
  scale_fill_gradientn(colours = brewer.pal(9,"Blues"))

p3 <- ggplot(BCAT1.BCAT2 %>% 
               dplyr::select(Xaxis,Cancer) %>% unique(),
             aes(y=1,x=Xaxis))+
  geom_col(aes(fill=Cancer),width = 1)+
  ggthemes::theme_few()+
  theme(legend.position = "bottom",
        #axis.text = element_text(color = "black"),
        #axis.text.x = element_text(size=9,angle = 90,vjust = 0.5,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  #ggsci::scale_color_nejm()+
  ggthemes::scale_fill_stata()+
  #scale_fill_manual(values = rev(c("#DC0000FF","#F39B7FFF","#00A087FF")))+
  theme(legend.position = "bottom")+
  #geom_vline(xintercept = accumulate(mid$Count,sum)+0.5)+
  labs(x="",y="") 

p4 <- ggplot(BCAT1.BCAT2 %>% 
               dplyr::select(Xaxis,Global_Cluster) %>% unique(),
             aes(y=1,x=Xaxis))+
  geom_col(aes(fill=Global_Cluster),width = 1)+
  ggthemes::theme_few()+
  theme(legend.position = "bottom",
        #axis.text = element_text(color = "black"),
        #axis.text.x = element_text(size=9,angle = 90,vjust = 0.5,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  #ggsci::scale_color_nejm()+
  #ggthemes::scale_fill_stata()+
  scale_fill_manual(values = color)+
  theme(legend.position = "bottom")+
  #geom_vline(xintercept = accumulate(mid$Count,sum)+0.5)+
  labs(x="",y="")

p5 <- ggplot(BCAT1.BCAT2 %>% 
               dplyr::select(Xaxis,Tissue) %>% unique(),
             aes(y=1,x=Xaxis))+
  geom_col(aes(fill=Tissue),width = 1)+
  ggthemes::theme_few()+
  theme(legend.position = "bottom",
        #axis.text = element_text(color = "black"),
        #axis.text.x = element_text(size=9,angle = 90,vjust = 0.5,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  #ggsci::scale_color_nejm()+
  #ggthemes::scale_fill_stata()+
  scale_fill_manual(values = rev(c("#DC0000FF","#F39B7FFF","#00A087FF")))+
  theme(legend.position = "bottom")+
  #geom_vline(xintercept = accumulate(mid$Count,sum)+0.5)+
  labs(x="",y="")

p12345 <- p1 %>% aplot::insert_bottom(p2,height = 1) %>%
  aplot::insert_top(p3,height = 0.1)%>%
  aplot::insert_top(p4,height = 0.1)%>%
  aplot::insert_top(p5,height = 0.1)
ggsave(p12345,filename = "Htero.SYMBOL.Myeloid.pdf",height = 8,width = 14)

#### Example => BCAT1 => LYM => cDC1 ####
Exa.BCAT1 <- BCAT1.BCAT2 %>% filter(SYMBOL == "BCAT1") %>%
  mutate(Xaxis=paste(Cancer,Tissue,sep = "_"))
p1 <- ggplot(Exa.BCAT1,aes(x=Xaxis,y=Global_Cluster,
                     color=AveExp,size=PerExp))+
  geom_point()+
  ggthemes::theme_few()+
  theme(legend.position = "bottom",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9,angle = 90,vjust = 0.5,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
        axis.text.y = element_text(size=9),
        #axis.ticks = element_blank(),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  scale_color_gradientn(colors = brewer.pal(9,"Reds"))+
  labs(x="",y="")
ggsave(p1,filename = "Htero.SYMBOL.Myeloid.BCAT1.pdf",height = 4,width = 6)

#### Percent Expressed + Average expression => Samples ####

library(data.table)
library(tidyverse)
library(Seurat)
library(SeuratObject)
library(magrittr)
files <- list.files(".",pattern = "csv$")
files <- files[!str_detect(files,"Hetero") & !str_detect(files,"Pancancer")]
BCAT1.BCAT2.All <- data.frame()
for (file in files) {
  project = file %>% str_remove_all("_.*") 
  print(project)
  BCAT1.BCAT2.All <- read.csv(file) %>%
    dplyr::select(Sub_Cluster,Global_Cluster,Sample,Tissue,Global_UMAP_1,Global_UMAP_2,ACAD8:MCCC2,CancerType) %>%
    gather(SYMBOL,Expression,ACAD8:MCCC2) %>%
    mutate(Cancer=project) %>%
    rbind.data.frame(BCAT1.BCAT2.All)
}
#BCAT1.BCAT2 <- read.csv("Hetero.Myeloid.BCAT1.BCAT2.csv")
write.csv(BCAT1.BCAT2.All,"Hetero.BCAT1.BCAT2.All.Samples.csv",row.names = F)


BCAT1.BCAT2.AveExp <- BCAT1.BCAT2.All %>%
  filter(Tissue != "L") %>%
  group_by(Cancer,Sample,Tissue,Global_Cluster,SYMBOL) %>%
  summarise(AveExp = mean(Expression)) %>% ungroup()

BCAT1.BCAT2.PerExp <- BCAT1.BCAT2.All %>%
  filter(Tissue != "L") %>%
  group_by(Cancer,Tissue,Sample,Global_Cluster,SYMBOL) %>%
  mutate(TotalCount=sum(Expression>=0)) %>%
  mutate(ExpCount = sum(Expression >0)) %>%
  dplyr::select(Cancer,Tissue,Global_Cluster,TotalCount,ExpCount) %>%
  unique() %>%
  mutate(PerExp = ExpCount*100/TotalCount)

BCAT1.BCAT2 <- merge(BCAT1.BCAT2.AveExp,BCAT1.BCAT2.PerExp,
                     by=c("Cancer","Tissue","Sample","Global_Cluster","SYMBOL")) %>%
  data.frame(check.names = F) %>%
  #filter(SYMBOL=="BCAT1") %>%
  mutate(mid = paste(Cancer,Tissue,sep = "_")) %>%
  mutate(Xaxis=paste(Global_Cluster,Tissue,Cancer,sep = "_")) %>%
  arrange(Global_Cluster,mid) %>%
  mutate(Xaxis=factor(Xaxis,levels = unique(.$Xaxis)))
write.csv(BCAT1.BCAT2,"Hetero.BCAT1.BCAT2.All.Samples.AveExp.PerExp.csv",row.names = F)

#### Example => BCAT1 => LYM => cDC1 ####
Exa.BCAT1 <- BCAT1.BCAT2 %>% filter(SYMBOL == "BCAT1") %>%
  mutate(Xaxis=paste(Cancer,Tissue,sep = "_")) %>%
  filter(!Global_Cluster %in% c("Monolike","Myeloid"))
p1 <- ggplot(Exa.BCAT1,aes(x=Cancer,y=AveExp))+
  geom_boxplot(aes(color=Tissue))+
  facet_wrap(~Global_Cluster,scales = "free_y",nrow = 2)+
  ggthemes::theme_few()+
  theme(legend.position = "bottom",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9,angle = 90,vjust = 0.5,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
        axis.text.y = element_text(size=9),
        #axis.ticks = element_blank(),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  #scale_color_gradientn(colors = brewer.pal(9,"Reds"))+
  ggsci::scale_fill_npg()+
  labs(x="")
ggsave(p1,filename = "Htero.SYMBOL.Myeloid.BCAT1.Boxplot.pdf",height = 5,width = 7)

#### BCAA score psudobulk => Figure ####
setwd("K:/TCGA/Anlysis/BCAA/Hetero/Myeloid/")
data <- read.csv("Pancancer.Myeloid.Bulk.BCAA.csv") %>%
  mutate(CellType=str_remove(cellType,"N_")) %>%
  mutate(CellType=str_remove(CellType,"P_")) %>%
  mutate(CellType=str_remove(CellType,"T_")) %>%
  mutate(Loc = str_remove(cellType,"_.*")) %>%
  mutate(Yaxis=paste(CellType,Loc,sep = "_"))

p1 <- ggplot(data,aes(y=Cancer,x=Yaxis))+
  geom_tile(aes(fill=BCAA),color="white")+
  ggthemes::theme_few()+
  theme(legend.position = "bottom",
        axis.text = element_text(color = "black"),
        axis.text.y = element_text(size=9), #,angle = 90,vjust = 0.3,hjust = 1
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="",y="")+
  scale_fill_gradientn(colors = brewer.pal(9,"Reds")) # YlOrRd Blues BuGn

library(ggsci)
color = c(pal_d3("category20")(20),pal_d3("category20b")(20),pal_d3("category20c")(20),pal_d3("category10")(10))
mid <- data %>% dplyr::select(Yaxis,CellType) %>% unique() %>%
  group_by(CellType) %>% summarise(Count=n())

p2 <- ggplot(data %>% dplyr::select(Yaxis,CellType) %>% unique(),
             aes(y=1,x=Yaxis))+
  geom_col(aes(fill=CellType),width = 1)+
  #ggthemes::theme_few()+
  #theme_void()+
  ggthemes::theme_few()+
  theme(legend.position = "bottom",
        #axis.text = element_text(color = "black"),
        #axis.text.x = element_text(size=9,angle = 90,vjust = 0.5,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  geom_vline(xintercept = accumulate(mid$Count,sum)+0.5)+
  scale_fill_manual(values = color)+
  labs(x="",y="")
#scale_fill_gradientn(colors = brewer.pal(9,"Reds")) # YlOrRd Blues BuGn
mid <- data %>% dplyr::select(Loc,CellType) %>% unique() %>%
  group_by(CellType) %>% summarise(Count=n())


p3 <- ggplot(data %>% dplyr::select(Yaxis,Loc) %>% unique(),
             aes(y=1,x=Yaxis))+
  geom_col(aes(fill=Loc),width = 1)+
  ggthemes::theme_few()+
  theme(legend.position = "bottom",
        #axis.text = element_text(color = "black"),
        #axis.text.x = element_text(size=9,angle = 90,vjust = 0.5,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  #ggsci::scale_color_nejm()+
  scale_fill_manual(values = rev(c("#DC0000FF","#F39B7FFF","#00A087FF")))+
  theme(legend.position = "bottom")+
  geom_vline(xintercept = accumulate(mid$Count,sum)+0.5)+
  labs(x="",y="")
#scale_fill_manual(values = color)
p123 <- p1 %>% aplot::insert_top(p2,height = 0.1) %>%
  aplot::insert_top(p3,height = 0.1)
ggsave(p123,filename = "Hetero.Myeloid.BCAAscore.ssGSEA.pdf",height = 3,width = 18)
ggsave(p123,filename = "Hetero.Myeloid.BCAAscore.ssGSEA.Legend.pdf",height = 10,width = 12)

#### BCAA score AUCell => Stat+Figure ####
setwd("K:/TCGA/Anlysis/BCAA/Hetero/Myeloid")
files <- list.files(".",pattern = "^AUC.*.csv")
files <- files[str_detect(files,"BCAAscore")]
for (file in files) {
  project <- (file %>% str_split("\\.") %>% unlist())[2]
  print(project)
  data <- read.csv(file) %>% 
    mutate(CellType=str_extract_all(MajorCluster,"_.*_")) %>%
    mutate(CellType=str_remove_all(CellType,"_"))
  p1 <- ggplot(data,aes(UMAP1,UMAP2,color=CellType))+
    geom_point(size=1)+
    facet_wrap(~tissue)+
    ggthemes::theme_few()+
    theme_test(12)+
    theme(#legend.position = "bottom",
          axis.text = element_text(color = "black"),
          axis.text.y = element_text(size=9), #,angle = 90,vjust = 0.3,hjust = 1
          axis.text.x = element_text(size=9),
          #axis.ticks.x = element_blank(),
          axis.title = element_text(size=12),
          panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
    #labs(x="",y="")
    ggsci::scale_color_npg()+
    guides(color = guide_legend(override.aes = list(size=4)))
    
  library(RColorBrewer)
  p2 <- ggplot(data %>% arrange(AUC),aes(UMAP1,UMAP2,color=AUC))+
    geom_point(size=1)+
    facet_wrap(~tissue)+
    ggthemes::theme_few()+
    theme_test(12)+
    theme(#legend.position = "bottom",
      axis.text = element_text(color = "black"),
      axis.text.y = element_text(size=9), #,angle = 90,vjust = 0.3,hjust = 1
      axis.text.x = element_text(size=9),
      #axis.ticks.x = element_blank(),
      axis.title = element_text(size=12),
      panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
    #labs(x="",y="")
    scale_color_gradientn(colours = brewer.pal(9,"Reds"))
  library(cowplot)
  p <- plot_grid(p1, p2, ncol = 1, align = "v")
  ggsave(p,filename = paste("../../Hetero.BCAA.AUCell.Score.Myeloid",project,"Figure.pdf",sep = "."),
         height = 4,width = 6)
}

data <- read.csv("AUCell.KIDNEY.BCAAscore.csv") %>% 
  mutate(CellType=str_extract_all(MajorCluster,"_.*_")) %>%
  mutate(CellType=str_remove_all(CellType,"_"))
library(ggpubr)
p123 <- ggplot(data%>% filter(CellType %in% c("cDC3","pDC")), 
       aes(tissue,AUC,color=tissue))+
  facet_wrap(~CellType)+
  geom_boxplot(linetype="dashed",color="black")+
  stat_boxplot(aes(ymin=..lower..,ymax=..upper..,
                   fill=tissue),
               color="black")+
  stat_boxplot(geom = "errorbar",aes(ymin=..ymax..,color=tissue),width=0.3,size=1)+
  stat_boxplot(geom = "errorbar",aes(ymax=..ymin..,color=tissue),width=0.3,size=1)+
  stat_compare_means(label = "p.format",label.x = 1.2,method = "t.test")+
  ggthemes::theme_few()+ theme_test(12)+
  #scale_size_continuous("-log10(Fisher's P)",range = c(1,4))+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        #axis.text.x = element_text(size=9),#,angle = 90,vjust = 0.5,hjust = 1
        #axis.text.y = element_text(size=9),
        #axis.ticks.x = element_blank(),
        #axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  theme(strip.background=element_rect(fill=c("#6F99AD7F",
                                             "#EE4C977F")))+
  scale_fill_manual(values = c("#E64B35FF","#4DBBD5FF"))+
  labs(x="")

ggsave(p123,filename = "AUCell.KIDNEY.Example.Boxplot.pdf",
       height = 3.5,width = 2.5)

files <- list.files(".",pattern = "^AUC.*.csv")
AUCell.Stat <- data.frame()
for (file in files) {
  project <- (file %>% str_split("\\.") %>% unlist())[2]
  print(project)
  data <- read.csv(file) %>% 
    mutate(CellType=str_extract_all(MajorCluster,"_.*_")) %>%
    mutate(CellType=str_remove_all(CellType,"_"))
  if (project %in% c("LYM")) {
    middata <- data %>% group_by(tissue,CellType) %>% summarise(Count=n()) %>%
      arrange(CellType) %>% ungroup() %>%
      spread(tissue,Count) %>% filter(P>=3) %>% filter(T>=3)
    AUCell.Stat <- data %>%
      filter(CellType %in% middata$CellType) %>%
      group_by(CellType) %>%
      rstatix::t_test(AUC~tissue,detailed = T) %>%
      rbind.data.frame(AUCell.Stat)
  }else if (project %in% c("PAAD","THCA")){
    middata <- data %>% group_by(tissue,CellType) %>% summarise(Count=n()) %>%
      arrange(CellType) %>% ungroup() %>%
      spread(tissue,Count) %>% filter(N>=3) %>% filter(T>=3)
    AUCell.Stat <- data %>%
      filter(CellType %in% middata$CellType) %>%
      group_by(CellType) %>%
      rstatix::t_test(AUC~tissue,detailed = T) %>%
      rbind.data.frame(AUCell.Stat)
  }else {
    AUCell.Stat <- data %>%
      #filter(CellType %in% middata$CellType) %>%
      group_by(CellType) %>%
      rstatix::t_test(AUC~tissue,detailed = T) %>%
      rbind.data.frame(AUCell.Stat)
  }

}
write.csv(AUCell.Stat,"AUCell.BCAA.Meyloid.Stat.csv",row.names = F)

#### BCAA SYMBOLs => T test ####
data <- read.csv("BCAA.SYMBOLS.Myeloid.Expression.csv") %>% 
  mutate(CellType=str_extract_all(MajorCluster,"_.*_")) %>%
  mutate(CellType=str_remove_all(CellType,"_"))

BCAA_SYMBOLS.Stat <- data.frame()
for (project in unique(data$Cancer)) {
  mid <- data %>% filter(Cancer==project) %>% 
    gather(SYMBOL,TPM,ACAD8:MCCC2)
  if (project %in% c("LYM")) {
    middata <- mid %>% group_by(tissue,CellType,SYMBOL) %>% summarise(Count=n()) %>%
      arrange(CellType) %>% ungroup() %>%
      spread(tissue,Count) %>% filter(P>=3) %>% filter(T>=3)
    BCAA_SYMBOLS.Stat <- mid %>%
      filter(CellType %in% middata$CellType) %>%
      group_by(SYMBOL,CellType) %>%
      rstatix::t_test(TPM~tissue,detailed = T) %>%
      mutate(Cancer=project) %>%
      rbind.data.frame(BCAA_SYMBOLS.Stat)
  }else if (project %in% c("PAAD","THCA")){
    middata <- mid %>% group_by(tissue,CellType,SYMBOL) %>% summarise(Count=n()) %>%
      arrange(CellType) %>% ungroup() %>%
      spread(tissue,Count) %>% filter(N>=3) %>% filter(T>=3)
    BCAA_SYMBOLS.Stat <- mid %>%
      filter(CellType %in% middata$CellType) %>%
      group_by(SYMBOL,CellType) %>%
      rstatix::t_test(TPM~tissue,detailed = T) %>%
      mutate(Cancer=project) %>%
      rbind.data.frame(BCAA_SYMBOLS.Stat)
  }else {
    BCAA_SYMBOLS.Stat <- mid %>%
      #filter(CellType %in% middata$CellType) %>%
      group_by(SYMBOL,CellType) %>%
      rstatix::t_test(TPM~tissue,detailed = T) %>%
      mutate(Cancer=project) %>%
      rbind.data.frame(BCAA_SYMBOLS.Stat)
  }
}
write.csv(BCAA_SYMBOLS.Stat,file = "BCAA_SYMBOLS.Meyloid.Stat.csv",row.names = F)

###### Heterogeneity => CCLE ##############
pkgs <- c(
  "Seurat", "SeuratWrappers", "ggplot2", "batchelor",
  "dplyr", "optparse", "reshape2", "data.table", "magrittr"
)
lapply(pkgs, function(x) require(package = x, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)) #nolint
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
suppressMessages(suppressWarnings(library(ggsci)))
suppressMessages(suppressWarnings(library(harmony)))
suppressMessages(suppressWarnings(library(scater)))
suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(readxl)))
suppressMessages(suppressWarnings(library(SingleCellExperiment)))
suppressMessages(suppressWarnings(library(scRNAseq)))
suppressMessages(suppressWarnings(library(scran)))
suppressMessages(suppressWarnings(library(Matrix)))
suppressMessages(suppressWarnings(library(scater)))
suppressMessages(suppressWarnings(library(ggthemes)))
suppressMessages(suppressWarnings(library(cowplot)))

data <- read.table("CPM_data.txt",header = T,row.names = 1) 
sce <- CreateSeuratObject(data)
sce[["percent.mt"]] <- PercentageFeatureSet(
  sce,
  pattern = "^MT-"
)
sce %<>% subset(subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20) %>% #nolint
  NormalizeData(normalization.method = "LogNormalize") %>%
  FindVariableFeatures(selection.method = "vst")
sce %<>% ScaleData()
BCAA <- c("ACAD8","ACADSB","ACAT1","ALDH6A1","AUH","BCAT1","BCAT2","BCKDHA","BCKDHB","DBT","DLD","ECHS1","HIBADH","HIBCH","HSD17B10","IVD","MCCC1","MCCC2")
gene <- as.list(BCAA)
sce <- AddModuleScore(sce,features = gene,name = 'BCAA')
mydata<- FetchData(sce,vars = c("BCAA1"))
write.csv(mydata,"FetchData.BCAAscore.csv")

BCAA.sce <- sce[BCAA,]
saveRDS(BCAA.sce,file = "CCLE.BCAA.Rds")
#### Figure ####
setwd("K:/TCGA/Pan-cancer cell line heterogeneity")
files <- list.files(".",pattern = "txt$")
files <- files[str_detect(files,"tSNE")]
Exp.Data <- read.csv("BCAA.CPM.csv",row.names = 1,check.names = F) %>%
  t()
Metadata <- read.table("Metadata.txt",sep="\t",header = T,check.names = F)
CCLE.SYMBOL=data.frame()
for (file in files) {
  project=file %>% str_remove_all("tSNE_") %>%
    str_remove_all("\\..*")
  CCLE.SYMBOL <- read.table(file,skip = 1,sep = "\t",header = T) %>%
    set_colnames(c("NAME","tSNE1","tSNE2")) %>%
    na.omit() %>%
    merge(Metadata,by="NAME") %>%
    dplyr::select(NAME,tSNE1,tSNE2,Cell_line) %>%
    mutate(Cancer=project) %>%
    mutate(Cell_line = str_remove_all(Cell_line,"_.*")) %>%
    merge(Exp.Data %>% data.frame(check.names = F) %>%
            rownames_to_column("NAME") %>%
            mutate(NAME=str_replace_all(NAME,"\\.","-")),by="NAME") %>%
    rbind.data.frame(CCLE.SYMBOL)
}
write.csv(CCLE.SYMBOL,"CCLE.SYMBOLS.BCAA.CPM.csv",row.names = F)

DesDir <- "../Anlysis/BCAA/Hetero/CCLE/"
for (project in unique(CCLE.SYMBOL$Cancer)) {
  print(project)
  middata <- CCLE.SYMBOL %>% #filter(SYMBOL=="BCAT1") %>%
    filter(Cancer==project)
  p1 <- ggplot(middata %>% arrange(BCAT1),aes(tSNE1,tSNE2,color=BCAT1))+
    geom_point(size=1)+
    #facet_wrap(~Tissue,ncol=3)+
    ggthemes::theme_few()+
    theme(#legend.position = "bottom",
          axis.text = element_text(color = "black"),
          axis.text.y = element_text(size=9), #,angle = 90,vjust = 0.3,hjust = 1
          axis.text.x = element_text(size=9),
          #axis.ticks.x = element_blank(),
          axis.title = element_text(size=12),
          panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
    #labs(x="",y="")+
    #scale_color_gradientn(colors = brewer.pal(9,"Reds"))
    scale_color_gradientn(colours = brewer.pal(9,"Reds")[1:9])
  p2 <- ggplot(middata %>% arrange(BCAT1),
               aes(tSNE1,tSNE2,color=Cell_line))+
    geom_point(size=0.5)+
    #facet_wrap(~Tissue,ncol=3)+
    ggthemes::theme_few()+
    theme(#legend.position = "bottom",
          axis.text = element_text(color = "black"),
          axis.text.y = element_text(size=9), #,angle = 90,vjust = 0.3,hjust = 1
          axis.text.x = element_text(size=9),
          #axis.ticks.x = element_blank(),
          axis.title = element_text(size=12),
          panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
    #labs(x="",y="")+
    #scale_color_gradientn(colors = brewer.pal(9,"Reds"))
    scale_color_manual(values = color)
  library(cowplot)
  
  ggsave(p1 %>% aplot::insert_right(p2,width = 1), #plot_grid(p1, p2, ncol = 2, align = "h"),
         filename = paste(DesDir,"Hetero.",project,"CellLine.BCAT1.pdf",sep = ""),
         height =3,width = 7)
}

#### Cell_Line => SYMBOLs ####
CCLE.SYMBOL <- read.csv("CCLE.SYMBOLS.BCAA.CPM.csv")
CCLE.SYMBOL2 <- CCLE.SYMBOL %>%
  gather(SYMBOL,Exp,ACAD8:MCCC2) %>%
  filter(!Cancer %in% c("Fibroblast","Gallbladder_Cancer")) %>%
  group_by(Cancer,Cell_line,SYMBOL) %>%
  summarise(Exp.M = mean(Exp)) %>%
  ungroup()
CCLE.SYMBOL2.Test <- CCLE.SYMBOL %>%
  gather(SYMBOL,Exp,ACAD8:MCCC2) %>%
  filter(!Cancer %in% c("Fibroblast","Gallbladder_Cancer")) %>%
  group_by(Cancer,SYMBOL) %>%
  rstatix::kruskal_test(Exp ~ Cell_line) %>%
  add_significance("p")

CCLE.SYMBOL2.Final <- merge(CCLE.SYMBOL2,CCLE.SYMBOL2.Test,by=c("Cancer","SYMBOL")) %>%
  arrange(Cancer,Cell_line) %>%
  mutate(Cell_line = factor(Cell_line,levels = unique(.$Cell_line))) 
p1 <- ggplot()+
  geom_tile(data = CCLE.SYMBOL2.Final %>%
              filter(!Cancer %in% c("Fibroblast","Gallbladder_Cancer")),
            aes(Cell_line,SYMBOL,fill=Exp.M),
            color="white")+
  geom_point(data=CCLE.SYMBOL2.Final %>%
               filter(p.signif!="ns"),
             aes(Cell_line,SYMBOL,fill=Exp.M),
             shape=1,color="grey90",size=0.1)+
  ggthemes::theme_few()+
  theme(legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_blank(), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    axis.ticks.x = element_blank(),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  #labs(x="",y="")+
  #scale_color_gradientn(colors = brewer.pal(9,"Reds"))
  scale_fill_gradientn(colours = brewer.pal(9,"Reds")[1:9])+
  labs(x="",y="",fill="Ave Exp")
library(rstatix)

p2 <- ggplot(data = CCLE.SYMBOL2.Final %>%
         filter(!Cancer %in% c("Fibroblast","Gallbladder_Cancer")) %>%
         dplyr::select(Cancer,Cell_line) %>% unique(),
       aes(Cell_line,1,fill=Cancer))+
  #geom_tile(color="white")+
  geom_col(width = 1)+
  #geom_text(aes(label=p.signif),color="white",angle=90)+
  ggthemes::theme_few()+
  theme(legend.position = "bottom",
    #axis.text = element_text(color = "black"),
    axis.text.x = element_blank(), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  #labs(x="",y="")+
  #scale_color_gradientn(colors = brewer.pal(9,"Reds"))
  #scale_fill_gradientn(colours = brewer.pal(9,"Reds")[1:9])+
  scale_fill_manual(values = color)+
  labs(x="",y="",color="Cancer")

p12.2 <- p1 %>% aplot::insert_top(p2,height = 0.1)
ggsave(p12.2,filename = paste(DesDir,"CCLE.SYMBOL.BCAA.pdf",sep = ""),
       height = 3,width = 18)


#### BCAA score => Figure ####
setwd("K:/TCGA/Anlysis/BCAA/Hetero/CCLE")
library(tidyverse)
library(data.table)
BCAA.ID <- read.csv("../../BCAA.ID.csv")
metadata <- fread("Metadata.txt")
data <- read.csv("CCLE.AUCell.SYMBOLS.csv") %>%
  dplyr::select(c("X","AUC")) %>%
  merge(metadata %>% mutate(NAME=str_replace_all(NAME,"-",".")),
        by.x="X",by.y="NAME") %>%
  mutate(celline=str_remove_all(Cell_line,"_.*"))

AUCell.Stat <- data.frame()
for (project in unique(data$Cancer_type)) {
  mid <- data %>% filter(Cancer_type == project) #%>%
    #gather(SYMBOL,TPM,ACAD8:MCCC2)
  if (length(unique(mid$celline)) > 1) {
    AUCell.Stat <- mid %>%
      #group_by(SYMBOL) %>%
      rstatix::t_test(AUC~celline,detailed = T) %>%
      mutate(Cancer=project) %>%
      plyr::rbind.fill(AUCell.Stat)
  }
}
write.csv(AUCell.Stat,file="AUCell.BCAA.CCLE.Stat.csv",row.names = F)

metadata <- fread("Metadata.txt")
data <- read.csv("CCLE.AUCell.SYMBOLS.csv") %>%
  dplyr::select(c("X","AUC")) %>%
  merge(metadata %>% mutate(NAME=str_replace_all(NAME,"-",".")),
        by.x="X",by.y="NAME") %>%
  mutate(celline=str_remove_all(Cell_line,"_.*"))
# unique(data$Cancer_type)
for (project in unique(data$Cancer_type)) {
  tsne.filename <- project %>%
    str_replace_all(" ","_") %>%
    str_replace_all("/","_") %>%
    paste("tSNE_",.,".txt",sep = "")
  figure.filename <- project %>% 
    str_replace_all(" ","_") %>%
    str_replace_all("/","_") %>%
    paste("../../Hetero.BCAA.AUCell.CCLE.",.,".pdf",sep = "")
  tsne <- fread(tsne.filename,skip = 2) %>%
    mutate(V1=str_replace_all(V1,"-",".")) %>%
    merge(data,by.x="V1",by.y="X")
  print(project)
  print(dim(tsne))
  print(length(unique(tsne$celline)))
  library(ggsci)
  color = c(pal_d3("category20")(20),pal_d3("category20b")(20),pal_d3("category20c")(20),pal_d3("category10")(10))
  if (length(unique(tsne$celline)) <= 70) {
    p1 <- ggplot(tsne %>% arrange(AUC),
                 aes(V2,V3,color=AUC))+
      geom_point(size=1)+
      ggthemes::theme_few()+ 
      theme_test(12)+
      #scale_size_continuous("-log10(Fisher's P)",range = c(1,4))+
      theme(#legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9),#,angle = 90,vjust = 0.5,hjust = 1
        axis.text.y = element_text(size=9),
        #axis.ticks.x = element_blank(),
        #axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
      scale_color_gradientn(colors = brewer.pal(9,"Reds"))+
      labs(x="tsne1",y="tsne2")
    p2 <- ggplot(tsne %>% arrange(AUC),
                 aes(V2,V3,color=celline))+
      geom_point(size=1)+
      ggthemes::theme_few()+ 
      theme_test(12)+
      #scale_size_continuous("-log10(Fisher's P)",range = c(1,4))+
      theme(#legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9),#,angle = 90,vjust = 0.5,hjust = 1
        axis.text.y = element_text(size=9),
        #axis.ticks.x = element_blank(),
        #axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
      #scale_color_gradientn(colors = brewer.pal(9,"Reds"))+
      #ggthemes::scale_color_stata()+
      scale_color_manual("Cell line",values = color)+
      labs(x="Tsne1",y="Tsne2",color="Cell line")+
      guides(color = guide_legend(override.aes = list(size=4)))
    p <- p1+p2
    ggsave(p,filename = figure.filename,height = 3.5,width = 9)
  }
}

#### BCAA SYMBOLs => T test ####
setwd("K:/TCGA/Anlysis/BCAA/Hetero/CCLE")
library(tidyverse)
library(data.table)
BCAA.ID <- read.csv("../../BCAA.ID.csv")
metadata <- fread("Metadata.txt")
data <- read.csv("CCLE.AUCell.SYMBOLS.csv") %>%
  dplyr::select(c("X",BCAA.ID$SYMBOL)) %>%
  merge(metadata %>% mutate(NAME=str_replace_all(NAME,"-",".")),
        by.x="X",by.y="NAME") %>%
  mutate(celline=str_remove_all(Cell_line,"_.*"))

SYMBOLS.Stat <- data.frame()
for (project in unique(data$Cancer_type)) {
  mid <- data %>% filter(Cancer_type == project) %>%
    gather(SYMBOL,TPM,ACAD8:MCCC2)
  if (length(unique(mid$celline)) > 1) {
    SYMBOLS.Stat <- mid %>%
      group_by(SYMBOL) %>%
      rstatix::t_test(TPM~celline,detailed = T) %>%
      mutate(Cancer=project) %>%
      plyr::rbind.fill(SYMBOLS.Stat)
  }
}
write.csv(SYMBOLS.Stat,file="BCAA_SYMBOLS.CCLE.Stat.csv",row.names = F)

######## Pancancer => scRNA #########
sce <- readRDS("sce.FindCluster.Rds")
library(ggsci)
color = c(pal_d3("category20")(20),pal_d3("category20b")(20),pal_d3("category20c")(20),pal_d3("category10")(10))

sce2 <- subset(sce,celltype != "Plasma" & celltype != "undefined" & celltype != "Erythrocyte")
#sce2$Axis <- paste(sce2$celltype,sce2$tissue,sce2$group,sep="_")
sce2$Axis <- paste(sce2$celltype,sce2$tissue,sce2$group,sep="_")
BCAA <- c("ACAD8","ACADSB","ACAT1","ALDH6A1","AUH","BCAT1","BCAT2","BCKDHA","BCKDHB","DBT","DLD","ECHS1","HIBADH","HIBCH","HSD17B10","IVD","MCCC1","MCCC2")
p1 <- DotPlot(sce2, features = BCAA,assay='RNA',group.by ="Axis")+ 
  coord_flip()+scale_color_gradientn(colors = brewer.pal(9,"Reds"))+
  scale_size_continuous(range = c(0.1,3))+
  theme(axis.text.x = element_text(size=6,angle=90,vjust = 0.3,hjust = 1,color = "black"),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="",y="")
ggsave(p1,filename = "Dotplot.BCAA.SYMBOLS.pdf",height = 6,width = 14)

sce3 <- subset(sce2,tissue == "PDAC")
P1<-FeaturePlot(object = sce3, features = c("BCAT1"),reduction='umap',
                split.by = "group",
                cols = c('grey','red'))
ggsave('Featureplot.PDAC.BCAT1.pdf',P1,height=4,width=11)

sce4 <- subset(sce2,celltype == "Epithelium")
sce4$Axis <- paste(sce4$cluster,sce4$tissue,sce4$group,sep="_")
p1 <- DotPlot(sce4, features = BCAA,assay='RNA',group.by ="Axis")+ 
  coord_flip()+scale_color_gradientn(colors = brewer.pal(9,"Reds"))+
  scale_size_continuous(range = c(0.1,3))+
  theme(axis.text.x = element_text(size=6,angle=90,vjust = 0.3,hjust = 1,color = "black"),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="",y="")
ggsave(p1,filename = "Dotplot.BCAA.SYMBOLS.Epithelium.pdf",height = 6,width = 14)

#### BCAA SYMBOLs => T test ####
setwd("K:/TCGA/Anlysis/BCAA/Hetero/GSE210347")
BCAA.ID <- read.csv("../../BCAA.ID.csv")
files <- list.files(".",pattern = "AUCell")
UMAP.D <- read.csv("embed_umap.csv")
data.T <- data.frame()
data.T2 <- data.frame()
for (file in files) {
  mid <- read.csv(file) %>% 
    dplyr::select(c("X",BCAA.ID$SYMBOL,cluster,celltype,group,tissue)) %>%
    gather(SYMBOL,TPM,ACAD8:MCCC2) %>% mutate(TPM=log2(TPM+1)) %>% 
    filter(!celltype %in% c("Erythrocyte","undefined")) 
  if (length(mid$group %>% unique()) > 1) {
    data.T <- mid %>%
      group_by(celltype,SYMBOL) %>%
      rstatix::t_test(TPM~group,detailed = T) %>%
      mutate(tissue=unique(mid$tissue)) %>%
      plyr::rbind.fill(data.T)
  }
  data.T2 <- mid %>% group_by(group,SYMBOL) %>% 
    rstatix::t_test(TPM~celltype,detailed = T) %>%
    mutate(tissue=unique(mid$tissue)) %>%
    plyr::rbind.fill(data.T2)
}
write.csv(data.T,file = "GSE210347.SYMBOLS.Ttest.CellType.Loction.csv",row.names = F)
write.csv(data.T2,file = "GSE210347.SYMBOLS.Ttest.Loction.CellType.csv",row.names = F)

data.T2 <- data.frame()
for (file in files) {
  mid <- read.csv(file) %>% 
    dplyr::select(c("X",BCAA.ID$SYMBOL,cluster,celltype,group,tissue)) %>%
    gather(SYMBOL,TPM,ACAD8:MCCC2) %>% mutate(TPM=log2(TPM+1)) %>% 
    filter(!celltype %in% c("Erythrocyte","undefined")) %>% 
    filter(celltype == "Epithelium")
  mid.count <- mid %>% group_by(group,SYMBOL,cluster) %>% 
    summarise(Count=n()) %>% 
    data.frame() %>% dplyr::select(-SYMBOL) %>% unique() %>% filter(Count >= 3)
  data.T2 <- mid %>% merge(mid.count,by=c("cluster","group")) %>%
    group_by(group,SYMBOL) %>% 
    rstatix::t_test(TPM~cluster,detailed = T) %>%
    mutate(tissue=unique(mid$tissue)) %>%
    plyr::rbind.fill(data.T2)
}
write.csv(data.T2,file = "GSE210347.SYMBOLS.Ttest.Epithelium.Loction.Cluster.csv",row.names = F)

data <- read.csv("GSE210347.PDAC.AUCell.SYMBOLS.csv") %>%
  merge(UMAP.D,by="X")
library(RColorBrewer)
p1 <- ggplot(data %>% filter(!celltype %in% c("Erythrocyte","undefined")) %>%
               arrange(BCAT1),
       aes(UMAP_1,UMAP_2,color=BCAT1))+
  geom_point(size=1)+
  facet_wrap(~group)+
  ggthemes::theme_few()+ 
  theme_test(12)+
  #scale_size_continuous("-log10(Fisher's P)",range = c(1,4))+
  theme(#legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9),#,angle = 90,vjust = 0.5,hjust = 1
        axis.text.y = element_text(size=9),
        #axis.ticks.x = element_blank(),
        #axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  scale_color_gradientn(colors = brewer.pal(9,"Reds"))

p2 <- ggplot(data %>% filter(!celltype %in% c("Erythrocyte","undefined")),# 
       aes(UMAP_1,UMAP_2,color=celltype))+
  geom_point(size=1)+
  facet_wrap(~group)+
  ggthemes::theme_few()+ 
  theme_test(12)+
  #scale_size_continuous("-log10(Fisher's P)",range = c(1,4))+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9),#,angle = 90,vjust = 0.5,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    #axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  ggsci::scale_color_npg()+
  guides(color = guide_legend(override.aes = list(size=5)))

library(patchwork)
p <- p1/p2
ggsave(p,filename = "../../Hetero.GSE210347.SYMBOLS.BCAT1.PDAC.pdf",
       height = 7.5,width = 10)

data <- read.csv("GSE210347.ICC.AUCell.SYMBOLS.csv")%>%
  merge(UMAP.D,by="X")
library(RColorBrewer)
p1 <- ggplot(data %>% filter(!celltype %in% c("Erythrocyte","undefined")) %>%
               arrange(BCAT1),
             aes(UMAP_1,UMAP_2,color=BCAT1))+
  geom_point(size=1)+
  facet_wrap(~group)+
  ggthemes::theme_few()+ 
  theme_test(12)+
  #scale_size_continuous("-log10(Fisher's P)",range = c(1,4))+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9),#,angle = 90,vjust = 0.5,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    #axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  scale_color_gradientn(colors = brewer.pal(9,"Reds"))

p2 <- ggplot(data %>% filter(!celltype %in% c("Erythrocyte","undefined")),# 
             aes(UMAP_1,UMAP_2,color=celltype))+
  geom_point(size=1)+
  facet_wrap(~group)+
  ggthemes::theme_few()+ 
  theme_test(12)+
  #scale_size_continuous("-log10(Fisher's P)",range = c(1,4))+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9),#,angle = 90,vjust = 0.5,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    #axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  ggsci::scale_color_npg()+
  guides(color = guide_legend(override.aes = list(size=5)))

library(patchwork)
p <- p1/p2
ggsave(p,filename = "../../Hetero.GSE210347.SYMBOLS.BCAT1.ICC.pdf",
       height = 7.5,width = 8.5)

data <- read.csv("GSE210347.colorectal.AUCell.SYMBOLS.csv")%>%
  merge(UMAP.D,by="X")
library(RColorBrewer)
p1 <- ggplot(data %>% filter(!celltype %in% c("Erythrocyte","undefined")) %>%
               arrange(BCAT1),
             aes(UMAP_1,UMAP_2,color=BCAT1))+
  geom_point(size=1)+
  facet_wrap(~group)+
  ggthemes::theme_few()+ 
  theme_test(12)+
  #scale_size_continuous("-log10(Fisher's P)",range = c(1,4))+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9),#,angle = 90,vjust = 0.5,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    #axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  scale_color_gradientn(colors = brewer.pal(9,"Reds"))

p2 <- ggplot(data %>% filter(!celltype %in% c("Erythrocyte","undefined")),# 
             aes(UMAP_1,UMAP_2,color=celltype))+
  geom_point(size=1)+
  facet_wrap(~group)+
  ggthemes::theme_few()+ 
  theme_test(12)+
  #scale_size_continuous("-log10(Fisher's P)",range = c(1,4))+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9),#,angle = 90,vjust = 0.5,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    #axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  ggsci::scale_color_npg()+
  guides(color = guide_legend(override.aes = list(size=5)))

library(patchwork)
p <- p1/p2
ggsave(p,filename = "../../Hetero.GSE210347.SYMBOLS.BCAT1.Colorectal.pdf",
       height = 7.5,width = 10)

#### BCAA score => AUCell => T test => Figure ####
setwd("K:/TCGA/Anlysis/BCAA/Hetero/GSE210347")
BCAA.ID <- read.csv("../../BCAA.ID.csv")
files <- list.files(".",pattern = "AUCell")
UMAP.D <- read.csv("embed_umap.csv")
data.T <- data.frame()
data.T2 <- data.frame()
for (file in files) {
  mid <- read.csv(file) %>% 
    dplyr::select(c("X",AUC,cluster,celltype,group,tissue)) %>%
    #gather(SYMBOL,TPM,ACAD8:MCCC2) %>% 
    #mutate(TPM=log2(TPM+1)) %>% 
    filter(!celltype %in% c("Erythrocyte","undefined")) 
  if (length(mid$group %>% unique()) > 1) {
    data.T <- mid %>%
      group_by(celltype) %>%
      rstatix::t_test(AUC~group,detailed = T) %>%
      mutate(tissue=unique(mid$tissue)) %>%
      plyr::rbind.fill(data.T)
  }
  data.T2 <- mid %>% group_by(group) %>% 
    rstatix::t_test(AUC~celltype,detailed = T) %>%
    mutate(tissue=unique(mid$tissue)) %>%
    plyr::rbind.fill(data.T2)
}
write.csv(data.T,file = "GSE210347.BCAA.AUCell.Ttest.CellType.Loction.csv",row.names = F)
write.csv(data.T2,file = "GSE210347.BCAA.AUCell.Ttest.Loction.CellType.csv",row.names = F)

data <- read.csv("GSE210347.PDAC.AUCell.SYMBOLS.csv") %>%
  merge(UMAP.D,by="X")
library(RColorBrewer)
p1 <- ggplot(data %>% filter(!celltype %in% c("Erythrocyte","undefined")) %>%
               arrange(AUC),
             aes(UMAP_1,UMAP_2,color=AUC))+
  geom_point(size=1)+
  facet_wrap(~group)+
  ggthemes::theme_few()+ 
  theme_test(12)+
  #scale_size_continuous("-log10(Fisher's P)",range = c(1,4))+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9),#,angle = 90,vjust = 0.5,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    #axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  scale_color_gradientn(colors = brewer.pal(9,"Reds"))

p2 <- ggplot(data %>% filter(!celltype %in% c("Erythrocyte","undefined")),# 
             aes(UMAP_1,UMAP_2,color=celltype))+
  geom_point(size=1)+
  facet_wrap(~group)+
  ggthemes::theme_few()+ 
  theme_test(12)+
  #scale_size_continuous("-log10(Fisher's P)",range = c(1,4))+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9),#,angle = 90,vjust = 0.5,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    #axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  ggsci::scale_color_npg()+
  guides(color = guide_legend(override.aes = list(size=5)))

library(patchwork)
p <- p1/p2
ggsave(p,filename = "../../Hetero.GSE210347.BCAA.AUCell.Score.PDAC.pdf",
       height = 7.5,width = 10)

data <- read.csv("GSE210347.ICC.AUCell.SYMBOLS.csv")%>%
  merge(UMAP.D,by="X")
library(RColorBrewer)
p1 <- ggplot(data %>% filter(!celltype %in% c("Erythrocyte","undefined")) %>%
               arrange(AUC),
             aes(UMAP_1,UMAP_2,color=AUC))+
  geom_point(size=1)+
  facet_wrap(~group)+
  ggthemes::theme_few()+ 
  theme_test(12)+
  #scale_size_continuous("-log10(Fisher's P)",range = c(1,4))+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9),#,angle = 90,vjust = 0.5,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    #axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  scale_color_gradientn(colors = brewer.pal(9,"Reds"))

p2 <- ggplot(data %>% filter(!celltype %in% c("Erythrocyte","undefined")),# 
             aes(UMAP_1,UMAP_2,color=celltype))+
  geom_point(size=1)+
  facet_wrap(~group)+
  ggthemes::theme_few()+ 
  theme_test(12)+
  #scale_size_continuous("-log10(Fisher's P)",range = c(1,4))+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9),#,angle = 90,vjust = 0.5,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    #axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  ggsci::scale_color_npg()+
  guides(color = guide_legend(override.aes = list(size=5)))

library(patchwork)
p <- p1/p2
ggsave(p,filename = "../../Hetero.GSE210347.BCAA.AUCell.Score.ICC.pdf",
       height = 7.5,width = 8.5)

data <- read.csv("GSE210347.colorectal.AUCell.SYMBOLS.csv")%>%
  merge(UMAP.D,by="X")
library(RColorBrewer)
p1 <- ggplot(data %>% filter(!celltype %in% c("Erythrocyte","undefined")) %>%
               arrange(AUC),
             aes(UMAP_1,UMAP_2,color=AUC))+
  geom_point(size=1)+
  facet_wrap(~group)+
  ggthemes::theme_few()+ 
  theme_test(12)+
  #scale_size_continuous("-log10(Fisher's P)",range = c(1,4))+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9),#,angle = 90,vjust = 0.5,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    #axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  scale_color_gradientn(colors = brewer.pal(9,"Reds"))

p2 <- ggplot(data %>% filter(!celltype %in% c("Erythrocyte","undefined")),# 
             aes(UMAP_1,UMAP_2,color=celltype))+
  geom_point(size=1)+
  facet_wrap(~group)+
  ggthemes::theme_few()+ 
  theme_test(12)+
  #scale_size_continuous("-log10(Fisher's P)",range = c(1,4))+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9),#,angle = 90,vjust = 0.5,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    #axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  ggsci::scale_color_npg()+
  guides(color = guide_legend(override.aes = list(size=5)))

library(patchwork)
p <- p1/p2
ggsave(p,filename = "../../Hetero.GSE210347.BCAA.AUCell.Score.Colorectal.pdf",
       height = 7.5,width = 10)

######## BCAT1 AML ########
setwd("K:/TCGA/Anlysis/BCAA/BCAT1.AML")
library(tidyverse)
library(GSVA)
library(GSEABase)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(limma)
library(tidyverse)
library(TCGAbiolinks)
library(DESeq2)
library(fgsea)
GO.KEGG.Enrichment <- function(Symbol,filename){
  GO.Enrich <- enrichGO(
    gene = Symbol,
    keyType="SYMBOL",
    OrgDb = "org.Hs.eg.db",
    ont = "ALL",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05
  )
  if (length(GO.Enrich) != 0) {
    write.csv(GO.Enrich@result,file = paste(filename,".Enrichment.GO.",length(Symbol),".csv",sep = ""))
    print("GO Done")
  }
  gene.id <- bitr(Symbol,fromType = 'SYMBOL',toType = c('ENTREZID'),OrgDb="org.Hs.eg.db")
  
  KEGG.Enrich <- enrichKEGG(gene = na.omit(gene.id$ENTREZID),
                            organism = "hsa",
                            keyType = "kegg",
                            pvalueCutoff = 1,qvalueCutoff = 1,
                            pAdjustMethod = "BH",
                            minGSSize = 5,use_internal_data = FALSE)
  if (length(KEGG.Enrich) != 0) {
    KEGG.Enrich2<-setReadable(KEGG.Enrich, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
    write.csv(KEGG.Enrich2@result,file = paste(filename,".Enrichment.KEGG.",length(na.omit(gene.id$ENTREZID)),".csv",sep = ""))
    print("KEGG Done")
  }
  #GO.KEGG=list()
  #GO.KEGG[["GO"]] = GO.Enrich
  #GO.KEGG[["KEGG"]] = KEGG.Enrich2
  #return(GO.KEGG)
}
Hallmark <- read.gmt("../../../msigdb/h.all.v7.5.1.symbols.gmt")
KEGG <- read.gmt("../../../msigdb/c2.cp.kegg.v7.5.1.symbols.gmt")
GO <- read.gmt("../../../msigdb/c5.go.v7.5.1.symbols.gmt")
#### GSE100778 => ####
setwd("K:/TCGA/Anlysis/BCAA/BCAT1.AML")
data1 <- read.csv("GSE100778.SKM1.top.table.tsv",sep = "\t") # BCAT1 Down
data1 <- data1 %>% filter(Gene.symbol != "")
ge = data1$logFC
names(ge) = data1$Gene.symbol
ge = sort(ge,decreasing = T)

GSEA.Res <- GSEA(ge, TERM2GENE = Hallmark,eps = 1e-100,pvalueCutoff =0.05)
saveRDS(GSEA.Res,file = "GSE100778.SKM1.GSEA.HALLMARK.RDS")
GSEA.Res <- GSEA(ge, TERM2GENE = KEGG,eps = 1e-100,pvalueCutoff =0.05)
saveRDS(GSEA.Res,file = "GSE100778.SKM1.GSEA.KEGG.RDS")
GSEA.Res <- GSEA(ge, TERM2GENE = GO,eps = 1e-100,pvalueCutoff =0.05)
saveRDS(GSEA.Res,file = "GSE100778.SKM1.GSEA.GO.RDS")

data1 <- read.csv("GSE100778.MOLM13.top.table.tsv",sep = "\t") # BCAT1 Down
data1 <- data1 %>% filter(Gene.symbol != "")
ge = data1$logFC
names(ge) = data1$Gene.symbol
ge = sort(ge,decreasing = T)

GSEA.Res <- GSEA(ge, TERM2GENE = Hallmark,eps = 1e-100,pvalueCutoff =0.05)
saveRDS(GSEA.Res,file = "GSE100778.MOLM13.GSEA.HALLMARK.RDS")
GSEA.Res <- GSEA(ge, TERM2GENE = KEGG,eps = 1e-100,pvalueCutoff =0.05)
saveRDS(GSEA.Res,file = "GSE100778.MOLM13.GSEA.KEGG.RDS")
GSEA.Res <- GSEA(ge, TERM2GENE = GO,eps = 1e-100,pvalueCutoff =0.05)
saveRDS(GSEA.Res,file = "GSE100778.MOLM13.GSEA.GO.RDS")

data1 <- read.csv("GSE100778.HL60.top.table.tsv",sep = "\t") # BCAT1 Down
data1 <- data1 %>% filter(Gene.symbol != "")
ge = data1$logFC
names(ge) = data1$Gene.symbol
ge = sort(ge,decreasing = T)

GSEA.Res <- GSEA(ge, TERM2GENE = Hallmark,eps = 1e-100,pvalueCutoff =0.05)
saveRDS(GSEA.Res,file = "GSE100778.HL60.GSEA.HALLMARK.RDS")
GSEA.Res <- GSEA(ge, TERM2GENE = KEGG,eps = 1e-100,pvalueCutoff =0.05)
saveRDS(GSEA.Res,file = "GSE100778.HL60.GSEA.KEGG.RDS")
GSEA.Res <- GSEA(ge, TERM2GENE = GO,eps = 1e-100,pvalueCutoff =0.05)
saveRDS(GSEA.Res,file = "GSE100778.HL60.GSEA.GO.RDS")

HALLMARK.R1 <- readRDS("GSE100778.SKM1.GSEA.HALLMARK.RDS")
HALLMARK.R2 <- readRDS("GSE100778.MOLM13.GSEA.HALLMARK.RDS")
HALLMARK.R3 <- readRDS("GSE100778.HL60.GSEA.HALLMARK.RDS")
TERMS <- c("HALLMARK_MYC_TARGETS_V1",
           "HALLMARK_MYC_TARGETS_V2",
           "HALLMARK_P53_PATHWAY",
           "HALLMARK_G2M_CHECKPOINT",
           "HALLMARK_E2F_TARGETS",
           "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
           "HALLMARK_MTORC1_SIGNALING",
           "HALLMARK_APOPTOSIS",
           "HALLMARK_FATTY_ACID_METABOLISM")

HALLMARK.D1 <- HALLMARK.R1@result %>% filter(ID %in% TERMS) %>%
  mutate(ID=str_remove_all(ID,"HALLMARK_")) %>%
  mutate(ID=factor(ID,levels = rev(.$ID))) %>%
  mutate(CellLine="SKM1") %>%
  mutate(NES=-NES)
HALLMARK.D2 <- HALLMARK.R2@result %>% filter(ID %in% TERMS) %>%
  mutate(ID=str_remove_all(ID,"HALLMARK_")) %>%
  mutate(ID=factor(ID,levels = rev(.$ID))) %>%
  mutate(CellLine="MOLM13")%>%
  mutate(NES=-NES)
HALLMARK.D3 <- HALLMARK.R2@result %>% filter(ID %in% TERMS) %>%
  mutate(ID=str_remove_all(ID,"HALLMARK_")) %>%
  mutate(ID=factor(ID,levels = rev(.$ID))) %>%
  mutate(CellLine="HL60")%>%
  mutate(NES=-NES)
HALLMARK.D <- rbind.data.frame(HALLMARK.D1,HALLMARK.D2) %>%
  rbind.data.frame(HALLMARK.D3)

library(RColorBrewer)
p1 <- ggplot(HALLMARK.D,aes(CellLine,ID))+
  #geom_segment(aes(x=0,xend=,y=ID,yend=ID),color="black")+
  geom_point(aes(fill=NES,size=-log10(p.adjust)),shape=21,color="black")+
  ggthemes::theme_few()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="",y="",title="")+
  theme(plot.title = element_text(vjust = 0.5))+
  scale_fill_gradient2(high = brewer.pal(11,"RdBu")[1],
                       low = brewer.pal(11,"RdBu")[11],
                       mid = "white")+ #,breaks=c(1.4,1.8,2.2)
  #ggsci::scale_fill_npg()+
  #scale_x_continuous(limits = c(0,3))+
  scale_size_continuous(range = c(2,4)) #,breaks = c(2,5,8)
ggsave(p1,filename = "GSE100778.GSEA.HALLMARK.pdf",height = 3,width = 4.3)

KEGG.R1 <- readRDS("GSE100778.SKM1.GSEA.KEGG.RDS")
KEGG.R2 <- readRDS("GSE100778.MOLM13.GSEA.KEGG.RDS")
KEGG.R3 <- readRDS("GSE100778.HL60.GSEA.KEGG.RDS")
TERMS <- c("KEGG_APOPTOSIS",
           "KEGG_DNA_REPLICATION",
           "KEGG_CELL_CYCLE",
           "KEGG_P53_SIGNALING_PATHWAY",
           "KEGG_FATTY_ACID_METABOLISM",
           "KEGG_FOCAL_ADHESION",
           "KEGG_HOMOLOGOUS_RECOMBINATION",
           "KEGG_MISMATCH_REPAIR",
           "KEGG_BASE_EXCISION_REPAIR",
           "KEGG_MTOR_SIGNALING_PATHWAY")
KEGG.D1 <- KEGG.R1@result %>% filter(ID %in% TERMS) %>%
  mutate(ID=str_remove_all(ID,"KEGG_")) %>%
  mutate(ID=factor(ID,levels = rev(.$ID))) %>%
  mutate(CellLine="SKM1") %>%
  mutate(NES=-NES)
KEGG.D2 <- KEGG.R2@result %>% filter(ID %in% TERMS) %>%
  mutate(ID=str_remove_all(ID,"KEGG_")) %>%
  mutate(ID=factor(ID,levels = rev(.$ID))) %>%
  mutate(CellLine="MOLM13")%>%
  mutate(NES=-NES)
KEGG.D3 <- KEGG.R2@result %>% filter(ID %in% TERMS) %>%
  mutate(ID=str_remove_all(ID,"KEGG_")) %>%
  mutate(ID=factor(ID,levels = rev(.$ID))) %>%
  mutate(CellLine="HL60")%>%
  mutate(NES=-NES)
KEGG.D <- rbind.data.frame(KEGG.D1,KEGG.D2) %>%
  rbind.data.frame(KEGG.D3)

p1 <- ggplot(KEGG.D,aes(CellLine,ID))+
  #geom_segment(aes(x=0,xend=,y=ID,yend=ID),color="black")+
  geom_point(aes(fill=NES,size=-log10(p.adjust)),shape=21,color="black")+
  ggthemes::theme_few()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="",y="",title="")+
  theme(plot.title = element_text(vjust = 0.5))+
  scale_fill_gradient2(high = brewer.pal(11,"RdBu")[1],
                       low = brewer.pal(11,"RdBu")[11],
                       mid = "white")+ #,breaks=c(1.4,1.8,2.2)
  #ggsci::scale_fill_npg()+
  #scale_x_continuous(limits = c(0,3))+
  scale_size_continuous(range = c(2,4),breaks = c(2,4,6)) #,breaks = c(2,5,8)
ggsave(p1,filename = "GSE100778.GSEA.KEGG.pdf",height = 3,width = 4.6)


library(GseaVis)
gseaNb(object = KEGG.R2,
       geneSetID = "KEGG_CELL_CYCLE",
       addPval = T,
       pvalX = 0.6,pvalY = 0.75,
       pCol = 'black',
       pHjust = 0,subPlot = 2)

#### GSE100779-82 => Methylation => Linux ####
files <- list.files(".",pattern="tsv$")
files <- files[str_detect(files,"Methylation")]

#### GSE103960 shBCAT1 ####
data1 <- read.csv("GSE103960.HD48.shBCAT1.top.table.tsv",sep = "\t") # Ctrl vs KD
data1 <- data1 %>% filter(Gene.symbol != "")
ge = -data1$logFC # BCAT1 Down
names(ge) = data1$Gene.symbol
ge = sort(ge,decreasing = T)

GSEA.Res <- GSEA(ge, TERM2GENE = Hallmark,eps = 1e-100,pvalueCutoff =0.05)
saveRDS(GSEA.Res,file = "GSE103960.HD48.GSEA.HALLMARK.RDS")
GSEA.Res <- GSEA(ge, TERM2GENE = KEGG,eps = 1e-100,pvalueCutoff =0.05)
saveRDS(GSEA.Res,file = "GSE103960.HD48.GSEA.KEGG.RDS")
GSEA.Res <- GSEA(ge, TERM2GENE = GO,eps = 1e-100,pvalueCutoff =0.05)
saveRDS(GSEA.Res,file = "GSE103960.HD48.GSEA.GO.RDS")

data1 <- read.csv("GSE103960.AML131.shBCAT1.top.table.tsv",sep = "\t") # Ctrl vs KD
data1 <- data1 %>% filter(Gene.symbol != "")
ge = -data1$logFC # BCAT1 Down
names(ge) = data1$Gene.symbol
ge = sort(ge,decreasing = T)

GSEA.Res <- GSEA(ge, TERM2GENE = Hallmark,eps = 1e-100,pvalueCutoff =0.05)
saveRDS(GSEA.Res,file = "GSE103960.AML131.GSEA.HALLMARK.RDS")
GSEA.Res <- GSEA(ge, TERM2GENE = KEGG,eps = 1e-100,pvalueCutoff =0.05)
saveRDS(GSEA.Res,file = "GSE103960.AML131.GSEA.KEGG.RDS")
GSEA.Res <- GSEA(ge, TERM2GENE = GO,eps = 1e-100,pvalueCutoff =0.05)
saveRDS(GSEA.Res,file = "GSE103960.AML131.GSEA.GO.RDS")

HALLMARK.R1 <- readRDS("GSE103960.AML131.GSEA.HALLMARK.RDS")
HALLMARK.R2 <- readRDS("GSE103960.HD48.GSEA.HALLMARK.RDS")
TERMS <- c("HALLMARK_MYC_TARGETS_V1",
           "HALLMARK_G2M_CHECKPOINT",
           "HALLMARK_E2F_TARGETS",
           "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
           "HALLMARK_MTORC1_SIGNALING",
           "HALLMARK_APOPTOSIS",
           "HALLMARK_FATTY_ACID_METABOLISM")
library(GseaVis)
HALLMARK.D1 <- HALLMARK.R1@result %>% filter(ID %in% TERMS) %>%
  mutate(ID=str_remove_all(ID,"HALLMARK_")) %>%
  mutate(ID=factor(ID,levels = rev(.$ID))) %>%
  mutate(CellLine="AML131")
HALLMARK.D2 <- HALLMARK.R2@result %>% filter(ID %in% TERMS) %>%
  mutate(ID=str_remove_all(ID,"HALLMARK_")) %>%
  mutate(ID=factor(ID,levels = rev(.$ID))) %>%
  mutate(CellLine="HD48")
HALLMARK.D <- rbind.data.frame(HALLMARK.D1,HALLMARK.D2)
library(RColorBrewer)
p1 <- ggplot(HALLMARK.D,aes(CellLine,ID))+
  #geom_segment(aes(x=0,xend=,y=ID,yend=ID),color="black")+
  geom_point(aes(fill=NES,size=-log10(p.adjust)),shape=21,color="black")+
  ggthemes::theme_few()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="",y="",title="")+
  theme(plot.title = element_text(vjust = 0.5))+
  scale_fill_gradientn(colours = brewer.pal(9,"Reds"),breaks=c(1.4,1.8,2.2))+
  #ggsci::scale_fill_npg()+
  #scale_x_continuous(limits = c(0,3))+
  scale_size_continuous(range = c(2,4),breaks = c(2,5,8))
ggsave(p1,filename = "GSE103960.AML131.GSEA.HALLMARK.pdf",height = 3,width = 5)


KEGG.R1 <- readRDS("GSE103960.AML131.GSEA.KEGG.RDS")
KEGG.R2 <- readRDS("GSE103960.HD48.GSEA.KEGG.RDS")
TERMS <- c("KEGG_APOPTOSIS",
           "KEGG_DNA_REPLICATION",
           "KEGG_CELL_CYCLE",
           "KEGG_P53_SIGNALING_PATHWAY",
           "KEGG_FATTY_ACID_METABOLISM",
           "KEGG_FOCAL_ADHESION",
           "KEGG_HOMOLOGOUS_RECOMBINATION",
           "KEGG_MISMATCH_REPAIR",
           "KEGG_BASE_EXCISION_REPAIR",
           "KEGG_MTOR_SIGNALING_PATHWAY")
library(GseaVis)
KEGG.D1 <- KEGG.R1@result %>% filter(ID %in% TERMS) %>%
  mutate(ID=str_remove_all(ID,"KEGG_")) %>%
  mutate(ID=factor(ID,levels = rev(.$ID))) %>%
  mutate(CellLine="AML131")
KEGG.D2 <- KEGG.R2@result %>% filter(ID %in% TERMS) %>%
  mutate(ID=str_remove_all(ID,"KEGG_")) %>%
  mutate(ID=factor(ID,levels = rev(.$ID))) %>%
  mutate(CellLine="HD48")
KEGG.D <- rbind.data.frame(KEGG.D1,KEGG.D2)
library(RColorBrewer)
p1 <- ggplot(KEGG.D,aes(CellLine,ID))+
  #geom_segment(aes(x=0,xend=,y=ID,yend=ID),color="black")+
  geom_point(aes(fill=NES,size=-log10(p.adjust)),shape=21,color="black")+
  ggthemes::theme_few()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="",y="",title="")+
  theme(plot.title = element_text(vjust = 0.5))+
  scale_fill_gradientn(colours = brewer.pal(9,"Reds"),breaks=c(1.6,2,2.4))+
  #ggsci::scale_fill_npg()+
  #scale_x_continuous(limits = c(0,3))+
  scale_size_continuous(range = c(2,4),breaks = c(2,6,10))
ggsave(p1,filename = "GSE103960.AML131.GSEA.KEGG.pdf",height = 3,width = 4.3)

#### GSE127181 => gabapentin => BCAT1 inhibitor ####
GPL <- read.csv("GPL21185-21174.txt",sep = "\t",header = T)
data1 <- read.csv("GSE127181.top.table.tsv",sep = "\t") %>%
  merge(GPL,by="ID") %>% filter(GENE_SYMBOL != "")

ge = data1$logFC #
names(ge) = data1$GENE_SYMBOL
ge = sort(ge,decreasing = T)

GSEA.Res <- GSEA(ge, TERM2GENE = Hallmark,eps = 1e-100,pvalueCutoff =0.05)
saveRDS(GSEA.Res,file = "GSE127181.AML.GSEA.HALLMARK.RDS")
GSEA.Res <- GSEA(ge, TERM2GENE = KEGG,eps = 1e-100,pvalueCutoff =0.05)
saveRDS(GSEA.Res,file = "GSE127181.AML.GSEA.KEGG.RDS")
GSEA.Res <- GSEA(ge, TERM2GENE = GO,eps = 1e-100,pvalueCutoff =0.05)
saveRDS(GSEA.Res,file = "GSE127181.AML.GSEA.GO.RDS")

HALLMARK.R <- readRDS("GSE127181.AML.GSEA.HALLMARK.RDS")
TERMS <- c("HALLMARK_MYC_TARGETS_V1",
           "HALLMARK_MYC_TARGETS_V2",
           "HALLMARK_G2M_CHECKPOINT",
           "HALLMARK_E2F_TARGETS",
           "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
           "HALLMARK_MTORC1_SIGNALING",
           "HALLMARK_APOPTOSIS",
           "HALLMARK_FATTY_ACID_METABOLISM")
library(GseaVis)
HALLMARK.D <- HALLMARK.R@result %>% filter(ID %in% TERMS) %>%
  mutate(ID=str_remove_all(ID,"HALLMARK_")) %>%
  mutate(ID=factor(ID,levels = rev(.$ID)))

p1 <- ggplot(HALLMARK.D,aes(NES,ID))+
  geom_segment(aes(x=0,xend=NES,y=ID,yend=ID),color="black")+
  geom_point(aes(fill=ID,size=-log10(p.adjust)),shape=21,color="black")+
  ggthemes::theme_few()+
  theme(#legend.position = "bottom",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
        axis.text.y = element_text(size=9),
        #axis.ticks.x = element_blank(),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="NES",y="",title="Ctrl vs Gabapentin")+
  theme(plot.title = element_text(vjust = 0.5))+
  ggsci::scale_fill_npg()+
  scale_x_continuous(limits = c(0,3))+
  scale_size_continuous(range = c(2,4))
ggsave(p1,filename = "GSE127181.AML.GSEA.HALLMARK.pdf",height = 3,width = 8)


KEGG.R <- readRDS("GSE127181.AML.GSEA.KEGG.RDS")
TERMS <- c("KEGG_APOPTOSIS",
           "KEGG_DNA_REPLICATION",
           "KEGG_FATTY_ACID_METABOLISM",
           "KEGG_FOCAL_ADHESION",
           "KEGG_HOMOLOGOUS_RECOMBINATION",
           "KEGG_MISMATCH_REPAIR",
           "KEGG_BASE_EXCISION_REPAIR",
           "KEGG_MTOR_SIGNALING_PATHWAY")
KEGG.D <- KEGG.R@result %>% filter(ID %in% TERMS) %>%
  mutate(ID=str_remove_all(ID,"KEGG_")) %>%
  arrange(desc(NES)) %>%
  mutate(ID=factor(ID,levels = rev(.$ID)))

p2 <- ggplot(KEGG.D,aes(NES,ID))+
  geom_segment(aes(x=0,xend=NES,y=ID,yend=ID),color="black")+
  geom_point(aes(fill=ID,size=-log10(p.adjust)),shape=21,color="black")+
  ggthemes::theme_few()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="NES",y="",title="Ctrl vs Gabapentin")+
  theme(plot.title = element_text(vjust = 0.5))+
  ggsci::scale_fill_npg()+
  scale_x_continuous(limits = c(0,3))+
  scale_size_continuous(range = c(2,4))
ggsave(p2,filename = "GSE127181.AML.GSEA.KEGG.pdf",height = 3,width = 7)


GO.R <- readRDS("GSE127181.AML.GSEA.GO.RDS")

#### GSE148701 => Macrophage ####
library(limma)
library(data.table)
files <- list.files("GSE148701_RAW/",pattern = "gz$")
data <- data.frame()
for (file in files) {
  data <- fread(file.path("GSE148701_RAW/",file)) %>% #t() %>%
    data.frame(check.names =F ) %>% remove_rownames() %>%
    column_to_rownames("V1") %>% t() %>%
    rbind.data.frame(data)
}
data <- data %>% t() %>% data.frame(check.names = F) %>%
  magrittr::set_colnames(c("Ctrl1","KD1","Ctrl2","KD2","Ctrl3","KD3"))
data <- data %>% magrittr::set_colnames(rev(c("Ctrl1","KD1","Ctrl2","KD2","Ctrl3","KD3")))
write.csv(data,file = file.path("GSE148701_RAW/","Count.csv"))
### DEG
library(DESeq2)
countData <- read.csv(file.path("GSE148701_RAW/","Count.csv"),row.names = 1)
countData <- countData[rowMeans(countData)>1,] 
condition <- factor(c(rep(c("Ctrl","KD"),length=6)))
colData <- data.frame(row.names=colnames(countData), condition)
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)
dds1 <- DESeq(dds) 
res <- results(dds1)
res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
res1 <- res1[order(res1$pvalue, res1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
write.csv(res1,file = "GSE148701.Macrophage.BCAT1.KD.csv",row.names = F)

GeneType <- read.table("../../../../../Metastasis/TissueMicrobiome/Human.GRC38.GeneType.txt",header = T,sep = ",")
GeneType.Unique <- GeneType %>% 
  dplyr::select(Gene.stable.ID,Gene.name,NCBI.gene..formerly.Entrezgene..ID) %>% unique() %>%
  magrittr::set_colnames(c("ENSEMBL","SYMBOL","EntreZ"))

countData2 <- countData %>% mutate(ENSEMBL=str_remove_all(rownames(.),"\\..*")) 

Res <- res1 %>% mutate(ENSEMBL=str_remove_all(rownames(.),"\\..*")) %>%
  merge(GeneType.Unique,by="ENSEMBL") %>%
  merge(countData2,by="ENSEMBL") %>%
  arrange(pvalue) %>%
  mutate(padj=p.adjust(pvalue,method = "BH",n = nrow(.))) %>%
  filter(SYMBOL != "")
write.csv(Res,file = "GSE148701.Macrophage.BCAT1.KD.csv",row.names = F)


ge = Res$log2FoldChange #
names(ge) = Res$SYMBOL
ge = sort(ge,decreasing = T)

GSEA.Res <- GSEA(ge, TERM2GENE = Hallmark,eps = 1e-100,pvalueCutoff =0.05)
saveRDS(GSEA.Res,file = "GSE148701.Macrophage.GSEA.HALLMARK.RDS")
GSEA.Res <- GSEA(ge, TERM2GENE = KEGG,eps = 1e-100,pvalueCutoff =0.05)
saveRDS(GSEA.Res,file = "GSE148701.Macrophage.GSEA.KEGG.RDS")
GSEA.Res <- GSEA(ge, TERM2GENE = GO,eps = 1e-100,pvalueCutoff =0.05)
saveRDS(GSEA.Res,file = "GSE148701.Macrophage.GSEA.GO.RDS")

HALLMARK.R <- readRDS("GSE148701.Macrophage.GSEA.HALLMARK.RDS")
library(GseaVis)
gseaNb(object = HALLMARK.R,
       geneSetID = "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
       addPval = T,
       pvalX = 0.6,pvalY = 0.75,
       pCol = 'black',
       pHjust = 0,subPlot = 2)

KEGG.R <- readRDS("GSE148701.Macrophage.GSEA.KEGG.RDS")
GO.R <- readRDS("GSE148701.Macrophage.GSEA.GO.RDS")

Res <- read.csv("GSE148701.Macrophage.BCAT1.KD.csv")
Up.BCAT1.KD <- Res %>% filter(log2FoldChange <= -1) %>% filter(padj <= 0.05)
GO.KEGG.Enrichment(Up.BCAT1.KD$SYMBOL,"GSE148701.Macrophage.Up.BCAT1.KD")

Up.BCAT1.Ctrl <- Res %>% filter(log2FoldChange >= 1) %>% filter(padj <= 0.05)
GO.KEGG.Enrichment(Up.BCAT1.Ctrl$SYMBOL,"GSE148701.Macrophage.Up.BCAT1.Ctrl")


#### Methylation => EPIC => ChAMP ####
library("ChAMP")
myLoad <- champ.load("GSE100779",arraytype="EPIC")
save(myLoad,file = "GSE100779.ChAMP.load.Rda")
#CpG.GUI(arraytype="EPIC")
champ.QC() # Alternatively QC.GUI(arraytype="EPIC")
myNorm <- champ.norm(arraytype="EPIC")
save(myNorm,file = "GSE100779.ChAMP.Norm.Rda")
#champ.SVD()
# If Batch detected, run champ.runCombat() here.This data is not suitable.
#myDMP <- champ.DMP(beta = myNorm,pheno=myLoad$pd$Sample_Group,arraytype="EPIC")
#DMP.GUI(DMP=myDMP[[1]],beta=myNorm,pheno=myLoad$pd$Sample_Group,arraytype="EPIC")
#save(myDMP,file="GSE100779.DMP.rda")

myDMR <- champ.DMR(beta=myNorm, pheno=myLoad$pd$Sample_Group, compare.group=c("Ctrl", "BCAT1OE"),
                   arraytype="EPIC",method="Bumphunter")
save(myDMR,file="GSE100779.DMR.rda")
DMR.GUI(arraytype="EPIC")

myGSEA <- champ.GSEA(arraytype="EPIC")
save(myGSEA,file="GSE100779.GSEA.rda")

#### GSE100782 => Methylation ####
setwd("K:/TCGA/Anlysis/BCAA/BCAT1.AML")
library(tidyverse)
library(data.table)
Methy.Data <- read.csv("GSE100782.MOLM13.Methylation.OE.20WEEK.top.table.tsv",sep = "\t") %>%
  filter(adj.P.Val<=0.05) #%>% filter(abs(logFC)>=0.5)
Anno <- read.csv("GPL21145_MethylationEPIC_15073387_v-1-0.csv",skip = 7,header = T) %>%
  dplyr::select(IlmnID,Name,UCSC_RefGene_Name,UCSC_RefGene_Accession,UCSC_RefGene_Group,UCSC_CpG_Islands_Name,Relation_to_UCSC_CpG_Island,Phantom4_Enhancers,Phantom5_Enhancers,DMR,HMM_Island,Regulatory_Feature_Name,Regulatory_Feature_Group) %>% 
  filter(UCSC_RefGene_Name != "")
Des <- Anno %>% filter(Name %in% Methy.Data$ID) %>%
  mutate(Gene.symbol=str_remove_all(UCSC_RefGene_Name,"\\;.*")) %>%
  filter(Gene.symbol %in% Gene.SYMBOL) %>%
  merge(Methy.Data,.,by.x="ID",by.y="Name") %>%
  filter(!str_detect(UCSC_RefGene_Group,"Body")) %>%
  arrange(Gene.symbol) %>%
  mutate(ID=factor(ID,levels = .$ID))

p1 <- ggplot(Des,aes(ID,logFC))+
  geom_segment(aes(x=ID,xend=ID,y=0,yend=logFC),color="black")+
  geom_point(aes(fill=Gene.symbol,size=-log10(adj.P.Val)),shape=21,color="black")+
  ggthemes::theme_few()+
  theme(legend.position = "bottom",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9,angle = 45,vjust = 0.5,hjust = 1),#element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
        axis.text.y = element_text(size=9),
        #axis.ticks.x = element_blank(),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(y="logFC",x="",title="Ctrl vs BCAT1 OE")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggsci::scale_fill_npg()+
  scale_y_continuous(limits = c(-0.2,0))+
  geom_hline(yintercept = c(-0.05,-0.1,-0.15),linetype=2)+
  scale_size_continuous(range = c(2,4))
p2 <- ggplot(Des %>% dplyr::select(ID,Gene.symbol) %>%
               unique(),
             aes(ID,1,fill=Gene.symbol))+
  geom_col(width = 1)+
  theme_void()+
  ggsci::scale_fill_npg()
p12 <- p1 %>% aplot::insert_bottom(p2,height = 0.2) 
ggsave(p12,filename = "SYMBOL.EXAMPLE.Methylation.pdf",height = 2,width = 9)
ggsave(p12,filename = "SYMBOL.EXAMPLE.Methylation.Legend.pdf",height = 8,width = 9)
write.csv(Des,file = "SYMBOL.EXAMPLE.Methylation.Figure.Data.csv",row.names = F)



data1 <- read.csv("GSE100778.SKM1.top.table.tsv",sep = "\t") %>%
  filter(adj.P.Val<=0.05) %>%
  filter(Gene.symbol %in%Gene.SYMBOL) %>%
  mutate(Project="GSE100778.SKM1")
data2 <- read.csv("GSE100778.MOLM13.top.table.tsv",sep = "\t")%>%
  filter(adj.P.Val<=0.05)%>%
  filter(Gene.symbol %in%Gene.SYMBOL) %>%
  mutate(Project="GSE100778.MOLM13")
data3 <- read.csv("GSE100778.HL60.top.table.tsv",sep = "\t")%>%
  filter(adj.P.Val<=0.05)%>%
  filter(Gene.symbol %in%Gene.SYMBOL) %>%
  mutate(Project="GSE100778.HL60")
  
data4 <- read.csv("GSE103960.HD48.shBCAT1.top.table.tsv",sep = "\t")%>%
  filter(adj.P.Val<=0.05) %>%
  filter(Gene.symbol %in%Gene.SYMBOL) %>%
  mutate(logFC=-logFC)%>%
  mutate(Project="GSE103960.HD48")# Ctrl vs KD
data <- rbind.data.frame(data1,data2) %>%
  rbind.data.frame(data3) %>%
  rbind.data.frame(data4)

Figure.Data <- data %>% filter(Gene.symbol != "BCAT1") %>%
  mutate(CellLine=str_remove_all(Project,".*\\.")) %>%
  mutate(Cancer=str_remove_all(Project,"\\..*")) %>%
  arrange(Gene.symbol) %>%
  mutate(Xaxis=paste(Project,ID,sep = "_")) %>%
  mutate(Xaxis=factor(Xaxis,levels = .$Xaxis))

p1 <- ggplot(Figure.Data,aes(Xaxis,logFC))+
  geom_segment(aes(x=Xaxis,xend=Xaxis,y=0,yend=logFC),color="black")+
  geom_point(aes(fill=Gene.symbol,size=-log10(adj.P.Val)),shape=21,color="black")+
  ggthemes::theme_few()+
  theme(legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_blank(),#element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    axis.ticks.x = element_blank(),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(y="logFC",x="",title="Ctrl vs shBCAT1")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggsci::scale_fill_npg()+
  scale_y_continuous(limits = c(-1.7,0))+
  geom_hline(yintercept = c(-0.5,-1,-1.5),linetype=2)+
  scale_size_continuous(range = c(2,4))

p2 <- ggplot(Figure.Data %>% dplyr::select(Xaxis,Gene.symbol) %>%
               unique(),
             aes(Xaxis,1,fill=Gene.symbol))+
  geom_col(width = 1)+
  theme_void()+
  ggsci::scale_fill_npg()
library(RColorBrewer)
p3 <- ggplot(Figure.Data %>% dplyr::select(Xaxis,Cancer) %>%
               unique(),
             aes(Xaxis,1,fill=Cancer))+
  geom_col(width = 1)+
  theme_void()+
  scale_fill_manual("",values = brewer.pal(8,"Paired"))

p123 <- p1 %>% aplot::insert_bottom(p2,height = 0.2) %>%
  aplot::insert_bottom(p3,height = 0.1)
ggsave(p123,filename = "SYMBOL.EXAMPLE.pdf",height = 2,width = 10)
ggsave(p123,filename = "SYMBOL.EXAMPLE.Legend.pdf",height = 8,width = 10)
write.csv(Figure.Data,file = "SYMBOL.EXAMPLE.Expression.Figure.Data.csv",row.names = F)



#GPL <- read.csv("GPL21185-21174.txt",sep = "\t",header = T)
#data6 <- read.csv("GSE127181.top.table.tsv",sep = "\t") %>%
#  merge(GPL,by="ID") %>% filter(GENE_SYMBOL != "") %>%
#  filter(adj.P.Val<=0.05)

#data5 <- read.csv("GSE103960.AML131.shBCAT1.top.table.tsv",sep = "\t")%>%
#  filter(adj.P.Val<=0.05)


Gene.SYMBOL <- intersect(data1$Gene.symbol,data2$Gene.symbol) %>%
  intersect(data3$Gene.symbol) %>%
  intersect(data4$Gene.symbol) %>%
  #intersect(data6$GENE_NAME) %>%
  intersect(Des$Gene.symbol) %>% sort()





######## Cancer cell States Signature ####
files <- list.files("../../CancerSEA",pattern = ".txt$")
#data <- read.csv("../../mRNA.PanCancer.Exp/TCGA-LIHC.mRNA.TPM.csv",check.names = F,row.names = 1)
library(data.table)
Score.Gene <- list()
for (file in files) {
  SigName = file %>% str_remove_all("\\..*")
  data2 <- fread(file.path("../../CancerSEA",file))
  Score.Gene[[SigName]] = data2$EnsembleID
}

files <- list.files(path = "../../mRNA.PanCancer.Exp/",pattern = "csv$")
Pancancer.CancerState <- data.frame()
for (file in files) {
  project <- file %>% str_remove_all("\\..*") %>% str_remove_all("TCGA-")
  print(project)
  phenotype <- read.csv(paste("../../Clinical.XENA/",project,"/","TCGA-",project,".GDC_phenotype.tsv",sep = ""),sep = "\t",header = T) %>% 
    dplyr::filter(sample_type.samples == "Primary Tumor" | 
                    sample_type.samples == "Primary Blood Derived Cancer - Peripheral Blood")
  TPM <- read.csv(paste("../../mRNA.PanCancer.Exp/TCGA-",project,".mRNA.TPM.csv",sep = ""),check.names = F,row.names = 1)
  TPM.log <- log2(TPM+1)
  
  GSEA.Res <- gsva(expr=as.matrix(TPM.log),
               gset.idx.list=Score.Gene, 
               method="ssgsea",
               kcdf="Gaussian" ,
               verbose=T)
  Pancancer.CancerState <- GSEA.Res %>% t() %>%
    data.frame(check.names = F) %>%
    mutate(Cancer=project) %>%
    rbind.data.frame(Pancancer.CancerState)
}
write.csv(Pancancer.CancerState,"../../CancerSEA/Pancancer.CancerState.csv")
#### RS-high vs RS-low => Stage ####
RiskScore.Cancer <- read.csv(file = "RiskGroup.COX.TCGA.RiskScore.csv") %>%
  mutate(ID=substring(SampleID,1,15))

PhenoType <- read.csv("../../mRNA.Phenotype.csv")
PhenoType <- PhenoType %>% filter(Sample.Type == "Primary Tumor" | 
                                    Sample.Type == "Solid Tissue Normal" | 
                                    Sample.Type == "Primary Blood Derived Cancer - Peripheral Blood") %>% 
  dplyr::select(Project.ID,Sample.ID,Sample.Type) %>%
  mutate(Type=ifelse(Sample.Type=="Solid Tissue Normal","Normal","Tumor")) %>%
  mutate(Cancer=str_remove(Project.ID,"TCGA-")) %>%
  filter(Type=="Tumor")

Risk.Cutoff <- read.csv(file = "RiskGroup.COX.TCGA.RiskScore.CutOff.csv")
Stage.Clinical <- read_xlsx("../../PanCancer.DiverseIndex.xlsx",sheet = 19)
data2<-Stage.Clinical[,-1]
# 18 Grade 19 Stage
Index.data2 <- merge(RiskScore.Cancer,data2,by="ID",by.y="SampleName") %>% 
  unique() %>%
  data.frame() %>%
  filter(SampleID %in% PhenoType$Sample.ID) %>%
  mutate(Stage=str_remove_all(Stage,"Stage ")) %>%
  #gather(Attribute,Score,StromalScore:ESTIMATEScore) %>%
  #dplyr::select(ID,SampleID,Attribute,Score,Cancer,total_risk_score) %>%
  group_by(Cancer) %>%
  mutate(cutpoint=median(total_risk_score)) %>%
  mutate(Group = if_else(total_risk_score >= cutpoint,"High","Low"))
write.csv(Index.data2,file = "RiskScore.Stage.csv",row.names = F)

Fisher.Stage <- data.frame()
for(cancer in unique(Index.data2$Cancer)){
  mid.data <- Index.data2 %>% filter(Cancer==cancer)
  table <- table(mid.data$Group,mid.data$Stage)
  test <- fisher.test(table)
  Fisher.Stage <- data.frame(Fisher.P = test$p.value,Cancer=cancer) %>%
    rbind.data.frame(Fisher.Stage)
}
write.csv(Fisher.Stage,file = "RiskScore.Stage.Fisher.P.csv")

col <- brewer.pal(9,"Greens")[c(2,4,5,6)] #c("#A6CEE3","#FDBF6F","#CAB2D6","#E31A1C")
names(col) <- sort(unique(Index.data2$Stage))
Figure.plot <- list()
for (cancer in unique(Index.data2$Cancer)) {
  Index.dat3 <- Index.data2 %>% filter(Cancer==cancer) %>%
    group_by(Cancer,Group,Stage) %>%
    summarise(Count=n()) %>%
    spread(Group,Count,fill=0)
  
  dat1 = aggregate(Index.dat3$High, by = list(Index.dat3$Stage), FUN = sum)
  dat1$per1 = dat1$x / sum(dat1$x)
  
  # for循环构建标签的相对位置
  for (i in seq(nrow(dat1), 1)) {
    if (i == nrow(dat1)) {
      dat1$per.y1[i] = dat1$per1[i] / 2
    }else{
      dat1$per.y1[i] = sum(dat1$per1[(i + 1):nrow(dat1)]) + dat1$per1[i] / 2
    }
  }
  dat1$label1 = paste(dat1$Group.1,'(',round(dat1$per1*100, 2),'%',')', sep = '')
  Index.dat3 = merge(Index.dat3, dat1[,c(1,3,4,5)], by.x = 'Stage', by.y = 'Group.1')
  
  dat2 = aggregate(Index.dat3$Low, by = list(Index.dat3$Stage), FUN = sum)
  dat2$per2 = dat2$x / sum(dat2$x)
  
  for (i in seq(nrow(dat2), 1)) {
    if (i == nrow(dat2)) {
      dat2$per.y2[i] = dat2$per2[i] / 2
    }else{
      dat2$per.y2[i] = sum(dat2$per2[(i + 1):nrow(dat2)]) + dat2$per2[i] / 2
    }
  }
  
  dat2$label2 = paste(dat2$Group.1,'(',round(dat2$per2*100, 2),'%',')', sep = '')
  Index.dat3 = merge(Index.dat3, dat2[,c(1,3,4,5)], by.x = 'Stage', by.y = 'Group.1')
  
  pvalue <- Fisher.Stage %>% filter(Cancer==cancer)
  p <- ggplot(Index.dat3) +
    # 绘制柱状图
    geom_bar(aes("a", 
                 per2, 
                 fill = Stage),
             stat = 'identity', width = 1.3) +
    # 添加标签
    #geom_text(aes(1.25, as.numeric(per.y2), 
    #              label = label2),
    #          size =2.5, color = 'black') +
    # 绘制柱状图
    geom_bar(aes("b", per1, fill = Stage), 
             stat = 'identity', width = .8, color = 'white') +
    # 添加标签
    #geom_text(aes(2, as.numeric(per.y1),label = label1),
    #          size = 2.5, color = 'black') +
    # 设置Y轴刻度
    scale_y_continuous(labels = scales::percent) +
    coord_polar(theta = "y") + # 转换坐标轴
    theme_void() +
    scale_fill_manual(values = col)+
    #ggsci::scale_fill_npg() + # 设置填充色
    theme(legend.position = 'none')+ # 隐藏图例
    ggtitle(paste(cancer,"\n","Fisher'P:",round(pvalue$Fisher.P,4),sep = ""))+
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(p,filename = paste("RS.Stage.",cancer,".PiePlot.pdf",sep = ""),
         height = 2.5,width = 2.5)
  Figure.plot[[cancer]] <- p
}

p321<-cowplot::plot_grid(plotlist = Figure.plot,
                         ncol = 4,align = 'hv')
ggsave(p321,filename = "RS.Stage.PanCancer.Fisher.pdf",height = 10,width = 8)

#### RS-high vs RS-low => Cancer State ####
Pancancer.CancerState <- read.csv("../../CancerSEA/Pancancer.CancerState.csv")
PhenoType <- read.csv("../../mRNA.Phenotype.csv")
PhenoType <- PhenoType %>% filter(Sample.Type == "Primary Tumor" | 
                                    Sample.Type == "Solid Tissue Normal" | 
                                    Sample.Type == "Primary Blood Derived Cancer - Peripheral Blood") %>% 
  dplyr::select(Project.ID,Sample.ID,Sample.Type) %>%
  mutate(Type=ifelse(Sample.Type=="Solid Tissue Normal","Normal","Tumor")) %>%
  mutate(Cancer=str_remove(Project.ID,"TCGA-")) %>%
  filter(Type=="Tumor")

setwd("K:/TCGA/Anlysis/BCAA")
RiskScore.Cancer <- read.csv(file = "RiskGroup.COX.TCGA.RiskScore.csv") %>%
  #mutate(ID=substring(SampleID,1,15))%>% 
  unique() %>%
  data.frame() %>%
  filter(SampleID %in% PhenoType$Sample.ID) %>%
  group_by(Cancer) %>%
  mutate(cutpoint=median(total_risk_score)) %>%
  mutate(Group = if_else(total_risk_score >= cutpoint,"High","Low")) %>%
  merge(Pancancer.CancerState %>% dplyr::select(-Cancer),by.x="SampleID",by.y="X") %>%
  dplyr::select(-cutpoint,-total_risk_score,-OS,-OS.time) %>%
  arrange(Cancer,Group)

Test.State <- RiskScore.Cancer %>%
  gather(State,Score,Angiogenesis:Quiescence) %>%
  group_by(Cancer,State) %>%
  rstatix::t_test(Score~Group,detailed = T)
write.csv(Test.State,"RiskScore.CancerState.Ttest.csv")



rt <- data.frame()
for(cancer in unique(RiskScore.Cancer$Cancer)){
  rt <- RiskScore.Cancer %>%
    filter(Cancer==cancer) %>%
    dplyr::select(-Cancer,-Group) %>%
    remove_rownames() %>%
    column_to_rownames("SampleID") %>%
    scale() %>%
    rbind(rt)
}

library(ComplexHeatmap)
library(circlize)
library(ggsci)
color = c(pal_d3("category20")(20),pal_d3("category20b")(20),pal_d3("category20c")(20),pal_d3("category10")(10))
col1 <- color[1:23]
names(col1) <- sort(unique(RiskScore.Cancer$Cancer))

col2 <- c("Black","Red")
names(col2) <- c("High","Low")
  
ha = HeatmapAnnotation(Cancer= RiskScore.Cancer$Cancer,
                       RS = RiskScore.Cancer$Group,
                       show_legend = rep(TRUE, 2),#是否要显示annotation legend
                       annotation_height = unit(rep(5,2), "mm"),#临床annotation的高度
                       col = list(Cancer=col1,
                                  RS=col2),
                       annotation_legend_param = list(
                         Cancer = list(title = "Cancers"),
                         RS = list(title = "RS")),
                       simple_anno_size = unit(0.3, "cm"),
                       na_col = "white",
                       annotation_name_side = "left"
)

ht_opt(legend_border='black',
       heatmap_border = TRUE,
       annotation_border = TRUE) 

col_fun <- colorRamp2(c(-2, 0,2), c("#377EB8", "white", "#E41A1C"))
subtype <- RiskScore.Cancer$Cancer
ht <- Heatmap(t(rt), col = col_fun, 
              name = "Z Score",
              cluster_rows = F, 
              cluster_columns = T,         
              show_row_names = TRUE, 
              row_names_side = "left", #行名位置
              show_column_names = FALSE,
              column_split = subtype, #根据数据修改
              #row_split = c(rep(1,3),rep(2,6),rep(3,15),rep(4,1),rep(5,6),rep(6,3),rep(7,6)),#根据数据修改
              #row_gap = unit(c(1,1,1,1,1,1) ,"mm"),
              column_gap = unit(c(2), "mm"),
              row_title = NULL,column_title = NULL,
              #right_annotation = harow,
              use_raster=F,
              top_annotation = ha)
pdf("RiskScore.CancerState.pdf",width = 16,height = 6)
draw(ht,padding = unit(c(0.1, 4, 0.1, 0.5), "cm"),
     annotation_legend_side = "right", 
     heatmap_legend_side = "right")
dev.off()

####################
####################
####################

library(TCGAbiolinks)
projects <- getGDCprojects()
projects <- projects$project_id

TCGA_dowload<-function(x,dirpath){
  query.rppa <- GDCquery(
    project = x, 
    data.category = "Proteome Profiling",
    data.type = "Protein Expression Quantification"
  )
  GDCdownload(query.rppa) 
  Proteins <- GDCprepare(query.rppa)
  saveRDS(Proteins,file = paste0(dirpath,x,"_protein.rds"))
}

for (j in projects) {
  print(j)
  try(TCGA_dowload(j,dirpath = "./CPTAC/"),silent = T)
}
