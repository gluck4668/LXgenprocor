
# The nine quad of genes and proteins

#---------------------------------------------------------------------

LXgenprocor <- function(gene_data,protein_data,FC,species){

  #list all the packages that have been installed
  all_packages <- data.frame(installed.packages())

  #To judge whether a package was installed. If not, it will be installed.
  pack <- data.frame(c("devtools","BiocManager","openxlsx","dplyr","psych",
                       "ggplot2","pheatmap", "ggrepel","pak","ggpubr") )

  # To judge whether a package was included in the all_packages: %in%
  pack$type <- pack[,1] %in% all_packages$Package

  for (i in 1:nrow(pack)){
    if(pack[i,2]==FALSE)
      install.packages(pack[i,1],update = F,ask = F)
  }
  rm(i)

  # 批量library
  packages <- as.character(pack[,1])

  for(i in packages){
    library(i, character.only = T)
  }
  rm(i)

  #-----------------
  if("tidyverse" %in% all_packages$Package==FALSE)
    pak::pak("tidyverse/tidyverse")
  library(tidyverse)

  #-----------------
  BiocManager_pack <- data.frame(c("clusterProfiler",
                                   "org.Hs.eg.db","org.Mm.eg.db","org.Rn.eg.db"))
  # human: "org.Hs.eg.db"
  # mouse: "org.Mm.eg.db"
  # rat: "org.Rn.eg.db"

  BiocManager_pack$type <- BiocManager_pack[,1] %in% all_packages$Package

  for (i in 1:nrow(BiocManager_pack)){
    if(BiocManager_pack[i,2]==FALSE)
      BiocManager::install(BiocManager_pack[i,1],update = F,ask = F)
  }

  # 批量library
  Bio_packages <- as.character(BiocManager_pack[,1])
  for(i in Bio_packages){
    library(i, character.only = T)
  }
  rm(i)

  #------处理gene data-------------------------------------------------------
  #读入数据；
  RNA <- read.xlsx(gene_data)
  colnames(RNA) <- c("gene_symbol","log2FC","pvalue")
  RNA$log2FC <- as.numeric(RNA$log2FC)
  RNA$pvalue <- as.numeric(RNA$pvalue)

  if(is.null(FC))
    RNA <- dplyr::filter(RNA,pvalue<0.05)  else
    RNA <- dplyr::filter(RNA,pvalue<0.05 & abs(log2FC)>=log2(FC))
  head(RNA)

  table(duplicated(RNA$gene_symbol))  # 查看gene_symbol是否有重复
  RNA <- distinct(RNA,gene_symbol,.keep_all = T) # 去重

  # hsa: all letters of the gene symbol must be UPPER（人类：全部大写）
  RNA_human <- RNA
  RNA_human$gene_symbol <-toupper(RNA$gene_symbol)

  # mouse/rat: the first letter is UPPER, and the rest of letters must be LOWER（鼠：首大写，其余小写）
  # gene_n <- dplyr::filter(RNA,tolower(str_sub(RNA$gene_symbol,1,1)) %in% 0:9) # 数字开头的基因
  # gene_L <- dplyr::filter(RNA,tolower(str_sub(RNA$gene_symbol,1,1)) %in% letters) # 字母开头的基因
  #gene_L$gene_symbol <- str_to_title(tolower(gene_L$gene_symbol)) # str_to_title()首字母大写
  #RNA_animal <-rbind(gene_n,gene_L)

  RNA_animal <- RNA

  #--------------------------------------------------------------------------
  spe <- trimws(species) %>% tolower()

  if (spe=="human"){
    keytypes(org.Hs.eg.db) # to check the key items
    geneID<-bitr(RNA_human$gene_symbol,
                 fromType = "SYMBOL",toType =c("SYMBOL","ENTREZID","ENSEMBL"),
                 OrgDb = org.Hs.eg.db)
  }


  if (spe=="rat"){
    keytypes(org.Rn.eg.db) # to check the key items
    geneID<-bitr(RNA_animal$gene_symbol,
                 fromType = "SYMBOL",toType =c("SYMBOL","ENTREZID","ENSEMBL"),
                 OrgDb = org.Rn.eg.db)
  }


  if (spe=="mouse"){
    keytypes(org.Mm.eg.db) # to check the key items
    geneID<-bitr(RNA_animal$gene_symbol,
                 fromType = "SYMBOL",toType =c("SYMBOL","ENTREZID","ENSEMBL"),
                 OrgDb = org.Mm.eg.db)
  }

  RNA$ID <- toupper(RNA$gene_symbol)
  geneID$ID <- toupper(geneID$SYMBOL)

  table(duplicated(geneID$SYMBOL))  # 查看symbol是否有重复
  geneID <- distinct(geneID,SYMBOL,.keep_all = T) # 去重

  RNA_data <- dplyr::inner_join(RNA,geneID,"ID")
  RNA_data <- RNA_data[,-c(4:6)]
  RNA_data <- RNA_data[,c(4,1:3)]
  colnames(RNA_data) <- c("ENSEMBL","gene_symbol","log2FC_gene","pvalue_gene")

  #-----------处理 protein data----------------------------------------

  protein <- read.xlsx(protein_data)
  protein[,2] <- as.numeric(protein[,2])
  protein[,3] <- as.numeric(protein[,3])
  protein[,2] <- log2(protein[,2])
  colnames(protein) <- c("UNIPROT","log2FC","pvalue")

  if(is.null(FC))
    protein <- dplyr::filter(protein,pvalue<0.05) else
    protein <- dplyr::filter(protein,pvalue<0.05 & abs(log2FC)>=log2(FC))

  head(protein)

  table(duplicated(protein$UNIPROT))  # 查看UNIPROT是否有重复
  protein <- distinct(protein,UNIPROT,.keep_all = T) # 去重

  if (spe=="human"){
    keytypes(org.Hs.eg.db) # to check the key items
    proteinID<-bitr(protein$UNIPROT,
                    fromType = "UNIPROT",toType =c("UNIPROT","SYMBOL","ENTREZID","ENSEMBL"),
                    OrgDb = org.Hs.eg.db)
  }

  if (spe=="rat"){
    keytypes(org.Rn.eg.db) # to check the key items
    proteinID<-bitr(protein$UNIPROT,
                    fromType = "UNIPROT",toType =c("UNIPROT","SYMBOL","ENTREZID","ENSEMBL"),
                    OrgDb = org.Rn.eg.db)
  }


  if (spe=="mouse"){
    keytypes(org.Mm.eg.db) # to check the key items
    proteinID <-bitr(protein$UNIPROT,
                     fromType = "UNIPROT",toType =c("UNIPROT","SYMBOL","ENTREZID","ENSEMBL"),
                     OrgDb = org.Mm.eg.db)
  }


  table(duplicated(proteinID$UNIPROT))  # 查看UNIPROT是否有重复
  proteinID <- distinct(proteinID,UNIPROT,.keep_all = T) # 去重

  protein_data <- dplyr::inner_join(protein,proteinID,"UNIPROT")
  protein_data <- protein_data[,c(6,1:3)]
  colnames(protein_data) <- c("ENSEMBL","UNIPROT","log2FC_protein","pvalue_protein")

  #合并gene_data和protein_data两个表格
  table(duplicated(RNA_data$ENSEMBL))
  RNA_data <- distinct(.data = RNA_data,ENSEMBL,.keep_all = T) #去重

  table(duplicated(protein_data$ENSEMBL))
  protein_data <- distinct(.data = protein_data,ENSEMBL,.keep_all = T) #去重

  data <- dplyr::inner_join(RNA_data,protein_data,"ENSEMBL")

  data <- na.omit(data)

  colnames(data) <- c("ENSEMBL","gene_symbol","log2FC_RNA","pvalue_gene","UNIPROT","log2FC_Protein","pvalue_protein")

  data$log2FC_RNA <- as.numeric(data$log2FC_RNA)
  data$log2FC_Protein <- as.numeric(data$log2FC_Protein)


  #对数据进行分组；
  #生成显著上下调数据标签；

  data_cor <- data

  data_cor$correlation <- case_when(data_cor$log2FC_RNA > 0 & data_cor$log2FC_Protein > 0 ~ "positive",
                                   data_cor$log2FC_RNA < 0 & data_cor$log2FC_Protein < 0 ~ "positive",
                                   data_cor$log2FC_RNA > 0 & data_cor$log2FC_Protein < 0 ~ "negative",
                                   data_cor$log2FC_RNA < 0 & data_cor$log2FC_Protein > 0 ~ "negative")

  postive_cor <- dplyr::filter(data_cor,correlation=="positive")
  postive_cor <- postive_cor[,-ncol(data_cor)]

  negative_cor <- dplyr::filter(data_cor,correlation=="negative")
  negative_cor <-negative_cor[,-ncol(data_cor)]


  if(dir.exists("analysis result")==FALSE)
    dir.create("analysis result")

  write.xlsx(data,"analysis result/data_cor.xlsx")
  write.xlsx(postive_cor,"analysis result/postive_cor.xlsx")
  write.xlsx(negative_cor,"analysis result/negative_cor.xlsx")

  #开始尝试绘图；

  p0 <-ggplot(data,aes(log2FC_RNA,log2FC_Protein,color="red"))

  p1 <- p0+geom_point(size = 2)+
           guides(color="none") # 不显示图例legend

  p1

  mytheme<-theme_bw()+
    theme(text=element_text(family = "sans",face="bold",colour ="black",size =14),
          #panel.grid =  element_line(linewidth = 0.3,colour = "gray60"),
          panel.grid.major = element_line(linewidth = 0.3,colour = "gray65"),
          panel.grid.minor = element_line(linewidth = 0.3,colour = "gray65"),
          panel.border = element_rect(linewidth = 0.8,colour = "black"),
          axis.line = element_blank(),
          axis.ticks = element_line(linewidth = 0.6,colour = "black"),
          axis.ticks.length = unit(1.5,units = "mm"))+
    theme(plot.title = element_text(hjust = 0.5))+
    theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))

  titile_text <- paste('Correlation Analysis')

  cor  <- WGCNA::cor(data$log2FC_RNA,data$log2FC_Protein,use = "complete.obs",method ="spearman")
  cor <- round(cor,4)

  lab = paste0("Coefficient=",cor)
  lab

  xmin <- min(data$log2FC_RNA)
  ymax <- max(data$log2FC_Protein)

  cor_value <- geom_text(x=xmin+0.9,y=ymax-0.2,label = lab, size=5,color="black")

  p2 <- p1+labs(title =titile_text)+mytheme+cor_value

  ggsave("analysis result/Correlation Analysis.png",
         p2,width=1000, height =800, dpi=150,units = "px")

  p2
  
  
  cor_line <- ggscatter(data = data,x = "log2FC_RNA",y = "log2FC_Protein", # ggpubr包
                        title="Correlation Analysis",
                        xlab="log2FC_genes",
                        ylab="log2FC_proteins",
                        add = "reg.line",
                        cor.method= "pearson",   # "pearson", "kendall", or "spearman".
                        conf.int = T,
                        cor.coef = T,
                        cor.coef.size=6,
                        color = "grey30", size = 2)
  cor_line
  
  cor_theme <- theme(plot.title = element_text(size = 20,hjust = 0.5))+
               theme(axis.title.x = element_text(size = 18),
                     axis.title.y = element_text(size = 18),
                     axis.text =  element_text(size = 14))+
               theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
  
  
  
  cor_line <- cor_line+cor_theme
  
  ggsave("analysis result/cor_line Analysis.png",
         cor_line,width=1000, height =800, dpi=150,units = "px")
  
  print("Please see the results in the folder of <analysis result>")
  
  cor_line

  
  }














