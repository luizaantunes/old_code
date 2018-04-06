# old_code

<p>Stat proteomica using file choose </p>
Estatistica, p-values and fdr-values </p>
install ROTS from http://www.btk.fi/research/research-groups/elo/software/rots/ </p>
library(ROTS)
library(limma)

#load data
filename <- file.choose()
Log2_data <- read.csv(filename, header=TRUE, row.names=1, stringsAsFactors=FALSE)
#DEP_vs <- cbind(subset(Input[ , grep(amostra1, names(Input))]), subset(Input[ , grep(amostra2, names(Input))]))

#function
stat.trial <- function(Log2_data, Raw_data) 
  {
  #teste t
  T.test.p <- apply(Log2_data, 1, function(x){t.test(x[1:2], x[3:4],var.equal = TRUE)$p.value})
  T.test.fdr <- p.adjust(T.test.p, method="BH")
  #limma (lebrando que limma só pode ter as colunas que contém os dados)
  design <- model.matrix(~0+factor(c(1,1,2,2)))
  colnames(design) <- c('Interest', 'Control')
  contrast.matrix <- makeContrasts (Interest-Control, levels=design)
  fit <- lmFit(Log2_data, design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  LIMMA.p <- fit2$p.value
  LIMMA.fdr <- as.data.frame(p.adjust(LIMMA.p,method="BH"))
  #ROTS (ROTs só pode ter as colunas que contém os dados)
  set.seed(123450)
  ROTS.out <- ROTS(data=Log2_data, groups=c(0,0,1,1), a1=0.4, a2=1)
  ROTS.summary <- summary(ROTS.out, fdr = 1)
  ROTS.summary <- data.frame(ROTS.summary[order(ROTS.summary[,1]),])
  #add pvalues
  Log2_data$T.test.p <- T.test.p
  Log2_data$T.test.fdr <- T.test.fdr
  Log2_data$LIMMA.p <- LIMMA.p
  Log2_data$LIMMA.fdr <- LIMMA.fdr [,1]
  Log2_data$ROTS.p <- ROTS.summary$pvalue
  Log2_data$ROTS.fdr <- ROTS.summary$FDR
  name <- basename(filename)
  write.csv(cbind(Log2_data, Raw_data[, -2]), file=paste(name, "p.csv", sep="")) 
}

#call function
stat.trial (Log2_data, Ratio_1vs3)
stat.trial (Log2_data, Ratio_1vs2)
stat.trial (Log2_data, Ratio_3vs4)
stat.trial (Log2_data, Ratio_4vs3)

######Vennn antigo

#Input files 
dv1vs3 <- read.csv("Stat_input_1vs3.csvp.csv", header=TRUE, row.names=1, stringsAsFactors=FALSE)
dv1vs2 <- read.csv("Stat_input_1vs2.csvp.csv", header=TRUE, row.names=1, stringsAsFactors=FALSE)
dv3vs4 <- read.csv("Stat_input_3vs4.csvp.csv", header=TRUE, row.names=1, stringsAsFactors=FALSE)
dv4vs3 <- read.csv("Stat_input_4vs3.csvp.csv", header=TRUE, row.names=1, stringsAsFactors=FALSE)
   
D_Venn <- function(dv)
   {
     dvenn <- dv[ , seq(5, 9, 2)]
     
     t.test.p <-row.names(Venn_select)[which(dvenn$T.test.p <= 0.05)]
     limma.p <- row.names(dvenn)[which(dvenn$LIMMA.p <= 0.05)]
     ROTS.p <- row.names(dvenn)[which(dvenn$ROTS.p <= 0.05)]
     
     f <- venn.diagram(list(t.test.p, limma.p, ROTS.p), 
                       filename = "venn4.tiff",
                       resolution = 300,
                       category.names = c("t.test", "limma", "ROTS"), 
                       lwd = 1.5,
                       lty = 1:3,
                       #col = "transparent",
                       fill = c("mediumseagreen", "cornflowerblue",  "sienna1"),
                       alpha = c(0.6, 0.75, 0.75),
                       cex = 3,
                       fontfamily = "3",
                       fontface = "bold",
                       #cat.col = c("darkorchid2", "cornflowerblue",  "yellow"),
                       cat.pos = -20,
                       cat.dist = c(0.1, 0.1, 0.025),
                       cat.col = c("darkgreen", "darkblue", "darkorange4"),
                       cat.cex = 3,
                       cat.fontface = "italic",
                       cat.fontfamily = "3",
                       cat.default.pos = "outer",
                       main="Random Gene Lists")
   }
     
   D_Venn(dv1vs3)
   D_Venn(dv1vs2)
   D_Venn(dv3vs4)
   D_Venn(dv4vs3)
}
#############################################################################################################
####Volcanos
library(plotly)
library(ggplot2)
library(gridExtra)
#install.packages('Rcpp')
   
##Input files (same as venn)
dv1vs3 <- read.csv("Stat_input_1vs3.csvp.csv", header=TRUE, row.names=1, stringsAsFactors=FALSE)
dv1vs2 <- read.csv("Stat_input_1vs2.csvp.csv", header=TRUE, row.names=1, stringsAsFactors=FALSE)
dv3vs4 <- read.csv("Stat_input_3vs4.csvp.csv", header=TRUE, row.names=1, stringsAsFactors=FALSE)
dv4vs3 <- read.csv("Stat_input_4vs3.csvp.csv", header=TRUE, row.names=1, stringsAsFactors=FALSE)
   
   
##Function
volcano_plot <- function(MATRIZFC, coluna_teste)
   {
   #keep only the fields needed for the plot
   MATRIZFC$Gene.names <- rownames(MATRIZFC)
   rownames(MATRIZFC) <- 1:nrow(MATRIZFC)
   MATRIZFC$coluna_teste <- MATRIZFC[, coluna_teste]
   MATRIZFCcut <- MATRIZFC[c("Gene.names", "log2FC", "coluna_teste")]
   # we want to highlight:
   MATRIZFCcut["group"] <- "NotSignificant"
   MATRIZFCcut[which((MATRIZFCcut["log2FC"]) < 1 & (MATRIZFCcut["log2FC"]) > -1 & MATRIZFCcut["coluna_teste"] <= 0.05), "group"] <- "Significant"
   MATRIZFCcut[which((MATRIZFCcut["log2FC"]) > 1 & MATRIZFCcut["coluna_teste"] < 0.05), "group"] <- "Upregulated" 
   MATRIZFCcut[which((MATRIZFCcut["log2FC"]) < -1 & MATRIZFCcut["coluna_teste"] < 0.05), "group"] <- "Downregulated"
   # Find and label the top peaks..
   top_peaks_d <- subset(MATRIZFCcut, group == "Downregulated", select = c(Gene.names, coluna_teste, log2FC)) 
   top_peaks_u <- subset(MATRIZFCcut, group == "Upregulated", select = c(Gene.names, coluna_teste, log2FC)) 
   top_peaks <- rbind(top_peaks_d,top_peaks_u)
   
   #Labels used in the plot
   #Add gene labels for all of the top genes we found
   a <- list()
   for (i in seq_len(nrow(top_peaks))) 
     {
     m <- top_peaks[i, ]
     a[[i]] <- list(x = m[["log2FC"]],
                    y = -log10(m[["coluna_teste"]]),
                    text = m[["coluna_teste"]],
                    yanchor = "top",
                    font = list(size = 13),
                    xref = "x",
                    yref = "y",
                    showarrow = FALSE) 
      }
   m <- list(l = 70, r = 50, b = 60, t = 60, pad = 4)      
   fx <- list(showline = TRUE, ticks = "outside", showgrid = FALSE, range = list(-4, 4), title = "Expressão diferencial (log2 ratio)", titlefont = list(size = 20), tickfont = list(size = 16), zeroline = FALSE)
   fy <- list(showline = TRUE, ticks = "outside", showgrid = FALSE, range = list(0, 4), title = "Log10 valor-p" , titlefont = list(size = 20), tickfont = list(size = 16), dtick = "L1")
   line1 <- list(type = "line", xref = "x", x0 = -4, x1 = 4, y0 = 1.30102999566398, y1 = 1.30102999566398,
                 line = list(color = 'gray', width = 1, dash = "dot"))
   line2<- list(type = "line", xref = "x", x0 = -1, x1 = -1, yref = "y", y0 = 0, y1= 4,
                line = list(color = 'gray', width = 1, dash = "dot"))
   line3<- list(type = "line", xref = "x", x0 = 1, x1 = 1, yref = "y", y0 = 0, y1= 4,
                line = list(color = 'gray', width = 1, dash = "dot"))
   # make the Plotly plot
   p <- plot_ly(data = MATRIZFCcut, x = ~log2FC, y = ~-log10(coluna_teste), 
                text = ~Gene.names, mode = "markers", size = ~abs(log2FC), 
                color = ~group, colors=c("firebrick3", "gray82", "gray50", "palegreen4")
                ) %>% 
        add_annotations(x= 3.2, y= 1.1, xref = "x", yref = "y", text = "p > 0.05", showarrow = F, font = list(color = "gray50", size = 14)
                ) %>%
        add_annotations(x= 3.2, y= 1.5, xref = "x", yref = "y", text = "p < 0.05", showarrow = F, font = list(color = "gray50", size = 14)
                ) %>%
        add_annotations(x= -3, y= 3.8, xref = "x", yref = "y", text = "Downregulated", showarrow = F, font = list(color = "#CD2626", size = 17)
                ) %>%
        add_annotations(x= 3, y= 3.8, xref = "x", yref = "y", text = "Upregulated", showarrow = F, font = list(color = "#548b54", size = 17)
                ) %>%
        layout(#annotations = a, #para bloquear que apareça
               title = "t test", titlefont = list(size = 22), width = 700, height = 500, 
               margin = m, xaxis = fx, yaxis = fy, shapes = list(line1, line2, line3), showlegend = FALSE)
   p #to load the graph at the viewer in R
   }
   
   volcano_plot(dv1vs3, "T.test.p")
   volcano_plot(dv1vs3, "LIMMA.p")
   volcano_plot(dv1vs3, "ROTS.p")
   
   volcano_plot(dv1vs2, "T.test.p")
   volcano_plot(dv1vs2, "LIMMA.p")
   volcano_plot(dv1vs2, "ROTS.p")
   
   volcano_plot(dv3vs4, "T.test.p")
   volcano_plot(dv3vs4, "LIMMA.p")
   volcano_plot(dv3vs4, "ROTS.p")
   
   volcano_plot(dv4vs3, "T.test.p")
   volcano_plot(dv4vs3, "LIMMA.p")
   volcano_plot(dv4vs3, "ROTS.p")
   
   
   ####Old Venn up
   #######Venn comparação de todos os grupos para o mesmo método estatístico
library(VennDiagram)
library(gplots)

##Input files
all_comp_p <- read.csv("Venn_input_all_p_downup.csv", header=TRUE, row.names=1, stringsAsFactors=FALSE)

all_t_test <- all_comp_p[,c(3,7,8,17)]
all_limma <- all_comp_p[,c(3,7,9,18)]
all_ROTS <- all_comp_p[,c(3,7,10,19)]

####up
Dp_Venn_up <- function(dvp)
{
  colnames(dvp) <- c("log2_1vs3", "log2_4vs3", "x1vs3", "x4vs3")
  c_1vs3 <- row.names(dvp)[which(dvp$x1vs3 <= 0.05 & dvp$log2_1vs3 > 0)]
  c_4vs3 <- row.names(dvp)[which(dvp$x4vs3 <= 0.05 & dvp$log2_4vs3 > 0)]
  
  f <- venn.diagram(list(c_1vs3, c_4vs3), 
                    filename = "venn_up.tiff",
                    resolution = 300,
                    category.names = c("1vs3", "4vs3"), 
                    lwd = 1.5,
                    lty = 1:2,
                    #col = "transparent",
                    fill = c("azure4",  "goldenrod2"),
                    alpha = c(0.6, 0.5),
                    cex = 3,
                    fontfamily = "3",
                    fontface = "bold",
                    cat.pos = c(0, 180),
                    cat.dist = c(0.05, 0.05),
                    cat.col = c("gray31", "darkorange4"),
                    cat.cex = 3,
                    cat.fontface = "italic",
                    cat.fontfamily = "3",
                    cat.default.pos = "outer")
}

Dp_Venn_up(all_t_test)
Dp_Venn_up(all_limma)
Dp_Venn_up(all_ROTS)

############Extrair a lista de genes do venn

###Extract UP
ext_Venn_up <- function(dvp)
{
  colnames(dvp) <- c("log2_1vs3", "log2_4vs3", "x1vs3", "x4vs3")
  c_1vs3 <- row.names(dvp)[which(dvp$x1vs3 <= 0.05 & dvp$log2_1vs3 > 0)]
  c_4vs3 <- row.names(dvp)[which(dvp$x4vs3 <= 0.05 & dvp$log2_4vs3 > 0)]
  inter <- get.venn.partitions(x=list(c_1vs3, c_4vs3))
  return(inter)
}

#rodando funcao para o teste t
inter_t_test_up <- ext_Venn_up(all_t_test)
inter_t_test_up2 <- inter_t_test_up[,4] #escolhendo a coluna 4 para ficar
inter_t_test_13 <- as.data.frame(inter_t_test_up2$`1`)
inter_t_test_13 <- gsub("Sep-07", "Sept7", inter_t_test_13[,1]) #correcao do nome para conversao para gene ID
inter_t_test_13 <- as.data.frame(inter_t_test_13)
#write.csv(inter_t_test_32, "t_32.csv")
inter_t_test_49 <- as.data.frame(inter_t_test_up2$`2`)
#write.csv(inter_t_test_99, "t_99.csv")
inter_t_test_25 <- as.data.frame(inter_t_test_up2$`3`)
inter_t_test_25 <- gsub("Sep-02", "Sept2", inter_t_test_25[,1])
inter_t_test_25 <- as.data.frame(inter_t_test_25)
#write.csv(inter_t_test_43, "t_43.csv")

#rodando funcao para o limma
inter_limma_up <- ext_Venn_up(all_limma)
inter_limma_up2 <- inter_limma_up[,4] #escolhendo a coluna 4 para ficar
inter_limma_22 <- as.data.frame(inter_limma_up2$`1`)
#write.csv(inter_limma_45, "l_45.csv")
inter_limma_92 <- as.data.frame(inter_limma_up2$`2`)
inter_limma_92 <- gsub("Sep-07", "Sept7", inter_limma_92[,1]) #correcao do nome para conversao para gene ID
inter_limma_92 <- as.data.frame(inter_limma_92)
#write.csv(inter_limma_196, "l_196.csv")
inter_limma_15 <- as.data.frame(inter_limma_up2$`3`)
#write.csv(inter_limma_31, "l_31.csv")

#rodando funcao para o ROTS
inter_ROTS_up <- ext_Venn_up(all_ROTS)
inter_ROTS_up2 <- inter_ROTS_up[,4] #escolhendo a coluna 4 para ficar
inter_ROTS_14 <- as.data.frame(inter_ROTS_up2$`1`)
#write.csv(inter_ROTS_26, "r_26.csv")
inter_ROTS_31 <- as.data.frame(inter_ROTS_up2$`2`)
#write.csv(inter_ROTS_77, "r_77.csv")
inter_ROTS_15 <- as.data.frame(inter_ROTS_up2$`3`)
#write.csv(inter_ROTS_28, "r_28.csv")


###Extract Down
ext_Venn_down <- function(dvp)
{
  colnames(dvp) <- c("log2_1vs3", "log2_4vs3", "x1vs3", "x4vs3")
  c_1vs3 <- row.names(dvp)[which(dvp$x1vs3 <= 0.05 & dvp$log2_1vs3 <= 0)]
  c_4vs3 <- row.names(dvp)[which(dvp$x4vs3 <= 0.05 & dvp$log2_4vs3 <= 0)]
  inter <- get.venn.partitions(x=list(c_1vs3, c_4vs3))
  return(inter)
}

#rodando funcao para o teste t
inter_t_test_down <- ext_Venn_down(all_t_test)
inter_t_test_down2 <- inter_t_test_down[,4] #escolhendo a coluna 4 para ficar
inter_t_test_17 <- as.data.frame(inter_t_test_down2$`1`)
inter_t_test_17 <- gsub("Sep-07", "Sept7", inter_t_test_17[,1]) #correcao do nome para conversao para gene ID
inter_t_test_17 <- as.data.frame(inter_t_test_17)
#write.csv(inter_t_test_32, "t_32.csv")
inter_t_test_52 <- as.data.frame(inter_t_test_down2$`2`)
p24 <- strtrim(inter_t_test_52[24, 1], 8)
inter_t_test_52 <- gsub(as.vector(inter_t_test_52[24, 1]), p24, inter_t_test_52[, 1])
inter_t_test_52 <- as.data.frame(inter_t_test_52) 
#write.csv(inter_t_test_99, "t_99.csv")
inter_t_test_20 <- as.data.frame(inter_t_test_down2$`3`)
#write.csv(inter_t_test_43, "t_43.csv")

#rodando funcao para o limma
inter_limma_down <- ext_Venn_down(all_limma)
inter_limma_down2 <- inter_limma_down[,4] #escolhendo a coluna 4 para ficar
inter_limma_22d <- as.data.frame(inter_limma_down2$`1`)
#write.csv(inter_limma_45, "l_45.csv")
inter_limma_105 <- as.data.frame(inter_limma_down2$`2`)
inter_limma_105 <- gsub("Sep-07", "Sept7", inter_limma_105[,1]) #correcao do nome para conversao para gene ID
inter_limma_105 <- as.data.frame(inter_limma_105)
inter_limma_105 <- gsub("Sep-02", "Sept2", inter_limma_105[,1]) #correcao do nome para conversao para gene ID
inter_limma_105 <- as.data.frame(inter_limma_105)
p52 <- strtrim(inter_limma_105[52, 1], 8)
inter_limma_105 <- gsub(as.vector(inter_limma_105[52, 1]), p52, inter_limma_105[, 1])
inter_limma_105 <- as.data.frame(inter_limma_105)
#write.csv(inter_limma_196, "l_196.csv")
inter_limma_17 <- as.data.frame(inter_limma_down2$`3`)
#write.csv(inter_limma_31, "l_31.csv")

#rodando funcao para o ROTS
inter_ROTS_down <- ext_Venn_down(all_ROTS)
inter_ROTS_down2 <- inter_ROTS_down[,4] #escolhendo a coluna 4 para ficar
inter_ROTS_12 <- as.data.frame(inter_ROTS_down2$`1`)
#write.csv(inter_ROTS_26, "r_26.csv")
inter_ROTS_46 <- as.data.frame(inter_ROTS_down2$`2`)
#write.csv(inter_ROTS_77, "r_77.csv")
inter_ROTS_13 <- as.data.frame(inter_ROTS_down2$`3`)
#write.csv(inter_ROTS_28, "r_28.csv")



##Venn para checar correlacao entre os testes estatisticos
f <- venn.diagram(list(inter_t_test_2$`1`, inter_limma_2$`1`, inter_ROTS_2$`1`),
                  filename = "vxx.tiff",
                  resolution = 300,
                  category.names = c("t_test", "limma", "ROTS"), 
                  lwd = 1.5,
                  lty = 1:3,
                  #col = "transparent",
                  fill = c("azure4", "dodgerblue3", "goldenrod2"),
                  alpha = c(0.6, 0.5, 0.5),
                  cex = 3,
                  fontfamily = "3",
                  fontface = "bold",
                  cat.pos = c(0, 0, 180),
                  cat.dist = c(0.05, 0.05, 0.05),
                  cat.col = c("gray31", "darkblue", "darkorange4"),
                  cat.cex = 3,
                  cat.fontface = "italic",
                  cat.fontfamily = "3",
                  cat.default.pos = "outer",
                  main="Stats")

#outra forma para ver a intersecao e quaise sao comuns (Adriene)
inter_ROTS_limma <- inter_ROTS[which(inter_ROTS[ , 1] %in% inter_limma[ , 1]), 1]
inter_ROTS_ttest <- inter_ROTS[which(inter_ROTS[ , 1] %in% inter_t_test[ , 1]), 1]

#################################################################################################
#########Gene ontology das intercecoes
library(clusterProfiler)
library(org.Mm.eg.db)
####New
##Chnging gene symbol to gene ID and preping input for GO
prep_GO <- function(inter_rots, inter_li, inter_t)
           {
           rots_ID <- bitr(inter_rots[ , 1], fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
           limma_ID <- bitr(inter_li[ , 1], fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
           t_test_ID <- bitr(inter_t[ , 1], fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
           return(list(rots_ID[, 2], limma_ID[, 2], t_test_ID[, 2]))
            }

e1 <- prep_GO(inter_ROTS_14,inter_limma_22, inter_t_test_13)
e2 <- prep_GO(inter_ROTS_31,inter_limma_92, inter_t_test_49)
e3 <- prep_GO(inter_ROTS_15,inter_limma_15, inter_t_test_25)

e4 <- prep_GO(inter_ROTS_12,inter_limma_22, inter_t_test_17)
e5 <- prep_GO(inter_ROTS_46,inter_limma_105, inter_t_test_52)
e6 <- prep_GO(inter_ROTS_13,inter_limma_17, inter_t_test_20)


ggo <- function(lista, o)
{
  ggo_r <- enrichGO(gene = lista[[1]],
                    OrgDb = org.Mm.eg.db,
                    keyType = "ENTREZID",
                    ont = o, 
                    #level = 3,
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    #qvalueCutoff  = 0.05,
                    minGSSize = 3,
                    readable = TRUE)
  
  head(ggo_r)
  dotplot(ggo_r, showCategory=30)
  
  ggo_li <- enrichGO(gene = lista[[2]],
                     OrgDb = org.Mm.eg.db,
                     keyType = "ENTREZID",
                     ont = o, 
                     #level = 3,
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     #qvalueCutoff  = 0.05,
                     minGSSize = 3,
                     readable = TRUE)
  
  head(ggo_li)
  dotplot(ggo_li, showCategory=30)
  
  ggo_t <- enrichGO(gene = lista[[3]],
                    OrgDb = org.Mm.eg.db,
                    keyType = "ENTREZID",
                    ont = o, 
                    #level = 3,
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    #qvalueCutoff  = 0.05,
                    minGSSize = 3,
                    readable = TRUE)
  
  head(ggo_t)
  dotplot(ggo_t, showCategory=30)
  
  mr <- merge_result(enrichResultList = list(ROTS=ggo_r, limma=ggo_li, ttest=ggo_t))
  dotplot(mr, showCategory = 10)
}

####UP
ggo(e1, "MF")
ggo(e1, "BP")
ggo(e1, "CC")

ggo(e2, "MF")
ggo(e2, "BP")
ggo(e2, "CC")

ggo(e3, "MF")
ggo(e3, "BP")
ggo(e3, "CC")


   
