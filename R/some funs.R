#hypertest():批量做mRNA和lncRNA的超几何分布检验
hypertest = function(lnc,pc,deMIR = NULL,lnctarget,pctarget){
  pcs = unique(pctarget[,1])
  pctarget = lapply(unique(pctarget[,1]), function(x){
    pctarget[,2][pctarget[,1]==x]
  })
  names(pctarget) = unique(pcs)

  lncs = unique(lnctarget[,1])
  lnctarget = lapply(unique(lnctarget[,1]), function(x){
    lnctarget[,2][lnctarget[,1]==x]
  })
  names(lnctarget) = lncs

  mir1 <- unique(unlist(lnctarget))
  mir2 <- unique(unlist(pctarget))

  mirs <- union(mir1,mir2) #合集
  popTotal <- length(mirs) #合集长度

  ceLNC <- lnc[lnc %in% names(lnctarget)] #有miRNA的lnc
  cePC <- pc[pc %in% names(pctarget)] #有miRNA的pc
  #ceMIR <- mir[mir %in% mirs]

  hyperOutput <- list()
  i <- 0
  for (lncID in ceLNC) {
    listTotal <- length(lnctarget[[lncID]]) #lnc 的miRNA数量
    for (gene in cePC) {
      i = i + 1
      ovlp <- intersect(lnctarget[[lncID]], pctarget[[gene]]) #交集

      popHits <- length(pctarget[[gene]])
      Counts <- length(ovlp)

      ovlpMIRs <- paste(ovlp, collapse = ',')
      foldEnrichment <- Counts/listTotal*popTotal/popHits
      pValue <- phyper(Counts-1, popHits, popTotal-popHits,
                       listTotal, lower.tail=FALSE, log.p=FALSE)

      ceMIR <- Reduce(intersect, list(ovlp, deMIR))
      deMIRs <- paste(ceMIR, collapse = ',')
      deMIRCounts <- length(ceMIR)

      hyperOutput[[i]] <- c(lncID, gene, Counts, listTotal,
                            popHits,popTotal,foldEnrichment,pValue,ovlpMIRs,
                            deMIRCounts, deMIRs)

    }
  }

  #hyperOutput <- Reduce(rbind, hyperOutput)  ## slower
  hyperOutput <- do.call(rbind, hyperOutput)
  #hyperOutput <- rbind_list(hyperOutput) ## not test

  colnames(hyperOutput) <- c('lncRNAs','Genes','Counts','listTotal',
                             'popHits','popTotal','foldEnrichment','hyperPValue','miRNAs',
                             'deMIRCounts','deMIRs')
  hyperOutput <- as.data.frame(as.matrix(hyperOutput),
                               stringsAsFactors=FALSE)
  hyperOutput <- hyperOutput[as.numeric(hyperOutput$Counts)>0,]

  #hyperOutput$FDR <- p.adjust(as.numeric(as.character(hyperOutput$pValue)),
  #method = 'fdr')
  #hyperOutput <- hyperOutput[hyperOutput$Counts>0,]
  #hyperOutput$lncRNAs <- ensembl2symbolFun(hyperOutput$lncRNAs)
  #hyperOutput$gene <- ensembl2symbolFun(hyperOutput$gene)

  if (is.null(deMIR)) {
    hyperOutput <- hyperOutput[,! colnames(hyperOutput) %in%
                                 c('deMIRCounts','deMIRs')]
  }
  hyperOutput = hyperOutput[order(hyperOutput$hyperPValue),]
  return (hyperOutput)
  message(paste0(sum(hyperOutput$hyperPValue<0.05)," pairs p<0.05"))
}

############################################################################################################
#plcortest():批量做mRNA和lncRNA的相关性检验
plcortest <- function(lnc_exp, mRNA_exp,cor_cutoff=0) {
  jp = list()
  for(i in 1:nrow(lnc_exp)){
    x = cbind(lnc_exp[i,],t(mRNA_exp))
    jp[[i]] = vector()
    for(j in 2: ncol(x)){
      k = cor.test(x[,1],x[,j])$p.value <0.05
      k2 = cor.test(x[,1],x[,j])$estimate >cor_cutoff
      if(k&k2){jp[[i]] = c(jp[[i]],colnames(x)[[j]])}
    }
  }
  names(jp) = rownames(lnc_exp)
  return(jp)
}



############################################################

library(sciplot)
library(reshape)
#载入绘图数据
err <- read.delim("clipboard")
head(err)
#将指定为x轴和分组的列设置为因子
err$Group <- as.factor(err$type)
err$Sample <- as.factor(err$type)
#对数据进行格式处理
err.m <- melt(err)
#绘制第一幅图
bargraph.CI(Sample, value, group = Group, data =err.m, 
            xlab = "Sample", ylab = "Number", 
            err.col = c(1), col = c(1))
#添加图例
legend("topright", legend = c("A","B"), bty = "n", 
       horiz = F, fill = c(1,2))

#除了绘制相邻条形图之外，bargraph.CI还可以使用split参数将两组相邻的条形改为对称形式。
bargraph.CI(Sample, value, group = Group, data =err.m, 
            xlab = "Sample", ylab = "Number", 
            err.col = c(1,2), col = c(1,2),split = TRUE)
legend("topright", legend = c("A","B"), bty = "n", 
       horiz = F, fill = c(1,2))





















