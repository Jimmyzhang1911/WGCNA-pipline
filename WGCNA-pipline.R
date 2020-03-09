rm(list=ls())
library(WGCNA)
library(dplyr)
library(ggplot2)
library(stringr)
library(pheatmap)
enableWGCNAThreads() #开启多核模式
stringsAsFactors = F
RNAseq_voom<-read.table("data", sep="\t",header=T, stringsAsFactors = F)

# 归一化函数，x为data，y=1按行，y=2按列 
nor<-function(x, y=1){
  a=nrow(x)
  b=ncol(x)
  x_min_tmp<-apply(x,y,min)
  x_extrem_tmp<-apply(x,y,max)-apply(x,y,min)
  if(y==1){
    x_min<-matrix(rep(x_min_tmp,b),byrow=F,ncol=b)
    x_extrem_tmp<-matrix(rep(x_extrem_tmp,b),byrow=F,ncol=b)}else{
      x_min<-matrix(rep(x_min_tmp,a),byrow=T,ncol=b)
      x_extrem_tmp<-matrix(rep(x_extrem_tmp,a),byrow=T,ncol=b)
    }
  tmp<-abs(x-x_min)/x_extrem_tmp
  return (tmp)
} 
#RNAseq_voom<-nor(RNAseq_voom)

## 因为WGCNA针对的是基因进行聚类，而一般我们的聚类是针对样本用hclust即可，所以这个时候需要转置。
WGCNA_matrix = t(RNAseq_voom[order(apply(RNAseq_voom,1,mad), decreasing = T)[1:10000],])
datExpr0 <- WGCNA_matrix  ## top 1000 mad genes
datExpr <- datExpr0 
head(datExpr)[1:5,1:5]

# 制作datTraits
if(T){
  trait<-strsplit(rownames(datExpr),"_")
  size<-factor(sapply(trait, "[", 2),levels = c("small","big"))
  day<-factor(sapply(trait, "[", 1),levels = c("X20d","X25d","X30d"))
  datTraits<-data.frame(day,size)
  rownames(datTraits)<-rownames(datExpr)
  glimpse(datTraits$size)}
save(datTraits,datExpr,file="input.RData")

#load("sft-net.RData")
#load("input.RData")

#（一）查看样品聚类情况
if(T){
  sampleTree = hclust(dist(datExpr), method = "average")
  plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")}


#（二）选定合适的软阈值
if(T){
  powers = c(c(1:10), seq(from = 12, to=20, by=2))
  # Call the network topology analysis function
  sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
  #设置网络构建参数选择范围，计算无尺度分布拓扑矩阵
  
  # Plot the results:
  ##sizeGrWindow(9, 5)
  par(mfrow = c(1,2));
  cex1 = 0.9;
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  # this line corresponds to using an R^2 cut-off of h
  abline(h=0.90,col="red")
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")}
glimpse(sft)


#（三）构建共表达矩阵
if(T){
  net = blockwiseModules(
    datExpr,
    power = sft$powerEstimate,
    maxBlockSize = 6000,
    TOMType = "unsigned", minModuleSize = 30,
    reassignThreshold = 0, mergeCutHeight = 0.25,
    numericLabels = TRUE, pamRespectsDendro = FALSE,
    saveTOMs = T, 
    verbose = 3)
    table(net$colors)}
glimpse(net)
load(net$TOMFiles[1], verbose=T) #提取TOM
save(sft, net, file="sft-net.RData")


#（四）模块的可视化
if(T){
  # Convert labels to colors for plotting
  mergedColors = labels2colors(net$colors)
  table(mergedColors)
  # Plot the dendrogram and the module colors underneath
  plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  ## assign all of the gene to their corresponding module 
  ## hclust for the genes.
}


#（五）关联模块和性状
head(datExpr)[1:6,1:6];glimpse(datExpr) #关联之前瞅一眼数据
head(datTraits);glimpse(datTraits)
if(T){
  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr)
  design=model.matrix(~0+ datTraits$size)
  colnames(design)=levels(datTraits$size)
  moduleColors <- labels2colors(net$colors)
  # Recalculate MEs with color labels
  MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
  MEs = orderMEs(MEs0); ##不同颜色的模块的ME值矩 (样本vs模块)
  moduleTraitCor = cor(MEs, design , use = "p");
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
  
  sizeGrWindow(10,6)
  # Will display correlations and their p-values
  textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")", sep = "");
  dim(textMatrix) = dim(moduleTraitCor)
  pdf("module-trait-relations.pdf")
  par(mar = c(6, 8.5, 3, 3));
  # Display the correlation values within a heatmap plot
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = colnames(design),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = greenWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
  dev.off()}


#（六）感兴趣性状的模块的具体基因分析
#基因的表达值与模块有相关性，同时与连续型性状也有相关性，
if(T){
  # names (colors) of the modules
  modNames = substring(names(MEs), 3) #从第三位开始取字符，原本MEbrown
  geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
  ## 算出每个模块跟基因的皮尔森相关系数矩阵
  ## MEs是每个模块在每个样本里面的值
  ## datExpr是每个基因在每个样本的表达量
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
  names(geneModuleMembership) = paste("MM", modNames, sep="");
  names(MMPvalue) = paste("p.MM", modNames, sep="");
  ## 只有连续型性状才能只有计算
  
  ## 这里把是否属于 Luminal 表型这个变量用0,1进行数值化。
  big = as.data.frame(design[,2]);
  names(big) = "big"
  geneTraitSignificance = as.data.frame(cor(datExpr, big, use = "p"));
  GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
  names(geneTraitSignificance) = paste("GS.", names(big), sep="");
  names(GSPvalue) = paste("p.GS.", names(big), sep="");
  
  module = "brown"
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
  pdf("scatter.pdf")
  sizeGrWindow(7, 7);
  par(mfrow = c(1,1));
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = "Gene significance for size",
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
  dev.off()}


#（七）网络的可视化，占用超大资源
if(T){
  load(net$TOMFiles[1], verbose=T)
  TOM <- as.matrix(TOM)
  dissTOM = 1-TOM
  # Transform dissTOM with a power to make moderately strong 
  # connections more visible in the heatmap
  pdf("Network heatmap plot, all genes.pdf")
  plotTOM = dissTOM^7
  # Set diagonal to NA for a nicer plot
  diag(plotTOM) = NA
  # Call the plot function
  
  # 这一部分特别耗时，行列同时做层级聚类
  TOMplot(plotTOM, net$dendrograms, moduleColors, 
          main = "Network heatmap plot, all genes")
  dev.off()}
if(T){
  # Recalculate module eigengenes
  MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
  ## 只有连续型性状才能只有计算
  ## 这里把是否属于 Luminal 表型这个变量用0,1进行数值化。
  big = as.data.frame(design[,2]);
  names(big) = "big"
  # Add the weight to existing module eigengenes
  MET = orderMEs(cbind(MEs, big))
  # Plot the relationships among the eigengenes and the trait
  sizeGrWindow(5,7.5);
  par(cex = 0.9)
  plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                        = 90)
  # Plot the dendrogram
  sizeGrWindow(6,6);
  par(cex = 1.0)
  ## 模块的聚类图
  plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                        plotHeatmaps = FALSE)
  # Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
  par(cex = 1.0)
  ## 性状与模块热图
  plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                        plotDendrograms = FALSE, xLabelsAngle = 90)
  }#最后画模块和性状的关系


#（八）模块中基因的提取和与性状的可视化
if(T){
  # Select module
  module = "brown";
  # Select module probes
  probes = colnames(datExpr) ## 我们例子里面的probe就是基因
  inModule = (moduleColors==module);
  modProbes = probes[inModule]; 
  head(modProbes)
  
  # 如果使用WGCNA包自带的热图就很丑。
  which.module="brown";
  dat=datExpr[,moduleColors==which.module ] 
  plotMat(t(scale(dat)),nrgcols=30,rlabels=T,
          clabels=T,rcols=which.module,
          title=which.module )
  datExpr[1:4,1:4]
  dat=t(datExpr[,moduleColors==which.module ] )
  library(pheatmap)
  pheatmap(dat ,show_colnames =F,show_rownames = F) #对那些提取出来的1000个基因所在的每一行取出，组合起来为一个新的表达矩阵
  n=t(scale(t(log(dat+1)))) # 'scale'可以对log-ratio数值进行归一化
  n[n>2]=2 
  n[n< -2]= -2
  n[1:4,1:4]
  pheatmap(n,show_colnames =F,show_rownames = F)
  group_list=datTraits$size
  ac=data.frame(g=group_list)
  rownames(ac)=colnames(n) 
  pheatmap(n,show_colnames =F,show_rownames = F,
           annotation_col=ac )
  # 可以很清晰的看到，所有的形状相关的模块基因
  # 其实未必就不是差异表达基因。
}#模块中基因和性状的可视化


#导出模块中的基因，输入颜色，vis，和cyt
getgene<-function(module = "color",vis=F,cyt=F){
  load(net$TOMFiles[1], verbose=T)
  if(vis==T){
    vis = exportNetworkToVisANT(modTOM,
                                file = paste("VisANTInput-", module, ".txt", sep=""),
                                weighted = TRUE,
                                threshold = 0)}
  else if(cyt==T){
    cyt = exportNetworkToCytoscape(
      modTOM,
      edgeFile = paste("CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
      nodeFile = paste("CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
      weighted = TRUE,
      threshold = 0.02,
      nodeNames = modProbes, 
      nodeAttr = moduleColors[inModule])}
  probes = colnames(datExpr) ## 我们例子里面的probe就是基因
  inModule = (moduleColors==module)
  modProbes = probes[inModule]
  return (modProbes)} #提取模块中的基因和导出

#进一步筛选
if(T){
  nTop = 30;
  IMConn = softConnectivity(datExpr[, modProbes]);
  top = (rank(-IMConn) <= nTop)
  filter <- modTOM[top, top]} #cyt太多可以过滤再导出

bule<-getgene(module="blue",vis=F,cyt=T)
