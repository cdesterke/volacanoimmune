devtools::install_github('dviraran/xCell')

library(xCell)
exprMatrix = read.table(expr_file,header=TRUE,row.names=1, as.is=TRUE)
xCellAnalysis(exprMatrix)


# usage : annot<-(data)
prepannot<-function(data)
	{
   	if(!require(tidyr)){install.packages("tidyr")}
    	library(tidyr)

		id<-colnames(data)
		id<-as.data.frame(id)
		row.names(id)<-id$id
		annot<-separate(data = id, col = id, into = c("group", "sample"), sep = "_")
		return(annot)
}



# usage : res<-limmaxcell(data,annot,control="HD")

limmaxcell<-function(data,annot,control="HD")
	{
		if(!require(limma)){
		
		if (!requireNamespace("BiocManager", quietly = TRUE))
			install.packages("BiocManager")
			BiocManager::install("limma")}
			
		library(limma)	

	# remove scores from xcell result before performing limma analysis
		data<-data[!grepl("*Score", rownames(data)), ]
	
	# prepare design matrix
	levels <- relevel(as.factor(annot$group),ref=control)
	design <- model.matrix(~levels)
	rownames(design) = colnames(data)
	
	# perform limma analysis
	tmp <- lmFit(data,design=design)
	fit <- eBayes(tmp)
	res = topTable(fit,number = nrow(data),coef=2)
	res
}






# usage : volxcell(res,thresp=0.5,thresfc=0.025,fold=45,title="Immune score post LIMMA")


volxcell<-function(res,thresp=0.5,thresfc=0.1, title="Volcano Plot on IMMUNE CELLS",fold=50)
	{
	# remove scores if not performed during limma analysis
	#res<-res[!grepl("*Score", rownames(res)), ]

	# set column of gene identifier starting of row names
		res$geneid<-rownames(res)
	
	# defined the x axis limits
	
		xaxis<-max(res$logFC)+ 0.5*max(res$logFC)
	
	# found up and down cells with fold change superior to theshold defined
		up <- rownames(res)[which(res$logFC>=thresfc)]
		
		down <- rownames(res)[which(res$logFC<=(1*-thresfc))]
		
	
	# conditional color and shape of the cells in links with their p-values and fold changes
		with(res, plot(logFC, -log10(P.Value), pch=8, col="#FFA500", main=title, xlim=c(-xaxis,xaxis),cex=1))
	
		with(subset(res , P.Value<thresp), points(logFC, -log10(P.Value), pch=20,cex=1.5,col="darkgray"))
	
		with(res[which(rownames(res) %in% up),], points(logFC, -log10(P.Value), pch=20,cex=AveExpr*fold, col="#FF00FF"))
	
		with(res[which(rownames(res) %in% down),], points(logFC, -log10(P.Value), pch=20, cex=AveExpr*fold, col="#00FA9A"))
	
	# add different threshold lines
		abline(h=-log10(thresp), col="darkgray",lty=2,lwd = 2)
	
		abline(v=thresfc, col="#FF00FF",lty=2,lwd = 2)
	
		abline(v=-thresfc, col="#00FA9A",lty=2,lwd = 2)

	# load calibrate library to write labels on plot
		if(!require(calibrate)){install.packages("calibrate")}
			library(calibrate)
	
	# add labels for regulated cells with significant p-values in limma 
		with(subset(res, P.Value<thresp & abs(logFC) >= thresfc), textxy(logFC, -log10(P.Value), labs=geneid, 
		cex=0.8,lwd=2,col="darkred"))
}
