library(EGAD)
library(WGCNA)
enableWGCNAThreads(nThreads = 8)
options(scipen=999) # disable scientifc format
options(stringsAsFactors=FALSE)
source("functionalConnectivity/FC.FUN.r")

# Import expression files
la<-read.table("Proteins/3.LA_normCounts_nolowcounts_prot.txt",head=TRUE,sep="\t")
lm<-read.table("Proteins/3.LM_normCounts_nolowcounts_prot.txt",head=TRUE,sep="\t")
ll<-read.table("Proteins/3.LL_normCounts_nolowcounts_prot.txt",head=TRUE,sep="\t")

# Import annotation
laAnno<-read.table("Proteins/la_proteins_functional_terms.txt",head=TRUE,sep="\t")
lmAnno<-read.table("Proteins/lm_proteins_functional_terms.txt",head=TRUE,sep="\t")
llAnno<-read.table("Proteins/lu_proteins_functional_terms.txt",head=TRUE,sep="\t")
anno2list<-function(anno,func=c("GO","KEGG")){
    if(func=="GO")
    {
        term<-as.character(anno[,3])
        clean<-sapply(term,function(x)grep("GO",unlist(strsplit(x,"`")),value=TRUE))
        names(clean)<-paste0(anno[,1],"N")
        clean<-unlist(clean)
        list<-data.frame(gene=gsub("N.*","",names(clean)),go=clean)

    }
    if(func=="KEGG")
    {
        term<-as.character(anno[,2])
        clean<-sapply(term,function(x)grep("KO",unlist(strsplit(x,"`")),value=TRUE))
        names(clean)<-anno[,1]
        clean<-unlist(clean)
        list<-data.frame(gene=names(clean),kegg=clean)
    }
    return(list)
}
kegg<-list( la=anno2list(laAnno,"KEGG"),ll=anno2list(llAnno,"KEGG"),lm=anno2list(lmAnno,"KEGG")  )
go  <-list( la=anno2list(laAnno,"GO"),ll=anno2list(llAnno,"GO"),lm=anno2list(lmAnno,"GO")  )


# Import symbiosis related genes
library(gdata)
symbiosis<-read.xls("Proteins/symbiotic_genes_lupines.xlsx",header=FALSE)
names(symbiosis)<-c("geneID","la","lu","lm","note")
symbiosis$term <-"symbiosis"
laSB <- symbiosis[,c("la","term","geneID")]
llSB <- symbiosis[,c("lu","term","geneID")]
lmSB <- symbiosis[,c("lm","term","geneID")]
sb<-list(la = laSB,lm=lmSB,ll=llSB)
 
# List sample information
info<-data.frame(species =c("L.albus","L.luteus","L. mariae-josephae"), proteinExpr=c("la","ll","lm"),proteinN=c(nrow(la),nrow(ll),nrow(lm)) )

# Test a series of network construction methods:
## binary - rank: 0.99 refers to top 1% of correlations as edges
## binary - Zscore: fisher's Z transformation applied to correlation r, Z=2 corresponding to r=0.964, Z=2.5 corresponding to r=0.987, Z=3 corresponding to r=0.995
## weighted - wgcna, previous wgcna analysis found 0.92-0.96 as hard threshold, and ~24 as soft threshold
netType<-data.frame(net=1:7,
edge=c(rep("binary",5),rep("weighted",2)),
metric=c(rep("rank",2),rep("Zscore",3),rep("wgcna",2)),
cutoff=c(0.99,0.995,2,2.5,3,12,24) )
netType$desc <- paste(netType$edge,netType$metric, netType$cutoff,sep=".")
netType
write.table(netType,"network.type.txt",row.names=FALSE,sep="\t")

# Loop three species, and each loop through network types
for(s in 1:3)
{
    message(paste0("Analyze functional connectivity for ",info$species[s]))
    
    # get expression datasets of exp and obs
    expr<-get(as.character(info$proteinExpr[s]))
    rownames(expr)<-expr$X
    expr<- expr[,-1]
    
    # Peason's correlation
    expr.r <- cor(t(expr),method = "pearson",nThread=8)
    
    # functional connectivity analysis
    for(i in 1:nrow(netType))
    {
        message(paste0("----------------------------------\n...network type ",netType$net[i],": ", netType$desc[i]))
        
        if(netType$edge[i] =="binary")
        {
            net <- build_binary_net(expr.r, metric = netType$metric[i], cutoff = netType$cutoff[i] )
        }
        if(netType$edge[i] =="weighted")
        {
            message(paste0("...build ",netType$edge[i]," network using soft threshold of ",netType$cutoff[i]))
            net <- adjacency(t(expr),power=netType$cutoff[i], type="signed")
        }
        
        
        for(func in c("kegg","go","sb")){
            message(paste0("......",func))
            annotation <- get(func)[[info$proteinExpr[s]]]
            labels <- prep_annot(expr, annotation, min = 10, max = 500)
            if(func=="sb")
            {
                fc = neighbor_voting(labels, net, output = 'AUROC')
                message(fc)
                assign(paste0(func,".fc"),fc)
            }else{
                gba = run_GBA(net, labels)
                message(gba[[3]])
                assign(paste0(func,".fc"),gba)
            }
        }
        save(kegg.fc,go.fc,sb.fc,file=paste0(info$proteinExpr[s],".",netType$desc[i],".rdata"))
    }
    gc()
}