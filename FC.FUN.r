##############
## Function ##
##############

## function to build binary network
library(psych)
build_binary_net <- function(gene.corr, metric = c("rank","Zscore","raw","none"), cutoff = 0.99 )
{
    message(paste0("...build binary network based on ",metric," of correlation coefficients using cutoff of ", cutoff))
    if (metric == "none")
    {
        net<-gene.corr
    }
    else{
        nGenes = ncol(gene.corr)
        if (metric == "rank") {
            # sample 10% of genes
            sub<-sample(1:nGenes,ceiling(nGenes/10))
            sub.corr<-gene.corr[sub,sub]
            threshold <- quantile(sub.corr[upper.tri(sub.corr)],cutoff)
        }
        else if (metric == "Zscore") {
            threshold <- fisherz2r(cutoff)
        } 
        else if (metric == "raw") {
            threshold <- cutoff
        }
        net<-ifelse(gene.corr<threshold,0,1)
        rownames(net)<-rownames(gene.corr)
        colnames(net)<-colnames(gene.corr)
        # message(paste0("...threshold=",threshold,", network containing ",nGenes," genes with ", (table(net)["1"]-nGenes)/2, " edges inferred."))
    }
    return(net)
}

## function to run functional connectivity by EGAD
run_egad<-function(network, annotation)
{
    terms<- unique(annotation[,2])
    annot <- make_annotations(annotation[,1:2],rownames(network),terms)
    # Run GBA
    res <- run_GBA(network, annot)
    return(res)
}
# res as a list of three elements: auc of each functional category, prediction of gene functions, average auc

## function to prepare annotation
prep_annot<-function(expr, annotation, min = 10, max = 500)
{
    # prep annotation
    terms<- unique(annotation[,2])
    annot <- make_annotations(annotation[,1:2],rownames(expr),terms)
    # restrict annotation to given size (min, max)
    labels <- filter_network_cols(annot, min, max)
    return(labels)
}

## function to run neigbour voting for subgenomes seperately
run_nv_homoeo<-function(network, annotation, min = 10, max = 500, nfold = 3)
{
    # prep annotation
    terms<- unique(annotation[,2])
    annot <- make_annotations(annotation[,1:2],rownames(network),terms)
    # restrict annotation to given size (min, max)
    labels <- filter_network_cols(annot, min, max)
    # Run neigbour voting
    ids<-rownames(network)
    netA<-network[grep("a$",ids),grep("a$",ids)]
    netD<-network[grep("d$",ids),grep("d$",ids)]
    message(paste0("......",ncol(netA)," At vs ",ncol(netD)," Dt genes."))
    res<-list(A =neighbor_voting(labels, netA, output = 'AUROC'), D=neighbor_voting(labels, netD, output = 'AUROC'))
    return(res)
}
# res as a list of three elements: auc of each functional category, prediction of gene functions, average auc

## function to compare neigbour voting results between subgenomes seperately
compare_nvs<-function(obs.fc, exp.fc, title="test", save.rdata=FALSE)
{
    # Wilcoxon rank-sum test:
    par(mfrow=c(2,2),oma = c(0, 0, 2, 0))
    plot_auc_compare(exp.fc$A[,1], exp.fc$D[,1], xlab ="AUC - At",ylab="AUC - Dt", title="Expected",pch=".")
    plot_auc_compare(obs.fc$A[,1], obs.fc$D[,1], xlab ="AUC - At",ylab="AUC - Dt", title="Observed",pch=".")
    plot_auc_compare(exp.fc$A[,1], obs.fc$A[,1], xlab ="AUC - exp",ylab="AUC - obs", title="At",pch=".")
    plot_auc_compare(exp.fc$D[,1], obs.fc$D[,1], xlab ="AUC - exp",ylab="AUC - obs", title="Dt",pch=".")
    mtext(title, outer = TRUE, cex = 2)
    if(save.rdata){save(obs.fc,exp.fc,file=paste0(title,".rdata"))}
}


## Functions: scatter plot of two AUCs
plot_auc_compare<-function (aucA, aucB, xlab = "AUROC", ylab = "AUROC", pch=20, xlim = c(0,
1), ylim = c(0, 1),title="comparison")
{
    a <- mean(aucA, na.rm = TRUE)
    b <- mean(aucB, na.rm = TRUE)
    plot(aucA, aucB, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab,
    bty = "n", cex.lab = 1.5, cex.axis = 1.2, pch = pch, main=title,cex.main=1.5)
    abline(0, 1, lwd = 3, col = "grey")
    abline(v = a, lty = 2, col = 2)
    abline(h = b, lty = 2, col = 2)
    text(0.8,1,paste0("r=",round(cor(aucA, aucB),3)))
}
