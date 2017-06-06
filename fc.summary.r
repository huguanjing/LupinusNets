
libray(ggplot2)
# place title in center
theme_update(plot.title = element_text(hjust = 0.5), legend.position="bottom", axis.text.x = element_text(angle = 90, hjust = 1))


## GO to folders with saved FC rdata
FCfiles<-grep(".rdata",list.files(),value=TRUE)
FCfiles
rm(sumFC)
rm(sumCorr)
for(file in FCfiles)
{
    print(file)
    fn<-load(file) #"kegg.fc" "go.fc"   "sb.fc"
    
    r<-data.frame(annotation=gsub(".fc","",fn), label=file,
    auc=c(kegg.fc[[3]],go.fc[[3]],sb.fc[1]) )

    if(exists("sumFC"))
    {
        sumFC<- rbind(sumFC,r)
    }else{sumFC<-r}
    
}

sumFC$species<-gsub("[.].*","",sumFC$label )
sumFC$netType<-gsub("l.[.]|.rdata","",sumFC$label )
sumFC$bw<-gsub("[.].*","",sumFC$netType)
write.table(sumFC,"summaryFC.txt",sep="\t",row.names=FALSE)

########################
## Examine FC results ##
########################
names(sumFC)
xtabs(auc~netType+species+annotation,data=sumFC)
summary(a1<-aov(auc~annotation+netType+species,data=sumFC))
summary(a2<-aov(auc~annotation+bw+species,data=sumFC))

TukeyHSD(x=a2, 'bw', conf.level=0.95)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# Fit: aov(formula = auc ~ annotation + bw + species, data = sumFC)
# $bw
#            diff        lwr        upr         p     adj
# weighted-binary 0.08871495 0.06231931 0.1151106     0

TukeyHSD(x=a2, 'species', conf.level=0.95)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# Fit: aov(formula = auc ~ annotation + bw + species, data = sumFC)
# $species
# diff         lwr        upr     p adj
# ll-la  0.07303901  0.03793831 0.10813971 0.0000167
# lm-la  0.05125895  0.01615825 0.08635965 0.0024754
# lm-ll -0.02178006 -0.05688076 0.01332064 0.3017599

TukeyHSD(x=a2, 'annotation', conf.level=0.95)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# Fit: aov(formula = auc ~ annotation + bw + species, data = sumFC)
# $annotation
# diff         lwr        upr     p adj
# kegg-go 0.02203047 -0.01307024 0.05713117 0.2937203
# sb-go   0.08073125  0.04563055 0.11583195 0.0000024
# sb-kegg 0.05870078  0.02360008 0.09380148 0.0004936


################
## Make plots ##
################
sumFC<-read.table("summaryFC.txt",sep="\t",header=TRUE)

pdf("plotFC.auc.pdf")
# wgcna outperforms others networks
ggplot(data=sumFC, aes(x=netType, y=auc, fill=bw)) + facet_grid(.~annotation) + geom_boxplot() +ggtitle("Functional connectivity: AUCs")
# wgcna.12
ggplot(data=sumFC[sumFC$bw=="weighted",], aes(x=annotation, fill=species, y=auc,)) + facet_grid(.~netType)+ geom_bar(stat="identity",position =position_dodge()) + ggtitle("Weighted Networks - Functional connectivity: AUCs")
dev.off()


####################################
## Highly connected GO categories ##
####################################
files <- grep("wgcna.24",FCfiles,value=TRUE)
for(file in files)
{
    print(file)
    fn<-load(file) #"kegg.fc" "go.fc"   "sb.fc"
    species <- gsub("[.].*","",file)
    
    r <- c(kegg.fc[[1]][,1], go.fc[[1]][,1])
    r<-data.frame(id=names(r),auc=as.numeric(r))
    names(r)[2]<-species
    
    if(exists("AUCs"))
    {
        AUCs<- merge(AUCs,r,by="id",all.x=TRUE,all.y=TRUE)
    }else{AUCs<-r}
    
}

# get GO terms
library(GO.db)
goterms = unlist(Term(GOTERM))
AUCs$Term <- goterms[AUCs$id]
write.table(AUCs,file="go.aucs.txt",row.names=FALSE,sep="\t")

######################################
## Highly connected KEGG categories ##
######################################

# get KEGG info
library(KEGG.db)
# I noticed that KO

KO:K09286

KEGG.db contains mappings based on older data because the original
resource was removed from the the public domain before the most
recent update was produced. This package should now be considered
deprecated and future versions of Bioconductor may not have it
available.  Users who want more current data are encouraged to look
at the KEGGREST or reactome.db packages

> xx <- as.list(KEGGPATHID2NAME)