# load data for bacteria DE
`B29.genes` <- read.delim("/home/m/maggi-brisbin/SiphComp/RSEM/B29.genes.results")
`B30.genes` <- read.delim("/home/m/maggi-brisbin/SiphComp/RSEM/B30.genes.results")
`B31.genes` <- read.delim("/home/m/maggi-brisbin/SiphComp/RSEM/B31.genes.results")
`B32.genes` <- read.delim("/home/m/maggi-brisbin/SiphComp/RSEM/B32.genes.results")
`LO2.genes` <- read.delim("/home/m/maggi-brisbin/SiphComp/RSEM/LO2.genes.results")
`LO3.genes` <- read.delim("/home/m/maggi-brisbin/SiphComp/RSEM/LO3.genes.results")
`LO4.genes` <- read.delim("/home/m/maggi-brisbin/SiphComp/RSEM/LO4.genes.results")
`LO5.genes` <- read.delim("/home/m/maggi-brisbin/SiphComp/RSEM/LO5.genes.results")

#make data frames
Genes<- data.frame(B29.genes$gene_id,B29.genes$expected_count,B30.genes$expected_count, B31.genes$expected_count,B32.genes$expected_count, LO2.genes$expected_count,LO3.genes$expected_count, LO4.genes$expected_count,LO5.genes$expected_count,LO5.genes$expected_count)
FPKM<- data.frame(B29.genes$gene_id,B29.genes$FPKM,B30.genes$FPKM,B31.genes$FPKM, B32.genes$FPKM, LO2.genes$FPKM, LO3.genes$FPKM,LO4.genes$FPKM,  LO5.genes$FPKM )
ercc_conc<-cms_095046 <- read.delim("~/Desktop/Siphamia/siph.genes.results/cms_095046.txt", stringsAsFactors=FALSE)
ercc_conc<-ercc_conc[,2:7]
ercc_conc<-ercc_conc[order(ercc_conc$ERCC.ID),]
ERCC<-FPKM[1:92,]

names(Genes)<-c("gene_id","B29","B30","B31", "B32","LO2","LO3","LO4", "LO5", "LO5b")
names(ERCC)<-c("gene_id","B29","B30","B31", "B32","LO2","LO3","LO4", "LO5" )

#Make data frame with just sample FPKM values for ERCC genes and concentration of ERCC genes in spike in mix
# make sure items in data frame are integers

#LO2 mix 1
LO2spkindf<-data.frame(ERCC$LO2,ercc_conc$concentration.in.Mix.1..attomoles.ul.)
names(LO2spkindf)<-c("FPKM", "conc")
LO2spkindf$FPKM<-as.integer(LO2spkindf$FPKM)
LO2spkindf$conc<-as.integer(LO2spkindf$conc)
row_sub = apply(LO2spkindf, 1, function(row) all(row !=0 ))
LO2spkindf=LO2spkindf[row_sub,]
LO2spkindf$FPKM<-log2(LO2spkindf$FPKM)
LO2spkindf$conc<-log2(LO2spkindf$conc)

#LO3 mix 2
LO3spkindf<-data.frame(ERCC$LO3,ercc_conc$concentration.in.Mix.2..attomoles.ul.)
names(LO3spkindf)<-c("FPKM", "conc")
LO3spkindf$FPKM<-as.integer(LO3spkindf$FPKM)
LO3spkindf$conc<-as.integer(LO3spkindf$conc)
row_sub = apply(LO3spkindf, 1, function(row) all(row !=0 ))
LO3spkindf=LO3spkindf[row_sub,]
LO3spkindf$FPKM<-log2(LO3spkindf$FPKM)
LO3spkindf$conc<-log2(LO3spkindf$conc)

#LO4 mix 2
LO4spkindf<-data.frame(ERCC$LO4,ercc_conc$concentration.in.Mix.2..attomoles.ul.)
names(LO4spkindf)<-c("FPKM", "conc")
LO4spkindf$FPKM<-as.integer(LO4spkindf$FPKM)
LO4spkindf$conc<-as.integer(LO4spkindf$conc)
row_sub = apply(LO4spkindf, 1, function(row) all(row !=0 ))
LO4spkindf=LO4spkindf[row_sub,]
LO4spkindf$FPKM<-log2(LO4spkindf$FPKM)
LO4spkindf$conc<-log2(LO4spkindf$conc)

#LO5 mix 1
LO5spkindf<-data.frame(ERCC$LO5,ercc_conc$concentration.in.Mix.1..attomoles.ul.)
names(LO5spkindf)<-c("FPKM", "conc")
LO5spkindf$FPKM<-as.integer(LO5spkindf$FPKM)
LO5spkindf$conc<-as.integer(LO5spkindf$conc)
row_sub = apply(LO5spkindf, 1, function(row) all(row !=0 ))
LO5spkindf=LO5spkindf[row_sub,]
LO5spkindf$FPKM<-log2(LO5spkindf$FPKM)
LO5spkindf$conc<-log2(LO5spkindf$conc)

#b29 mix 1
B29spkindf<-data.frame(ERCC$B29,ercc_conc$concentration.in.Mix.1..attomoles.ul.)
names(B29spkindf)<-c("FPKM", "conc")
B29spkindf$FPKM<-as.integer(B29spkindf$FPKM)
B29spkindf$conc<-as.integer(B29spkindf$conc)
row_sub = apply(B29spkindf, 1, function(row) all(row !=0 ))
B29spkindf=B29spkindf[row_sub,]
B29spkindf$FPKM<-log2(B29spkindf$FPKM)
B29spkindf$conc<-log2(B29spkindf$conc)

#b30 mix 2
B30spkindf<-data.frame(ERCC$B30,ercc_conc$concentration.in.Mix.2..attomoles.ul.)
names(B30spkindf)<-c("FPKM", "conc")
B30spkindf$FPKM<-as.integer(B30spkindf$FPKM)
B30spkindf$conc<-as.integer(B30spkindf$conc)
row_sub = apply(B30spkindf, 1, function(row) all(row !=0 ))
B30spkindf=B30spkindf[row_sub,]
B30spkindf$FPKM<-log2(B30spkindf$FPKM)
B30spkindf$conc<-log2(B30spkindf$conc)

#b31 mix 2
B31spkindf<-data.frame(ERCC$B31,ercc_conc$concentration.in.Mix.2..attomoles.ul.)
names(B31spkindf)<-c("FPKM", "conc")
B31spkindf$FPKM<-as.integer(B31spkindf$FPKM)
B31spkindf$conc<-as.integer(B31spkindf$conc)
row_sub = apply(B31spkindf, 1, function(row) all(row !=0 ))
B31spkindf=B31spkindf[row_sub,]
B31spkindf$FPKM<-log2(B31spkindf$FPKM)
B31spkindf$conc<-log2(B31spkindf$conc)

#b32 mix 1
B32spkindf<-data.frame(ERCC$B32,ercc_conc$concentration.in.Mix.1..attomoles.ul.)
names(B32spkindf)<-c("FPKM", "conc")
B32spkindf$FPKM<-as.integer(B32spkindf$FPKM)
B32spkindf$conc<-as.integer(B32spkindf$conc)
row_sub = apply(B32spkindf, 1, function(row) all(row !=0 ))
B32spkindf=B32spkindf[row_sub,]
B32spkindf$FPKM<-log2(B32spkindf$FPKM)
B32spkindf$conc<-log2(B32spkindf$conc)

#function creating text for line equation and r squared value
#enter name of data frame
lm_eqn <- function(df){
  y<-df$FPKM
  x<-df$conc
  m <- lm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(coef(m)[1], digits = 2), 
                        b = format(coef(m)[2], digits = 2), 
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));                 
}

# plot ERCC genes FPKM values v. the concentration of genes from ERCC spikein
#use Log2 x and y axis
#add line equation and rsquared to the plot
## tried to add linear regression line and could not get past the error messages
library("ggplot2", lib.loc="~/Desktop/R-3.2.2/library")
#LO2
spikeinPlotLO2<-ggplot(LO2spkindf, aes(x=conc,y=FPKM)) +geom_point()+ggtitle("LO2/spike-in 1")+ geom_smooth(method=lm, se=FALSE, color="black") +
  theme(panel.background =element_blank(), panel.grid.major = element_blank(), panel.grid.minor =element_blank(), panel.border =element_blank()) +
  theme(axis.line = element_line(color = 'black'))+ geom_text(x = 4.5, y = 15, label = lm_eqn(LO2spkindf), parse = TRUE) +
  labs(x=expression(paste("Log"[2]," concentration", sep=""))) +
  labs(y=expression(paste("Log"[2]," FPKM", sep="")))

#LO3
spikeinPlotLO3<-ggplot(LO3spkindf, aes(x=conc,y=FPKM)) +geom_point()+ggtitle("LO3/spike-in 2")+
  theme(panel.background =element_blank(), panel.grid.major = element_blank(), panel.grid.minor =element_blank(), panel.border =element_blank()) +
  theme(axis.line = element_line(color = 'black')) +
  geom_text(x = 4.5, y = 15, label = lm_eqn(LO3spkindf), parse = TRUE) +
  geom_smooth(method=lm, se=FALSE, color="black") +
  labs(x=expression(paste("Log"[2]," concentration", sep=""))) +
  labs(y=expression(paste("Log"[2]," FPKM", sep="")))

# LO4
spikeinPlotLO4<-ggplot(LO4spkindf, aes(x=conc,y=FPKM)) +geom_point()+ggtitle("LO4/spike-in 2")+
  theme(panel.background =element_blank(), panel.grid.major = element_blank(), panel.grid.minor =element_blank(), panel.border =element_blank()) +
  theme(axis.line = element_line(color = 'black')) +
  geom_text(x = 4.5, y = 15, label = lm_eqn(LO4spkindf), parse = TRUE) +
  geom_smooth(method=lm, se=FALSE, color="black") +
  labs(x=expression(paste("Log"[2]," concentration", sep=""))) +
  labs(y=expression(paste("Log"[2]," FPKM", sep="")))

#LO5
spikeinPlotLO5<-ggplot(LO5spkindf, aes(x=conc,y=FPKM)) +geom_point()+ggtitle("LO5/spike-in 1")+
  theme(panel.background =element_blank(), panel.grid.major = element_blank(), panel.grid.minor =element_blank(), panel.border =element_blank()) +
  theme(axis.line = element_line(color = 'black')) +
  geom_text(x = 4.5, y = 15, label = lm_eqn(LO5spkindf), parse = TRUE) +
  geom_smooth(method=lm, se=FALSE, color="black") +
  labs(x=expression(paste("Log"[2]," concentration", sep=""))) +
  labs(y=expression(paste("Log"[2]," FPKM", sep="")))

#B29
spikeinPlotB29<-ggplot(B29spkindf, aes(x=conc,y=FPKM)) +geom_point()+ggtitle("B29/spike-in 1")+
  theme(panel.background =element_blank(), panel.grid.major = element_blank(), panel.grid.minor =element_blank(), panel.border =element_blank()) +
  theme(axis.line = element_line(color = 'black')) +
  geom_text(x = 5, y = 11, label = lm_eqn(B29spkindf), parse = TRUE)+
  geom_smooth(method=lm, se=FALSE, color="black") +
  labs(x=expression(paste("Log"[2]," concentration", sep=""))) +
  labs(y=expression(paste("Log"[2]," FPKM", sep="")))

#B30
spikeinPlotB30<-ggplot(B30spkindf, aes(x=conc,y=FPKM)) +geom_point()+ggtitle("B30/spike-in 2")+
  theme(panel.background =element_blank(), panel.grid.major = element_blank(), panel.grid.minor =element_blank(), panel.border =element_blank()) +
  theme(axis.line = element_line(color = 'black')) +
  geom_text(x = 5, y = 11, label = lm_eqn(B30spkindf), parse = TRUE)+
  geom_smooth(method=lm, se=FALSE, color="black") +
  labs(x=expression(paste("Log"[2]," concentration", sep=""))) +
  labs(y=expression(paste("Log"[2]," FPKM", sep="")))

#B31
spikeinPlotB31<-ggplot(B31spkindf, aes(x=conc,y=FPKM)) +geom_point()+ggtitle("B31/spike-in 2")+
  theme(panel.background =element_blank(), panel.grid.major = element_blank(), panel.grid.minor =element_blank(), panel.border =element_blank()) +
  theme(axis.line = element_line(color = 'black')) +
  geom_text(x = 5, y = 11, label = lm_eqn(B31spkindf), parse = TRUE)+
  geom_smooth(method=lm, se=FALSE, color="black") +
  labs(x=expression(paste("Log"[2]," concentration", sep=""))) +
  labs(y=expression(paste("Log"[2]," FPKM", sep="")))

#B32
spikeinPlotB32<-ggplot(B32spkindf, aes(x=conc,y=FPKM)) +geom_point()+ggtitle("B32/spike-in 1")+
  theme(panel.background =element_blank(), panel.grid.major = element_blank(), panel.grid.minor =element_blank(), panel.border =element_blank()) +
  theme(axis.line = element_line(color = 'black')) +
  geom_text(x = 5, y = 11, label = lm_eqn(B32spkindf), parse = TRUE)+
  geom_smooth(method=lm, se=FALSE, color="black") +
  labs(x=expression(paste("Log"[2]," concentration", sep=""))) +
  labs(y=expression(paste("Log"[2]," FPKM", sep="")))

# put together concentration v. FPKM plots
library("gridExtra", lib.loc="~/Desktop/R-3.2.2/library")
grid.arrange(spikeinPlotLO2, spikeinPlotLO3, spikeinPlotLO4, spikeinPlotLO5, ncol=2, nrow=2)
grid.arrange(spikeinPlotB29, spikeinPlotB30, spikeinPlotB31, spikeinPlotB32, ncol=2, nrow=2)



#remove ERCC genes from dataframe Genes
Genes<-Genes[93:4061,]

View (Genes)

## in case you want to write data frames to csv to use in other programs
write.csv(Genes, file="/Users/brisbin/desktop/Siphamia/siph.genes.results/Genes.csv")
write.csv(FPKM, file="/Users/brisbin/desktop/Siphamia/siph.genes.results/FPKM.csv")

#Load edgeR
library("edgeR", lib.loc="~/Desktop/R-3.2.2/library")
#set working directory
setwd("/Users/brisbin/desktop/Siphamia/edgeR_wd")

#read in data
counts <- Genes[ , -c(1,ncol(Genes)) ]
rownames(counts)<-Genes[ , 1] #gene names
colnames(counts)<-paste(c("B29","B30","B31","B32","LO2","LO3","LO4", "LO5"))
View(counts)

##Basic info about the data
#size of data frame
dim(counts)
#column sums - library size
colSums(counts)
#library size in millions
colSums(counts)/1e06
# Number of genes with low counts
table( rowSums( counts ) )[ 1:30 ] 

#convert count matrix to edgeR Object
#create a group variable that tells edgeR which samples belong to which group
group<- c(rep("B",4), rep("LO",4))
cds<-DGEList(counts, group=group)
names(cds)
#original count data
head(cds$counts)
#contains sample information
head(cds$samples)
# How many genes have 0 counts across all samples
sum( cds$all.zeros )


##need to filter out low count reads bc impossible to detect differential expression
#keep only genes with atleast 1 read per million reads in at least 3 samples
#once this is done, can calculate normalisation factors which correct for different compositions of samples
# effective library size = product of actual library size and these factors

cds <- cds[rowSums(1e+06 * cds$counts/expandAsMatrix(cds$samples$lib.size, dim(cds)) > 1) >= 3, ]
dim( cds )
cds <- calcNormFactors( cds )
cds$samples

# effective library sizes
cds$samples$lib.size * cds$samples$norm.factors
#Normalize gene counts to library size --- done?

##from edger docs: an MD plot can show the performance of the TMM normalization, visualizes the library size-adjusted log-fold change between two
#libraries (the difference) against the average log-expression across those libraries (the mean)
#The following MD plot is generated by comparing sample 1 against an artificial library
#constructed from the average of all other samples.
plotMD(cpm(cds, log=TRUE), column=8)
abline(h=0, col="red", lty=2, lwd=2)
#Ideally, the bulk of genes should be centred at a log-fold change of zero. This indicates
#that any composition bias between libraries has been successfully removed. This quality
#check should be repeated by constructing a MD plot for each sample.




#Multi-Dimensional Scaling Plot
#measures similarity of samples and projects measure into 2 dimensions
plotMDS( cds , main = "MDS Plot for Bacteria Count Data", labels = colnames( cds$counts ) )

#Estimating Dispersions
#1st: calculate common dispersion
#each gene get assigned the same dispersion estimate
#output of the estimation includes estimate andother elements added to object cds
cds <- estimateCommonDisp( cds )
names( cds )
#The estimate
cds$common.dispersion


#with common dispersion, can estimate tagwise dispersions
#each gene will get its own dispersion estimate
#tagwise dispersions are squeezed towards a common value
#amt of squeezing is governed by parameter prior.n
#higher prior.n --> closer the estimates will be to common dispersion
#default value is nearest interger to 50/(#samples-#groups)
# 50/(6-2)= 12.5 -->13

cds <- estimateTagwiseDisp( cds)
names( cds )
summary( cds$tagwise.dispersion )

#Mean-Variance Plot
#see how well the dispersion factors fit the data
# see raw variance of counts (grey dots), 
#the variance using tagwise dispersions(light blue dots)
#variances using common dispersion (solid blue line)
#variance=mean aka poisson variance (solid black line)

meanVarPlot <- plotMeanVar( cds , show.raw.vars=TRUE ,
                            show.tagwise.vars=TRUE ,
                            show.binned.common.disp.vars=FALSE ,
                            show.ave.raw.vars=FALSE ,
                            dispersion.method = "qcml" , NBline = TRUE ,
                            nbins = 100 ,
                            pch = 16 ,
                            xlab ="Mean Expression (Log10 Scale)" ,
                            ylab = "Variance (Log10 Scale)" ,
                            main = "Bacteria Mean-Variance Plot" )

#Testing
# exactTest() performs pairwise tests for diff exp between 2 groups
# pair indicates which groups should be compared
# output is a list of elements, one of which is a table of results
# de.poi sets dispersion parameter to zero - poisson model - can compare negative binomial results to this

de.cmn <- exactTest( cds , dispersion="common", pair = c( "B" , "LO" ) )
de.tgw <- exactTest( cds , dispersion="tagwise", pair = c( "B" , "LO" ) )
de.poi <- exactTest( cds , dispersion = 1e-06 , pair = c( "B" , "LO" ) )

names( de.tgw )
de.tgw$comparison # which groups have been compared
head( de.tgw$table ) # the results table is in same order of count matrix.

# this does not contain p-values adjusted for multiple testing
# topTags() takes output from exactTest(), adjusts raw p-values using False Discover Rate (FDR) correction
# returns top differentially expressed genes
# sort.by argument allws sorting by p-value, concentration, or fold change
# topTags() returns original counts of the top differentially expressed genes
# set n parameter to number of genes, save the entire topTags() results table

# Top tags for tagwise analysis
options( digits = 3 ) # print only 3 digits
topTags( de.tgw , n = 20 , sort.by = "p.value" ) # top 20 DE genes
# Back to count matrix for tagwise analysis
cds$counts[ rownames( topTags( de.tgw , n = 15 )$table ) , ]
# Sort tagwise results by Fold-Change instead of p-value
resultsByFC.tgw <- topTags( de.tgw , n = nrow( de.tgw$table ) , sort.by = "logFC" )$table
head( resultsByFC.tgw )
# Store full topTags results table
resultsTbl.cmn <- topTags( de.cmn , n = nrow( de.cmn$table ) )$table
resultsTbl.tgw <- topTags( de.tgw , n = nrow( de.tgw$table ) )$table
resultsTbl.poi <- topTags( de.poi , n = nrow( de.poi$table ) )$table
head( resultsTbl.tgw )

#compare p-values to significance level and determine # of diff exp genes
# sig level = ).05
# decideTestsDGE() to determine how many diff exp genes are up or down regulated compared to control

# Names/IDs of DE genes
de.genes.cmn <- rownames( resultsTbl.cmn )[ resultsTbl.cmn$PValue <= 0.05 ]
de.genes.tgw <- rownames( resultsTbl.tgw )[ resultsTbl.tgw$PValue <= 0.05 ]
de.genes.poi <- rownames( resultsTbl.poi )[ resultsTbl.poi$PValue <= 0.05 ]
# Amount significant
length( de.genes.cmn )
length( de.genes.tgw )
length( de.genes.poi )
# Percentage of total genes
length( de.genes.cmn ) / nrow( resultsTbl.cmn ) * 100
length( de.genes.tgw ) / nrow( resultsTbl.tgw ) * 100
length( de.genes.poi ) / nrow( resultsTbl.poi ) * 100
# Up/Down regulated summary for tagwise results
summary( decideTestsDGE( de.tgw , p.value = 0.05 ) ) # the adjusted p-values are used here



#Visualize results
#spread of expression levels for DE genes
#plot concentrations for top 100 DE genes for each analysis

#MA plot shows relationship between concentration and fold-change across genes
# DE genes are red
# non DE genes are black
# organge dots = counts were zero in all samples in one group
# blue line is at log-FC of 2 to represent a level for biological significance
##An MA plot is an application of a Bland???Altman plot for visual representation of two channel DNA microarray gene expression data which has been transformed onto the M (log ratios) and A (mean average) scale.
par( mfrow=c(2,1) )
plotSmear( cds , de.tags=de.genes.poi , main="Bacteria Poisson MA plot" ,
           pair = c("B","LO") ,
           cex = .35 ,
           xlab="Log Concentration" , ylab="Log Fold-Change" )
abline( h = c(-2, 2) , col = "dodgerblue" )
plotSmear( cds , de.tags=de.genes.tgw , main="Bacteria Tagwise MA plot" ,
           pair = c("B","LO") ,
           cex = .35 ,
           xlab="Log Concentration" , ylab="Log Fold-Change" )
abline( h=c(-2,2) , col="dodgerblue" )
par( mfrow=c(1,1) )

# log difference between DE genes under negative binomial/poisson model and tagwise in MA plot
# plot top 500 DE genes, rest are black

par( mfrow = c(2,1) )
plotSmear( cds , de.tags=de.genes.poi[1:500] , main=" Bacteria Poisson MA plot" ,
           pair=c("B","LO") ,
           cex=.5 ,
           xlab="Log Concentration" , ylab="Log Fold-Change" )
abline( h=c(-2,2) , col="dodgerblue" )
plotSmear( cds , de.tags=de.genes.tgw[1:500] , main="Bacteria Tagwise MA plot" ,
           pair=c("B","LO") ,
           cex = .5 ,
           xlab="Log Concentration" , ylab="Log Fold-Change" )
abline( h=c(-2,2) , col="dodgerblue" )
par( mfrow=c(1,1) )
