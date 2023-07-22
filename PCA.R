library(mixOmics) 

fecal <- read.csv("fecalMB2.csv", row.names=1, check.names=T, header=T)
serum <- read.csv("serum.csv", row.names=1, check.names=T, header=T)
set1n <- read.csv("set1n.csv", row.names=1, check.names=T, header=T)
stratify <- read.csv("stratefy.csv", row.names=1, check.names=T, header=T)
unstratify <- read.csv("unstratify.csv", row.names=1, check.names=T, header=T) #raw_unstratify <- read.csv("rawl.csv", row.names=1, check.names=T, header=T)
dim(stratify); dim(unstratify); dim(fecal)
mapping = read.table('../Mapping.txt', header=T, sep='\t')
row.names(mapping) <- mapping$sample
group=mapping$group

lumen=t(fecal) #as.matrix convert the data into matrix format before performing CCA
microb = t(set1n)
strat=t(stratify)
unstrat=t(unstratify)
sera=t(serum)
#write.csv(file="set1n.csv",otu_table(set1n))
dim(lumen);dim(microb);dim(strat); dim(unstrat);dim(sera)

data=strat
##PCA
MyResult.pca <- pca(data, ncomp=10, center = TRUE, scale = TRUE)  # 1 Run the method
MyResult.pca
pdf(file="PC_strat.pdf")
plot(MyResult.pca) # Amount of variance explained and choice of number of components => choose the number of components
explainedVariance <- tune.pca(data, ncomp = 10, center = TRUE, scale = FALSE) 
plot(explainedVariance) # plot these values
dev.off()

MyResult.pca <- pca(data, ncomp=2) 
pdf(file="PCA_strat.pdf")
plotIndiv(MyResult.pca,ncomp=2, ind.names = FALSE, group=group, legend=TRUE, title = "PCA", legend.title = "Group", scale=TRUE) # 2 Plot the samples
dev.off()

selectVar(MyResult.pca, comp=2)$value[1:10,] #extracted coefficient values of variables. loading weights. A large absolute value indicate the importance of the variable in this PC.
pdf(file="PCA_PC1loading_strat.pdf")
plotLoadings(MyResult.pca, group=group, comp = 1, ndisplay = 10, size.name = rel(0.8), contrib = 'max', method = 'mean')
dev.off()
pdf(file="PCA_PC2loading_strat.pdf")
plotLoadings(MyResult.pca, group=group, comp = 2, ndisplay = 10, size.name = rel(0.8), contrib = 'max', method = 'mean')
dev.off()
plotVar(MyResult.pca,cex=3,cutoff = 0.8) 
final.pca <- pca(data, ncomp = 2, scale = TRUE, center = TRUE)
biplot(final.pca, cex = 0.7, 
       xlabs = paste(group, 1:nrow(data)), #simplify names
       group = group,  # colour by sample class
       title = 'sPCA comp 1 - 2')
background <- background.predict(final.pca, comp.predicted=2,
                                 dist = "max.dist") 
plotIndiv(MyResult.spca, comp = 1:2, group = group,
          ind.names = FALSE, title = "Maximum distance",
          legend = TRUE,  background = background3)

MyResult.spca <- spca(data, ncomp=8, center = TRUE, scale = TRUE)  # 1 Run the method
MyResult.spca
pdf(file="SPC_strat.pdf")
plot(MyResult.spca) # Amount of variance explained and choice of number of components => choose the number of components
dev.off()
explainedVariance <- tune.pca(data, ncomp = 10, center = TRUE, scale = FALSE) 
plot(explainedVariance) # plot these values

MyResult.spca <- spca(data, ncomp=2)
#MyResult.spca <- spca(data, ncomp=2, keepX=c(25,25))
#pdf(file="sPCA_strat.pdf")
plotIndiv(MyResult.spca, ind.names = FALSE, ncomp=2, group=group, legend=TRUE, title = "Sparse PCA", legend.title = "Group") # 2 Plot the samples
#dev.off()

final.spca <- spca(data, ncomp = 2, keepX = c(50,50))
biplot(final.spca, cex = 0.7, 
       xlabs = paste(group, 1:nrow(data)), #simplify names
       group = group,  # colour by sample class
       title = 'sPCA comp 1 - 2')

data=unstrat
##PCA
MyResult.pca <- pca(data, ncomp=10, center = TRUE, scale = TRUE)  # 1 Run the method
MyResult.pca
pdf(file="PC_unstrat.pdf")
plot(MyResult.pca) # Amount of variance explained and choice of number of components => choose the number of components
explainedVariance <- tune.pca(data, ncomp = 10, center = TRUE, scale = FALSE) 
plot(explainedVariance) # plot these values
dev.off()

MyResult.pca <- pca(data, ncomp=2) 
pdf(file="PCA_unstrat.pdf")
plotIndiv(MyResult.pca,ncomp=2, ind.names = FALSE, group=group, legend=TRUE, title = "PCA", legend.title = "Group", scale=TRUE) # 2 Plot the samples
dev.off()

selectVar(MyResult.pca, comp=2)$value[1:10,] #extracted coefficient values of variables. loading weights. A large absolute value indicate the importance of the variable in this PC.
pdf(file="PCA_PC1loading_unstrat.pdf")
plotLoadings(MyResult.pca, group=group, comp = 1, ndisplay = 20, size.name = rel(0.8), contrib = 'max', method = 'mean')
dev.off()
pdf(file="PCA_PC2loading_unstrat.pdf")
plotLoadings(MyResult.pca, group=group, comp = 2, ndisplay = 20, size.name = rel(0.8), contrib = 'max', method = 'mean')
dev.off()
plotVar(MyResult.pca,cex=3,cutoff = 0.8) 
#final.pca <- pca(data, ncomp = 2, scale = TRUE, center = TRUE)

MyResult.spca <- spca(data, ncomp=8, center = TRUE, scale = TRUE)  # 1 Run the method
MyResult.spca
pdf(file="sPC_unstrat.pdf")
plot(MyResult.spca) # Amount of variance explained and choice of number of components => choose the number of components
explainedVariance <- tune.spca(data, ncomp = 8, center = TRUE, scale = FALSE) 
plot(explainedVariance) # plot these values
dev.off()

MyResult.spca <- spca(data, ncomp=4)
pdf(file="sPCA_unstrat.pdf")
plotIndiv(MyResult.spca, ind.names = FALSE, ncomp=2, group=group, legend=TRUE, title = "Sparse PCA", legend.title = "Group") # 2 Plot the samples
dev.off()

final.spca <- spca(data, ncomp = 2, keepX = c(50,50))
biplot(final.spca, cex = 0.7, 
       xlabs = paste(group, 1:nrow(data)), #simplify names
       group = group,  # colour by sample class
       title = 'sPCA comp 1 - 2')

data=lumen
name="lumen"
##PCA
MyResult.pca <- pca(data, ncomp=10, center = TRUE, scale = TRUE)  # 1 Run the method
MyResult.pca
pdf(file="PC_lumen.pdf") #paste("./PLSDA/","weight",".pdf",sep="")
plot(MyResult.pca) # Amount of variance explained and choice of number of components => choose the number of components
explainedVariance <- tune.pca(data, ncomp = 10, center = TRUE, scale = FALSE) 
plot(explainedVariance) # plot these values
dev.off()

MyResult.pca <- pca(data, ncomp=4) 
pdf(file="PCA_lumen.pdf")
plotIndiv(MyResult.pca,ncomp=4, ind.names = FALSE, group=group, legend=TRUE, title = "PCA", legend.title = "Group", scale=TRUE) # 2 Plot the samples
dev.off()

selectVar(MyResult.pca, comp=1)$value[1:30,] #extracted coefficient values of variables. loading weights. A large absolute value indicate the importance of the variable in this PC.
selectVar(MyResult.pca, comp=2)$value[1:21,] #extracted coefficient values of variables. loading weights. A large absolute value indicate the importance of the variable in this PC.
selectVar(MyResult.pca, comp=3)$value[1:27,] #extracted coefficient values of variables. loading weights. A large absolute value indicate the importance of the variable in this PC.
selectVar(MyResult.pca, comp=4)$value[1:33,] #extracted coefficient values of variables. loading weights. A large absolute value indicate the importance of the variable in this PC.
pdf(file="PCA_PC1loading_lumen.pdf")
plotLoadings(MyResult.pca, group=group, comp = 1, ndisplay = 30, size.name = rel(0.8), contrib = 'max', method = 'mean')
dev.off()
pdf(file="PCA_PC2loading_lumen.pdf")
plotLoadings(MyResult.pca, group=group, comp = 2, ndisplay = 21, size.name = rel(0.8), contrib = 'max', method = 'mean')
dev.off()
pdf(file="PCA_PC3loading_lumen.pdf")
plotLoadings(MyResult.pca, group=group, comp = 1, ndisplay = 27, size.name = rel(0.8), contrib = 'max', method = 'mean')
dev.off()
pdf(file="PCA_PC4loading_lumen.pdf")
plotLoadings(MyResult.pca, group=group, comp = 2, ndisplay = 33, size.name = rel(0.8), contrib = 'max', method = 'mean')
dev.off()
pdf(file="PCA_variable_lumen.pdf")
plotVar(MyResult.pca,cex=3,cutoff = 0.8) 
dev.off()
MyResult.spca <- spca(data, ncomp=8, center = TRUE, scale = TRUE)  # 1 Run the method
MyResult.spca
pdf(file="sPC_lumen.pdf")
plot(MyResult.spca) # Amount of variance explained and choice of number of components => choose the number of components
explainedVariance <- tune.spca(data, ncomp = 10, center = TRUE, scale = FALSE) 
plot(explainedVariance) # plot these values
dev.off()
MyResult.spca <- spca(data, ncomp=4)
pdf(file="sPCA_lumen.pdf")
plotIndiv(MyResult.spca, ind.names = FALSE, ncomp=4, group=group, legend=TRUE, title = "Sparse PCA", legend.title = "Group") # 2 Plot the samples
dev.off()

set.seed(333) # for reproducibility with this case study, remove otherwise
test.keepX <- c(seq(5, 25, 5)) # set the number of variable values to be tested
tune.spca.res <- tune.spca(data, ncomp = 4, # generate the first 3 components
                           nrepeat = 5, # repeat the cross-validation 5 times
                           folds = 3, # use 3 folds for the cross-validation
                           test.keepX = test.keepX)
plot(tune.spca.res) # plot the optimisation output
tune.spca.res$choice.keepX # how many variables per component is optimal
final.spca <- spca(data, ncomp = 4, # based off figure 1, three components is best
                   keepX = tune.spca.res$choice.keepX)
pdf(file="biPlot_lumen.pdf")
biplot(final.spca, cex = 0.7, 
       xlabs = paste(group, 1:nrow(data)), #simplify names
       group = group,  # colour by sample class
       title = 'sPCA comp 1 - 2')
dev.off()

data=microb
##PCA
MyResult.pca <- pca(data, ncomp=10, center = TRUE, scale = TRUE)  # 1 Run the method
MyResult.pca
pdf(file="PC_microb.pdf") #paste("./PLSDA/","weight",".pdf",sep="")
plot(MyResult.pca) # Amount of variance explained and choice of number of components => choose the number of components
explainedVariance <- tune.pca(data, ncomp = 10, center = TRUE, scale = FALSE) 
plot(explainedVariance) # plot these values
dev.off()

MyResult.pca <- pca(data, ncomp=2) 
pdf(file="PCA_microb.pdf")
plotIndiv(MyResult.pca,ncomp=2, ind.names = FALSE, group=group, legend=TRUE, title = "PCA", legend.title = "Group", scale=TRUE) # 2 Plot the samples
dev.off()

asvs=selectVar(MyResult.pca, comp=1)$value
write.csv(file="asvs_comp1.csv",asvs)

selectVar(MyResult.pca, comp=1)$value[1:10,] #extracted coefficient values of variables. loading weights. A large absolute value indicate the importance of the variable in this PC.
selectVar(MyResult.pca, comp=2)$value[1:10,] #extracted coefficient values of variables. loading weights. A large absolute value indicate the importance of the variable in this PC.
pdf(file="PCA_PC1loading_microb.pdf")
plotLoadings(MyResult.pca, group=group, comp = 1, ndisplay = 20, size.name = rel(0.8), contrib = 'max', method = 'mean')
dev.off()
pdf(file="PCA_PC2loading_microb.pdf")
plotLoadings(MyResult.pca, group=group, comp = 2, ndisplay = 10, size.name = rel(0.8), contrib = 'max', method = 'mean')
dev.off()
pdf(file="PCA_variable_microb.pdf")
plotVar(MyResult.pca,cex=3,cutoff = 0.8) 
dev.off()
#final.pca <- pca(data, ncomp = 2, scale = TRUE, center = TRUE)
background <- background.predict(MyResult.pca, comp.predicted=2,
                                 dist = "max.dist") 
plotIndiv(MyResult.pca, comp = 1:2, group = group,
          ind.names = FALSE, title = "Maximum distance",
          legend = TRUE,  background = background3)

MyResult.spca <- spca(data, ncomp=8, center = TRUE, scale = TRUE)  # 1 Run the method
MyResult.spca
pdf(file="sPC_microb.pdf")
plot(MyResult.spca) # Amount of variance explained and choice of number of components => choose the number of components
explainedVariance <- tune.spca(data, ncomp = 8, center = TRUE, scale = FALSE) 
plot(explainedVariance) # plot these values
dev.off()

MyResult.spca <- spca(data, ncomp=2)
pdf(file="sPCA_microb.pdf")
plotIndiv(MyResult.spca, ind.names = FALSE, ncomp=2, group=group, legend=TRUE, title = "Sparse PCA", legend.title = "Group") # 2 Plot the samples
dev.off()

set.seed(333) # for reproducibility with this case study, remove otherwise
test.keepX <- c(seq(5, 25, 5)) # set the number of variable values to be tested

tune.spca.res <- tune.spca(data, ncomp = 2, # generate the first 3 components
                           nrepeat = 5, # repeat the cross-validation 5 times
                           folds = 3, # use 3 folds for the cross-validation
                           test.keepX = test.keepX)
plot(tune.spca.res) # plot the optimisation output
tune.spca.res$choice.keepX # how many variables per component is optimal
final.spca <- spca(data, ncomp = 2, # based off figure 1, three components is best
                   keepX = tune.spca.res$choice.keepX)
final.spca <- spca(data, ncomp = 2, # based off figure 1, three components is best
                   keepX = c(10,10))
pdf(file="biPlot_microbe.pdf")
biplot(final.spca, cex = 0.7, 
       xlabs = paste(group, 1:nrow(data)), #simplify names
       group = group,  # colour by sample class
       title = 'sPCA comp 1 - 2')
dev.off()

data=sera
##PCA
MyResult.pca <- pca(data, ncomp=10, center = TRUE, scale = TRUE)  # 1 Run the method
MyResult.pca
pdf(file="PC_serum.pdf") #paste("./PLSDA/","weight",".pdf",sep="")
plot(MyResult.pca) # Amount of variance explained and choice of number of components => choose the number of components
explainedVariance <- tune.pca(data, ncomp = 10, center = TRUE, scale = FALSE) 
plot(explainedVariance) # plot these values
dev.off()

MyResult.pca <- pca(data, ncomp=4) 
pdf(file="PCA_serum.pdf")
plotIndiv(MyResult.pca,ncomp=2, ind.names = FALSE, group=group, legend=TRUE, title = "PCA", legend.title = "Group", scale=TRUE) # 2 Plot the samples
dev.off()

selectVar(MyResult.pca, comp=1)$value[1:21,] #extracted coefficient values of variables. loading weights. A large absolute value indicate the importance of the variable in this PC.
selectVar(MyResult.pca, comp=2)$value[1:21,] #extracted coefficient values of variables. loading weights. A large absolute value indicate the importance of the variable in this PC.
selectVar(MyResult.pca, comp=3)$value[1:34,] #extracted coefficient values of variables. loading weights. A large absolute value indicate the importance of the variable in this PC.
selectVar(MyResult.pca, comp=4)$value[1:26,] #extracted coefficient values of variables. loading weights. A large absolute value indicate the importance of the variable in this PC.

pdf(file="PCA_PC1loading_serum.pdf")
plotLoadings(MyResult.pca, group=group, comp = 1, ndisplay = 21, size.name = rel(0.8), contrib = 'max', method = 'mean')
dev.off()
pdf(file="PCA_PC2loading_serum.pdf")
plotLoadings(MyResult.pca, group=group, comp = 2, ndisplay = 21, size.name = rel(0.8), contrib = 'max', method = 'mean')
dev.off()
pdf(file="PCA_variable_serum.pdf")
plotVar(MyResult.pca,cex=3,cutoff = 0.8) 
dev.off()

MyResult.spca <- spca(data, ncomp=8, center = TRUE, scale = TRUE)  # 1 Run the method
MyResult.spca
pdf(file="sPC_serum.pdf")
plot(MyResult.spca) # Amount of variance explained and choice of number of components => choose the number of components
explainedVariance <- tune.pca(data, ncomp = 10, center = TRUE, scale = FALSE) 
plot(explainedVariance) # plot these values
dev.off()

MyResult.spca <- spca(data, ncomp=4)

pdf(file="sPCA_serum.pdf")
plotIndiv(MyResult.spca, ind.names = FALSE, ncomp=4, group=group, legend=TRUE, title = "Sparse PCA", legend.title = "Group") # 2 Plot the samples
dev.off()

set.seed(333) # for reproducibility with this case study, remove otherwise
test.keepX <- c(seq(5, 25, 5)) # set the number of variable values to be tested

tune.spca.res <- tune.spca(data, ncomp = 4, # generate the first 3 components
                           nrepeat = 5, # repeat the cross-validation 5 times
                           folds = 3, # use 3 folds for the cross-validation
                           test.keepX = test.keepX)
plot(tune.spca.res) # plot the optimisation output
tune.spca.res$choice.keepX # how many variables per component is optimal
final.spca <- spca(data, ncomp = 4, # based off figure 1, three components is best
                   keepX = tune.spca.res$choice.keepX)
pdf(file="biPlot_serum.pdf")
biplot(final.spca, cex = 0.7, 
       xlabs = paste(group, 1:nrow(data)), #simplify names
       group = group,  # colour by sample class
       title = 'sPCA comp 1 - 2')
dev.off()


