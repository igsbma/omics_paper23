library(mixOmics); packageVersion('mixOmics')

fecal <- read.csv("fecalMB2.csv", row.names=1, check.names=T, header=T)
serum <- read.csv("serum2.csv", row.names=1, check.names=F, header=T)
genus <- read.csv("genus.csv", row.names=1, check.names=T, header=T)
family <- read.csv("family.csv", row.names=1, check.names=T, header=T)
set1n <- read.csv("set1n2.csv", row.names=1, check.names=T, header=T)
stratify <- read.csv("stratefy.csv", row.names=1, check.names=T, header=T)
unstratify <- read.csv("unstratify.csv", row.names=1, check.names=T, header=T) #raw_unstratify <- read.csv("rawl.csv", row.names=1, check.names=T, header=T)
dim(stratify); dim(unstratify); dim(fecal)
mapping = read.table('Mapping.txt', header=T, sep='\t')
row.names(mapping) <- mapping$sample
group=mapping$group

lumen=t(fecal) #as.matrix convert the data into matrix format before performing CCA
microb = t(set1n)
strat=t(stratify)
unstrat=t(unstratify)
sera=t(serum)
genera=t(genus)
families=t(family)
dim(lumen);dim(microb);dim(strat); dim(unstrat);dim(sera);dim(genera); dim(families)
stim.col <- c("darkblue", "purple", "green4","red3")

#sparse PLSDA
#lumen
data=lumen
file <- paste("lumen",".pdf",sep="")
pdf(file, width=11, height=8.5)
lumen.splsda <- splsda(data, group, ncomp = 3)  # 1 Run the method
plotIndiv(lumen.splsda, ncomp = 3, ind.names = FALSE, ellipse=TRUE, group=group, star=TRUE, legend=TRUE, title = "sPLS-DA-all", legend.title = "Group", X.label = 'PLS-DA 1', Y.label = 'PLS-DA 2') # 2 Plot the samples
cim <- cim(lumen.splsda, clust.method = c("ward", "ward"), comp = c(1,2), row.sideColors = stim.col[factor(group)], legend=list(legend = unique(group), col=stim.col, title = "Group", cex=0.9))
background <- background.predict(lumen.splsda, comp.predicted=2, dist = "max.dist") 
plotIndiv(lumen.splsda, comp = 1:2, group = group,
          ind.names = FALSE, title = "Maximum distance-all",
          legend = TRUE,  background = background)
auc.plsda <- auroc(lumen.splsda)
plotLoadings(lumen.splsda, contrib = 'max', method = 'mean',legend.title = "group", comp = 1, ndisplay = 50) #70
plotLoadings(lumen.splsda, contrib = 'max', method = 'mean',legend.title = "group", comp = 2, ndisplay = 50)
plotLoadings(lumen.splsda, contrib = 'max', method = 'mean',legend.title = "group", comp = 3, ndisplay = 30)

#selected
lumen.splsda2 <- splsda(data, group, keepX=c(50,40,30), ncomp=3)  # 1 Run the method
plotIndiv(lumen.splsda2, ind.names = FALSE, ellipse=TRUE, group=group, star=TRUE, legend=TRUE, title = "sPLS-DA", legend.title = "Group", X.label = 'PLS-DA 1', Y.label = 'PLS-DA 2') # 2 Plot the samples
background2 <- background.predict(lumen.splsda2, comp.predicted=2,
                                  dist = "max.dist") 
plotIndiv(lumen.splsda2, comp = 1:2, group = group,
          ind.names = FALSE, title = "Maximum distance",
          legend = TRUE,  background = background2)
cim <- cim(lumen.splsda2,clust.method = c("ward", "ward"), comp = c(1,2),scale=TRUE, row.sideColors = stim.col[factor(group)], legend=list(legend = unique(group), col=stim.col, title = "Group", cex=0.8),margins = c(8, 8))
auc.plsda2 <- auroc(lumen.splsda2)
list.keepX <- c(5:10,  seq(20, 100, 10))
list.keepX # to output the grid of values tested
set.seed(30) # for reproducbility in this vignette, otherwise increase nrepeat
tune.splsda.lumen <- tune.splsda(data, group, ncomp = 4, # we suggest to push ncomp a bit more, e.g. 4
                                 validation = 'Mfold',
                                 folds = 5, dist = 'max.dist', progressBar = FALSE,
                                 measure = "BER", test.keepX = list.keepX,
                                 nrepeat = 50)   # we suggest nrepeat = 50
error <- tune.splsda.lumen$error.rate
ncomp <- tune.splsda.lumen$choice.ncomp$ncomp # optimal number of components based on t-tests on the error rate
ncomp
select.keepX <- tune.splsda.lumen$choice.keepX[1:ncomp]  # optimal number of variables to select
select.keepX
plot(tune.splsda.lumen)   #, col = color.jet(ncomp))
lumen.splsda.final <- splsda(data, group, ncomp = ncomp, keepX = select.keepX)
legend=list(legend = levels(group), # set of classescol = unique(color.mixo(group)), # set of colours
            title = "Group", # legend title
            cex = 0.7) # legend size
plotIndiv(lumen.splsda.final, ind.names = FALSE, ellipse=TRUE, group=group, star=TRUE, legend=TRUE, title = "sPLS-DA final result", legend.title = "Group", X.label = 'PLS-DA 1', Y.label = 'PLS-DA 2') # 2 Plot the samples
cim <- cim(lumen.splsda.final, clust.method = c("ward", "ward"), comp = c(1,2,3), row.sideColors = stim.col[factor(group)], legend=list(legend = unique(group), col=stim.col, title = "Group", cex=0.9),margins = c(8, 8))
background3 <- background.predict(lumen.splsda.final, comp.predicted=2,
                                  dist = "max.dist") 
plotIndiv(lumen.splsda.final, ncomp = 3, group = group,
          ind.names = FALSE, title = "Maximum distance",
          legend = TRUE,  background = background3)
plotLoadings(lumen.splsda.final, contrib = 'max', method = 'mean', legend.title = "group", comp = 1, ndisplay = 35)
plotLoadings(lumen.splsda.final, contrib = 'max', method = 'mean', legend.title = "group", comp = 2, ndisplay = 30)
plotLoadings(lumen.splsda.final, contrib = 'max', method = 'mean', legend.title = "group", comp = 3, ndisplay = 50)
# form new perf() object which utilises the final model
perf.splsda.lumen <- perf(lumen.splsda.final, 
                          folds = 5, nrepeat = 50, # use repeated cross-validation
                          validation = "Mfold", dist = "max.dist",  # use max.dist measure
                          progressBar = FALSE)
# plot the stability of each feature for the first three components, 'h' type refers to histogram
par(mfrow=c(1,3))
plot(perf.splsda.lumen$features$stable[[1]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features', 
     main = '(a) Comp 1', las =2)
plot(perf.splsda.lumen$features$stable[[2]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features', 
     main = '(b) Comp 2', las =2)
plot(perf.splsda.lumen$features$stable[[3]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features',
     main = '(c) Comp 3', las =2)
dev.off()

#microb
data=microb
file <- paste("microb2",".pdf",sep="")
pdf(file, width=11, height=8.5)
microb.splsda <- splsda(data, group, ncomp = 3)  # 1 Run the method
plotIndiv(microb.splsda, ncomp = 3, ind.names = FALSE, ellipse=TRUE, group=group, star=TRUE, legend=TRUE, title = "sPLS-DA-all", legend.title = "Group", X.label = 'PLS-DA 1', Y.label = 'PLS-DA 2') # 2 Plot the samples
cim <- cim(microb.splsda, clust.method = c("ward", "ward"), comp = c(1,2,3), row.sideColors = stim.col[factor(group)], legend=list(legend = unique(group), col=stim.col, title = "Group", cex=0.9),margins = c(8, 8))
background <- background.predict(microb.splsda, comp.predicted=2, dist = "max.dist") 
plotIndiv(microb.splsda, comp = 1:2, group = group,
          ind.names = FALSE, title = "Maximum distance-all",
          legend = TRUE,  background = background)
auc.plsda <- auroc(microb.splsda)
plotLoadings(microb.splsda, contrib = 'max', method = 'mean',legend.title = "group", comp = 1, ndisplay = 75)
plotLoadings(microb.splsda, contrib = 'max', method = 'mean',legend.title = "group", comp = 2, ndisplay = 60)
plotLoadings(microb.splsda, contrib = 'max', method = 'mean',legend.title = "group", comp = 3, ndisplay = 25)
#selected
microb.splsda2 <- splsda(data, group, keepX=c(50,40,25), ncomp=3)  # 1 Run the method
plotIndiv(microb.splsda2, ind.names = FALSE, ellipse=TRUE, group=group, star=TRUE, legend=TRUE, title = "sPLS-DA", legend.title = "Group", X.label = 'PLS-DA 1', Y.label = 'PLS-DA 2') # 2 Plot the samples
background2 <- background.predict(microb.splsda2, comp.predicted=2,
                                  dist = "max.dist") 
plotIndiv(microb.splsda2, comp = 1:2, group = group,
          ind.names = FALSE, title = "Maximum distance",
          legend = TRUE,  background = background2)
cim <- cim(microb.splsda2,clust.method = c("ward", "ward"), comp = c(1,2,3), row.sideColors = stim.col[factor(group)], legend=list(legend = unique(group), col=stim.col, title = "Group", cex=0.9),margins = c(8, 8))
auc.plsda2 <- auroc(microb.splsda2)
list.keepX <- c(5:10,  seq(20, 100, 10))
list.keepX # to output the grid of values tested
set.seed(30) # for reproducbility in this vignette, otherwise increase nrepeat

tune.splsda.microb <- tune.splsda(data, group, ncomp = 4, # we suggest to push ncomp a bit more, e.g. 4
                                 validation = 'Mfold',
                                 folds = 5, dist = 'max.dist', progressBar = FALSE,
                                 measure = "BER", test.keepX = list.keepX,
                                 nrepeat = 50)   # we suggest nrepeat = 50
error <- tune.splsda.microb$error.rate
ncomp <- tune.splsda.microb$choice.ncomp$ncomp # optimal number of components based on t-tests on the error rate
ncomp
select.keepX <- tune.splsda.microb$choice.keepX[1:ncomp]  # optimal number of variables to select
select.keepX
select.keepX <- c(30,25,20)
names(select.keepX) <- c("comp1","comp2","comp3")
plot(tune.splsda.microb)   #, col = color.jet(ncomp))
microb.splsda.final <- splsda(data, group, ncomp = ncomp, keepX = select.keepX)
legend=list(legend = levels(group), # set of classescol = unique(color.mixo(group)), # set of colours
            title = "Group", # legend title
            cex = 0.7) # legend size
plotIndiv(microb.splsda.final, ind.names = FALSE, ellipse=TRUE, group=group, star=TRUE, legend=TRUE, title = "sPLS-DA final result", legend.title = "Group", X.label = 'PLS-DA 1', Y.label = 'PLS-DA 2') # 2 Plot the samples
#plotIndiv(microb.splsda.final, comp = c(2,3), ind.names = FALSE, ellipse=TRUE, group=group, star=TRUE, legend=TRUE, title = "sPLS-DA final result", legend.title = "Group", X.label = 'PLS-DA 1', Y.label = 'PLS-DA 2') # 2 Plot the samples
cim <- cim(microb.splsda.final, clust.method = c("ward", "ward"), comp = c(1,2,3), row.sideColors = stim.col[factor(group)], legend=list(legend = unique(group), col=stim.col, title = "Group", cex=0.9),margins = c(8, 8))
background3 <- background.predict(microb.splsda.final, comp.predicted=2,
                                  dist = "max.dist") 
plotIndiv(microb.splsda.final, ncomp = 3, group = group,
          ind.names = FALSE, title = "Maximum distance",
          legend = TRUE,  background = background3)
plotLoadings(microb.splsda.final, contrib = 'max', method = 'mean', legend.title = "group", comp = 1, ndisplay = 30)
plotLoadings(microb.splsda.final, contrib = 'max', method = 'mean', legend.title = "group", comp = 2, ndisplay = 25)
plotLoadings(microb.splsda.final, contrib = 'max', method = 'mean', legend.title = "group", comp = 3, ndisplay = 20)
# form new perf() object which utilises the final model
perf.splsda.microb <- perf(microb.splsda.final, 
                          folds = 5, nrepeat = 50, # use repeated cross-validation
                          validation = "Mfold", dist = "max.dist",  # use max.dist measure
                          progressBar = FALSE)
# plot the stability of each feature for the first three components, 'h' type refers to histogram
par(mfrow=c(1,3))
plot(perf.splsda.microb$features$stable[[1]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features', 
     main = '(a) Comp 1', las =2)
plot(perf.splsda.microb$features$stable[[2]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features', 
     main = '(b) Comp 2', las =2)
plot(perf.splsda.microb$features$stable[[3]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features',
     main = '(c) Comp 3', las =2)
dev.off()

#genus
data=genera
file <- paste("genus",".pdf",sep="")
pdf(file, width=11, height=8.5)
genus.splsda <- splsda(data, group, ncomp = 3)  # 1 Run the method
plotIndiv(genus.splsda, ncomp = 3, ind.names = FALSE, ellipse=TRUE, group=group, star=TRUE, legend=TRUE, title = "sPLS-DA-all", legend.title = "Group", X.label = 'PLS-DA 1', Y.label = 'PLS-DA 2') # 2 Plot the samples
cim <- cim(genus.splsda, clust.method = c("ward", "ward"), comp = c(1,2), row.sideColors = stim.col[factor(group)], legend=list(legend = unique(group), col=stim.col, title = "Group", cex=0.9),margins = c(8, 8))
background <- background.predict(genus.splsda, comp.predicted=2, dist = "max.dist") 
plotIndiv(genus.splsda, comp = 1:2, group = group,
          ind.names = FALSE, title = "Maximum distance-all",
          legend = TRUE,  background = background)
auc.plsda <- auroc(genus.splsda)
plotLoadings(genus.splsda, contrib = 'max', method = 'mean',legend.title = "group", comp = 1, ndisplay = 50)
plotLoadings(genus.splsda, contrib = 'max', method = 'mean',legend.title = "group", comp = 2, ndisplay = 35)
plotLoadings(genus.splsda, contrib = 'max', method = 'mean',legend.title = "group", comp = 3, ndisplay = 15)
#selected
genus.splsda2 <- splsda(data, group, keepX=c(50,35,15), ncomp=3)  # 1 Run the method
plotIndiv(genus.splsda2, ind.names = FALSE, ellipse=TRUE, group=group, star=TRUE, legend=TRUE, title = "sPLS-DA", legend.title = "Group", X.label = 'PLS-DA 1', Y.label = 'PLS-DA 2') # 2 Plot the samples
background2 <- background.predict(genus.splsda2, comp.predicted=2,
                                  dist = "max.dist") 
plotIndiv(genus.splsda2, comp = 1:2, group = group,
          ind.names = FALSE, title = "Maximum distance",
          legend = TRUE,  background = background2)
cim <- cim(genus.splsda2,clust.method = c("ward", "ward"), comp = c(1,2,3), row.sideColors = stim.col[factor(group)], legend=list(legend = unique(group), col=stim.col, title = "Group", cex=0.9),margins = c(8, 8))
auc.plsda2 <- auroc(genus.splsda2)
list.keepX <- c(5:10,  seq(20, 100, 10))
list.keepX # to output the grid of values tested
set.seed(30) # for reproducbility in this vignette, otherwise increase nrepeat
tune.splsda.genus <- tune.splsda(data, group, ncomp = 4, # we suggest to push ncomp a bit more, e.g. 4
                                  validation = 'Mfold',
                                  folds = 5, dist = 'max.dist', progressBar = FALSE,
                                  measure = "BER", test.keepX = list.keepX,
                                  nrepeat = 50)   # we suggest nrepeat = 50
error <- tune.splsda.genus$error.rate
ncomp <- tune.splsda.genus$choice.ncomp$ncomp # optimal number of components based on t-tests on the error rate
ncomp
select.keepX <- tune.splsda.genus$choice.keepX[1:ncomp]  # optimal number of variables to select
select.keepX
plot(tune.splsda.genus)   #, col = color.jet(ncomp))
genus.splsda.final <- splsda(data, group, ncomp = ncomp, keepX = select.keepX)
legend=list(legend = levels(group), # set of classescol = unique(color.mixo(group)), # set of colours
            title = "Group", # legend title
            cex = 0.7) # legend size
plotIndiv(genus.splsda.final, ind.names = FALSE, ellipse=TRUE, group=group, star=TRUE, legend=TRUE, title = "sPLS-DA final result", legend.title = "Group", X.label = 'PLS-DA 1', Y.label = 'PLS-DA 2') # 2 Plot the samples
cim <- cim(genus.splsda.final, clust.method = c("ward", "ward"), comp = c(1,2,3), row.sideColors = stim.col[factor(group)], legend=list(legend = unique(group), col=stim.col, title = "Group", cex=0.9),margins = c(8, 8))
background3 <- background.predict(genus.splsda.final, comp.predicted=2,
                                  dist = "max.dist") 
plotIndiv(genus.splsda.final, ncomp = 3, group = group,
          ind.names = FALSE, title = "Maximum distance",
          legend = TRUE,  background = background3)
plotLoadings(genus.splsda.final, contrib = 'max', method = 'mean', legend.title = "group", comp = 1, ndisplay = 7)
plotLoadings(genus.splsda.final, contrib = 'max', method = 'mean', legend.title = "group", comp = 2, ndisplay = 40)
plotLoadings(genus.splsda.final, contrib = 'max', method = 'mean', legend.title = "group", comp = 3, ndisplay = 7)
# form new perf() object which utilises the final model
perf.splsda.microb <- perf(genus.splsda.final, 
                           folds = 5, nrepeat = 50, # use repeated cross-validation
                           validation = "Mfold", dist = "max.dist",  # use max.dist measure
                           progressBar = FALSE)
par(mfrow=c(1,3))
plot(perf.splsda.microb$features$stable[[1]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features', 
     main = '(a) Comp 1', las =2)
plot(perf.splsda.microb$features$stable[[2]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features', 
     main = '(b) Comp 2', las =2)
plot(perf.splsda.microb$features$stable[[3]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features',
     main = '(c) Comp 3', las =2)
dev.off()


#strat
data=strat
file <- paste("strat",".pdf",sep="")
pdf(file, width=11, height=8.5)
strat.splsda <- splsda(data, group, ncomp = 3)  # 1 Run the method
plotIndiv(strat.splsda, ncomp = 3, ind.names = FALSE, ellipse=TRUE, group=group, star=TRUE, legend=TRUE, title = "sPLS-DA-all", legend.title = "Group", X.label = 'PLS-DA 1', Y.label = 'PLS-DA 2') # 2 Plot the samples
cim <- cim(strat.splsda, clust.method = c("ward", "ward"), comp = c(1,2), row.sideColors = stim.col[factor(group)], legend=list(legend = unique(group), col=stim.col, title = "Group", cex=0.9))
background <- background.predict(strat.splsda, comp.predicted=2, dist = "max.dist") 
plotIndiv(strat.splsda, comp = 1:2, group = group,
          ind.names = FALSE, title = "Maximum distance-all",
          legend = TRUE,  background = background)
auc.plsda <- auroc(strat.splsda)
plotLoadings(strat.splsda, contrib = 'max', method = 'mean',legend.title = "group", comp = 1, ndisplay = 50)
plotLoadings(strat.splsda, contrib = 'max', method = 'mean',legend.title = "group", comp = 2, ndisplay = 50)
plotLoadings(strat.splsda, contrib = 'max', method = 'mean',legend.title = "group", comp = 3, ndisplay = 30)
strat.splsda2 <- splsda(data, group, keepX=c(50,30,10), ncomp=3)  # 1 Run the method
plotIndiv(strat.splsda2, ind.names = FALSE, ellipse=TRUE, group=group, star=TRUE, legend=TRUE, title = "sPLS-DA", legend.title = "Group", X.label = 'PLS-DA 1', Y.label = 'PLS-DA 2') # 2 Plot the samples
background2 <- background.predict(strat.splsda2, comp.predicted=2,
                                  dist = "max.dist") 
plotIndiv(strat.splsda2, comp = 1:2, group = group,
          ind.names = FALSE, title = "Maximum distance",
          legend = TRUE,  background = background2)
cim <- cim(strat.splsda2,clust.method = c("ward", "ward"), comp = c(1,2,3), row.sideColors = stim.col[factor(group)], legend=list(legend = unique(group), col=stim.col, title = "Group", cex=0.9))
auc.plsda2 <- auroc(strat.splsda2)
#The number of components to retain ncomp. The rule of thumb is usually K−1 where K is the number of classes, but it is worth testing a few extra components.
list.keepX <- c(5:10,  seq(20, 100, 10))
list.keepX # to output the grid of values tested
set.seed(30) # for reproducbility in this vignette, otherwise increase nrepeat
tune.splsda.strat <- tune.splsda(data, group, ncomp = 4, # we suggest to push ncomp a bit more, e.g. 4
                                  validation = 'Mfold',
                                  folds = 5, dist = 'max.dist', progressBar = FALSE,
                                  measure = "BER", test.keepX = list.keepX,
                                  nrepeat = 50)   # we suggest nrepeat = 50
error <- tune.splsda.strat$error.rate
ncomp <- tune.splsda.strat$choice.ncomp$ncomp # optimal number of components based on t-tests on the error rate
ncomp
select.keepX <- tune.splsda.strat$choice.keepX[1:ncomp]  # optimal number of variables to select
select.keepX
plot(tune.splsda.strat)   #, col = color.jet(ncomp))
strat.splsda.final <- splsda(data, group, ncomp = ncomp, keepX = select.keepX)
legend=list(legend = levels(group), # set of classescol = unique(color.mixo(group)), # set of colours
            title = "Group", # legend title
            cex = 0.7) # legend size
plotIndiv(strat.splsda.final, ind.names = FALSE, ellipse=TRUE, group=group, star=TRUE, legend=TRUE, title = "sPLS-DA final result", legend.title = "Group", X.label = 'PLS-DA 1', Y.label = 'PLS-DA 2') # 2 Plot the samples
cim <- cim(strat.splsda.final, clust.method = c("ward", "ward"), comp = c(1,2,3), row.sideColors = stim.col[factor(group)], legend=list(legend = unique(group), col=stim.col, title = "Group", cex=0.9))
background3 <- background.predict(strat.splsda.final, comp.predicted=2,
                                  dist = "max.dist") 
plotIndiv(strat.splsda.final, ncomp = 3, group = group,
          ind.names = FALSE, title = "Maximum distance",
          legend = TRUE,  background = background3)
plotLoadings(strat.splsda.final, contrib = 'max', method = 'mean', legend.title = "group", comp = 1, ndisplay = 30)
plotLoadings(strat.splsda.final, contrib = 'max', method = 'mean', legend.title = "group", comp = 2, ndisplay = 20)
plotLoadings(strat.splsda.final, contrib = 'max', method = 'mean', legend.title = "group", comp = 3, ndisplay = 50)
# form new perf() object which utilises the final model
perf.splsda.strat <- perf(strat.splsda.final, 
                           folds = 5, nrepeat = 50, # use repeated cross-validation
                           validation = "Mfold", dist = "max.dist",  # use max.dist measure
                           progressBar = FALSE)
par(mfrow=c(1,3))
plot(perf.splsda.strat$features$stable[[1]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features', 
     main = '(a) Comp 1', las =2)
plot(perf.splsda.strat$features$stable[[2]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features', 
     main = '(b) Comp 2', las =2)
plot(perf.splsda.strat$features$stable[[3]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features',
     main = '(c) Comp 3', las =2)
dev.off()

data=unstrat
file <- paste("unstrat",".pdf",sep="")
pdf(file, width=11, height=8.5)
unstrat.splsda <- splsda(data, group, ncomp = 3)  # 1 Run the method
plotIndiv(unstrat.splsda, ncomp = 3, ind.names = FALSE, ellipse=TRUE, group=group, star=TRUE, legend=TRUE, title = "sPLS-DA-all", legend.title = "Group", X.label = 'PLS-DA 1', Y.label = 'PLS-DA 2') # 2 Plot the samples
cim <- cim(unstrat.splsda, clust.method = c("ward", "ward"), comp = c(1,2), row.sideColors = stim.col[factor(group)], legend=list(legend = unique(group), col=stim.col, title = "Group", cex=0.9))
background <- background.predict(unstrat.splsda, comp.predicted=2, dist = "max.dist") 
plotIndiv(unstrat.splsda, comp = 1:2, group = group,
          ind.names = FALSE, title = "Maximum distance-all",
          legend = TRUE,  background = background)
auc.plsda <- auroc(unstrat.splsda)
plotLoadings(unstrat.splsda, contrib = 'max', method = 'mean',legend.title = "group", comp = 1, ndisplay = 50)
plotLoadings(unstrat.splsda, contrib = 'max', method = 'mean',legend.title = "group", comp = 2, ndisplay = 50)
plotLoadings(unstrat.splsda, contrib = 'max', method = 'mean',legend.title = "group", comp = 3, ndisplay = 50)
#selected
unstrat.splsda2 <- splsda(data, group, keepX=c(50,30,10), ncomp=3)  # 1 Run the method
plotIndiv(unstrat.splsda2, ind.names = FALSE, ellipse=TRUE, group=group, star=TRUE, legend=TRUE, title = "sPLS-DA", legend.title = "Group", X.label = 'PLS-DA 1', Y.label = 'PLS-DA 2') # 2 Plot the samples
background2 <- background.predict(unstrat.splsda2, comp.predicted=2,
                                  dist = "max.dist") 
plotIndiv(unstrat.splsda2, comp = 1:2, group = group,
          ind.names = FALSE, title = "Maximum distance",
          legend = TRUE,  background = background2)
cim <- cim(unstrat.splsda2,clust.method = c("ward", "ward"), comp = c(1,2,3), row.sideColors = stim.col[factor(group)], legend=list(legend = unique(group), col=stim.col, title = "Group", cex=0.9))
auc.plsda2 <- auroc(unstrat.splsda2)
list.keepX <- c(5:10,  seq(20, 100, 10))
list.keepX # to output the grid of values tested
set.seed(30) # for reproducbility in this vignette, otherwise increase nrepeat
tune.splsda.unstrat <- tune.splsda(data, group, ncomp = 4, # we suggest to push ncomp a bit more, e.g. 4
                                  validation = 'Mfold',
                                  folds = 5, dist = 'max.dist', progressBar = FALSE,
                                  measure = "BER", test.keepX = list.keepX,
                                  nrepeat = 50)   # we suggest nrepeat = 50
error <- tune.splsda.unstrat$error.rate
ncomp <- tune.splsda.unstrat$choice.ncomp$ncomp # optimal number of components based on t-tests on the error rate
ncomp
select.keepX <- tune.splsda.unstrat$choice.keepX[1:ncomp]  # optimal number of variables to select
select.keepX
plot(tune.splsda.unstrat)   #, col = color.jet(ncomp))
unstrat.splsda.final <- splsda(data, group, ncomp = ncomp, keepX = select.keepX)
legend=list(legend = levels(group), # set of classescol = unique(color.mixo(group)), # set of colours
            title = "Group", # legend title
            cex = 0.7) # legend size
plotIndiv(unstrat.splsda.final, ind.names = FALSE, ellipse=TRUE, group=group, star=TRUE, legend=TRUE, title = "sPLS-DA final result", legend.title = "Group", X.label = 'PLS-DA 1', Y.label = 'PLS-DA 2') # 2 Plot the samples
cim <- cim(unstrat.splsda.final, clust.method = c("ward", "ward"), comp = c(1,2,3), row.sideColors = stim.col[factor(group)], legend=list(legend = unique(group), col=stim.col, title = "Group", cex=0.9))
background3 <- background.predict(unstrat.splsda.final, comp.predicted=2,
                                  dist = "max.dist") 
plotIndiv(unstrat.splsda.final, ncomp = 3, group = group,
          ind.names = FALSE, title = "Maximum distance",
          legend = TRUE,  background = background3)
plotLoadings(unstrat.splsda.final, contrib = 'max', method = 'mean', legend.title = "group", comp = 1, ndisplay = 20)
plotLoadings(unstrat.splsda.final, contrib = 'max', method = 'mean', legend.title = "group", comp = 2, ndisplay = 50)
plotLoadings(unstrat.splsda.final, contrib = 'max', method = 'mean', legend.title = "group", comp = 3, ndisplay = 50)
perf.splsda.unstrat <- perf(unstrat.splsda.final, 
                           folds = 5, nrepeat = 50, # use repeated cross-validation
                           validation = "Mfold", dist = "max.dist",  # use max.dist measure
                           progressBar = FALSE)
par(mfrow=c(1,3))
plot(perf.splsda.unstrat$features$stable[[1]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features', 
     main = '(a) Comp 1', las =2)
plot(perf.splsda.unstrat$features$stable[[2]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features', 
     main = '(b) Comp 2', las =2)
plot(perf.splsda.unstrat$features$stable[[3]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features',
     main = '(c) Comp 3', las =2)
dev.off()

#sera
data=sera
file <- paste("sera",".pdf",sep="")
pdf(file, width=11, height=8.5)
sera.splsda <- splsda(data, group, ncomp = 3)  # 1 Run the method
plotIndiv(sera.splsda, ncomp = 3, ind.names = FALSE, ellipse=TRUE, group=group, star=TRUE, legend=TRUE, title = "sPLS-DA-all", legend.title = "Group", X.label = 'PLS-DA 1', Y.label = 'PLS-DA 2') # 2 Plot the samples
cim <- cim(sera.splsda, clust.method = c("ward", "ward"), comp = c(1,2), row.sideColors = stim.col[factor(group)], legend=list(legend = unique(group), col=stim.col, title = "Group", cex=0.9))
background <- background.predict(sera.splsda, comp.predicted=2, dist = "max.dist") 
plotIndiv(sera.splsda, comp = 1:2, group = group,
          ind.names = FALSE, title = "Maximum distance-all",
          legend = TRUE,  background = background)
auc.plsda <- auroc(sera.splsda)
plotLoadings(sera.splsda, contrib = 'max', method = 'mean',legend.title = "group", comp = 1, ndisplay = 50)
plotLoadings(sera.splsda, contrib = 'max', method = 'mean',legend.title = "group", comp = 2, ndisplay = 50)
plotLoadings(sera.splsda, contrib = 'max', method = 'mean',legend.title = "group", comp = 3, ndisplay = 50)
#selected
sera.splsda2 <- splsda(data, group, keepX=c(50,30,10), ncomp=3)  # 1 Run the method
plotIndiv(sera.splsda2, ind.names = FALSE, ellipse=TRUE, group=group, star=TRUE, legend=TRUE, title = "sPLS-DA", legend.title = "Group", X.label = 'PLS-DA 1', Y.label = 'PLS-DA 2') # 2 Plot the samples
background2 <- background.predict(sera.splsda2, comp.predicted=2,
                                  dist = "max.dist") 
plotIndiv(sera.splsda2, comp = 1:2, group = group,
          ind.names = FALSE, title = "Maximum distance",
          legend = TRUE,  background = background2)
cim <- cim(sera.splsda2,clust.method = c("ward", "ward"), comp = c(1,2,3), row.sideColors = stim.col[factor(group)], legend=list(legend = unique(group), col=stim.col, title = "Group", cex=0.9))
auc.plsda2 <- auroc(sera.splsda2)
#The number of components to retain ncomp. The rule of thumb is usually K−1 where K is the number of classes, but it is worth testing a few extra components.
list.keepX <- c(5:10,  seq(20, 100, 10))
list.keepX # to output the grid of values tested
set.seed(30) # for reproducbility in this vignette, otherwise increase nrepeat
tune.splsda.sera <- tune.splsda(data, group, ncomp = 4, # we suggest to push ncomp a bit more, e.g. 4
                                  validation = 'Mfold',
                                  folds = 5, dist = 'max.dist', progressBar = FALSE,
                                  measure = "BER", test.keepX = list.keepX,
                                  nrepeat = 50)   # we suggest nrepeat = 50
error <- tune.splsda.sera$error.rate
ncomp <- tune.splsda.sera$choice.ncomp$ncomp # optimal number of components based on t-tests on the error rate
ncomp
select.keepX <- tune.splsda.sera$choice.keepX[1:ncomp]  # optimal number of variables to select
select.keepX
plot(tune.splsda.sera)   #, col = color.jet(ncomp))
sera.splsda.final <- splsda(data, group, ncomp = ncomp, keepX = select.keepX)
legend=list(legend = levels(group), # set of classescol = unique(color.mixo(group)), # set of colours
            title = "Group", # legend title
            cex = 0.7) # legend size
plotIndiv(sera.splsda.final, ind.names = FALSE, ellipse=TRUE, group=group, star=TRUE, legend=TRUE, title = "sPLS-DA final result", legend.title = "Group", X.label = 'PLS-DA 1', Y.label = 'PLS-DA 2') # 2 Plot the samples
cim <- cim(sera.splsda.final, clust.method = c("ward", "ward"), comp = c(1,2,3), row.sideColors = stim.col[factor(group)], legend=list(legend = unique(group), col=stim.col, title = "Group", cex=0.9))
background3 <- background.predict(sera.splsda.final, comp.predicted=2,
                                  dist = "max.dist") 
plotIndiv(sera.splsda.final, ncomp = 3, group = group,
          ind.names = FALSE, title = "Maximum distance",
          legend = TRUE,  background = background3)
plotLoadings(sera.splsda.final, contrib = 'max', method = 'mean', legend.title = "group", comp = 1, ndisplay = 40)
plotLoadings(sera.splsda.final, contrib = 'max', method = 'mean', legend.title = "group", comp = 2, ndisplay = 40)
plotLoadings(sera.splsda.final, contrib = 'max', method = 'mean', legend.title = "group", comp = 3, ndisplay = 45)
perf.splsda.sera <- perf(sera.splsda.final, 
                           folds = 5, nrepeat = 50, # use repeated cross-validation
                           validation = "Mfold", dist = "max.dist",  # use max.dist measure
                           progressBar = FALSE)
# plot the stability of each feature for the first three components, 'h' type refers to histogram
par(mfrow=c(1,3))
plot(perf.splsda.sera$features$stable[[1]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features', 
     main = '(a) Comp 1', las =2)
plot(perf.splsda.sera$features$stable[[2]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features', 
     main = '(b) Comp 2', las =2)
plot(perf.splsda.sera$features$stable[[3]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features',
     main = '(c) Comp 3', las =2)
dev.off()



