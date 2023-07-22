library(mixOmics) 

fecal <- read.csv("fecalMB2.csv", row.names=1, check.names=T, header=T)
serum <- read.csv("serum2.csv", row.names=1, check.names=T, header=T)
set1n <- read.csv("set1n2.csv", row.names=1, check.names=T, header=T)
genus <- read.csv("genus.csv", row.names=1, check.names=T, header=T)
family <- read.csv("family.csv", row.names=1, check.names=T, header=T)
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
dim(lumen);dim(microb);dim(strat); dim(unstrat);dim(sera); dim(genera); dim(families)
stim.col <- c("darkblue", "purple", "green4","red3")

postscript("pca_family.eps",hor=T)
pca.lumen <- pca(lumen, ncomp = 10, center = TRUE, scale = TRUE)
plot(pca.lumen)
pca.species <- pca(microb, ncomp = 10, center = TRUE, scale = TRUE)
plot(pca.species)
pca.genus <- pca(genera, ncomp = 10, center = TRUE, scale = TRUE)
plot(pca.genus)
pca.family <- pca(families, ncomp = 10, center = TRUE, scale = TRUE)
plot(pca.family)
dev.off()

optimal.ncomp = 3
lumenSpecies.spls <- spls(microb,lumen,scale=TRUE, ncomp=5, mode = "regression")  #default is 2 components; data are scaled (variance = 1, strongly advised here); by default a PLS regression mode should be used 
plotIndiv(lumenSpecies.spls,group = group,legend = TRUE, legend.title = 'Group', ind.names = FALSE, title = "sPLS-all-Species/Lumen", scale = TRUE)       
perf.lumenSpecies.spls<- perf(lumenSpecies.spls, validation = 'Mfold',
                             folds = 10, nrepeat = 5) 
plot(perf.lumenSpecies.spls, criterion = 'Q2.total')
list.keepX <- c(seq(20, 50, 5))
list.keepY <- c(3:10) 
tune.lumenSpecies.spls <- tune.spls(microb,lumen, ncomp = 3,
                                  test.keepX = list.keepX,
                                  test.keepY = list.keepY,
                                  nrepeat = 1, folds = 10, # use 10 folds
                                  mode = 'regression', measure = 'cor') 
postscript("tune.eps",hor=T)
plot(tune.lumenSpecies.spls)         # use the correlation measure for tuning
dev.off()
tune.lumenSpecies.spls$choice.keepX
optimal.keepX <- tune.lumenSerum.spls$choice.keepX 
optimal.keepY <- tune.lumenSerum.spls$choice.keepY
optimal.keepY <- c(25,20,10)
names(optimal.keepY) <- c("comp1","comp2","comp3")
final.speciesLumen.spls <- spls(microb,lumen, ncomp = optimal.ncomp, 
                              keepX = optimal.keepX,
                              keepY = optimal.keepY,
                              mode = "regression") # explanitory approach being used, # hence use regression mode
pdf(file="variables.pdf")
plotVar(final.speciesLumen.spls, cutoff=0.7, legend = c("Species", "Lumen"), legend.title="Type", cex=c(3,3))        
dev.off()

pdf(file="PLS_final_species.pdf")
plotIndiv(final.speciesLumen.spls, ind.names = FALSE, 
          rep.space = "X-variate", # plot in X-variate subspace
          group = group, # colour by time group
          pch = as.factor(group), 
          col.per.group = color.mixo(1:4), 
          legend = TRUE, legend.title = 'Group')
plotIndiv(final.speciesLumen.spls, ind.names = FALSE,
          rep.space = "Y-variate", # plot in Y-variate subspace
          group = group, # colour by time group
          pch = as.factor(group), 
          col.per.group = color.mixo(1:4), 
          legend = TRUE, legend.title = 'Group')
#Samples are projected into the space spanned by the averaged components of both datasets.
plotIndiv(final.speciesLumen.spls, ind.names = FALSE, 
          rep.space = "XY-variate", # plot in averaged subspace
          group = group, # colour by time group
          pch = as.factor(group), # select symbol
          col.per.group = color.mixo(1:4),                      # by dose group
          legend = TRUE, legend.title = 'Group')
dev.off()
col.tox <- color.mixo(1:4) # create set of colours
plotIndiv(final.speciesLumen.spls, ind.names = FALSE, 
          rep.space = "XY-variate", # plot in averaged subspace
          axes.box = "both", col = col.tox, style = '3d')

#visualise the level of agreement between data sets
postscript("arrow_species_lumen_final.eps",hor=T)
plotArrow(final.speciesLumen.spls, 
          group=group, legend = TRUE,
          arrow.alpha=1,
          ind.names.position = c("start", "end"),
          ind.names.size = 4,
          pch.size = 3, #circule size
          arrow.size = 1, #arrow thickness
          arrow.length = 0.3, #arrow head
          X.label = 'PLS comp 1', Y.label = 'PLS comp 2', legend.title = 'Group', title = 'Lumen -> Serum')
dev.off()

pdf(file="network_species.pdf")
network(final.speciesLumen.spls, comp = 1:3, cutoff = 0.85, name.save = 'PLSnetwork', color.node = c("cyan", "pink"),
        shape.node = c("rectangle", "circle"), color.edge = color.spectral(100), lty.edge = "solid", lwd.edge =  1, show.edge.labels = FALSE, interactive = FALSE)
network(final.speciesLumen.spls, comp = 2, cutoff = 0.75, name.save = 'PLSnetwork', color.node = c("cyan", "pink"),
        shape.node = c("rectangle", "circle"), color.edge = color.spectral(100), lty.edge = "solid", lwd.edge =  1, show.edge.labels = FALSE, interactive = FALSE)
network(final.speciesLumen.spls, comp = 3, cutoff = 0.6, name.save = 'PLSnetwork', color.node = c("cyan", "pink"),
        shape.node = c("rectangle", "circle"), color.edge = color.spectral(100), lty.edge = "solid", lwd.edge =  1, show.edge.labels = FALSE, interactive = FALSE)
dev.off()

pdf(file="cim_species.pdf")
cim(final.speciesLumen.spls, xlab = "Lumen", ylab = "Species", clust.method = c("ward", "ward"), comp = 1:3, margins = c(8, 8))
cim(final.speciesLumen.spls, xlab = "Lumen", ylab = "Species", clust.method = c("ward", "ward"), comp = 2, margins = c(8, 8))
cim(final.speciesLumen.spls, xlab = "Lumen", ylab = "Species", clust.method = c("ward", "ward"), comp = 3, margins = c(8, 8))
dev.off()

####Genus to lumen

optimal.ncomp = 3
##initial model
genusLumen.spls <- spls(genera,lumen,scale=TRUE, ncomp=5, mode = "regression")  #default is 2 components; data are scaled (variance = 1, strongly advised here); by default a PLS regression mode should be used 
plotIndiv(genusLumen.spls,group = group,legend = TRUE, legend.title = 'Group', ind.names = FALSE, title = "sPLS-all-Genus/Lumen", scale = TRUE)       
perf.genusLumen.spls<- perf(genusLumen.spls, validation = 'Mfold',
                              folds = 10, nrepeat = 5) 
plot(perf.genusLumen.spls, criterion = 'Q2.total')
list.keepX <- c(seq(20, 50, 5))
list.keepY <- c(3:10) 
tune.genusLumen.spls <- tune.spls(genera,lumen, ncomp = 3,
                                    test.keepX = list.keepX,
                                    test.keepY = list.keepY,
                                    nrepeat = 1, folds = 10, # use 10 folds
                                    mode = 'regression', measure = 'cor') 
postscript("tune.eps",hor=T)
plot(tune.genusLumen.spls)         # use the correlation measure for tuning
dev.off()
tune.genusLumen.spls$choice.keepX
optimal.keepX <- tune.genusLumen.spls$choice.keepX 
optimal.keepX
optimal.keepY <- tune.genusLumen.spls$choice.keepY
optimal.keepY #    3     4     4 
optimal.keepY <- c(45,25,30)
names(optimal.keepY) <- c("comp1","comp2","comp3")
final.genusLumen.spls <- spls(genera,lumen, ncomp = optimal.ncomp, 
                                keepX = optimal.keepX,
                                keepY = optimal.keepY,
                                mode = "regression") # explanitory approach being used, # hence use regression mode
pdf(file="variables.pdf")
plotVar(final.genusLumen.spls, cutoff=0.6, legend = c("Genus", "Lumen"), legend.title="Type", cex=c(3,3))        
dev.off()

pdf(file="PLS_final_genus.pdf")
plotIndiv(final.genusLumen.spls, ind.names = FALSE, 
          rep.space = "X-variate", # plot in X-variate subspace
          group = group, # colour by time group
          pch = as.factor(group), 
          col.per.group = color.mixo(1:4), 
          legend = TRUE, legend.title = 'Group')
plotIndiv(final.genusLumen.spls, ind.names = FALSE,
          rep.space = "Y-variate", # plot in Y-variate subspace
          group = group, # colour by time group
          pch = as.factor(group), 
          col.per.group = color.mixo(1:4), 
          legend = TRUE, legend.title = 'Group')
#Samples are projected into the space spanned by the averaged components of both datasets.
plotIndiv(final.genusLumen.spls, ind.names = FALSE, 
          rep.space = "XY-variate", # plot in averaged subspace
          group = group, # colour by time group
          pch = as.factor(group), # select symbol
          col.per.group = color.mixo(1:4),                      # by dose group
          legend = TRUE, legend.title = 'Group')
dev.off()
col.tox <- color.mixo(1:4) # create set of colours
plotIndiv(final.genusLumen.spls, ind.names = FALSE, 
          rep.space = "XY-variate", # plot in averaged subspace
          axes.box = "both", col = col.tox, style = '3d')

#visualise the level of agreement between data sets
postscript("arrow_genus_lumen_final2.eps",hor=T)
plotArrow(final.genusLumen.spls, 
          group=group, legend = TRUE,
          arrow.alpha=1,
          ind.names.position = c("start", "end"),
          ind.names.size = 4,
          pch.size = 3, #circule size
          arrow.size = 1, #arrow thickness
          arrow.length = 0.3, #arrow head
          X.label = 'PLS comp 1', Y.label = 'PLS comp 2', legend.title = 'Group', title = 'Lumen -> Serum')
dev.off()

pdf(file="network_genus.pdf")
network(final.genusLumen.spls, comp = 1, cutoff = 0.75, name.save = 'PLSnetwork', color.node = c("cyan", "pink"),
        shape.node = c("rectangle", "circle"), color.edge = color.spectral(100), lty.edge = "solid", lwd.edge =  1, show.edge.labels = FALSE, interactive = FALSE)
network(final.genusLumen.spls, comp = 2, cutoff = 0.6, name.save = 'PLSnetwork', color.node = c("cyan", "pink"),
        shape.node = c("rectangle", "circle"), color.edge = color.spectral(100), lty.edge = "solid", lwd.edge =  1, show.edge.labels = FALSE, interactive = FALSE)
network(final.genusLumen.spls, comp = 3, cutoff = 0.6, name.save = 'PLSnetwork', color.node = c("cyan", "pink"),
        shape.node = c("rectangle", "circle"), color.edge = color.spectral(100), lty.edge = "solid", lwd.edge =  1, show.edge.labels = FALSE, interactive = FALSE)
dev.off()

pdf(file="cim_genus_set.pdf")
cim(final.genusLumen.spls, xlab = "Lumen", ylab = "Genus", clust.method = c("ward", "ward"), comp = 1:2, margins = c(8, 8))
cim(final.genusLumen.spls, xlab = "Lumen", ylab = "Genus", clust.method = c("ward", "ward"), comp = 1, margins = c(8, 8))
cim(final.genusLumen.spls, xlab = "Lumen", ylab = "Genus", clust.method = c("ward", "ward"), comp = 2, margins = c(8, 8))
cim(final.genusLumen.spls, xlab = "Lumen", ylab = "Genus", clust.method = c("ward", "ward"), comp = 3, margins = c(8, 8))
dev.off()


postscript("pca_serum.eps",hor=T)
pca.lumen <- pca(lumen, ncomp = 10, center = TRUE, scale = TRUE)
pca.serum <- pca(sera, ncomp = 10, center = TRUE, scale = TRUE)
plot(pca.lumen)
plot(pca.serum)
dev.off()
optimal.ncomp = 3

lumenSerum.spls <- spls(lumen,sera,scale=TRUE, ncomp=3, mode = "regression")  #default is 2 components; data are scaled (variance = 1, strongly advised here); by default a PLS regression mode should be used 
plotIndiv(lumenSerum.spls,group = group,legend = TRUE, legend.title = 'Group', ind.names = FALSE, title = "sPLS-all-Lumen/Serum", scale = TRUE)       
perf.lumenSerum.spls <- perf(lumenSerum.spls, validation = 'Mfold',
                        folds = 10, nrepeat = 5) 
plot(perf.lumenSerum.spls, criterion = 'Q2.total')
list.keepX <- c(seq(20, 50, 5))
list.keepY <- c(3:10) 
tune.lumenSerum.spls <- tune.spls(lumen,sera, ncomp = 3,
                             test.keepX = list.keepX,
                             test.keepY = list.keepY,
                             nrepeat = 1, folds = 10, # use 10 folds
                             mode = 'regression', measure = 'cor') 
postscript("tune.eps",hor=T)
plot(tune.lumenSerum.spls)         # use the correlation measure for tuning
dev.off()
tune.lumenSerum.spls$choice.keepX
optimal.keepX <- tune.lumenSerum.spls$choice.keepX 
optimal.keepY <- tune.lumenSerum.spls$choice.keepY
optimal.keepX <- c(25,35,30)
names(optimal.keepX) <- c("comp1","comp2","comp3")

optimal.keepY <- c(45,35,35)
names(optimal.keepY) <- c("comp1","comp2","comp3")
final.lumenSerum.spls <- spls(lumen,sera, ncomp = optimal.ncomp, 
                         keepX = optimal.keepX,
                         keepY = optimal.keepY,
                         mode = "regression") # explanitory approach being used, # hence use regression mode
pdf(file="heatmap_correlation.pdf")
cim(final.lumenSerum.spls, xlab = "serum", ylab = "lumen", clust.method = c("ward", "ward"), comp = 1, margins = c(8, 8))
cim(final.lumenSerum.spls, xlab = "serum", ylab = "lumen", clust.method = c("ward", "ward"), comp = 2, margins = c(8, 8))
cim(final.lumenSerum.spls, xlab = "serum", ylab = "lumen", clust.method = c("ward", "ward"), comp = 3, margins = c(8, 8))
dev.off()

plotVar(final.lumenSerum.spls, cutoff=0.7, legend = c("Lumen", "Serum"), legend.title="Type", cex=c(4,4))        
plotVar(final.lumenSerum.spls, cutoff=0.75, legend = c("Lumen", "Serum"), legend.title="Type", cex=c(4,4))    
plotVar(final.lumenSerum.spls, cutoff=0.75, legend = c("Lumen", "Serum"), var.names = c(FALSE, TRUE))      
plotVar(final.lumenSerum.spls, cutoff=0.75, legend = c("Lumen", "Serum"), var.names = c(TRUE, FALSE))        
dev.off()

pdf(file="PLS_final.pdf")
#Samples are projected into the space spanned by the components associated to each data set (or block).
plotIndiv(final.lumenSerum.spls, ind.names = FALSE, 
          rep.space = "X-variate", # plot in X-variate subspace
          group = group, # colour by time group
          pch = as.factor(group), 
          col.per.group = color.mixo(1:4), 
          legend = TRUE, legend.title = 'Group')
plotIndiv(final.lumenSerum.spls, ind.names = FALSE,
          rep.space = "Y-variate", # plot in Y-variate subspace
          group = group, # colour by time group
          pch = as.factor(group), 
          col.per.group = color.mixo(1:4), 
          legend = TRUE, legend.title = 'Group')
#Samples are projected into the space spanned by the averaged components of both datasets.
plotIndiv(final.lumenSerum.spls, ind.names = FALSE, 
          rep.space = "XY-variate", # plot in averaged subspace
          group = group, # colour by time group
          pch = as.factor(group), # select symbol
          col.per.group = color.mixo(1:4),                      # by dose group
          legend = TRUE, legend.title = 'Group')
dev.off()
col.tox <- color.mixo(1:4) # create set of colours
plotIndiv(final.lumenSerum.spls, ind.names = FALSE, 
          rep.space = "XY-variate", # plot in averaged subspace
          axes.box = "both", col = col.tox, style = '3d')

#visualise the level of agreement between data sets
postscript("arrow_lumen_serum_final2.eps",hor=T)
plotArrow(final.lumenSerum.spls, 
          group=group, legend = TRUE,
          arrow.alpha=1,
          ind.names.position = c("start", "end"),
          ind.names.size = 4,
          pch.size = 3, #circule size
          arrow.size = 1, #arrow thickness
          arrow.length = 0.3, #arrow head
          X.label = 'PLS comp 1', Y.label = 'PLS comp 2', legend.title = 'Group', title = 'Lumen -> Serum')
dev.off()

pdf(file="network.pdf")
network(final.lumenSerum.spls, comp = 1, cutoff = 0.8, name.save = 'PLSnetwork', color.node = c("cyan", "pink"),
        shape.node = c("rectangle", "circle"), color.edge = color.spectral(100), lty.edge = "solid", lwd.edge =  1, show.edge.labels = FALSE, interactive = FALSE)
network(final.lumenSerum.spls, comp = 2, cutoff = 0.7, name.save = 'PLSnetwork', color.node = c("cyan", "pink"),
        shape.node = c("rectangle", "circle"), color.edge = color.spectral(100), lty.edge = "solid", lwd.edge =  1, show.edge.labels = FALSE, interactive = FALSE)
network(final.lumenSerum.spls, comp = 3, cutoff = 0.6, name.save = 'PLSnetwork', color.node = c("cyan", "pink"),
        shape.node = c("rectangle", "circle"), color.edge = color.spectral(100), lty.edge = "solid", lwd.edge =  1, show.edge.labels = FALSE, interactive = FALSE)
dev.off()

#coordinates <- plotVar(final.lumenSerum.spls, plot = FALSE)
#cim(lumenSerum.spls, xlab = "lumen", ylab = "serum", clust.method = c("ward", "ward"), comp = 2, margins = c(5, 5), save = 'pdf')
cim(final.lumenSerum.spls, xlab = "lumen", ylab = "serum", clust.method = c("ward", "ward"), comp = 1:3, margins = c(8, 8))
cim(final.lumenSerum.spls, xlab = "lumen", ylab = "serum", clust.method = c("ward", "ward"), comp = 2, margins = c(8, 8))
cim(final.lumenSerum.spls, xlab = "lumen", ylab = "serum", clust.method = c("ward", "ward"), comp = 3, margins = c(8, 8))
plotVar(lumenSerum.spls, cutoff=0.7, legend = c("Lumen", "Serum"), var.names = c(FALSE, TRUE))      
plotVar(lumenSerum.spls, cutoff=0.7, legend = c("Lumen", "Serum"), var.names = c(TRUE, FALSE))        

plotLoadings(lumenSerum.spls, contrib = 'max', method = 'median', legend.title = "group", comp = 1, ndisplay = 50)

lumenSerum.spls2 <- spls(lumen,sera, keepX = c(50,50,50), keepY = c(50,50,50), ncomp=3, scale=TRUE, mode = "regression")  #default is 2 components; data are scaled (variance = 1, strongly advised here); by default a PLS regression mode should be used 
plotIndiv(lumenSerum.spls2,group = group,legend = TRUE, legend.title = 'Group', ind.names = FALSE, title = "sPLS-subset-Lumen/Serum", scale = TRUE)       
network(final.lumenSerum.spls, comp = 1, cutoff = 0.8, name.save = 'PLSnetwork', color.node = c("cyan", "pink"),
        shape.node = c("rectangle", "circle"), color.edge = color.spectral(100), lty.edge = "solid", lwd.edge =  1, show.edge.labels = FALSE, interactive = FALSE)
network(lumenSerum.spls2, comp = 2, cutoff = 0.8, name.save = 'PLSnetwork', color.node = c("mistyrose", "lightcyan"),
        shape.node = c("rectangle", "circle"), color.edge = color.spectral(100), lty.edge = "solid", lwd.edge =  1, show.edge.labels = FALSE, interactive = FALSE)
cim(lumenSerum.spls2, xlab = "lumen", ylab = "serum", clust.method = c("ward", "ward"), comp = 1, margins = c(8, 8))
cim(lumenSerum.spls2, xlab = "lumen", ylab = "serum", clust.method = c("ward", "ward"), comp = 2, margins = c(8, 8))
cim(lumenSerum.spls2, xlab = "lumen", ylab = "serum", clust.method = c("ward", "ward"), comp = 3, margins = c(8, 8))

plotArrow(lumenSerum.spls2,group=group, legend = TRUE,
          X.label = 'PLS comp 1', Y.label = 'PLS comp 2', legend.title = 'Group')


