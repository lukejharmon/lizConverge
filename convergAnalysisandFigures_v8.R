library(geometry)
library(geiger)
library(pcaMethods)
library(MASS)
library(TreeSim)

setwd("~/Documents/lizardMorphospace/analyses/")

# get trait data
sqMorphA<-read.csv("squamorph_v2.csv", na.strings=c("", " ", "."))

# pull out just continuously distributed characters
sqMorphContinuous<-sqMorphA[1:384, 4:26]

species<-as.character(sqMorphA[1:384, 1])

#manual translation of outdated names
altSpecies<-species
altSpecies[species=="Cosymbotus"]<-"Hemidactylus craspedotus" # former monotypic genus
altSpecies[species=="Gecko"]<-"Gekko albofasciolatus" # arbitrary
altSpecies[species=="Geckonia"]<-"Tarentola chazaliae" # former monotypic genus
altSpecies[species=="Microscalabo"]<-"Lygodactylus bivittis" # former monotypic genus
altSpecies[species=="Palmatogecko"]<-"Pachydactylus vanzyli" # former genus w 2 species
altSpecies[species=="Tenuidactylu"]<-"Cyrtopodion caspium" # not sure what happened here
altSpecies[species=="Teratolepis"]<-"Hemidactylus albofasciatus" # small genus absorbed into hemidac?
altSpecies[species=="Chamaeleojax"]<-"Trioceros jacksonii" # reclassified
altSpecies[species=="Chamaeleorhi"]<-"Furcifer rhinoceratus" # has to be? rhi?
altSpecies[species=="amphibbarbat"]<-"Pogona barbata" # I am pretty sure these three are correct
altSpecies[species=="amphibisolep"]<-"Ctenophorus isolepis" #
altSpecies[species=="amphibnuchal"]<-"Ctenophorus nuchalis" #
# I think "Stellio" is correct - the listed species has a synonym of Stellio stellio
altSpecies[species=="Sator"]<-"Sceloporus angustus" #
altSpecies[species=="Scelopmagist"]<-"Sceloporus magister" #
altSpecies[species=="Scelopoccide"]<-"Sceloporus occidentalis" #
altSpecies[species=="Sceloporussi"]<-"Sceloporus rufidorsum" # I have no idea what this might be, but it's in the right genus at least
altSpecies[species=="Aperopristis"]<-"Leiosaurus catamarcensis" # old name
altSpecies[species=="Aptycholaemu"]<-"Anisolepis longicauda" # old name
altSpecies[species=="Chamaeleolis"]<-"Anolis chamaeleonides" # old name
altSpecies[species=="Chamaelinoro"]<-"Anolis barbouri" # old name
altSpecies[species=="Phenacosauru"]<-"Anolis bellipeniculus" # old name
altSpecies[species=="Ophryoessoid"]<-"Stenocercus caducus" # old name
altSpecies[species=="Proctotrectu"]<-"Stenocercus pectinatus" # old name - might be something else
altSpecies[species=="Tropidurusto"]<-"Tropidurus torquatus" # 
altSpecies[species=="Sauresia"]<-"Celestus agasepsoides" # old name
altSpecies[species=="Wetmorena"]<-"Celestus haetianus" # i think this is exact
altSpecies[species=="Pelamis"]<-"Hydrophis curtus" # 
altSpecies[species=="Python"]<-"Python molurus" # 
altSpecies[species=="Zaocys"]<-"Ptyas fusca" # 
altSpecies[species=="Aporosaura"]<-"Meroles anchietae" # exact
altSpecies[species=="Meroles"]<-"Meroles ctenodactylus" # choosing different species
altSpecies[species=="Cabrita"]<-"Ophisops leschenaultii" # old genus
altSpecies[species=="Lacertalep"]<-"Timon lepidus" # must be correct, a real trick
altSpecies[species=="Platyplacopu"]<-"Takydromus intermedius" # 
altSpecies[species=="Anops"]<-"Amphisbaena acrobeles" # 
altSpecies[species=="Leposternon"]<-"Amphisbaena cerradensis" # 
altSpecies[species=="Pantodactylu"]<-"Cercosaura quadrilineata" # 
altSpecies[species=="Prionodactyl"]<-"Cercosaura argulus" # 
altSpecies[species=="Anotis"]<-"Nannoscincus gracilis" # could be other things but all lygosomines
altSpecies[species=="Apterygodon"]<-"Dasia vittata" # only hit
altSpecies[species=="Riopa"]<-"Ablepharus pannonicus" # honestly no idea - tons of things
altSpecies[species=="Siaphos"]<-"Saiphos equalis" # 
altSpecies[species=="Trachydosaur"]<-"Tiliqua rugosa" # 
altSpecies[species=="Fitzsimonsia"]<-"Typhlacontias brevipes" # I think - locality matches, seems right
altSpecies[species=="Typhlacontia"]<-"Typhlacontias gracilis" # Random different species in genus
altSpecies[species=="Neoseps"]<-"Plestiodon reynoldsi" # 
altSpecies[381]<-"Platysaurus capensis" # just choosing a different species
altSpecies[384]<-"Lygosoma albopunctata" # just choosing a different species

tnrs<-read.csv("squam_TNRS.csv")
family<-character(length=length(altSpecies))
subfam<-character(length=length(altSpecies))
spMatchForFamily<-character(length=length(altSpecies))
for(i in 1:length(altSpecies)) {
	mm<-agrep(altSpecies[i], tnrs[,2])
	if(length(mm)>0) {
		aa<-adist(altSpecies[i],tnrs[mm,2], partial=T)
		md<-which(aa[1,]==min(aa[1,]))
		spMatchForFamily[i]<-as.character(tnrs[mm,2][md[1]])
		family[i]<-as.character(tnrs[mm,3][md[1]])
		subfam[i]<-as.character(tnrs[mm,4][md[1]])
	}
}

cbind(species, altSpecies, spMatchForFamily, subfam, family)


pyronTree<-read.tree("squam_shl_dates (1).tre")

# genus tree - incomplete - DO NOT USE THIS FOR ANALYSES
genus<-unlist(strsplit(spMatchForFamily, " "))[1:length(spMatchForFamily)*2-1]

treeGenera<-unlist(strsplit(pyronTree$tip.label, "_"))[1:length(pyronTree$tip.label)*2-1]
toDrop<-rep(T, length(treeGenera))
for(i in 1:length(species)) {
	ww<-which(treeGenera==genus[i])
	toDrop[ww[1]]<-F
}
genusTree<-drop.tip(pyronTree, pyronTree$tip.label[toDrop])
genusTree$tip.label<-unlist(strsplit(genusTree$tip.label, "_"))[1:length(genusTree$tip.label)*2-1]
length(unique(genus))
length(genusTree$tip.label)

# genus tree is missing the following 16 genera:
genus[!(genus %in% genusTree$tip.label)]

plot(genusTree)

#subfamily tree
sfl<-levels(as.factor(subfam))
translation<-matrix(nrow=length(sfl), ncol=2)
toDrop<-rep(T, length(treeGenera))
for(i in 1:length(sfl)) {
	# get species in subfamily
	ww<-which(subfam==sfl[i])

	goodGenera<-genus[ww] %in% genusTree$tip.label
	ww<-ww[goodGenera]
	
	# choose one
	if(length(ww)==1) {
		ss<-spMatchForFamily[ww]
	} else {
		ss<-spMatchForFamily[sample(ww, 1)]
	}
	# pull genus name (to make matching work)
	sg<-strsplit(ss, " ")[[1]][1]
	
	if(!(sg %in% treeGenera))
		cat("ERROR!", i, sg)
		
	ww2<-which(treeGenera==sg)
	
	if(length(ww2)==1) {
		tipToSave<-ww2
	} else {	
		tipToSave<-sample(ww2, 1)
	}
	
	toDrop[tipToSave]<-F
	translation[i,]<-c(pyronTree$tip.label[tipToSave], sfl[i])
}

subfamTree<-drop.tip(pyronTree, pyronTree$tip.label[toDrop])

mm<-match(subfamTree$tip.label, translation[,1])
subfamTree$tip.label<-translation[mm,2]

pdf("subfamTree.pdf", width=10, height=20)
plot(subfamTree)
dev.off()


#family tree
fl<-levels(as.factor(family))
translation2<-matrix(nrow=length(fl), ncol=2)
toDrop<-rep(T, length(treeGenera))
for(i in 1:length(fl)) {
	# get species in subfamily
	ww<-which(family==fl[i])

	goodGenera<-genus[ww] %in% genusTree$tip.label
	ww<-ww[goodGenera]
	
	# choose one
	if(length(ww)==1) {
		ss<-spMatchForFamily[ww]
	} else {
		ss<-spMatchForFamily[sample(ww, 1)]
	}
	# pull genus name (to make matching work)
	sg<-strsplit(ss, " ")[[1]][1]
	
	if(!(sg %in% treeGenera))
		cat("ERROR!", i, sg)
		
	ww2<-which(treeGenera==sg)
	if(length(ww2)==1) {
		tipToSave<-ww2
	} else {	
		tipToSave<-sample(ww2, 1)
	}
	toDrop[tipToSave]<-F
	translation2[i,]<-c(pyronTree$tip.label[tipToSave], fl[i])
}

famTree<-drop.tip(pyronTree, pyronTree$tip.label[toDrop])

mm<-match(famTree$tip.label, translation2[,1])
famTree$tip.label<-translation2[mm,2]

pdf("famTree.pdf", width=10, height=20)
plot(famTree)
dev.off()

write.csv(cbind(1:59, levels(as.factor(family))), "familyTable.csv")

# finally, make the "clade tree"
# these are the groups I think we should use

cladeTrans<-read.csv("subfamtrans.csv")
clade<-character(length=length(subfam))
for(i in 1:length(subfam)) {
	rr<-which(cladeTrans[,1]==subfam[i])
	clade[i]<-as.character(cladeTrans[rr,2])
}

table(clade)

#clade tree
cl<-levels(as.factor(clade))
translation3<-matrix(nrow=length(cl), ncol=2)
toDrop<-rep(T, length(treeGenera))

for(i in 1:length(cl)) {
	# get species in clade
	ww<-which(clade==cl[i])

	goodGenera<-genus[ww] %in% genusTree$tip.label
	ww<-ww[goodGenera]
	
	# choose one
	if(length(ww)==1) {
		ss<-spMatchForFamily[ww]
	} else {
		ss<-spMatchForFamily[sample(ww, 1)]
	}
	# pull genus name (to make matching work)
	sg<-strsplit(ss, " ")[[1]][1]
	
	if(!(sg %in% treeGenera))
		cat("ERROR!", i, sg)
		
	ww2<-which(treeGenera==sg)
	if(length(ww2)==1) {
		tipToSave<-ww2
	} else {	
		tipToSave<-sample(ww2, 1)
	}
	toDrop[tipToSave]<-F
	translation3[i,]<-c(pyronTree$tip.label[tipToSave], cl[i])
}

cladeTree<-drop.tip(pyronTree, pyronTree$tip.label[toDrop])

mm<-match(cladeTree$tip.label, translation3[,1])
cladeTree$tip.label<-translation3[mm,2]

pdf("cladeTree.pdf", width=10, height=20)
plot(cladeTree)
dev.off()


# Table for SI
df<-data.frame(species, genus, clade, family, subfam)
oo<-order(df[,3], df[,2], df[,1])
colnames(df)<-c("Jonathan name", "Genus match (may not be exact)", "Clade", "Family", "Subfamily")
write.csv(df[oo,], "SortedSpeciesAndClade_forSI.csv")


# analyses by clade
## family mean trait values
cladeMeans<-aggregate(sqMorphContinuous, by=list(clade), FUN=mean, na.rm=T)
# family sample size
cladeN<-table(clade)

# get diversities from Pyron
nn<-dim(tnrs)[1]
allspeciesclade<-character(length=nn)
for(i in 1:nn) {
	rr<-which(as.character(cladeTrans[,1])==as.character(tnrs[i,4]))
	if(length(rr)>0) {
		allspeciesclade[i]<-as.character(cladeTrans[rr,2])
	} else {
		cat(as.character(tnrs[i,4]), "\n")
	}	
}
cladeDiv<-table(allspeciesclade)[-1]

# Create Table 1: Clade, diversity, age, sample size
# diversity
df<-data.frame(cladeDiv)


tipNodes<-cladeTree$edge[,2]<=length(cladeTree$tip.label)
tipAge<-cladeTree$edge.length[tipNodes]
names(tipAge)<-cladeTree$tip.label[cladeTree$edge[tipNodes,2]]

mm<-match(names(cladeDiv), names(tipAge))
tipAge[mm]
df[,2]<-tipAge[mm]

plot(df[,1:2], log="x")

df[,3]<-cladeN

colnames(df)<-c("Diversity", "Age", "SampleSize")

write.csv(df, "Table1_cladeAgeDiversity.csv")

plot(df[,2], df[,1], log="y", xlab="Age (my)", ylab="Diversity", pch=19)
summary(lm(log(df[,1])~df[,2]))
summary(lm(df[,3]~df[,1]))


pdf(file="cladeTreeTri.pdf", width=8, height=20)

nt<-read.tree(text=write.tree(ladderize(cladeTree)))
lct<-nt
tips<-nt$edge[,2]<=length(nt$tip.label)
nt$edge.length[tips]<-0.005

layout(matrix(c(1,1, 2,2, 2), nrow=1))
par(xpd=NA)
palette(rainbow(100, start=0.7, end=1/3))

plot(nt, show.tip.label=F)
sl<-c(1, 10, 100, 1000)
legend(0.02,7,fill=log(sl+1)*10, legend=sl, title="Taxa", cex=1.3)

totalDepth<-max(branching.times(cladeTree))
#for(xx in c(0, 0.05, 0.1, 0.15)) {
#	xc<-totalDepth-xx
#	lines(c(xc, xc), c(0, -0.2))
#	text(xc, -1, as.character(xx), cex=1.5)
#}
#lines(c(totalDepth-0, totalDepth-1.5), c(0, 0))


#text(0.07, -2, "MYA", cex=1.5)

ms<-max(cladeDiv)
for(i in 1:length(nt$tip.label))
{
	tt<-nt$tip.label[i]
	ns<-cladeDiv[tt]
	bn<-which(nt$edge[,2]==i)
	oldBranch<-lct$edge.length[bn]
	xb<-totalDepth-oldBranch+0.005
	polygon(c(xb, totalDepth, totalDepth, xb), c(i, i+0.2, i-0.2, i), lwd=2, col=log(ns+1)*10)
	#lines(c(xb, totalDepth), c(i, i+0.2), lwd=3)
	#lines(c(xb, totalDepth), c(i, i-0.2), lwd=3)
	#lines(c(totalDepth, totalDepth), c(i-0.2, i+0.2), lwd=3)

}
plot.new()

for(i in 1:length(nt$tip.label))
{
	tt<-nt$tip.label[i]
	ns<-cladeDiv[tt]
	bn<-which(nt$edge[,2]==i)
	yc<-(i-1)/(length(nt$tip.label)-1)
	text(0.23,yc, lct$tip.label[i], adj=c(0,0.5), cex=1.3)

}
dev.off()

zeros<-sqMorphContinuous[,]==0

varsWithZeros<-which(apply(zeros, 2, sum, na.rm=T)!=0)

# All species, limbs excluded, size included
# we can take out size by just chopping off PC1

sqMorph2<-log(sqMorphContinuous[,])
sqMorph2[,varsWithZeros]<-log(sqMorphContinuous[,varsWithZeros]+1)


sm2<-scale(sqMorph2[,1:12])
pcls<-pca(sm2, nPcs=5, method="nipals")
summary(pcls)
#slplot(pcls)

ps3<-scores(pcls)

write.csv(loadings(pcls), "loadings_with_size.csv")


makeHullPlot<-function(x, y, clade, cladesToHighlight, pch=19, cex=1, lineColor=rgb(0, 0, 0, .99), ...) {
plot(x, y, pch=pch, cex=cex, col=rgb(190/255, 190/255, 190/255, .99), ...)
ll<-levels(as.factor(clade))
for(i in cladesToHighlight) {
	ok<-clade==ll[i]
	pp<-cbind(x[ok], y[ok])
	points(pp, cex=cex, pch=pch, col=rgb(0, 0, 0, .99))
	if(sum(ok)>2) {
		hh<-convhulln(pp)
		for(j in 1:nrow(hh))
			lines(x=c(pp[hh[j,1],1], pp[hh[j,2],1]), y=c(pp[hh[j,1],2], pp[hh[j,2],2]), col=lineColor)
		
		}
	if(sum(ok)==2)
		lines(x=pp[1:2,1], y=pp[1:2,2], col=lineColor)	
		
	if(sum(ok)==1)
		points(x=pp[1], y=pp[2],col=lineColor, cex=3)
	}
}	


#3 x 3 = 9 + 43 = 52
#13 x 4

row1<-c(1, 1, 1, 2:11)
row2<-c(1, 1, 1, 12:21)
row3<-c(1, 1, 1, 22:31)
row4<-32:44

plotMatrix<-rbind(row1, row2, row3, row4)

pdf("morphoRaw_withSize.pdf", width=16, height=8)

layout(plotMatrix)
plot(ps3[,1:2], pch=19, cex=2, col=rgb(0, 0, 0, .99), axes=T, frame=F)
par(mar=c(0, 0, 0, 0))
for(i in 1:43) {
	makeHullPlot(ps3[,1], ps3[,2], clade, cladesToHighlight=i, axes=F, frame=F)
	text(5, -2, as.character(i), cex=2)
	}
dev.off()	

pdf("morphoRaw_noSize.pdf", width=16, height=8)

layout(plotMatrix)
plot(ps3[,2:3], pch=19, cex=2, col=rgb(0, 0, 0, .99), axes=T, frame=F)
par(mar=c(0, 0, 0, 0))
for(i in 1:43) {
	makeHullPlot(ps3[,2], ps3[,3], clade, cladesToHighlight=i, axes=F, frame=F)
	text(2, -2.3, as.character(i), cex=2)
	}
dev.off()	
	
	

# justify OU models?

cm<-aggregate(ps3, by=list(clade), FUN=mean)
cladeMeans<-cm[,-1]
rownames(cladeMeans)<-cm[,1]

fcbm<-fitContinuous(cladeTree, cladeMeans)
fcou<-fitContinuous(cladeTree, cladeMeans, model="OU")

aicDiff<-numeric(5)
alpha<-numeric(5)
for(i in 1:5) {
	aicDiff[i]<-fcbm[[i]]$opt$aic-fcou[[i]]$opt$aic
	alpha[i]<-fcou[[i]]$opt$alpha
}

aicDiff
alpha
max(alpha)
# OU is never supported, and maxes out at 0.01
# lets try reanalysing with alpha = 0.01, 0.1
	
	
# are clades in distinct regions of morphospace? MANOVA

# default analysis


mm3<-manova(ps3~as.factor(clade))
summary(mm3, test="Wilks")
# significant - and...

rndWilks<-numeric(1000)
for(i in 1:1000) {
	xx<-createRndData()
	
	ss<-strsplit(rownames(xx),split="_")
	rcl<-character(length=length(ss))
	for(j in 1:length(ss)) rcl[j]<-ss[[j]][1]
	mr<-manova(xx~as.factor(rcl))
	sm<-summary(mr, test="Wilks")
	rndWilks[i]<-sm[[4]][1,2]
	
	cat(i, " done\n")
	}
	
hist(rndWilks)	

actSumm<-summary(mm3, test="Wilks")
realWilks<-actSumm[[4]][1,2]


(sum(rndWilks <= realWilks)+1)/(length(rndWilks)+1)


# MORE clustered than Brownian motion

hist(rndWilks)
arrows(realWilks, 10, realWilks, 0)



rndWilksOU1<-numeric(1000)
for(i in 1:1000) {
	xx<-createRndData(alpha=0.01)
	
	ss<-strsplit(rownames(xx),split="_")
	rcl<-character(length=length(ss))
	for(j in 1:length(ss)) rcl[j]<-ss[[j]][1]
	mr<-manova(xx~as.factor(rcl))
	sm<-summary(mr, test="Wilks")
	rndWilksOU1[i]<-sm[[4]][1,2]
	
	cat(i, " done\n")
	}
	
hist(rndWilksOU1)	

(sum(rndWilksOU1 <= realWilks)+1)/(length(rndWilksOU1)+1)



rndWilksOU2<-numeric(1000)
for(i in 1:1000) {
	xx<-createRndData(alpha=0.1)
	
	ss<-strsplit(rownames(xx),split="_")
	rcl<-character(length=length(ss))
	for(j in 1:length(ss)) rcl[j]<-ss[[j]][1]
	mr<-manova(xx~as.factor(rcl))
	sm<-summary(mr, test="Wilks")
	rndWilksOU2[i]<-sm[[4]][1,2]
	
	cat(i, " done\n")
	}
	
hist(rndWilksOU2)	

(sum(rndWilksOU2 <= realWilks)+1)/(length(rndWilksOU2)+1)

# significant in all cases

# Without size.



mm3noSize<-manova(ps3[,-1]~as.factor(clade))
summary(mm3noSize, test="Wilks")
# significant - and...

rndWilks<-numeric(1000)
for(i in 1:1000) {
	xx<-createRndData(isSize=F)
	
	ss<-strsplit(rownames(xx),split="_")
	rcl<-character(length=length(ss))
	for(j in 1:length(ss)) rcl[j]<-ss[[j]][1]
	mr<-manova(xx~as.factor(rcl))
	sm<-summary(mr, test="Wilks")
	rndWilks[i]<-sm[[4]][1,2]
	
	cat(i, " done\n")
	}
	
hist(rndWilks)	

actSumm<-summary(mm3noSize, test="Wilks")
realWilks<-actSumm[[4]][1,2]


(sum(rndWilks <= realWilks)+1)/(length(rndWilks)+1)


# MORE clustered than Brownian motion

hist(rndWilks)
arrows(realWilks, 10, realWilks, 0)



rndWilksOU1<-numeric(1000)
for(i in 1:1000) {
	xx<-createRndData(alpha=0.01, isSize=F)
	
	ss<-strsplit(rownames(xx),split="_")
	rcl<-character(length=length(ss))
	for(j in 1:length(ss)) rcl[j]<-ss[[j]][1]
	mr<-manova(xx~as.factor(rcl))
	sm<-summary(mr, test="Wilks")
	rndWilksOU1[i]<-sm[[4]][1,2]
	
	cat(i, " done\n")
	}
	
hist(rndWilksOU1)	

(sum(rndWilksOU1 <= realWilks)+1)/(length(rndWilksOU1)+1)



rndWilksOU2<-numeric(1000)
for(i in 1:1000) {
	xx<-createRndData(alpha=0.1, isSize=F)
	
	ss<-strsplit(rownames(xx),split="_")
	rcl<-character(length=length(ss))
	for(j in 1:length(ss)) rcl[j]<-ss[[j]][1]
	mr<-manova(xx~as.factor(rcl))
	sm<-summary(mr, test="Wilks")
	rndWilksOU2[i]<-sm[[4]][1,2]
	
	cat(i, " done\n")
	}
	
hist(rndWilksOU2)	

(sum(rndWilksOU2 <= realWilks)+1)/(length(rndWilksOU2)+1)





# discriminant analysis



singlets<-which(table(clade)==1)
ok<-!(clade %in% names(singlets))

nClades<-nlevels(factor(clade[ok]))

prior<-rep(1/nClades, nClades)
names(prior)<-levels(factor(clade[ok]))

lll<-lda(ps3[ok,], factor(clade[ok]), CV=T, prior=prior)
ttt<-table(lll$class, clade[ok])

sum(diag(ttt))
sum(diag(ttt))/sum(ok)

mm<-matrix(rep(colSums(ttt), nClades), nrow= nClades, ncol= nClades, byrow=T)
propTTT<-ttt/mm

lll2<-lda(ps3[ok,-1], factor(clade[ok]), CV=T, prior=prior)
ttt2<-table(lll2$class, clade[ok])

sum(diag(ttt2))
sum(diag(ttt2))/sum(ok)

mm2<-matrix(rep(colSums(ttt2), nClades), nrow= nClades, ncol= nClades, byrow=T)
propTTT2<-ttt2/mm2


wrongOrRight<-lll$class== clade[ok]
wrongOrRight2<-lll2$class== clade[ok]

# get an order that "looks nice"
pp<-t(100-(propTTT)*100+1)
tt<-nj(pp)
tt2<-reorder(tt)
nn<-length(tt2$tip.label)
oo<-tt2$edge[,2][tt$edge[,2]<=nn]
newOrder<-tt2$tip.label[oo]

np<-propTTT[newOrder, newOrder]

# output results for Jonathan

write.csv(ttt, file="correctByClade_size.csv")
write.csv(ttt[newOrder, newOrder], file="correctByClade_size_reorder.csv")

pdf("correctByCladeHeatMap_size.pdf", width=5, height=5)
palette(gray(1:101/101))


plot("", xlim=c(1, nClades), ylim=c(1, nClades), axes=F, xlab="Actual", ylab="Inferred")

proportions<-t(100-(np)*100+1)

for(i in 1: nClades)
	for(j in 1: nClades)
		points(i, nClades-j+1, col= proportions[i,j], pch=15, cex=0.8)
		
for(i in 1: nClades) points(i, nClades-i+1, pch=22, cex=1, col=1)

axis(1, at=1: nClades, labels= newOrder, cex.axis=0.3, las=2)
axis(2, at=1: nClades, labels= rev(newOrder), cex.axis=0.3, las=2)
legend(17,5, fill=100-c(0, 25, 50, 75, 100)+1, legend=c(0, 0.25, 0.5, 0.75, 1), cex=0.4)

dev.off()


# get an order that "looks nice"
pp<-t(100-(propTTT2)*100+1)
tt<-nj(pp)
tt2<-reorder(tt)
nn<-length(tt2$tip.label)
oo<-tt2$edge[,2][tt$edge[,2]<=nn]
newOrder2<-tt2$tip.label[oo]

np2<-propTTT2[newOrder2, newOrder2]

pdf("correctByCladeHeatMap_nosize.pdf", width=5, height=5)
palette(gray(1:101/101))


plot("", xlim=c(1, nClades), ylim=c(1, nClades), axes=F, xlab="Actual", ylab="Inferred")

proportions<-t(100-(np2)*100+1)

for(i in 1: nClades)
	for(j in 1: nClades)
		points(i, nClades-j+1, col= proportions[i,j], pch=15, cex=0.8)
		
for(i in 1: nClades) points(i, nClades-i+1, pch=22, cex=1, col=1)

axis(1, at=1: nClades, labels= newOrder2, cex.axis=0.3, las=2)
axis(2, at=1: nClades, labels=rev(newOrder2), cex.axis=0.3, las=2)
legend(17,5, fill=100-c(0, 25, 50, 75, 100)+1, legend=c(0, 0.25, 0.5, 0.75, 1), cex=0.4)

dev.off()

write.csv(ttt2, file="correctByClade_nosize.csv")
write.csv(ttt2[newOrder2, newOrder2], file="correctByClade_nosize_reorder.csv")



# with size
# how many were misclassified into each clade
x1a<-rowSums(ttt)-diag(ttt)
sort(x1a)


# how many from each clade were misclassified 
x2<-colSums(ttt)-diag(ttt)
# all wrong
which(x2==colSums(ttt))
# all right
which(x2==0)

sort(x2)

temp<-ttt2
diag(temp)<-0
max(temp)
which(temp==8)

# no size
# how many were misclassified into each clade
x1b<-rowSums(ttt2)-diag(ttt2)
sort(x1b)

# how many from each clade were misclassified 
x2<-colSums(ttt2)-diag(ttt2)
# all wrong
which(x2==colSums(ttt2))
# all right
which(x2==0)

sort(x2)


# k-means analysis

kk3<-list()
ss3<-numeric(20)
for(i in 1:20) {
	kk3[[i]]<-kmeans(ps3, centers=i, nstart=5)
	ss3[i]<-kk3[[i]]$betweenss/kk3[[i]]$totss
}

plot(ss3)



kk4<-list()
ss4<-numeric(20)
for(i in 1:20) {
	kk4[[i]]<-kmeans(ps3[,-1], centers=i, nstart=5)
	ss4[i]<-kk4[[i]]$betweenss/kk4[[i]]$totss
}

plot(ss4)

# choose k = 6 by elbow criterion


pdf("k-means_elbow.pdf")
layout(matrix(1:2, nrow=1))

plot(ss3, type="l", xlab="k", ylab="Variance explained by groups")
points(ss3, pch=19)
arrows(6, ss3[6]+0.1,6, ss3[6]+0.02, length=0.1)

plot(ss4, type="l", xlab="k", ylab="Variance explained by groups")
points(ss4, pch=19)
arrows(6, ss4[6]+0.1,6, ss4[6]+0.02, length=0.1)

dev.off()


kmt3<-table(kk3[[6]]$cluster, clade)
prop3<-kmt3/rowSums(kmt3)

kmt4<-table(kk4[[6]]$cluster, clade)
prop4<-kmt4/rowSums(kmt4)

limbs<-sqMorphContinuous[,13]!=0

kleg<-as.data.frame(table(clade, limbs, kk3[[6]]$cluster))
write.csv(kleg[kleg[,4]!=0,], file="kmeans_clade_leg_size.csv")

kleg<-as.data.frame(table(clade, limbs, kk4[[6]]$cluster))
write.csv(kleg[kleg[,4]!=0,], file="kmeans_clade_leg_nosize.csv")

cor(t(kmt3), t(kmt4))

cor(t(prop3), t(prop4))

pdf("k-means_byClade_size.pdf", width=10, height=10)

x<-rep(1:6, each=43)
y<-rep(1:43, times=6)
circleSize3<-matrix(t(kmt3), nrow=1)[1,]*0.1
plot(x, y, cex= circleSize3, axes=F, xlab="K-means cluster", ylab="Clade", xlim=c(1, 7.5))
for(i in 1:43) lines(c(0.8, 6.2), c(i, i), lwd=0.5, col="grey")
points(x, y, cex= circleSize3, pch=21, bg="white", col="black")
axis(1, lwd=0, at=1:6)
axis(2, at=1:43, cex.axis=1, las=1, lwd=0)

leg<-c(1, 5, 10, 15, 20, 25)
points(rep(7, 6), 1:6*1.4, cex= leg*0.1, pch=21, bg="white", col="black")
text(rep(7.5, 6), 1:6*1.4,labels=as.character(leg))

dev.off()


pdf("k-means_byClade_nosize.pdf", width=10, height=10)

x<-rep(1:6, each=43)
y<-rep(1:43, times=6)
circleSize4<-matrix(t(kmt4), nrow=1)[1,]*0.1

plot(x, y, cex= circleSize4, axes=F, xlab="K-means cluster", ylab="Clade", xlim=c(1, 7.5))
for(i in 1:43) lines(c(0.8, 6.2), c(i, i), lwd=0.5, col="grey")
points(x, y, cex= circleSize4, pch=21, bg="white", col="black")
axis(1, lwd=0, at=1:6)
axis(2, at=1:43, cex.axis=1, las=1, lwd=0)
leg<-c(1, 5, 10, 15, 20, 25)
points(rep(7, 6), 1:6*1.4, cex= leg*0.1, pch=21, bg="white", col="black")
text(rep(7.5, 6), 1:6*1.4,labels=as.character(leg))

dev.off()


makeHulls<-function(x, y, clade, cladesToHighlight, pch=19, cex=1, lineColor=rgb(0, 0, 0, .99), ...) {
ll<-levels(as.factor(clade))
for(i in cladesToHighlight) {
	ok<-clade==ll[i]
	pp<-cbind(x[ok], y[ok])
	if(sum(ok)>2) {
		hh<-convhulln(pp)
		for(j in 1:nrow(hh))
			lines(x=c(pp[hh[j,1],1], pp[hh[j,2],1]), y=c(pp[hh[j,1],2], pp[hh[j,2],2]), col=lineColor)
		
		}
	if(sum(ok)==2)
		lines(x=pp[1:2,1], y=pp[1:2,2], col=lineColor)	
		
	if(sum(ok)==1)
		points(x=pp[1], y=pp[2],col=lineColor, cex=3)
	}
}	



row1<-c(1, 1, 2, 3)
row2<-c(1, 1, 4, 5)
row3<-c(1, 1, 6, 7)

plotMatrix<-rbind(row1, row2, row3)

cn<-as.numeric(as.factor(clade))

pdf("morphoKmeans_nosize.pdf", width=10, height=10)
plot(ps3[,2], ps3[,3],pch="",xlab="PC2", ylab="PC3")
text(ps3[,2], ps3[,3],cn)
makeHulls(ps3[,2], ps3[,3], kk4[[6]]$cluster, cladesToHighlight=1:6, axes=T, frame=T, xlab="PC2", ylab="PC3")
dev.off()	

#pdf("morphoKmeans_nosize.pdf", width=10, height=10)
#makeHullPlot(ps3[,2], ps3[,3], kk4[[6]]$cluster, cladesToHighlight=1:6, axes=T, frame=T, xlab="PC2", ylab="PC3")
#dev.off()	

pdf("morphoKmeans_size.pdf", width=10, height=10)
plot(ps3[,1], ps3[,2],pch="",xlab="PC1", ylab="PC2")
text(ps3[,1], ps3[,2],cn)
makeHulls(ps3[,1], ps3[,2], kk3[[6]]$cluster, cladesToHighlight=1:6, axes=T, frame=T, xlab="PC1", ylab="PC2")
dev.off()		
	

mmddSize<-manova(ps3~as.factor(kk3[[6]]$cluster))
summary(mmddSize, test="Wilks")
# significant - and...

rndWilks<-numeric(1000)
for(i in 1:1000) {
	xx<-createRndData(isSize=T)
	
	ss<-strsplit(rownames(xx),split="_")
	rcl<-character(length=length(ss))
	for(j in 1:length(ss)) rcl[j]<-ss[[j]][1]
	mr<-manova(xx~as.factor(kk3[[6]]$cluster))
	sm<-summary(mr, test="Wilks")
	rndWilks[i]<-sm[[4]][1,2]
	
	cat(i, " done\n")
	}
	
hist(rndWilks)	

actSumm<-summary(mmddSize, test="Wilks")
realWilks<-actSumm[[4]][1,2]


(sum(rndWilks <= realWilks)+1)/(length(rndWilks)+1)


# MORE clustered than Brownian motion

hist(rndWilks)
arrows(realWilks, 10, realWilks, 0)



rndWilksOU1<-numeric(1000)
for(i in 1:1000) {
	xx<-createRndData(alpha=0.01, isSize=T)
	
	ss<-strsplit(rownames(xx),split="_")
	rcl<-character(length=length(ss))
	for(j in 1:length(ss)) rcl[j]<-ss[[j]][1]
	mr<-manova(xx~as.factor(rcl))
	mr<-manova(xx~as.factor(kk3[[6]]$cluster))
	rndWilksOU1[i]<-sm[[4]][1,2]
	
	cat(i, " done\n")
	}
	
hist(rndWilksOU1)	

(sum(rndWilksOU1 <= realWilks)+1)/(length(rndWilksOU1)+1)



rndWilksOU2<-numeric(1000)
for(i in 1:1000) {
	xx<-createRndData(alpha=0.1, isSize=T)
	
	ss<-strsplit(rownames(xx),split="_")
	rcl<-character(length=length(ss))
	for(j in 1:length(ss)) rcl[j]<-ss[[j]][1]
	mr<-manova(xx~as.factor(kk3[[6]]$cluster))
	sm<-summary(mr, test="Wilks")
	rndWilksOU2[i]<-sm[[4]][1,2]
	
	cat(i, " done\n")
	}
	
hist(rndWilksOU2)	

(sum(rndWilksOU2 <= realWilks)+1)/(length(rndWilksOU2)+1)



# repeat without size
mmddNoSize<-manova(ps3[,-1]~as.factor(kk4[[6]]$cluster))
summary(mmddNoSize, test="Wilks")
# significant - and...

rndWilks<-numeric(1000)
for(i in 1:1000) {
	xx<-createRndData(isSize=F)
	
	ss<-strsplit(rownames(xx),split="_")
	rcl<-character(length=length(ss))
	for(j in 1:length(ss)) rcl[j]<-ss[[j]][1]
	mr<-manova(xx~as.factor(kk4[[6]]$cluster))
	sm<-summary(mr, test="Wilks")
	rndWilks[i]<-sm[[4]][1,2]
	
	cat(i, " done\n")
	}
	
hist(rndWilks)	

actSumm<-summary(mmddNoSize, test="Wilks")
realWilks<-actSumm[[4]][1,2]


(sum(rndWilks <= realWilks)+1)/(length(rndWilks)+1)



rndWilksOU1<-numeric(1000)
for(i in 1:1000) {
	xx<-createRndData(alpha=0.01, isSize=F)
	
	ss<-strsplit(rownames(xx),split="_")
	rcl<-character(length=length(ss))
	for(j in 1:length(ss)) rcl[j]<-ss[[j]][1]
	mr<-manova(xx~as.factor(kk4[[6]]$cluster))
	sm<-summary(mr, test="Wilks")

	rndWilksOU1[i]<-sm[[4]][1,2]
	
	cat(i, " done\n")
	}
	
hist(rndWilksOU1)	

(sum(rndWilksOU1 <= realWilks)+1)/(length(rndWilksOU1)+1)



rndWilksOU2<-numeric(1000)
for(i in 1:1000) {
	xx<-createRndData(alpha=0.1, isSize=F)
	
	ss<-strsplit(rownames(xx),split="_")
	rcl<-character(length=length(ss))
	for(j in 1:length(ss)) rcl[j]<-ss[[j]][1]
	mr<-manova(xx~as.factor(kk4[[6]]$cluster))
	sm<-summary(mr, test="Wilks")
	rndWilksOU2[i]<-sm[[4]][1,2]
	
	cat(i, " done\n")
	}
	
hist(rndWilksOU2)	

(sum(rndWilksOU2 <= realWilks)+1)/(length(rndWilksOU2)+1)


# but what happens if we do the df on the sims?


mmddSize<-manova(ps3~as.factor(kk3[[6]]$cluster))
summary(mmddSize, test="Wilks")
# significant - and...

rndWilks<-numeric(1000)
for(i in 1:1000) {
	xx<-createRndData(isSize=T)
	
	rkm<-kmeans(xx, centers=6, nstart=5)

	mr<-manova(xx~as.factor(rkm$cluster))
	sm<-summary(mr, test="Wilks")
	rndWilks[i]<-sm[[4]][1,2]
	
	cat(i, " done\n")
	}
	

actSumm<-summary(mmddSize, test="Wilks")
realWilks<-actSumm[[4]][1,2]


(sum(rndWilks <= realWilks)+1)/(length(rndWilks)+1)

	
mmddSize<-manova(ps3~as.factor(kk3[[6]]$cluster))
summary(mmddSize, test="Wilks")
# significant - and...

rndWilks<-numeric(1000)
for(i in 1:1000) {
	xx<-createRndData(isSize=T, alpha=0.01)
	
	rkm<-kmeans(xx, centers=6, nstart=5)

	mr<-manova(xx~as.factor(rkm$cluster))
	sm<-summary(mr, test="Wilks")
	rndWilks[i]<-sm[[4]][1,2]
	
	cat(i, " done\n")
	}
	

actSumm<-summary(mmddSize, test="Wilks")
realWilks<-actSumm[[4]][1,2]


(sum(rndWilks <= realWilks)+1)/(length(rndWilks)+1)


mmddSize<-manova(ps3~as.factor(kk3[[6]]$cluster))
summary(mmddSize, test="Wilks")
# significant - and...

rndWilks<-numeric(1000)
for(i in 1:1000) {
	xx<-createRndData(isSize=T, alpha=0.1)
	
	rkm<-kmeans(xx, centers=6, nstart=5)

	mr<-manova(xx~as.factor(rkm$cluster))
	sm<-summary(mr, test="Wilks")
	rndWilks[i]<-sm[[4]][1,2]
	
	cat(i, " done\n")
	}
	

actSumm<-summary(mmddSize, test="Wilks")
realWilks<-actSumm[[4]][1,2]


(sum(rndWilks <= realWilks)+1)/(length(rndWilks)+1)
	
	


mmddnoSize<-manova(ps3[,-1]~as.factor(kk4[[6]]$cluster))
summary(mmddnoSize, test="Wilks")
# significant - and...

rndWilks<-numeric(1000)
for(i in 1:1000) {
	xx<-createRndData(isSize=F)
	
	rkm<-kmeans(xx, centers=6, nstart=5)

	mr<-manova(xx~as.factor(rkm$cluster))
	sm<-summary(mr, test="Wilks")
	rndWilks[i]<-sm[[4]][1,2]
	
	cat(i, " done\n")
	}
	

actSumm<-summary(mmddnoSize, test="Wilks")
realWilks<-actSumm[[4]][1,2]


(sum(rndWilks <= realWilks)+1)/(length(rndWilks)+1)

	
rndWilks<-numeric(1000)
for(i in 1:1000) {
	xx<-createRndData(isSize=F, alpha=0.01)
	
	rkm<-kmeans(xx, centers=6, nstart=5)

	mr<-manova(xx~as.factor(rkm$cluster))
	sm<-summary(mr, test="Wilks")
	rndWilks[i]<-sm[[4]][1,2]
	
	cat(i, " done\n")
	}
	

actSumm<-summary(mmddnoSize, test="Wilks")
realWilks<-actSumm[[4]][1,2]


(sum(rndWilks <= realWilks)+1)/(length(rndWilks)+1)



rndWilks<-numeric(1000)
for(i in 1:1000) {
	xx<-createRndData(isSize=F, alpha=0.1)
	
	rkm<-kmeans(xx, centers=6, nstart=5)

	mr<-manova(xx~as.factor(rkm$cluster))
	sm<-summary(mr, test="Wilks")
	rndWilks[i]<-sm[[4]][1,2]
	
	cat(i, " done\n")
	}
	

actSumm<-summary(mmddnoSize, test="Wilks")
realWilks<-actSumm[[4]][1,2]


(sum(rndWilks <= realWilks)+1)/(length(rndWilks)+1)
		
		
		
chisq.test(kk3[[6]]$cluster, clade, simulate.p.value=T)
chisq.test(kk4[[6]]$cluster, clade, simulate.p.value=T)


# significant - and...

chts<-chisq.test(kk3[[6]]$cluster, clade)[[1]]


rndChi<-numeric(1000)
for(i in 1:1000) {
	xx<-createRndData(isSize=T)
	
	rkm<-kmeans(xx, centers=6, nstart=5)

	rr<-chisq.test(rkm$cluster, clade)[[1]]


	rndChi[i]<-rr
	
	cat(i, " done\n")
	}

(sum(rndChi >= chts)+1)/(length(rndChi)+1)

	

rndChi<-numeric(1000)
for(i in 1:1000) {
	xx<-createRndData(isSize=T, alpha=0.01)
	
	rkm<-kmeans(xx, centers=6, nstart=5)

	rr<-chisq.test(rkm$cluster, clade)[[1]]


	rndChi[i]<-rr
	
	cat(i, " done\n")
	}

(sum(rndChi >= chts)+1)/(length(rndChi)+1)



rndChi<-numeric(1000)
for(i in 1:1000) {
	xx<-createRndData(isSize=T, alpha=0.1)
	
	rkm<-kmeans(xx, centers=6, nstart=5)

	rr<-chisq.test(rkm$cluster, clade)[[1]]


	rndChi[i]<-rr
	
	cat(i, " done\n")
	}

(sum(rndChi >= chts)+1)/(length(rndChi)+1)



# significant - and...

chts<-chisq.test(kk4[[6]]$cluster, clade)[[1]]


rndChi<-numeric(1000)
for(i in 1:1000) {
	xx<-createRndData(isSize=F)
	
	rkm<-kmeans(xx, centers=6, nstart=5)

	rr<-chisq.test(rkm$cluster, clade)[[1]]


	rndChi[i]<-rr
	
	cat(i, " done\n")
	}

(sum(rndChi >= chts)+1)/(length(rndChi)+1)

	

rndChi<-numeric(1000)
for(i in 1:1000) {
	xx<-createRndData(isSize=F, alpha=0.01)
	
	rkm<-kmeans(xx, centers=6, nstart=5)

	rr<-chisq.test(rkm$cluster, clade)[[1]]


	rndChi[i]<-rr
	
	cat(i, " done\n")
	}

(sum(rndChi >= chts)+1)/(length(rndChi)+1)



rndChi<-numeric(1000)
for(i in 1:1000) {
	xx<-createRndData(isSize=F, alpha=0.1)
	
	rkm<-kmeans(xx, centers=6, nstart=5)

	rr<-chisq.test(rkm$cluster, clade)[[1]]


	rndChi[i]<-rr
	
	cat(i, " done\n")
	}

(sum(rndChi >= chts)+1)/(length(rndChi)+1)




	
limbs<-sqMorphContinuous[,13]!=0

table(kk3[[6]]$cluster, limbs)
table(kk4[[6]]$cluster, limbs)

plot(as.factor(kk3[[6]]$cluster), as.factor(limbs))
plot(as.factor(kk4[[6]]$cluster), as.factor(limbs))


kleg<-as.data.frame(table(clade, limbs, kk3[[6]]$cluster))
write.csv(kleg[kleg[,4]!=0,], file="kmeans_clade_leg_size.csv")

kleg<-as.data.frame(table(clade, limbs, kk4[[6]]$cluster))
write.csv(kleg[kleg[,4]!=0,], file="kmeans_clade_leg_nosize.csv")


# effect of limblessness

# Can we differentiate between limbed and limbless forms, even without looking at the limbs?
limbs<-sqMorphContinuous[,13]!=0



mm3<-manova(ps3~as.factor(clade)*as.factor(limbs))
summary(mm3, test="Wilks")

mm4<-manova(ps3[,-1]~as.factor(clade)*as.factor(limbs))
summary(mm4, test="Wilks")


newdata<-data.frame(as.factor(rep(levels(as.factor(clade)), 2)), as.factor(rep(c(T, F), each=43)))
colnames(newdata)<-c("clade", "limbs")

predict(mm3, newdata)->mmp3
predict(mm4, newdata)->mmp4

tl<-table(clade, limbs)



tl<-table(clade, limbs)
ww<-which(tl[,1]!=0 & tl[,2]!=0)

pdf("limbLoss_size.pdf", width=10, height=10)

pointColor<-character(length=length(limbs))
pointColor[limbs]<-"grey"
pointColor[!limbs]<-"white"

plot(ps3[,1:2], pch=21, col="black", bg= pointColor, cex=1)

ty=c(1, 5, 6, 1, 5, 6)
wd=c(1, 1, 1, 3, 3, 3)
for(j in 1:6) {
	i<-ww[j]
	rr<-which(newdata[,1]==names(i))
	arrows(x0=mmp3[rr,1][1], x1=mmp3[rr,1][2], y0=mmp3[rr,2][1], y1=mmp3[rr,2][2], lwd=wd[j], lty=ty[j],length=0.1)
	}
legend("bottomright", lwd=wd, lty=ty, legend=c("Clade 2", "Clade 3", "Cordylidae", "Gymnophthalmidae", "Lygosominae", "Scincinae"), cex=1.5)
dev.off()

pdf("limbLoss_nosize.pdf", width=10, height=10)

plot(ps3[,2:3], pch=21, col="black", bg= pointColor, cex=1)
for(j in 1:6) {
	i<-ww[j]
	rr<-which(newdata[,1]==names(i))
	arrows(x0=mmp4[rr,1][1], x1=mmp4[rr,1][2], y0=mmp4[rr,2][1], y1=mmp4[rr,2][2], , lwd=wd[j], lty=ty[j],length=0.1)
	}
legend("bottomleft", lwd=wd, lty=ty, legend=c("Clade 2", "Clade 3", "Cordylidae", "Gymnophthalmidae", "Lygosominae", "Scincinae"), cex=1.5)

dev.off()


# we get a certain amount of clustering with k=6
# is this better than chance?






cladeMeans<-aggregate(ps4, by=list(clade), FUN=mean, na.rm=T)

distOrigin<-function(x) sqrt(sum(x^2))
cladeDist<-apply(cladeMeans[,-1], 1, distOrigin)

summary(lm(log(diversity)~cladeDist))

summary(lm(cladeDisparity~cladeDist))
summary(lm(age~cladeDist))



pdf("nature/explainDiversity.pdf", width=12, height=4)
layout(matrix(1:3, nrow=1))
plot(diversity~age, log="y", pch=19, xlab="Age (mya)", ylab="Clade diversity")
plot(cladeDisparity~age, pch=19, xlab="Age (mya)", ylab="Clade disparity")
plot(cladeDist, diversity, log="y", pch=19, xlab="Clade Distance", ylab="Clade diversity")
fm<-lm(log(cladeDiv[,3])~cladeDist)
ny<-exp(predict(fm))
nx<-sort(cladeDist)
ny<-ny[order(cladeDist)]
lines(nx, ny, lwd=2)
dev.off()




cladeMeans2<-aggregate(ps4, by=list(clade), FUN=mean, na.rm=T)

plot(cladeMeans2[,2], cladeDiv[,3], log="y")

ss<-smooth.spline(x=cladeMeans2[,2], y=log(cladeDiv[,3]), spar=0.9)
lines(ss$x, exp(ss$y))

plot(cladeMeans2[,3], cladeDiv[,3], log="y")

ss<-smooth.spline(x=cladeMeans2[,3], y=log(cladeDiv[,3]), spar=0.9)
lines(ss$x, exp(ss$y))

sqSlope<-numeric(11)
sqP<-numeric(11)

for(i in 2:12) {

	fm<-lm(log(cladeDiv[,3])~abs(cladeMeans2[,i]))
	sqSlope[i-1]<-fm[[1]][2]
	sqP[i-1]<-summary(fm)[[4]][2,4]
	
	}
	


cbind(sqSlope, sqP)

pdf("nature/centerPCA.pdf")
layout(matrix(1:12, nrow=4, byrow=T))
par(mar=c(0, 0, 0, 0))

for(i in 1:11) {
	plot(cladeMeans2[,i+1], cladeDiv[,3], log="y", pch=19, axes=F, frame=T)

ss<-smooth.spline(x=cladeMeans2[,i+1], y=log(cladeDiv[,3]), spar=1)
lines(ss$x, exp(ss$y), col="red")
	lines(x=c(0,0), y=c(0.01, 10000), lty=2)
	
	}

dev.off()

sqQR1<-matrix(nrow=5, ncol=4)
sqQR2<-matrix(nrow=5, ncol=4)


for(i in 2:6) {
	fm<-lm(log(cladeDiv[,3])~ cladeMeans[,i]+I(cladeMeans[,i]^2))
	sqQR1[i-1,1:2]<-fm[[1]][2:3]
	sqQR1[i-1,3:4]<-summary(fm)[[4]][2:3,4]
	
	fm<-lm(log(cladeDiv[,3])~ cladeMeans2[,i]+I(cladeMeans2[,i]^2))
	sqQR2[i-1,1:2]<-fm[[1]][2:3]
	sqQR2[i-1,3:4]<-summary(fm)[[4]][2:3,4]
	}
	




