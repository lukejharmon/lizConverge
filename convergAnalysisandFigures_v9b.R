library(geometry)
library(geiger)
library(pcaMethods)
library(MASS)
library(TreeSim)

setwd("~/Documents/lizardMorphospace/analyses/")

# get trait data
sqMorphA<-read.csv("squamorph_v3.csv", na.strings=c("", " ", "."))





# Taxonomic name resolution

fullName<-as.character(sqMorphA[1:384, 4])
# The one Jonathan couldn't pin down:
noFullName<-which(is.na(fullName))
fullName[noFullName]<-as.character(sqMorphA[noFullName,1])

# delete the one that is repeated
repeated<-which(fullName=="Fitzsimonsia brevipes")
fullName<-fullName[-repeated]
sqMorphA <- sqMorphA[-repeated,]

# pull out just continuously distributed characters
sqMorphContinuous<-sqMorphA[1:383, 5:27]

#manual clean-up of a few names
# Subspecies designated; remove subspecies
# If the subspecies name matches a full name then elevate
# otherwise drop subspecies name

# drop
fullName[fullName=="Leiocephlus carinatus aquarius"]<-"Leiocephlus carinatus"
fullName[fullName=="Petrosaurus thalassinus thalassinus"]<-"Petrosaurus thalassinus"
fullName[fullName=="Sceloporus occidentalis bocourtii"]<-"Sceloporus occidentalis"
fullName[fullName=="Cylindrophis ruffus ruffus"]<-"Cylindrophis ruffus"
fullName[fullName=="Xenodon rabdocephalus rabdocephalus"]<-"Xenodon rabdocephalus"
fullName[fullName=="Agamodon anguliceps anguliceps"]<-"Agamodon anguliceps"
fullName[fullName=="Blanus strauchi aporus"]<-"Blanus strauchi"
fullName[fullName=="Monopeltis sphenorhynchus mauricei"]<-"Monopeltis sphenorhynchus"
fullName[fullName=="Tretioscincus bifasciatus kugleri"]<-"Tretioscincus bifasciatus"
fullName[fullName=="Ateuchosaurus pellopleurus browni"]<-"Ateuchosaurus pellopleurus"

# elevate
fullName[fullName=="Chamaesaura anguina tenuior"]<-"Chamaesaura tenuior"
fullName[fullName=="Pseudocordylus subviridis transvaalensis"]<-"Pseudocordylus transvaalensis"
fullName[fullName=="Goniurosaurus kuroiwae orientalis"]<-"Goniurosaurus orientalis" # ref: Burbrink table





# Species missing from taxonomy; substitute different species in same genus
fullName[fullName=="Bunopus blanfordi"]<-"Bunopus crassicauda"

# update taxonomy
fullName[fullName=="Coleodactylus amazonicus"]<-"Chatogekko amazonicus" # ref: reptile database
fullName[fullName=="Cosymbotus platyurus"]<-"Hemidactylus platyurus" # ref: reptile database
fullName[fullName=="Geckonia chazaliae"]<-"Tarentola chazaliae" # ref: reptile database
fullName[fullName=="Homonota horrida"]<-"Homonota fasciata" # ref: reptile database
fullName[fullName=="Microscalabotes vivittis"]<-"Lygodactylus bivittis" # ref: reptile database and Jonathan is a bad speller
fullName[fullName=="Pachydactylus bibronii"]<-"Chondrodactylus bibronii" # ref: reptile database
fullName[fullName=="Palmatogecko rangei"]<-"Pachydactylus rangei" # ref: reptile database
fullName[fullName=="Pseudogonatodes amazonicus"]<-"Pseudogonatodes guianensis" # ref: reptile database
fullName[fullName=="Teratolepis fasciata"]<-"Hemidactylus imbricatus" # ref: reptile database
fullName[fullName=="Hoplodactylus maculatus"]<-"Mokopirirakau nebulosus" # ref: reptile database
fullName[fullName=="Oedura rhombifer"]<-"Amalosia rhombifer" # ref: reptile database
fullName[fullName=="Bradypodion tenue"]<-"Kinyongia tenue" # ref: reptile database
fullName[fullName=="Chamaeleo jacksonii"]<-"Trioceros jacksonii" # ref: reptile database
fullName[fullName=="Chamaeleo rhinoceratus"]<-"Furcifer rhinoceratus" # ref: reptile database
fullName[fullName=="Agma bibronii"]<-"Agama impalearis" # ref: reptile database and spelling
fullName[fullName=="Physignathus lesueurii"]<-"Intellagama lesueurii" # ref: reptile database
fullName[fullName=="Aptycholaemus longicauda"]<-"Anisolepis longicauda" # ref: reptile database
fullName[fullName=="Phenacosaurus heterodermus"]<-"Anolis heterodermus" # ref: reptile database
fullName[fullName=="Chamaeleolis chamaeleonides"]<-"Anolis chamaeleonides" # ref: reptile database
fullName[fullName=="Chamaelinorops barbouri"]<-"Anolis barbouri" # ref: reptile database
fullName[fullName=="Ophryoessoides iridescens"]<-"Stenocercus iridescens" # ref: reptile database
fullName[fullName=="Anniella stebbinsi"]<-"Anniella pulchra" # ref: reptile database; I think this is a regression, actually, to an older name
fullName[fullName=="Ophisaurus apodus"]<-"Pseudopus apodus" # ref: reptile database
fullName[fullName=="Sauresia sepsoides"]<-"Celestus sepsoides" # ref: reptile database
fullName[fullName=="Wemorena haetiana haetiana"]<-"Celestus haetianus" # ref: reptile database
fullName[fullName=="Python reticulatus"]<-"Broghammerus reticulatus" # ref: reptile database, where this 
species is listed as Malayopython reticulatus; named here to match tnrs
fullName[fullName=="Zaocys luzonensis"]<-"Ptyas luzonensis" # ref: reptile database
fullName[fullName=="Tenuidactylus caspius"]<-"Cyrtopodion caspium" # ref: reptile database (regress to older name here)
fullName[fullName=="Tenuidactylus rohtasfortai"]<-"Cyrtopodion rohtasfortai" # ref: reptile database
fullName[fullName=="Vipera russelii"]<-"Daboia russelii" # ref: reptile database
fullName[fullName=="Oligodon analepticos"]<-"Oligodon ocellatus" # ref: reptile database
fullName[fullName=="Pelamis platura"]<-"Hydrophis platurus" # ref: reptile database
fullName[fullName=="Cabrita leschenaulti"]<-"Ophisops leschenaulti" # ref: reptile database
fullName[fullName=="Lacerta lepida"]<-"Timon lepidus" # ref: reptile database
fullName[fullName=="Platyplacopus dorsalis"]<-"Takydromus dorsalis" # ref: This is not a synonym in reptile database, but many other Takydromus species have Platypoacopus as a synonym so I am sure this is correct
fullName[fullName=="Anops kingii"]<-"Amphisbaena kingii" # ref: reptile database
fullName[fullName=="Leposternon scutigerum"]<-"Amphisbaena scutigerum" # ref: reptile database
fullName[fullName=="Crocodilurus lacertinus"]<-"Crocodilurus amazonicus" # ref: reptile database
fullName[fullName=="Gymnodactylus antillensis"]<-"Gonatodes antillensis" # ref: reptile database
fullName[fullName=="Gymnophthalmus rubricauda"]<-"Vanzosaura rubricauda" # ref: reptile database
fullName[fullName=="Neusticurus strangulatus"]<-"Potamites strangulatus" # ref: reptile database
fullName[fullName=="Pantodactylus schreibersii"]<-"Cercosaura schreibersii" # ref: reptile database
fullName[fullName=="Prionodactylus oshaughnessyi"]<-"Cercosaura argulus" # ref: reptile database
fullName[fullName=="Apterygodon vittatus"]<-"Dasia vittata" # ref: reptile database
fullName[fullName=="Cophoscincus durum"]<-"Cophoscincopus durus" # ref: reptile database
fullName[fullName=="Hemiergis maccoyi"]<-"Nannoscincus maccoyi" # ref: reptile database as Anepischetosia maccoyi, with nannoscincus as synonym
fullName[fullName=="Lygosoma afer"]<-"Mochlus afer" # ref: reptile database
fullName[fullName=="Marmorosphax euryotis"]<-"Celatiscincus euryotis" # ref: reptile database
fullName[fullName=="Riopa punctata"]<-"Lygosoma punctata" # ref: reptile database
fullName[fullName=="Sphenomorphus incertus"]<-"Scincella incerta" # ref: reptile database
fullName[fullName=="Trachydosaurus rugosus"]<-"Tiliqua rugosa" # ref: reptile database
fullName[fullName=="Androngo elongatus"]<-"Amphiglossus elongatus" # ref: reptile database
fullName[fullName=="Fitzsimonsia brevipes"]<-"Typhlacontias brevipes" # ref: reptile database
fullName[fullName=="Neoseps reynoldsi"]<-"Plestiodon reynoldsi" # ref: reptile database
fullName[fullName=="Sphenops delislei"]<-"Chalcides delislei" # ref: reptile database
fullName[fullName=="Typhlosaurus lineatus"]<-"Acontias rieppeli" # ref: reptile database; there are other possible mathces in this genus too
fullName[fullName=="Sphenops delislei"]<-"Chalcides delislei" # ref: reptile database
fullName[fullName=="Cordylus giganteus"]<-"Smaug giganteus" # ref: reptile database

# Resolve uncertainty
fullName[fullName=="Ptyas korros maybe??? I think probable"]<-"Ptyas korros" # OK I guess we can keep that
fullName[fullName=="probabily 207644 Thamnophis sirtalis"]<-"Thamnophis sirtalis" # OK I guess we can keep that
fullName[fullName=="Diplometopon"]<-"Diplometopon zarudnyi" # only species in genus


# species uncertain: random choice of common species in genus
fullName[fullName=="Crotaphythus"]<-"Crotaphytus collaris"
fullName[fullName=="Rhinophis"]<-"Rhinophis blythii"
fullName[fullName=="Asymblepharus"]<-"Asymblepharus alaicus"
fullName[fullName=="Egernia sp."]<-"Egernia epsisolus"
fullName[fullName=="Lygisaurus"]<-"Lygisaurus abscondita"



tnrs<-read.csv("squam_TNRS.csv")
family<-character(length=length(fullName))
subfam<-character(length=length(fullName))
spMatchForFamily<-character(length=length(fullName))
for(i in 1:length(fullName)) {
	mm<-agrep(fullName[i], tnrs[,2])
	if(length(mm)>0) {
		aa<-adist(fullName[i],tnrs[mm,2], partial=T)
		md<-which(aa[1,]==min(aa[1,]))
		spMatchForFamily[i]<-as.character(tnrs[mm,2][md[1]])
		family[i]<-as.character(tnrs[mm,3][md[1]])
		subfam[i]<-as.character(tnrs[mm,4][md[1]])
	}
}

newMatches<-cbind(fullName, spMatchForFamily, subfam, family)
badMatches<-which(fullName != spMatchForFamily)

# These are fuzzy matches that I have checked against reptile database and so on
# and I believe to be corredt
approvedSubstitutionList<-c(20, 48, 77, 96, 106, 109, 113, 127, 137, 138, 164, 172, 191, 192, 195, 216, 236, 239, 249, 252, 254, 257, 272, 300, 301, 331, 336, 341, 346, 363)

# filter these out
badMatches <- badMatches[!(badMatches %in% approvedSubstitutionList)]

length(badMatches)
# this should be 0 - then we're good.

length(fullName) # sample size



pyronTree<-read.tree("squam_shl_dates (1).tre")

# species matches
alteredPyronNames<-sub("_", " ", pyronTree$tip.label)
numExactMatch<-length(which(fullName %in% alteredPyronNames))


# genus tree - incomplete - DO NOT USE THIS FOR ANALYSES
genus<-unlist(strsplit(spMatchForFamily, " "))[1:length(spMatchForFamily)*2-1]

treeGenera<-unlist(strsplit(pyronTree$tip.label, "_"))[1:length(pyronTree$tip.label)*2-1]
toDrop<-rep(T, length(treeGenera))
for(i in 1:length(fullName)) {
	ww<-which(treeGenera==genus[i])
	toDrop[ww[1]]<-F
}
genusTree<-drop.tip(pyronTree, pyronTree$tip.label[toDrop])
genusTree$tip.label<-unlist(strsplit(genusTree$tip.label, "_"))[1:length(genusTree$tip.label)*2-1]
length(unique(genus))
length(genusTree$tip.label)

# genus tree is missing the following 20 genera:
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
	
	# one species is named differently in the tree
	# and in the taxonomy: Epicta albifrons = Leptotyphlops albifrons
	
	if(i==54){
		ss<-"Leptotyphlops albifrons"
	} else if(length(ww)==1) {
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

# subfam tree is complete
subfam[!(subfam %in% subfamTree$tip.label)]


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
	if(i == 35) { # see subfam tree for notes
		ss<-"Leptotyphlops albifrons"
	} else if(length(ww)==1) {
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
	if(i == 29) {
		ss<-"Leptotyphlops albifrons"
	} else if(length(ww)==1) {
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


# Table for SI - this was heavily edited by hand to relate TNRS code above
df<-data.frame(fullName, genus, clade, family, subfam)
oo<-order(df[,3], df[,2], df[,1])
colnames(df)<-c("Full name in data", "Genus match (may not be exact)", "Clade", "Family", "Subfamily")
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


# sm2<-scale(sqMorph2[,1:12]) # uncomment for correlation matrix analysis
sm2<-sqMorph2[,1:12] # uncomment for covariance matrix analysis


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
#sort(x1b)

# how many from each clade were misclassified 
x2<-colSums(ttt2)-diag(ttt2)
# all wrong
which(x2==colSums(ttt2))
# all right
which(x2==0)

sort(x2)

# confusing pairs with size
offdiag<-ttt
diag(offdiag)<-0
sort(offdiag, decreasing=T)[1:3]

ww<-which(offdiag==max(offdiag), arr.ind=T)
c(rownames(ttt)[ww[1]], colnames(ttt)[ww[2]])

ww2<-which(offdiag==8, arr.ind=T)
c(rownames(ttt)[ww2[1]], colnames(ttt)[ww2[2]])

colSums(ttt)-diag(ttt)


# confusing pairs no size
offdiag<-ttt2
diag(offdiag)<-0
sort(offdiag, decreasing=T)[1:3]

ww<-which(offdiag==max(offdiag), arr.ind=T)
c(rownames(ttt)[ww[1]], colnames(ttt)[ww[2]])

ww2<-which(offdiag==7, arr.ind=T)
c(rownames(ttt)[ww2[1]], colnames(ttt)[ww2[2]])

colSums(ttt)-diag(ttt)

# attractors
sort(rowSums(ttt)-diag(ttt))
sort(rowSums(ttt2)-diag(ttt2))


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
	




