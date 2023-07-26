# QTL exercise using r/qtl & r/qtl2

# 1. Setup and load packages

Clean up a workplace first
```
rm(list=ls())
```

Install the r/qtl and r/qtl2 package. It takes ca. 10 min. to complete.
```
xxxxxxxxxxxxxxxx("qtl",repos="http://cran.us.r-project.org") # install the qtl package
xxxxxxxxxxxxxxxx("qtl2",repos="http://cran.us.r-project.org") # install the qtl2 package
```

# 2. qtl exercise with RILs

Load the r/qtl package
```
xxxxxxxxx(qtl) # load the qtl package
```

Load genotypes and phenotypes as a "cross" object
```
colkas <- read.cross(format="csvs",dir="./",
                     genfile="ColKasFloweringGeno.csv",
                     phefile = "ColKasFloweringPheno.csv",
                     na.strings = c("-"), estimate.map=FALSE,
                     crosstype = "riself")

xxxxxx(colkas) # see summary of the cross object
totmar(colkas) # total no. of genetic markers

plotMissing(colkas) # check missing genotypes

par(mfcol=c(2,1)); par(mai=c(1,1,0.25,0.25)) # change the plot parameters
plotPheno(colkas, pheno.col=2) # Sweden days-to-bolting (SWDTF) # plot the results of Sweden days-to-bolting (SWDTF)
plotPheno(colkas, pheno.col=3) # Spain days-to-bolting (SPDTF) # plot the results of Spain days-to-bolting (SPDTF)
```

Estimate and plot genetic map
```
newmap <- est.map(colkas, error.prob=0.01)
colkas <- replace.map(colkas, newmap)
plotMap(colkas)
```

See the genotype and points of recombination
```
colkas <- calc.errorlod(colkas, error.prob=0.01)
plotGeno(colkas, chr=4)
```

Calculate genotype probabilities
```
colkas_genoprob <- calc.genoprob(colkas, step=1)
colkas_genoprob$geno$"4"$prob
```

QTL mapping for flowering time under Swedish climates
```
scanSWDTF <- scanone(colkas_genoprob, pheno.col=colkas$pheno$SWDTF,
                     method="hk") # hk = Haley-Knott regression
xxx(mfcol=c(1,2)); xxx(mai=c(1,1,0.25,0.25)) # change the plot parameters
xxxx(scanSWDTF); xxxx(scanSWDTF, chr=4) # plot the results
```

QTL mapping for flowering time under Spanish climates
```
scanSPDTF <- scanone(colkas_genoprob, pheno.col=colkas$pheno$SPDTF,
                     method="hk") # hk = Haley-Knott regression
xxx(mfcol=c(1,2)); xxx(mai=c(1,1,0.25,0.25)) # change the plot parameters
xxxx(scanSPDTF); xxxx(scanSPDTF, chr=4) # plot the results
```

# 3. qtl2 exercise with MAGIC lines

Load the qtl2 package
```
xxxxxxxx(qtl2) # load the qtl2 package
```

Download dataset from the r/qtl2 website
```
file <- paste0("https://raw.githubusercontent.com/rqtl/",
               "qtl2data/master/ArabMAGIC/arabmagic_tair9.zip")
magic <- read_cross2(file)
head(magic$pheno) # see phenotypes

xxxxxxx(magic) # see summary of the cross2 object
magic$gmap$"4"[1:50] # markers at the top of chr. 4
```

Visualize recombination
```
xxx(mfcol=c(1,2)); xxx(mai=c(1,1,0.25,0.25)) # change the plot parameters
plot_onegeno(magic$geno, map=magic$gmap, ind=3) # Position in Mbp
plot_onegeno(magic$geno, map=magic$gmap, ind=4) # use "ind=" to change individuals
```

Calculate genotype probabilities
```
map2 <- insert_pseudomarkers(magic$gmap, step=1)
magic_p <- calc_genoprob(magic, map=map2)
plot_genoprob(magic_p, map=map2, chr=4, ind=3) # use "ind=" to change individuals
```

QTL mapping for flowering time in MAGIC lines
```
res_scan2 <- scan1(magic_p, pheno=magic$pheno[,1])
xxx(mfcol=c(1,2)); xxx(mai=c(1,1,0.25,0.25)) # change the plot parameters
xxxx(res_scan2, map=map2); xxxx(res_scan2, map=map2, chr=4) # plot the results
```
