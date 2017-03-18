# By Elisabeth Meyers 2016

source("myFunctions.R") # custom functions
library(SNPediaR) # get the snp study summaries text
library(rsnps)    # get the user data
library(stringr)  # used for some string splitting
library(tm)       # for parsing
library(wordcloud2) # visualizing words

# snps of interest
rs.vec = list(rs1 = "rs4680", rs2 = "rs6323")

#------------------------------Get data

#Acquiring genotype data from Rsnps (only run once)
#rs.set1 = allgensnp(snp = rs.vec$rs1, df = TRUE)
#rs.set2 = allgensnp(snp = rs.vec$rs2, df = TRUE)
#write.table(rs.set1, file = "rs.set1.csv", row.names = TRUE, col.names = TRUE, sep=",")
#write.table(rs.set2, file = "rs.set2.csv", row.names = TRUE, col.names = TRUE, sep=",")

# read in user snp data that we wrote out above
set1.4680 = read.csv("rs.set1.csv")
set2.6323 = read.csv("rs.set2.csv")

#Getting phenotypic data from Rsnps (only run once)
#set1.phenotypes.dat = getPhenoData(rs.set1)
#set2.phenotypes.dat = getPhenoData(rs.set1)
#write.csv(set1.phenotypes.dat, "phenotypesRs1.csv")
#write.csv(set1.phenotypes.dat, "phenotypesRs2.csv")

# read in the phenotype data we wrote from above
set1.phenotypes.dat = read.csv("phenotypesRs1.csv",colClasses = c("NULL",NA))
set2.phenotypes.dat = read.csv("phenotypesRs2.csv",colClasses = c("NULL",NA))

## get counts and names of valid genotypes
set1.genotype.counts = (getGenoCounts(set1.4680))
set2.genotype.counts = (getGenoCounts(set2.6323))

## get text data from snpediaR and pass in genotype names to query the db
set1.snpedia.dat = getData(rs.vec$rs1, set1.genotype.counts)
set2.snpedia.dat = getData(rs.vec$rs2, set2.genotype.counts)

#-----------------------------Parsing & Cleanup

set1.snpedia.parsed = parseSnpediar(set1.snpedia.dat)
set2.snpedia.parsed = parseSnpediar(set2.snpedia.dat)

# snpedia cleanup
set1.snpedia.parsed = dataCleanup(set1.snpedia.parsed)
set2.snpedia.parsed = dataCleanup(set2.snpedia.parsed)

# phenotypes parsing
set1.rsnps.parsed = parsePheno(set1.phenotypes.dat)
set2.rsnps.parsed = parsePheno(set2.phenotypes.dat)

# phenotypes cleanup 
set1.rsnps.parsed = dataCleanup(set1.rsnps.parsed)
set2.rsnps.parsed = dataCleanup(set2.rsnps.parsed)

# change the column names
set1.rsnps.parsed = setColNames(set1.rsnps.parsed,"text")
set2.rsnps.parsed = setColNames(set2.rsnps.parsed,"text")
set1.snpedia.parsed = setColNames(set1.snpedia.parsed, "text")
set2.snpedia.parsed = setColNames(set2.snpedia.parsed, "text")

# ----------------------Frequency tables

# word counts for phenotype traits text
set1.rsnps.freq = table(set1.rsnps.parsed)
set2.rsnps.freq = table(set2.rsnps.parsed)

#sort the word counts tables for phenotypes traits
set1.rsnps.freq = as.data.frame(sort(set1.rsnps.freq, decreasing=TRUE))
set2.rsnps.freq = as.data.frame(sort(set2.rsnps.freq, decreasing=TRUE))

# word counts for snpediar study text
set1.snpedia.freq = table(set1.snpedia.parsed)
set2.snpedia.freq = table(set2.snpedia.parsed)

# sort the word frequency tables for snpedia study text
set1.snpedia.freq = as.data.frame(sort(set1.snpedia.freq, decreasing=TRUE))
set2.snpedia.freq = as.data.frame(sort(set2.snpedia.freq, decreasing=TRUE))

# column names changed for merging purposes
colnames(set1.rsnps.freq)[colnames(set1.rsnps.freq)=="set1.rsnps.parsed"] = "text"
colnames(set2.rsnps.freq)[colnames(set2.rsnps.freq)=="set2.rsnps.parsed"] = "text"
colnames(set1.snpedia.freq)[colnames(set1.snpedia.freq)=="set1.snpedia.parsed"] = "text"
colnames(set2.snpedia.freq)[colnames(set2.snpedia.freq)=="set2.snpedia.parsed"] = "text"

# -----------------------------Analysis and visualization
# 1)  Run sections seperately to view graphs properly

set1.rsnps.most.freq = set1.rsnps.freq[set1.rsnps.freq$Freq>150,]

# creating a dot chart, by word frequency
dotchart(as.numeric(set1.rsnps.most.freq$Freq),labels=set1.rsnps.most.freq$text,cex=.7,groups= -set1.rsnps.most.freq$Freq,
         main="frequency counts for rs4880",
         xlab="Frequency of word appearance", gcolor="black", color=c("blue","red","purple","brown"),
)

# merge
crossover1 = merge(set1.snpedia.freq, set1.rsnps.freq, by.x="text", by.y="text", sort=F)
crossover2 = merge(set2.snpedia.freq, set2.rsnps.freq, by.x="text", by.y="text", sort=F)

# 2) -------------------Word Clouds

# we need one column with word frequencies to make a word cloud
# so we take the average of the two datasets and then display the word clouds respectively
crossoverw1 = meanFreq(crossover1)
crossoverw2= meanFreq(crossover2)
crossoverw1$Freq = as.numeric(as.character(crossoverw1$Freq))
crossoverw2$Freq = as.numeric(as.character(crossoverw2$Freq))

# creating word clouds
wordcloud2(crossoverw1, size=1,color='random-dark', backgroundColor="white",fontFamily = "Arial")
wordcloud2(crossoverw2, size=1,color='random-light', backgroundColor="white", fontFamily = "Arial")

# 3) ------------------- Chi square
name.list = levels(droplevels(crossover1$text))
name.list
tbl = as.data.frame(rbind(crossover1$Freq.x,crossover1$Freq.y))
names(tbl) = c(name.list)
row.names(tbl) = c("study","phenotypes")
test = chisq.test(tbl, simulate.p.value = TRUE)
t = data.frame(test$expected)

qqnorm(test$expected)
qqline(test$expected, col = 2,lwd=2,lty=2)

