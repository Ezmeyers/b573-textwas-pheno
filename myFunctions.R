biocLite("biomaRt")
biocLite("SNPediaR")
# By Elisabeth Meyers 2016 

#--------Functions

# this function takes a snp rsid, such as "rs5832" and names
# is a list of genotypes(AA,TT..)
# then it returns a data frame with wiki json text from snpediaR database
# that is related to the rsid+genotype (if any data can found) 
# that is related to the rsid+genotype (if any data can found)

#------------getData() from snpediar
getData = function(x,names){
  list.names = sapply(names,"[")
  nl = names(list.names)
  # get types from snp data and associate with subsets
  n = length(nl)
  indices = 1:n
  # set columns which is the number of genotypes passed in
  cols.needed = length(nl)
  # create a data frame with the specified cols
  # first have to create a matrix to create the cols
  snp.text <- matrix(0, ncol = cols.needed)
  snp.text = data.frame(snp.text)
  for(i in indices){
    # change the column names to match snps passed in, one per snp and genotype
    # will generate something like "rs5832AA | rs5832TT |
    colnames(snp.text)[colnames(snp.text)==paste("X",i,sep="")] = paste(x,nl[i],sep="")
    # need a specific format for querying the DB which is "rs5832(A;A)"
    # seperate the genotypes and add needed chars (A;A) 
    # then concatenate with the rsid
    g1 = substr(nl[i],start=1,stop=1)
    g2 = substr(nl[i],start=2,stop=2)
    text.titles = paste(x,"(",g1,";",g2,")",sep="")
    # can use to specify tags to query if we want this later
    tags = x
    # query the data base for snpedia wiki text returned as json data
    pgs = getPages(titles=c(text.titles))
    # if not null, push to the corresponding column
    if(!is.null(pgs[[1]])==TRUE){
      snp.text[i] = pgs
    }
  }
  # remove zero columns
  snp.text = snp.text[!sapply(snp.text, function(x) all(x == 0))]
  return(snp.text)
}

#------------getGenoCounts() Count genotypes in the dataset 
getGenoCounts = function(x){
  # this will hold our counts 
  genotype.list = list()
  # get counts by creating a table
  counts = table(x$genotype)
  # how many categories of elements we counted
  n = length(counts)
  indices = 1:n
  # loop through elements found
  for(i in indices){
    # look at genotype data where >10 cases
    if(!counts[[i]]<10){
      # push names to list (AA,TT..)
      genotype.list[i] = names(counts[i])
      # push counts to list
      genotype.list[[i]] = counts[i]
    }
  }
  # Remove nulls before returning
  genotype.list = genotype.list[!sapply(genotype.list,is.null)]
  return(genotype.list)
}



# -------------------------parse phenotype text data
parsePheno = function(x){
  options(stringsAsFactors=FALSE)
  n = dim(x)[1]
  indices = 1:n
  text=matrix(ncol=1)
  text = data.frame()
  for(i in indices){
    currRow = str_extract_all((x[i,1]),"\\b(?!list\\b|[0-9])\\w{3,}+")
    text= rbind(text, currRow)
    text[1:2,] = text[1:2,] = NA
  }
  text = na.omit(text)
  return(text)
}

# ----------- rename columns
# simple helper function to change column names
setColNames = function(x,newName){
  colnames(x)[colnames(x)==names(x)] = newName
  return(x)
}

# ---------------data cleanup
# remove stop words and set to lower case
dataCleanup = function(text){
  # convert to lowercase
  text[,1] = tolower(text[,1])
  # remove uneeded common junk words
  junk = c("variation","yes","Yes","doi","pmid","www","can","see")
  text[,1] = tm::removeWords(x = text[,1], c(stopwords("en"),junk))
  # remove blank lines where junk words removed
  text[text==""] = NA
  text = na.omit(text)
  return(text)
}

# ----------- parsing for snpediar data
# only difference is we parse by column here
parseSnpediar = function(x){
  n = ncol(x)
  indices = 1:n
  text = data.frame()
  for(i in indices){
    currRow = str_extract_all((x[1,i]),"\\b(?!list\\b|[0-9])\\w{2,}+")
    text = rbind(text,currRow)
  }
  text = na.omit(text)
  return(text)
}


# -----------------gets phenotype traits data from rsnps
# takes in a data frame of users from rsnps and returns their phenotype data (if any available)
getPhenoData = function (x){
  options(stringsAsFactors=FALSE)
  numrows = nrow(x)
  # create new data frame.
  phenotype.data = as.data.frame(matrix(nrow = numrows, ncol=2, byrow=TRUE))
  # change row names to the user_id's passed in
  rownames(phenotype.data) = paste(x[,4])
  # order by user_id
  phenotype.data = phenotype.data[ order(as.numeric(row.names(phenotype.data))),]
  # change column name
  colnames(phenotype.data)[colnames(phenotype.data)=="V1"] = "Phenotype data"
  # this line is here because I had a strange error where I could not set
  # column names on a data frame with <2 columns. 
  phenotype.data[,2]= NULL
  # rsnps DB can only be queried by consequtive ID's if you are doing more than one user_id query
  # so this line allows you to have less overall queries (takes less time)
  consecutive.ids = as.matrix(split(x$user_id, cumsum(c(1, diff(x$user_id) != 1))))
  colnames(consecutive.ids)[1] = 'V1'
  # length of consecutive ids
  n = length(consecutive.ids)
  # will use to loop through ID
  indices = 1:n
  # store the current row
  currentRow = list()
  for(i in indices){
    # current row
    currentRow = phenotypes(consecutive.ids[i,1])
    # lump all the sublists into a big row
    row = unique(rapply(currentRow, function(x) head(x, 1)))
    # unlist the large row and basically collapse all
    row.concat = paste(rle(unlist(row)),sep="")
    # this will occur if the user doesn't have data, so we are cleaning these values out
    if(length(row)<=2){
      phenotype.data[i,1] = NA
    }
    else{
      phenotype.data = rbind(phenotype.data,row.concat[[2]])
    }
  }
  # remove NA values
  phenotype.data = na.omit(phenotype.data)
  return(phenotype.data)
}

# get average number of frequencies
# data set passed in has a column of words and two frequency columns per word
meanFreq = function(x){
  n=nrow(x)
  indices = 1:n
  x$Freq = ""
  for(i in indices){
    mean = c(x$Freq.x[i], x$Freq.y[i])
    x$Freq[i] = as.numeric(mean(mean))
    x$Freq[i] = round(as.numeric(x$Freq[i]))
  }
  x$Freq.x = NULL
  x$Freq.y = NULL
  return(x)
}
