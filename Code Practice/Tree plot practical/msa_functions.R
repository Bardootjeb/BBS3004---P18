## fn should be a file name that specifies a fasta formatted file
read.fasta <- function(fn){
  data <- readLines(fn)
  id.lines <- grep("^>", data)
  seq.beg <- id.lines + 1
  seq.end <- c(id.lines - 1, length(data))[-1]
  seq <- mapply( function(b,e){
    paste( data[b:e], collapse="" )},
    seq.beg, seq.end )
  names(seq) <- sub("^>", "", data[id.lines])
  seq
}

pep <- read.fasta("Human_SOX17_orthologues.fa")

pep.al <- read.fasta("Human_SOX17_orthologues_al.fa")

cdna <- read.fasta("Human_SOX17_orthologues_cds.fa")

## this will tell me the class and length of the pep object
summary(pep)
## this returned:
##   Length     Class      Mode 
##      195 character character
## indicating that I have 195 sequences

## This will give me the first 6 sections:
head(pep)

## This gives me the number of characters for each name/section
nchar(pep)