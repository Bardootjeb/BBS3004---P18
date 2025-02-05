## This is the scoring matrix
NWscore <- function(scores, seq, i, j, match, mismatch, gap){
  alt.scores <- vector(mode='numeric', length=3)
  M <- ifelse( seq[[1]][i-1] == seq[[2]][j-1], match, mismatch )
  alt.scores[1] <- scores[i-1, j] + gap
  alt.scores[2] <- scores[i, j-1] + gap
  alt.scores[3] <- scores[i-1, j-1] + M
  max.i <- which.max(alt.scores)
  c(max.i, alt.scores[max.i])
}

## This is the function created to align NWscores and NWtables?
## seq should be a character vector
NWalign <- function( seq, match, mismatch, gap.o, gap.e) {
  ## insert code to split character vector to list that can be used
  seq.n <- strsplit(seq, "")
  nw.tables <- NWinit(length(seq.n[[1]]), length(seq.n[[2]]), gap)
  for(i in 2:nrow(nw.tables$scores)){
    for(j in 2:ncol(nw.tables$scores)){
      tmp <- NWscore(nw.tables$scores, seq.n, i, j, match, mismatch, gap)
      ## assign the new values into the two tables:
      nw.tables$scores[i, j] <- tmp[2]
      nw.tables$ptr[i, j] <- tmp[1]
    }
  } 
  nw.tables
}

## This is the initialization matrix
NWinit <- function(l1, l2, gap){
  scores <<- matrix(nrow=l1+1, ncol=l2+1)
  ptr <<- scores
  scores[1,1] <- 0
  ptr[1,1] <- 0
  for(i in 2:nrow(scores)){
    scores[i,1] <- scores[i-1,1] + gap
    ptr[i,1] <- 1
  }
  for(i in 2:ncol(scores)){
    scores[1,i] <- scores[1,i-1] + gap
    ptr[1,i] <- 2
  }
  list('scores'=scores, 'ptr'=ptr)
}

## this function uses the pointer information in the ptr matrix 
## to reconstruct the aligned sequences by following the path 
## of diagonal and horizontal movements back from the end of 
## the alignment.
extract.alignment <- function(ptr, seq1, seq2,
                              row=nrow(ptr), colum=ncol(ptr)){
  seq1.al <- c()
  seq2.al <- c()
  while(ptr[row, column] != 0){
    p <- ptr[row, column]
    if( bitwAnd(1, p) ){
      seq1.al <- c( seq1[row-1], seq1.al )
      row <- row - 1
    }else{
      seq1.al <- c( "-", seq1.al )
    }
    if( bitwAnd(2, p) ){
      seq2.al <- c( seq2[column-1], seq2.al )
      column <- column - 1
    }else{
      seq2.al <- c( "-", seq2.al )
    }
  }
  list(seq1.al, seq2.al)
}
