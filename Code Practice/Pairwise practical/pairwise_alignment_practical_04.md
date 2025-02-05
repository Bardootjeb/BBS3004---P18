---
title: Using affine gap penalties
output: pdf_document
---

# Affine and normal gap penalties

When scoring an alignment we usually assign a score at each position of the
alignment and take the sum. In the simplest case we assign three different
scores:

1. match: the residues are identical
2. mismatch: the residues are different
3. gap: there is a gap inserted into one of the sequences

That is if we have the alignment:

```
ACAT--TGAGC
| |   ||| 
AGAGCTTGA--
```

and we assign score, $g=-8$, $m=4$, $d=-4$ for the gaps, matches and
mismatches ($d$ for different) we get an alignment score of:

$$
S = 4g + 5m + 2d = (4\times -8) + (5\times 4) + (2\times -4) = -20
$$

However, we generally don't use these kind of penalties because insertion and
deletion mutations do not necessarily involve a single nucleotide and because
aligned subsequences are more likely to be informative. Instead we use a
separate penalty for opening and extending gaps. That is that we have a
penalty associated with each run of gaps and a separate penalty associated
with how long that run of gaps is. The gap extension penalty is usually much
smaller than the gap opening (also known as insertion) penalty. Instead of
having, $g=-8$, we might define $g_o=-8$ and $g_e=-1$ for opening and
extending a gap respectively. For the alignment given above we would then get:

$$
S = 2g_o + 2g_e + 5m + 2d = (2\times -8) + (2\times-1) + (5\times 4) + (2\times -4) = -6
$$

## Affine penalties and dynamic programming
When aligning sequences using dynamic programming (Needleman-Wunsch or
Smith-Waterman), we determine the score at each position of the score matrix
by:

$$ 
S_{i,j} = \max
  \begin{cases}
    S_{i-1,j} + P_g \\
	S_{i,j-1} + P_g \\
	S_{i-1,j-1} + M
\end{cases}
$$	

where

$$ M = \begin{cases}
	m, & \text{if}\ A_{i-1} = B_{j-1} \\
	mm, & \text{otherwise}
	\end{cases}
$$

and $P_g$ is a single penalty associated with a gap (i.e. not an affine gap penalty).

If we are using gap penalties we can no longer use a single gap penalty
because we have to determine whether the operation involves opening a new, or
extending a previous run of gaps. To do this we can modify the equation to:

$$ 
S_{i,j} = \max
  \begin{cases}
    S_{i-1,j} + G_u \\
	S_{i,j-1} + G_l \\
	S_{i-1,j-1} + M
\end{cases}
$$	

And define the variables $G_u$ and $G_l$ depending on the route used to obtain
the maximal score in the cell above ($G_u$) and to left ($G_l$) of the current
cell ($i,j$). Remember that a leftwards movement in the table is associated
with adding a gap in the 'vertical' sequence and a downwards
movement is associated with adding gap in the 'horizontal' sequence.

In order to determine the route leading up the prior cells we need to
interrogate the pointer table we set along with the scores. I previously
referred to this as $T$ (in pairwise\_alignment\_practical.pdf), where I
assigned the values 1-3 to indicate whether the maximal score at a given cell
was obtained by moving from the cell above (1), the cell to the left (2) or
the cell above and to the left (3). We can then define $G_u$ and $G_l$ as:

$$
G_u = \begin{cases}
	g_e & \text{if}\ T_{i-1,j} = 1 \\
	g_o & \text{otherwise}
	\end{cases}
$$

$$
G_l = \begin{cases}
	g_e & \text{if}\ T_{i,j-1} = 2 \\
	g_o & \text{otherwise}
	\end{cases}
$$

## Implementing the affine gap penalty in R

We can modify the function found in pairwise\_alignment\_practical.pdf to use
affine gap penalties. That function was defined as follows:

```{r, nmscore_1, eval=FALSE}
NWscore <- function(scores, seq, i, j, match, mismatch, gap){
    ...
}
```

In order to determine whether a gap should be considered as a gap extension or
insertion we need to also pass the pointer table (contaning all the arrows or
paths) and the penalties associated with gap opening (`gap.o`) and extension
(`gap.e).

```{r, nmscore_2, eval=FALSE}
NWscore <- function(scores, pointers, seq, i, j, match, mismatch, gap.o, gap.e){
    ...
}
```

And then modify the body of the function to use these values:

```{r, nmscore_3, eval=FALSE}
NWscore <- function(scores, pointers, seq, i, j, match, mismatch, gap.o, gap.e){
    ## make a vector of scores in order to allow us to use
    ## which.max()
    ## construct a vector of three elements that holds numeric values
    alt.scores <- vector(mode='numeric', length=3)
    M <- ifelse( seq[[1]][i-1] == seq[[2]][j-1], match, mismatch )
    gap.u <- ifelse( pointers[i-1,j] == 1, gap.e, gap.o )
    gap.l <- ifelse( pointers[i,j-1] == 2, gap.e, gap.o )
    alt.scores[1] <- scores[i-1, j] + gap.u
    alt.scores[2] <- scores[i, j-1] + gap.l
    alt.scores[3] <- scores[i-1, j-1] + M
    ## which.max is a convenient function which does what you might expect
    max.i <- which.max(alt.scores)
    ## in R functions return the evalution of the last statement
    ## Hence the following statement will return a vector of two values:
    ## 1. The index associated with the maximal score.
    ## 2. The maximum score itself.
    c(max.i, alt.scores[max.i])
}		
```

We also need to modify the `NWinit` function to handle affine gap
penalties. In this case we can simply specify that the argument `gap` should
be a vector of length two, containing the gap opening and extension penalties:

```{r, nmscore_4, eval=FALSE}
## gap should contain: gap opening, followed by gap extension
NWinit <- function(l1, l2, gap){
    scores <- matrix(nrow=l1+1, ncol=l2+1)
    ## to make the pointer table we can simply copy the scores
    ## table as it has exactly the same dimensions
    ptr <- scores
    ## set the first value to 0
    scores[1,1] <- 0
    ptr[1,1] <- 0
    ## set the values of the first column
    for(i in 2:nrow(scores)){
        scores[i,1] <- scores[i-1,1] + ifelse(ptr[i-1,1] == 1, gap[2], gap[1])
        ptr[i,1] <- 1
    }
    ## and then the values of the first row
    for(i in 2:ncol(scores)){
        scores[1,i] <- scores[1,i-1] + ifelse(ptr[1,i-1] == 2, gap[2], gap[1])
        ptr[1,i] <- 2
    }
    list('scores'=scores, 'ptr'=ptr)
}
```

To use the function you can try:

```{r, nmscore_5, dev='pdf', fig.cap="Needleman Wunsch matrix", warning=FALSE, fig.width=7}
match <- 4
mismatch <- -4
gap <- c(-6, -1)
seq <- c("ACTAGACGAT", "TAGAGACGTTA")
seq.n <- strsplit(seq, "")

## remember this will return a list of two tables
nw.tables <- NWinit(length(seq.n[[1]]), length(seq.n[[2]]), gap)
## check that these are OK before you go any further.

## loop through all cells not yet calculated
for(i in 2:nrow(nw.tables$scores)){
    for(j in 2:ncol(nw.tables$scores)){
        tmp <- NWscore(nw.tables$scores, nw.tables$ptr, seq.n,
                       i, j, match, mismatch, gap[1], gap[2])
        ## assign the new values into the two tables:
        nw.tables$scores[i, j] <- tmp[2]
        nw.tables$ptr[i, j] <- tmp[1]
    }
}

## we can draw the tables..
s1.l <- length(seq.n[[1]])
s2.l <- length(seq.n[[2]])
## 
plot.new()
plot.window( xlim=c(-1, s2.l+2), ylim=c(-s1.l-2, 1), asp=1 )

text( 0, -2:(-s1.l-1), seq.n[[1]] )
text( 2:(s2.l+1), 0, seq.n[[2]] )

segments(x0=0.5, y0=-0.5:(-s1.l-1.5), x1=s2.l+1.5, y1=-0.5:(-s1.l-1.5))
segments(x0=0.5:(1.5+s2.l), y0=-0.5, y1=-s1.l-1.5)

x <- matrix(data=1:ncol(nw.tables$scores),
            nrow=nrow(nw.tables$scores), ncol=ncol(nw.tables$scores), byrow=TRUE)

y <- matrix(data=-1:-nrow(nw.tables$scores),
            nrow=nrow(nw.tables$scores), ncol=ncol(nw.tables$scores), byrow=FALSE)
text(x, y, nw.tables$scores)

## and add some arrows
## define two tables with 0 as the default values
x.delta <- matrix(0, nrow=nrow(nw.tables$ptr), ncol=ncol(nw.tables$ptr))
y.delta <- x.delta

## Here let's use a bit of bitwise magic..
## explained in more detail later..
y.delta[ as.logical(bitwAnd(1, nw.tables$ptr)) ] <- 0.75
x.delta[ as.logical(bitwAnd(2, nw.tables$ptr)) ] <- -0.75

## draw some arrows
arrows( x+x.delta/2, y+y.delta/2, x+x.delta, y+y.delta, length=0.05, col='red' )
```

The score and pointer matrices obtained from this (see figure) are reasonable
and can be used to extract the following alignment.

```
ACTAGA--CGAT-
  ||||  || |
--TAGAGACGTTA
```

\newpage

## Converting the Needleman-Wunsch to the Smith-Waterman local alignment

In order to turn perform a local alignment we simply add an additional
parameter to the equation defining the score table:

$$ 
S_{i,j} = \max
  \begin{cases}
    S_{i-1,j} + P_g \\
	S_{i,j-1} + P_g \\
	S_{i-1,j-1} + M \\
	0
\end{cases}
$$	

That is, only non-negative scores are allowed. This means that the score
matrix usually contains mostly 0s, and can be recognised by the fact that the
first column and the first row will be filled with 0. The alignment is started
from the position of the maximum score instead of from the bottom right-hand
corner and stops when a 0 is reached.



