source("nw_functions.R")

##Here are some variables
match <- 4
mismatch <- -4
gap <- -8
gap.o <- -8
gap.e <- -4
seq <- c("ACTAGACGAT", "TAGAGACGTTA")
seq.n <- strsplit(seq, "")

##Return list of two tables
nw.tables <- NWalign(seq, match, mismatch, gap.o, gap.e)

## Here are two variables for the plot
s1.l <- length(seq.n[[1]])
s2.l <- length(seq.n[[2]])

## first create a plot window
plot.new()

## Set up a plotting surface with specific coordinate ranges
plot.window( xlim=c(-1, s2.l+2), ylim=c(-s1.l-2, 1) )

## This gives me the text for the X and Y-axis
text( 0, -2:(-s1.l-1), seq.n[[1]] )
text( 2:(s2.l+1), 0, seq.n[[2]] )

## To draw the lines for the plot we have to take half of the 
## coordinates of the points to get to the mid position.
# For the X-axis:
segments(x0=0.5, y0=-0.5:(-s1.l-1.5), x1=s2.l+1.5, y1=-0.5:(-s1.l-1.5))

## For the Y-axis:
segments(x0=0.5:(1.5+s2.l), y0=-0.5, y1=-s1.l-1.5)

## We set up two matrices of x and y coordinates of the same dimensions 
## as the score matrix for the scores
x <- matrix(data=1:ncol(nw.tables$scores),
            nrow=nrow(nw.tables$scores), ncol=ncol(nw.tables$scores), byrow=TRUE)

y <- matrix(data=-1:-nrow(nw.tables$scores),
            nrow=nrow(nw.tables$scores), ncol=ncol(nw.tables$scores), byrow=FALSE)

## This draws in the text into the cells
text(x, y, nw.tables$scores)

## With this line of code we set two new tables with a base value of 0
x.delta <- matrix(0, nrow=nrow(nw.tables$ptr), ncol=ncol(nw.tables$ptr))
y.delta <- x.delta

##Here the guy uses bitwise magic??? From what I understand, it checks
## if a value is either 1 or 2 and then assigns the y.delta or x.delta
y.delta[ as.logical(bitwAnd(1, nw.tables$ptr)) ] <- 0.75
x.delta[ as.logical(bitwAnd(2, nw.tables$ptr)) ] <- -0.75

## This will give us the red arrows
arrows( x+x.delta/2, y+y.delta/2, x+x.delta, y+y.delta, length=0.05, col='red' )

## NWalign calls both NWinit and NW scores.
nw.tables <- NWalign(seq, match, mismatch, gap.o, gap.e)

## This SHOULD visualize the matrix, but it doesn't...
visualise.matrix(nw.tables)

## This line of code takes the extract alignment function and...
seq.al <- extract.alignment(nw.tables$ptr, nw.tables$seq[[1]], nw.tables$seq[[2]])

