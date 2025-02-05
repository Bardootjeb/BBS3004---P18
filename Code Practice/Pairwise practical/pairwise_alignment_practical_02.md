---
title: Visualising the score and pointer matrices
---

# Preamble

If you have followed the instructions in part 1 of this practical you will hopefully have implemented some code that generates score and pointer tables (as matrices). At this point you should ask yourself if your code does what it is supposed to do. One way of doing this is to run your code on small examples where it is possible to visualise the result and then to manually confirm that the code does what it's supposed to do.

# Visualising the matrices

If you have followed the instructions in the previous sections you may have obtained a data structure called `nw.tables` (or similar) that contains two matrices (`nw.tables$scores` and `nw.tables$ptr`) containing the scores and pointers. If the sequences you specified are small enough, then you can of course simply look at these and see if they are what they are supposed to be. You could also export the matrices as text files and read them into a spreadsheet program of some sort, but that would be troublesome. Instead I will explain how you can visualise the tables using `R` low-level drawing commands. To be honest, this isn't likely to be that useful, but it's a resonable introduction to how to draw stuff in `R`; and it really *should* allow you to check if your code does what it's supposed to do.

## Drawing outline

Although `R` has lots of built in `plot` commands and many more can be obtained from various packages, it is useful to be able to create your own from the built in functions. To draw arbitrary things in `R` we:

1.  Clear any pre-existing plot using `plot.new()`.
2.  Set up a plotting surface with specific coordinate ranges using `plot.window()`.
3.  Draw lines, polygons and text based on some data that we wish to visualise.

The actual drawing commands are rather simple; often most of the time is spent making suitable data structures that contain the coordinates of the lines or other objects that we wish to draw.

## Drawing the score table

If we have a score table created from two sequences following the previous instructions, then this will have `length(seq1) + 1` rows and `length(seq2) + 1` columns. We have to map the positions of entries in the table to plotting coordinates. Since we haven't yet set up the coordinate surface that we will draw in, we can use any system we like. The simplest mapping you can do, is to simply map rows and columns to positions in x and y. Then the plotting window will need to have limits defined by the table. This can easily be set up.

### Setting up a plotting surface

``` r
## Here seq.n is a list of character vectors defined in the previous
## practical
## nw.tables contains the score and ptr matrix, and should exist in your
## workspace / session

## first create a plot window
plot.new()

## then define the coordinate system using plot.window()
## plot.window takes two arguments, xlim and ylim
## xlim and ylim should both be vectors containing two values:
## The smallest and largest coordinate in the x and y dimensions
## respectively. Here I simply set both ranges to -1, to the length
## of the sequences + 2. This should give enough space to draw everything
## we need.
plot.window( xlim=c(-1, length(seq.n[[2]])+2), ylim=c(-1, length(seq.n[[1]])+2) )

## we can check this by drawing a point at each corner of the drawing area:
## using the points() function. This takes a vector of x and y coordinates
## as arguments:
points( c(-1, -1, rep(length(seq.n[[2]])+2, 2)),
       c(-1, length(seq.n[[1]])+2, -1, length(seq.n[[1]])+2))

## you may note that we keep repeating the length(seq.n[[1]]). We could maybe
## create a couple of variable for later use:
s1.l <- length(seq.n[[1]])
s2.l <- length(seq.n[[2]])

## but of course, we should have created those before the plot.window command
## hmm. Oh, well, never mind...
```

### Drawing the sequence at reasonable positions

Since we decided that the position of cell `(r,c)` will be drawn at the x-y position `(c,r)` (`r` is the row and `c` is the column, so `c` maps to `x` and `r` maps to `y`), it should be clear that the first sequence should be drawn at an x-position of 0 and with y-positions starting from 2 to `s1.l + 1`, and that the positions of the second sequence should be determined in the same manner.

``` r
## draw the first sequence using the text command.
## the first two arguments are the x and y coordinates
text( 0, 2:(s1.l+1), seq.n[[1]] )
text( 2:(s2.l+1), 0, seq.n[[2]] )
```

### Drawing some lines in the table

To draw the borders in the table, we want to draw at the mid positions between the cells. Since the cell positions are at (1, 2, 3, 4, ..., n), the border positions should be at (0.5, 1.5, 2.5, 3.5, ..., n+0.5). We can use the `segments` function. It's arguments are the (left, bottom, right, top) positions of the lines that you want to draw.

``` r
## horizontal lines (y position varies)
segments(x0=0.5, y0=0.5:(1.5+s1.l), x1=s2.l+1.5, y1=0.5:(1.5+s1.l))
## vertical lines (here using a more efficient way of writing
## the equation
segments(x0=0.5:(1.5+s2.l), y0=0.5, y1=s1.l+1.5)
```

### Adding the scores to the table

To draw the scores we can set up two matrices of x and y coordinates of the same dimensions as the score matrix. This can be done using the `matrix` function paying attention to whether the values are added `byrow` or not.

``` r
x <- matrix(data=1:ncol(nw.tables$scores),
            nrow=nrow(nw.tables$scores), ncol=ncol(nw.tables$scores), byrow=TRUE)

y <- matrix(data=1:nrow(nw.tables$scores),
            nrow=nrow(nw.tables$scores), ncol=ncol(nw.tables$scores), byrow=FALSE)

## Look at x and y and make sure what the above commands did.
## now we can simply use the text command and draw the scores:
text(x, y, nw.tables$scores)
```

If you've got this far, you may have noted that we have drawn the table kind of upside down compared to how it's normally drawn. That is, that we move from the bottom left to the top right. It should not be difficult for you to work out how you can reverse the table so that it looks like would expect it to. I will recap the above code modified so that it draws things in the more common manner.

### Putting it together and reversing y without all the comments!

``` r
## more reasonable order
s1.l <- length(seq.n[[1]])
s2.l <- length(seq.n[[2]])
## 
plot.new()
plot.window( xlim=c(-1, s2.l+2), ylim=c(-s1.l-2, 1) )

text( 0, -2:(-s1.l-1), seq.n[[1]] )
text( 2:(s2.l+1), 0, seq.n[[2]] )

segments(x0=0.5, y0=-0.5:(-s1.l-1.5), x1=s2.l+1.5, y1=-0.5:(-s1.l-1.5))
segments(x0=0.5:(1.5+s2.l), y0=-0.5, y1=-s1.l-1.5)

x <- matrix(data=1:ncol(nw.tables$scores),
            nrow=nrow(nw.tables$scores), ncol=ncol(nw.tables$scores), byrow=TRUE)

y <- matrix(data=-1:-nrow(nw.tables$scores),
            nrow=nrow(nw.tables$scores), ncol=ncol(nw.tables$scores), byrow=FALSE)
text(x, y, nw.tables$scores)
```

### Drawing some arrows in the table

Drawing the numbers in the table was reasonably easy. Visualising the pointer matrix is a little bit more complicated; however we can use the `arrows()` function. This works very much like the `segments()` function, except that it puts an arrowhead at the and of each segment (or beginning, or any combination you like). To use it we simply have to work out what direction the arrow should point in, and thus where it should end. This can be done reasonably simply by using the `ptr` table.

``` r
## We want to define whether the arrows point left, up or diagonally.
## I don't quite remember, but looking at the NWscore function it seems
## that pointer values indicate:
## 1. up
## 2. left
## 3. diagonally up and left.

## This means that we can define two matrices x.delta and y.delta
## that embody this:
## ptr  direction   xdelta    ydelta
## 1    up          0         0.75
## 2    left        -0.75     0
## 3    diagonally  -0.75     0.75
##
## Note that up is a positive direction.

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

Does it look like your code implements the Needleman-Wunsch correctly? In my case it looks like it creates a pretty bad alignment, but it does seem to follow the rules. The quality of the alignment seems to be more related to the high cost of gap insertion and extension. It's possible that using an affine scoring system would give much better alignment. The next session discusses modifiying the code to use affine gap penalties.

### Wrapping it all up as a function

You don't want to write or run all of the above code every time you create another alignment. Instead you want to create a function that performs the above steps. Note that you should probably make a function that does all of the steps of the NWA shown in the previous section as well. I will leave you to consider how that can be done.
