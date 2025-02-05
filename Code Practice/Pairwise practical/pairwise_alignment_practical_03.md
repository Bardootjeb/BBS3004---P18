---
title : Extracting the alignment
---

# Preamble

In the previous sections you have seen how you can calculate the score and
pointer matrices using both fixed and affine gap penalties. To complete the
alignment we need to extract the aligned sequences from the pointer
matrix. This is conceptually simple, but dealing with strings and characters
in `R` is a bit clumsy, and so the resulting code may not look particularly
*elegant*[^elegant].

[^elegant]: Elegant code is code that is concise *and* easy to understand;
    this makes it easier to write without making mistakes and makes it likely
    to be more efficient.
	
# Extraction of the aligned sequences

Aligned sequences are simply the original sequences with gap characters added
the appropriate locations. The pointer matrix tells us where the gaps
should go. To extract the alignments we start either at;

1. the bottom right hand
   corner of the pointer matrix (Needleman-Wunsch, global alignment) or,
2. at the position of the highest scoring cell in the score matrix
   (Smith-Waterman, local alignment) and then,
   
follow the pointers (arrows) until we come across a cell that doesn't point to
any other cell. Since we do not know beforehand how long the alignment will
be, we can't easily use a `for` loop; instead we will use a `while` loop. 
The code in the loop body (enclosed by `{` and `}`) is repeatedly executed
until some condition is met. We will write a separate function that performs
the extraction; if we specify the starting cell as part of the arguments then
we should be able to use the same function for both Needleman-Wunsch and
Smith-Waterman.

```R
## ptr should be a pointer matrix
## row and column specify the position from where we will
##   start the extraction. In the case of Needleman-Wunsch
##   this will always be at the bottom right hand corner, so
##   we can specify default arguments for these two arguments
## seq1 and seq2 are the sequences which have been aligned
##   seq1 is associated with rows and seq2 with columns
## ptr values:
## 1: up
## 2: left
## 3: up and left
extract.alignment <- function(ptr, seq1, seq2,
                              row=nrow(ptr), colum=ncol(ptr)){
    seq1.al <- c()
    seq2.al <- c()
    while(ptr[row, column] != 0){
        ## assign ptr[row, column] to p since we
        ## modify row and column variables below.
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
    ## return a list of the sequences with gaps inserted
    list(seq1.al, seq2.al)
}

```

## Bitwise AND ??

There are two lines of code in the `extract.alignment()` function that
probably confuse you:

```R
if( bitwAnd(1, p) )
```

The `bitwAnd(a,b)` function performs the bitwise *AND* operation on `a` and
`b`. That is it performs an `AND` operation on each pair of bits in its
arguments. Remember that the values in the pointer table can be between 0
to 3. The binary representation of integers zero to three are, 00, 01, 10,
and 11:

decimal  bit 2  bit 1
-------  -----  -----
      0      0      0
	  1      0      1
	  2      1      0
	  3      1      1
	  
The `bitwAnd(a,b)` operation combines the first bit of `a` with the first
bit of `b`, and then the second bit of `a` and with the second bit of `b` and
returns the resulting value:

arg   decimal   binary
---  --------  -------
  a        2      10
  b        3      11
AND        2      10

Hence the bitwise combination of 3 (11) and 2 (10) will give 2 (10) since the
combination of the first bits (1 and 0 for 3 and 2 respectively) gives 0 and
the combination of the second bits (1 and 1) gives 1 resulting in 10.

You may be confused by the fact that the first bit in `10` is the `0` and not
the `1`. This is because we count the bits from the least significant to the
most significant; you can also remember that the least significant bit
multiplies 2^0^ (1), the second bit 2^1^, and so on. But when we write down a
number we write the *most* signficant digit first (eg. the 2 in 239 signifies
200 and is thus the *most* significant).

All of this means, that as we chose 1, 2, and 3 to indicate up, left and up
$+$ left respectively, we can simply ask whether the first and second bits are
set to 0 or 1. If the first bit is set to 1 then we need to incorporate a
letter from the first sequence and move up by one; if the second bit is set to
1 then we need to incorporate a letter from the second sequence and move left
by 1. If the first or second bit is not set then we should incorporate a gap
character ("-") into the respective aligned sequence itself. We use the
`bitwAnd()` function to interrogate the bits of the pointer value:

 ptr  binary   ptr AND 1   ptr AND 2   movement
---- -------- ----------- ----------- ----------
0          00        00           00   stop
1          01        01           00   up
2          10        00           10   left
3          11        01           10   left $+$ up

We then use the return value of the `bitwAnd()` operation as a conditional
statement; a return value of `0` will be coerced to `FALSE`, all non-zero
values will be interpreted as `TRUE`. Hence the code following the `if`
statement will be executed if the `bitwAnd()` call returns a non-0 value.

## Extending the aligned sequence

The lines following the `if(bitwAnd( ... ))` statements extends the
aligned sequence by adding a letter in front of the aligned sequences (these
start out as empty vectors, defined by `seq1.al <- c()` and `seq2.al <- c()`)
using the concatenate function (`c()`) and modifies the `row` and `column`
values appropriately. The first conditional block checks whether the first bit
is set (indicating an up or up + left pointer):

```R
## if the first bit is set (ptr value equals 1 or 3)
if( bitwAnd(1, p) ){
    ## add a letter to the beginning of the seq1.al vector
    seq1.al <- c( seq1[row-1], seq1.al )
    ## and decrease the row counter
    row <- row - 1
}else{
    ## we do not move upwards, and so should add a gap to seq1.al
    seq1.al <- c( "-", seq1.al )
}
```

At the end of the function we create a list containing the two aligned
sequences which gets returned by the function.

## Putting it together

I have put all the functions that were introduced in previous parts and here
together in a file (`functions.R`) that can be sourced from R in order to run
the alignments:

```R
## Note that for this to work, the correct funtions.R must be present
## in the current working directory. If the file is not here, then you'll
## need to specify the path (location) to the file.
source('functions.R')

## first let us define some variables:
match <- 4
mismatch <- -4
gap.o <- -6
gap.e <- -1
seq <- c("ACTAGACGAT", "TAGAGACGTTA")

## NWalign is a function that includes the loop previously present
## in the main script. It calls NWinit and NWscore
nw.tables <- NWalign(seq, match, mismatch, gap.o, gap.e)
visualise.matrix(nw.tables)
seq.al <- extract.alignment(nw.tables$ptr, nw.tables$seq[[1]], nw.tables$seq[[2]])
```

## Fixing your mistakes

Usually when you write a piece of code it will have some mistakes in it. These
mistakes can be logical (eg. you forget about the effect of the order in which
variables are changed) or simple typos (eg. you typed i instead of j). This
means that you must always test your code in some way or other. In this case
this is easiest done using short sequences where you can manually check if the
code does what it's supposed to do. But what do you do when you find that the
code runs but returns something that somehow isn't correct. Simply looking at
your code carefully is unlikely to work; you usually need to use some form of
debugging to work out what is going on. There are two main ways of doing this:

### Add some print statements

This is the simplest way to work out what your code is going on. You can
simply print out the values of variables from strategic positions in your
code. That then allows you to see where things start to go wrong. 

### Run the function through the debugger

`R` comes with a debugger; if you call `debug()` on a function then `R` allows
you to step through each line of code of that function. At each point you can
inspect the values of variables using all standard `R` functions (eg. you can
use plotting functions to inspect complex data sets). This allows you a lot
more control than simply printing out stuff to the console. However, it can
take time, and you ideally should learn all the ways in which you can use the
debugger and this involves reading the manual. Sometimes, printing out stuff
can actually be simpler, though you will need to modify the function to remove
the `print` statements afterwards.
