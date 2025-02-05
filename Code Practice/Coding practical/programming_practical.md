% General programming concepts and R

# Preamble

This practical is designed to introduce you to general concepts in
programming. However, all example codes are in `R` which is a very non-general
programming language. This means that I will give many examples of code
(eg. loops) which you may *not* want to use in `R` scripts. I will try to
also give the more `R`-ish way of expressing the same intent for most of
these.

# Basic data types

Computers only handle numbers, and to the computer these numbers are always
encoded in binary as long strings of 0s and 1s. However, the central
processing unit (CPU, the chip where calculations are performed) does not
handle one bit at a time, but instead handles words containing specific
numbers of bits. The shortest word length handled natively is 8 bits; a word
of 8 bits is referred to as a *byte* (or more reasonably in French as an
*octet*). A byte can encode 256 different values and can thus be used to
either encode values from 0 to 255 or -128 to +127[^two_complement]. To
represent larger values we use word lengths of 16, 32 or 64 bits which can
represent $2^{16}$, $2^{32}$ and $2^{64}$ different values respectively.

[^two_complement]: This is true when signed integers are represented using
	two's complement. There are alternatives, but two's complement is the most
	common form.

## Text

To a computer, text is also only a sequence of numbers; the only difference is
in how the numbers are handled. For a given encoding there is a mapping
between the numerical values of the bit-words and characters. The most basic[^basic]
encoding is known as `ASCII`[^ascii] text. `ASCII` uses 7 bits of individual
bytes to represent characters. The first 32 values (0-31) are used for control
characters[^ctl], followed by printable characters and finally a value to
represent the *delete* action. The letters of the English alphabet fall within
65-90 and 97-123 for upper and lower case letters respectively. This may
seem like an arbitrary encoding, but in fact if you consider the binary representation
of the numbers the encoding starts to make a lot of sense.


[^basic]: There may well be more basic encodings, but I am not aware of them;
    but ASCII is pretty simple.
	
[^ascii]: American Standard Code for Information Interchange.

[^ctl]: Control characters change how the text is displayed and include things
    like tabs, new lines, carriage returns and even one that makes your
    computer go 'ping'.

\pagebreak
	
```

          2 3 4 5 6 7       30 40 50 60 70 80 90 100 110 120
        -------------      ---------------------------------
       0:   0 @ P ` p     0:    (  2  <  F  P  Z  d   n   x
       1: ! 1 A Q a q     1:    )  3  =  G  Q  [  e   o   y
       2: " 2 B R b r     2:    *  4  >  H  R  \  f   p   z
       3: # 3 C S c s     3: !  +  5  ?  I  S  ]  g   q   {
       4: $ 4 D T d t     4: "  ,  6  @  J  T  ^  h   r   |
       5: % 5 E U e u     5: #  -  7  A  K  U  _  i   s   }
       6: & 6 F V f v     6: $  .  8  B  L  V  `  j   t   ~
       7: ' 7 G W g w     7: %  /  9  C  M  W  a  k   u  DEL
       8: ( 8 H X h x     8: &  0  :  D  N  X  b  l   v
       9: ) 9 I Y i y     9: '  1  ;  E  O  Y  c  m   w
       A: * : J Z j z
       B: + ; K [ k {
       C: , < L \ l |
       D: - = M ] m }
       E: . > N ^ n ~
       F: / ? O _ o DEL
```

The above tables show the `ASCII` encoding in both hexadecimal (left) and
decimal. In hexadecimal the first printable character is `0x20` (space, top left),
`A` and `a` are `0x41` and `0x61` respectively (i.e. they are exactly 32 spaces
apart). `ASCII` encoding makes it possible to convert between lower case and
upper case by flipping the 6^th^ bit representing 32 (i.e. `0x20`).

If you don't happen to only use English you may find the `ASCII` encoding to
be somewhat unsatisfactory and there are a number of other
encodings that use more bits to encode larger sets of characters. However in
order to represent all characters from all languages requires more than 16
bits[^encod]. Remember that the CPU only handles[^word_length] words of 8, 16,
32 and 64 bits natively. This means that we would need a 32 bit word length to
be able to handle all types of character types. That would require 4 times as
much memory, storage and network capacity than `ASCII` and would be rather
wasteful. Instead, today we generally use `UTF8` encoding. `UTF8` is a
variable length encoding that uses between 1 and 4 bytes where the number of bytes
used for a character is encoded in the first byte of a set.
This makes it impossible to directly index a set of `UTF8` bytes
(i.e. you can't know which of a sequence of bytes represents a given
character; instead you have to parse the bytes from the beginning to the end),
but nevertheless the compactness of the representation makes it the best
option for general representation of text.

[^encod]: Well, I've read this in various places. But I haven't actually gone
    out and counted or tried to find a reference where someone else has done.
	
[^word_length]: This isn't strictly speaking true; there are CPU specific
	extensions that allow you to handle longer words (at least 256 bits, but
	probably 512 as well); but there are no native data types of this
	length. Usually you would use such longer words to perform an operation on
	multiple words simultaneously (eg. SIMD).

You don't need to know the details of the different ways of encoding
text. However you should understand `ASCII` encoding and the fact that dealing
with multilingual text data is complex[^complex]. Understanding `ASCII` is
necessary in order to understand the `fastq` file format and an understanding
of binary encoding in general is needed to understand the use of bitwise flags
in the `sam` format.

[^complex]: Complex; it's difficult to handle correctly and in general a pain
    to deal with. For an example you can read about the differences between
    Python2 and Python3.

### Text data in R

In `R` text is handled by 'character' vectors. Each element of a character
vector holds[^holds] a sequence of bytes representing a set of characters. `R`
treats these characters as `UTF8`, but whether this will work or not may
depend on the setup of your computer. To access single characters (usually
single bytes) you have to use functions like `substring`; however `R` does not
have a data type that handles single characters, so these functions actually
return sequences of characters that can be one character long[^string].

[^string]: You do not need to understand how these sequences of characters are
    encoded, but I am happy to explain if you would like to know.

[^holds]: It is more accurate to say that each element *points* at a sequence
    of bytes. In fact each element contains the location in memory where the
    bytes are held. On a 64 bit computer such a memory address takes up 64
    bits (8 bytes); this means that even an empty string will require at least
    9 bytes (one for the address, and one for the terminal 0 value). But for
    the purpose of this course you can think of the elements as holding or
    containing the strings.

In general, when handling data for bioinformatics we usually care only
about sequences of amino acids or nucleotides, so we are not usually affected
by different forms of text encoding. However, sample descriptions and other
metadata may contain non-ASCII characters and this can potentially be an issue
if your system does not support those character sets as it is not enough for R
to do handle the strings correctly; your system must also know how to draw the
characters.

## Numbers

There two broad categories of numbers; whole numbers (integers), that have no
fractional parts (i.e. 2, but not 2.5), and continuous numbers that can represent
quantities lying between the integers (eg. 2.5 lies between 2 and 3).

### Integers

Integers are relatively easy to represent and deal with. However, even
integers can be divided up into signed and unsigned integers. Unsigned
integers only represent positive numbers whereas signed integers can represent
both negative and positive numbers.

Very small numbers (-128 $\rightarrow$ +127, or 0 $\rightarrow$ 255) can be
represented using single bytes. To represent larger numbers we use either 2, 4
or 8 bytes. Of these, the 4 byte representation is the most common and it can
hold values either between -2e9 to +2e9 (signed) or between 0 to 4e9. Although
these may seem like large numbers it is relatively easy to encounter
situations where they do not suffice. For example, mammalian genomes usually
contain about 3e9 bases; that is to say that you cannot define a position in
the genome using a single signed 32 bit integer. You can also consider trying
to represent the value of a company like Apple or Tesla using a 32 bit
integer.

### Floating point or real numbers

It should be fairly obvious that we often want to handle non-integral
quantities as these pop up all over the place. The problem with representing
such numbers is that there are an infinite number of values between every pair
of integers and we have to somehow choose which ones will be represented by
the encoding. To handle this problem the numbers are split into an exponent
and a significand[^float] and the values calculated as:

$$ S \times 2^{E} $$

where $S$ is the significand and $E$ is the exponent. You do not need to know
these details; however, you should understand that this means that more
quantities can be represented between 0 and 1 than between 100 and 101, and so on,
until you cannot represent adjacent integers.

[^float]: I am just taking these terms from the Wikipedia page
    (<https://en.wikipedia.org/wiki/Single-precision_floating-point_format>),
	which has a rather nice explanation of the actual encoding of a 32 bit
    floating point value. There is obviously a bit more than what I've
    included above.
	
Fractional numbers are usually encoded using either 32 or 64 bits. 64 bits
encoding results both in a larger number and a larger range of quantites that
can be represented.


### Numeric values in R

`R` has only two types of numeric values. 32 bit integers and 64 bit floating
point. Both are signed so the range of integral values that can be represented
in `R` is fairly limited. Most of the time you will not be aware of how
numeric values are represented in `R` and in practice it usually doesn't
matter unless you have very specific requirement.

64 bit floating point values are often referred to as `double` precision
values. You can check whether a value in `R` is encoded as a 64 bit floating
point or 32 bit integer using the functions `is.integer` and `is.double`.

As for text data (strings), in R everything is a vector; there is no way to
have variables that can only hold a single value. A single value is held by a
single element vector.

#### Integer or real

Try the following:

```R
a <- 10
is.integer(a)
is.double(a)

b <- 10L
is.integer(b)
is.double(b)

c <- a + b
is.integer(c)
is.double(c)

d <- b
is.integer(d)

e <- b + 1
is.integer(e)

b <- b + 1
is.integer(b)
## oops

## and lets try ..
a <- c(10L, 10)
is.integer(a)
is.double(a)

b <- c(10L, 15L)
is.integer(b)

is.integer( b + 1 )
is.integer( b + 1L )
```

From doing this you should be able to infer that `R` defaults to using
floating point values; however, you can force R to treat numbers as integers
by appending an `L` to the end of the number. Any operation that combines an
integer with a floating point will result in a floating point number.

You should be able to guess from the above that any reading in of data from
external sources (eg. `read.table`) that looks like numbers will result in `R`
converting these to floating point rather than integers. To force numbers to
be encoded as integers or floating point we can use `as.integer` and
`as.double`:

```R
a <- as.integer(10)
is.integer(a)

a <- as.integer(10.1)
is.integer(a)
a
## oops..

as.double(a)
is.double( as.double(a) )

## vector format
a <- 1:10
is.integer(a)
is.integer( as.double(a) )
## as expected.

is.integer( seq(1, 10, 1))
## so there is a difference!
```

#### Precision and limitations

We can try to test how large a number can be held by an integer by
continuously doubling it and seeing what happens. This can be done using a
for loop (which will be explained in more detail later on), an sapply
statement or using a vectorised expression:

```R
## double the value of x for a given number of times
x <- 1L
for(i in 1:34){
    print(paste(i, x))
    x <- x * 2L
}
print(paste(i, x))

## In R we can try to use sapply (will be explained later)
## to do the above in a bit more compact form

sapply( 1:34, function(i){ 2L^i })

## what happened there?
is.integer( sapply( 1:34, function(i){ 2L^i }) )

## ahah. It got coerced to floats.. How about:
sapply( 1:34, function(i){ as.integer(2L^i) })
## OK, that's more consistent.

## try to understand why this is slightly different:
sapply( 1:34, function(i){ as.integer(2L^i - 1) })

## but this is R, so we can do the same by just writing:
as.integer( 2^(1:34) )

## we can also try..
a <- 2000000000L
is.integer(a)

a + a
## oops
```

From the above you should see that we can run out of integer numbers fairly
easily if we are performing operations where all the operands are integers;
`R` will actively coerce integers to floats in many cases, but will not do so
for overflows or underflows (when the value becomes too large or small to be
represented by the data type used).

To test the precision of a floating point numbers we can add small numbers to
increasingly larger numbers. Again we can do this using `sapply` and see what
happens:

```R
t(sapply(10^(1:21), function(x){ x == x + 10^(-4:4) }))
```

Here I have wrapped the `sapply` statement in `t` merely to make the results
display more nicely. I will cover the use of `t` in detail later on. In short,
it *transposes* the data; that is it swaps columns and rows around.

To work out what the above statement does (if you are having trouble to
understand) you can try to break it down into it's constituent parts:

```R
## what are we applying your function to:
10^(1:21)

## ok, that's just a bunch of numbers.

## what does the function do?
## simulate the function by assigning a value to x
## for convenience we can make a vector and choose values
## from that:

a <- 10^(1:21)

## try with the first element of a:
x <- a[1]
x == x + 10^(-4:4)

## how about the 15th element
x <- a[15]
x == x + 10^(-4:4)

## you should recognise something in the results;

## so what does
x == x + 10^(-4:4)

## actually do. Again break down each element and work it out:
10^(-4:4)

## and lets use a smaller x value
x <- a[3]
x

x + 10^(-4:4)

x == x + 10^(-4:4)

## to work out why I've used t() in the above statement. Just try it without:
sapply(10^(1:21), function(x){ x == x + 10^(-4:4) })
## I think it is easier to look at using t(), but there is no real difference.
```

The above should demonstrate to you that 64 bit floats have a precision of
about 1 up to values about 10^15^. If you don't understand that, either look
at the code above and think about it for a bit, or ask your neighbour or me.

## Other data types

In addition to text and numeric data, lower level programs also handle things
referred to as *references* and *pointers*. These are data types that refer or
point to other variables. In general, *pointers* are usually used in lower
level programming languages where they are equivalent to memory
addresses. *References* are usually some kind of abstraction of *pointers*,
and are used in similar ways. `R` does not have any explicit reference types
and variables can be thought of as containing data that gets copied during
assignments and when used as arguments to functions.[^pointers]

[^pointers]: Obviously this isn't actually true, but for the purpose of
    using R you don't need to know what really happens unless you handle large
    data sets.

Some computing languages (like `R`) also have an explicit logical (Boolean)
data type; this is used to hold values which can only be `TRUE` or
`FALSE`. Generally you can use such values interchangably with numeric values,
but this will depend a bit on the language.

`R` also has a number of other data types, most of which you will not need to
care about. However, it is useful to know that functions are can be treated as
objects of the class `function`. Functions can not really be considered as a
*basic* data type, but I have mentioned them here since they are built into
the `R` language.

# Variable types

When we program we assign values to variables which we then manipulate in
order to do something useful. The set of variables that we use can be said to
represent a model of the thing we are interested in. Frequently this data
model tries to represent physical or conceptual objects. In bioinformatics the
data model may provide information about genes, genomes, sequences, gene
expression gene modification and many other things. The first challenge when
designing a program is to come up with a way in which to represent the data
that we want to manipulate or analyse. To do this, we must understand the set
of variable types that we can use.

Most programming languages have the following variable types:

1. scalars: hold a single value
2. vectors: hold a block of contiguous values that can be accessed by their
   position in the sequence
3. maps: key-value stores where specific values can be accessed and set by a
   some sort of key. This key is often a string (a short bit of text), but can
   conceptually be of any data type.
4. composite types: variables that contain a set of variables. These are often
   referred to as *structs* or *classes* depending on the language. They are
   not defined by the language but can be defined by combining sets of other
   variables.
   
`R` is a bit of an outlier in that it does not have a scalar variable type;
the simplest data type is instead a vector that holds a sequence of
values. `R` also does not have a native map type, but elements of vectors can
be accessed using *names* and this can be used in a similar manner to maps in
other languages.

`R` has at least two different composite data types: `S3` and `S4`
classes. `S3` classes are not formally defined, and are usually implemented by
setting the *class* attribute of a list object. `S4` classes need to be
formally defined and are more complex to implement and to use. They will not
be covered here. `R` also has something referred to as an `R` (for reference)
class; however, this is very rarely used and will also not be covered here.

## `R` variable types

Apart from `S4` classes, almost everything in `R` can be considered as a
vector that has some attributes set:

1. `vector`: a sequence of objects all of which are of the same data type.
2. `matrix`: a two dimensional table; actually implemented as a vector with the
   attribute `dim` set to a vector of length two. The first value gives the
   number of rows and the second the number of columns.
3. `array`: an n-dimensional sequence; implemented as a vector with the
   attribute `dim` set to a vector of length `n`.
4. `list`: a vector that can hold elements of different data types.
5. `data.frame`: a list of lists of equal length. Used as a table that can
   contain different data types.
   
Although all can sort of be considered as vectors we will mostly treat them as
different data types.

### `R` vectors

Generally sequences of objects of same type; vectors can be created in many
different ways and are returned by many functions. Vectors
can be used directly in mathematical operations, and many built-in functions
will perform operations on all elements of a vector.

Consider the following code:

```R
## vectors can be created using the concatenation function:
a <- c(3, 2, 1)
a
class(a)
is.integer(a)
is.vector(a)

## R allows you to do:
a <- c(3, 2, 'b')
a
class(a)
## coerced to a character class
## this kind of conversion can happen when you read in data.

as.integer(a)
## works, but introduces NA (missing value)

## c() allows you to set names:
a <- c('a'=1, 'b'=2, 'c'=3)
a

## we can access elements of the vector using numeric indices or the names:
a[1]

a['a']

## we can access several elements in one go by creating a vector of indices
## or names:

a[c(1,3)]

a[c('a', 'c')]

## note that you can do:
a[c(1,1,2,2,1,3,3,1)]

## and
c <- c('a', 'a', 'c', 'b', 'a', 'c')
a[c]

## this can often be useful. If you wanted to visualise a sequence:
nuc.cols <- c(A='red', C='green', G='blue', T='black')

## we can make a random sequence with the sample command:
seq <- sample(c('A', 'C', 'T', 'G'), size=20, replace=TRUE)

## we can then draw that by doing:
plot.new()
plot.window(xlim=c(0, 1+length(seq)), ylim=c(0,50))
x <- 1:length(seq)
rect(x, 20, x+1, 30, col=nuc.cols)
text(x+0.5, 25, seq, col='white', cex=2)

## and we can use a logical vector to select elements:
b <- a == 2
b

a[b]

b <- a > 1
a[b]
## note how the names are preserved

## vectors can be created as integer ranges
a <- 1:10
a
b <- 10:1
b
class(a)
is.integer(a)

## which makes integers

## the seq function makes sequences of vectors:
c <- seq(1, 10, 1)
d <- seq(1, 10, 0.5)
c
d

## We can also use the vector function
a <- vector(mode='numeric', length=10)
a

## ahah, defaults to 0

## how about a logical vector?
b <- vector(mode='logical', length=10)
b

## you can even do:
b <- vector(mode='list', length=10)
b
is.list(b)

## you can extend a vector by indexing
## an out of range element:
a <- c(1, 2, 3)
a[4] <- 5
a

a[7] <- 8
a

### and c() can be used to concatenate many vectors in one go
c(a, 100, 32:40, c(5, 4, 3, 2, 1))

## which flattens out into a single vector
```

Note that the above is all pretty `R` specific; in other languages you would
create vectors in different ways.

#### Arithmetic with vectors

One of the strengths of `R` is that it is easy to combine vectors using
mathematical operations. Consider the following:

```R
## make a vector:
a <- 1:10

a + 1

## add another vector to a:
a + c(1, 10)

## multiply
a * 2

## against another vector
a * c(0, 1)

## how about:
a * c(0, 1, 0)
## works but with a warning:

## what do the following do:
a %/% 2
a %% 2

## try:
b <- a %/% 2
r <- a %% 2

b * 2

## and finally; the magic happens!
b * 2 + r
```

### `R` matrices

Again, matrices are pretty `R` specific. However, in order to use `R`
effectively you need to be aware of how they work. Try out the following:

```R
m1 <- matrix(1:12, nrow=4)
m1

m2 <- matrix(1:12, ncol=4, byrow=TRUE)
m2

## the magic function t() transposes matrices:
t(m2)

## you can compare two matrices:
m1 == m2

## if they have the same dimensions:
m1 == t(m2)

## to access a single value specify the row and column:
m1[2,3]

## to specify a complete row specify only the row:
m2[2,]

## you can also specify several rows:
m1[2:3,]

m2[2:3,]

## note the results of the following:
dim(m1[2:3,])
dim(m1[2,]) ## especially this one
dim(m1[2,,drop=FALSE])

## the above is explained by:
class(m1[2:3,])
class(m1[2,]) ## especially this one
class(m1[2,,drop=FALSE])

## similarly you can specify columns and combinations
## of columns and rows:

m1[,2]
m2[,2:3]

m1[1:2,2:3]

## and you can even make a bigger matrix:
## through subsetting
m1[c(1,1,1,1:4),]

## You can use all of these to assign values as well:

m3 <- m1
m3[1,] <- m3[1,] * 3

m1
m3

## you can specify a single value:
m3[1,] <- 0
m3

## you can combine vectors using cbind and rbind
## these combine by column and row respectively:

## remind yourself of what m1 and m2 look like:
cbind( m1, m2 )  ## not ok..

cbind( m1, m2[1,] )
cbind( m1, m2[1:2,] ) ## not OK again..

cbind( m1, t(m2) )  ## ok..

rbind( m1, t(m2) ) ## try to understand whats happened here

## cbind and rbind can also make matrices from vectors:
cbind(1:10, 10:1)

rbind(1:10, 10:1)

### Just like with vectors we can have names. But in this case we have one
### set of names for rows and one for the columns:

## first for convenience, check out the letters and Letters objects defined by R
letters
LETTERS

rownames(m1) <- letters[ 1:nrow(m1) ]
colnames(m1) <- LETTERS[ 1:ncol(m1) ]

m1

## you can now do as before, substituting labels instead:

m1['a',]
m1['b','A']
m1[,'A']
```

#### Matrix maths

R is pretty much made for doing linear algebra with matrices. As you might
expect this means that you can do matrix multiplication as well as elementwise
arithemetic.

```R
m1 <- matrix(1:12, nrow=4)
m1

m2 <- matrix(1:12, ncol=4, byrow=TRUE)
m2

## arithmetic with a single value does what you
## would expect:

m1 + 2

## The following is IMPORTANT, and not obvious:
m3 <- m1 + 1:4

## you can compare by element wise subtraction:
m3 - m1

## compare to:
m4 <- m2 + 1:4
m4

m4 - m2

## you should notice that the vector is filled by column;
## If I wanted to add c(1,2,3) to every row of m1 I would have to
## first transpose it, and then transpose it back:

t(m1) + 1:3
t( t(m1) + 1:3 )
## learning to think about using transpositions like this is quite useful.

## you can also do matrix multiply:
m1 %*% m2

## as long as the dimensions of the matrics are good:
m2 %*% m1
m1 %*% m1  ## oops

## there are also a wealth of linear algebra functions for
## matrix maths. But that is quite a specific field.
```

#### Vectors to matrices

I did say something about matrices being vectors with the dimensions set. Well
you can do the following:

```R
a <- 1:12
class(a)
a

dim(a) <- c(3,4)
class(a)
a

dim(a) <- c(4,3)
a

dim(a) <- NULL
a
```

### `R` arrays

Arrays are kind of unusual objects and you probably won't come across them too
often. However, one use is to encode pixel based pictures that you can view
with the `rasterImage` function. Colour in computer images are usually
represented as three values representing the primary colours red, green and
blue. Hence a colour picture can be represented as an array of dimensions
height times width times three (one value for each colour), and that is what
`rasterImage` expects for colour pictures.

```R
## make a cube of 3x3x3 dimensions and fill with
## 1-27

a1 <- array(1:27, dim=c(3,3,3))
a1
## note the order of the filling

## make an array with dimensions 4,3,3
a2 <- array(1:36, dim=c(4,3,3))
a2

## modify the dimensions:
dim(a2) <- c(3,3,4)
a2

## or since we can divide 4 by 2 we can even do
dim(a2) <- c(2,3,3,2)
a2

## to access slices or other elements as usual
a2[1,,,]
a2[,1,,]
a2[,,,1]
```

### `R` lists

`R` lists are just special kind of vectors that can contain elements of
different types. This means that `R` elements can also be lists and this
allows you to make more complex data structures. Like vectors, lists can be
named.

There are really only two peculiarities of the `R` `list` data structure:

1. To access a range of elements you use `[]`. A range of elements (even if of
   length 1) is still a `list` and so cannot be used as the actual value. To
   access the value of an element as its native type you have to use
   `[[]]`.
2. Elements of named lists can also be accessed using the `$` notation (see
   below). This can sometimes be convenient when accessing elements.

```R
## to make a list, either use
## 1. list( )
## 2. vector( mode='list' )
## 3. or append an existing list

l1 <- list('a'=1, 'b'=matrix(1:12, nrow=4), 'd'=4:10, 'e'="hello")

## note that the following also works
## and produces exactly the same thing.
l1 <- list(a=1, b=matrix(1:12, nrow=4), d=4:10, e="hello")

## to access an a range use []
l1[1:2]

## note how that was printed to screen
## esp the $ signs
class(l1[1:2])
## should be a list..

l1[1]

class(l1[1])
class(l1[2])

## note the difference:
l1[[1]]

class(l1[[1]])

l1[[2]]

class(l1[[2]])

## as usual you can also use names. The same rules apply:
l1['a']
class(l1['a'])

l1[c('a', 'b')]

l1[['b']]

class(l1[['b']])

## and you can use the $ notation:
## in which case you don't need quotes..
l1$b

## note:
class(l1$b)

## Unfortunately you cannot do:
el.name <- 'b'
class$el.name

## You can of course do:
names(l1)

## and
names(l1) <- letters[5 + 1:4]

names(l1)

## if you use c to concatate a list:
l2 <- c(l1, vector(mode='numeric', length=10))
## you will get something unexpected
l2

## And doing this also gets an unexpected result
l2 <- list(l1, vector(mode='numeric', length=10))
l2
length(l2)

## how about:
l2 <- c(l1, list(vector(mode='numeric', length=10)))
l2

length(l2)
## maybe that was what we wanted in the first place. But
## now look at the names attribute:

names(l2)

## to set the 5th name we can do..(that is empty at the moment)
names(l2)[5] <- 'K'

names(l2)

l2
```

### `R` Data frames

The `data.frame` structure is essentially a list of lists of the same
length. Normally you obtain dataframes from functions that read files
(eg. `read.table`), but you can also make a `data.frame` using the
`data.frame` function.

```R
df1 <- data.frame( letters[1:10], 1:10 )
df1

## but it's niceer to give column names:
df1 <- data.frame( 'al'=letters[10:1], 'num'=1:10 )
df1

## that looks nicer.
## what's the length of the data frame?
length(df1)

## the length is the number of columns..
## what happens if I select the first element?

df1[1]

class(df1[1])

## what about
df1[[1]]
## eh??

class(df1[[1]])
## The result of this operation depends on whether you
## are using R version 4 or later. In earlier versions
## this would give you a factor; not a character vector.

## how about ..
class(df1[[2]])
## more reasonable

## the class of df1[[1]] is a factor because R likes to convert strings
## to factors. This behaviour can cause weirdness and in bioinformatices it's
## often better to suppress it. This can be done by:
df1 <- data.frame( 'al'=letters[1:10], 'num'=1:10, stringsAsFactors=FALSE )

## Though you probably don't have to do this as you are probably using
## up to date versions of R.

class(df1[[1]])
## possibly more in line with what we expect.

### Having seen the above you may guess that we can also access elements using
### the $ operator:

df1$al
df1$num

#### You can also treat the data frame as a table, specifying rows and columns:
df1[1,1]
df1[1,2]

## or
df1[1:3,]

## note the difference between:
class( df1[1,1] )
class( df1[1,] )
class( df1[,2] )
```

## Control flow (conditionals, loops and functions)

The simplest of all programs consist of a set of statements executed in linear
order as in the following short script:

```R
source <- readLines( "programming_practical.md" )
print("The first line of the source:")
print(source[1])
print("The last line of the sourece:")
print(source[ length(source) ])
```

Here the execution of the program follows in a strictly linear order
proceeding from the top (line 1) to the bottom with every line executed
exactly once. This is easy to understand, but is not sufficient to do
interesting or useful things.

The script simply reads in the contents of a file (this file as it happens)
and prints out the first and the last line of the file along with a couple of
lines informing the user what happens. You can try to run the code if you
like; it should work as long as you have downloaded this file or have any file
of the same name in current working directory. You can of course change the
name of the file, but note that you may get some unexpected result if the file
is not in a text format[^text].

[^text]: There is often a confusion about what it means for a file to be a
	text format, and I have seen things like Microsoft Word files referred to
	as 'text files'. Word files are definitely not text files, but it is
	actually not that easy to define a rule that can determine whether a file
	is in a 'text format' or not. However, for a file to be text format, then
	every byte or set of bytes within the file should map to a character of
	some sort. It also means that numbers are not represented as integers or
	floating point numbers but by strings of characters of some sort.

You may also note that this script cannot do very much, and that it might
result in an error (for example if the file specified does not exist or is not
readable). The script also cannot go through every line of the file because
the author cannot possibly know how many lines there are. If there are no
lines it will also fail, and if there is only one line it will print the same
thing out twice.

To deal with these issues and to be able to make the script do something
useful, we make use of *conditionals*, *loops* and more structured statements
like the *apply* family that allow you to run functions on the elements of
data collections.

### Conditionals

The basic conditional expression, which you can find in more or less any
programming language is the `if-else` statement. We can modify the script
above to first check if the file exists. Note that if we are going to do this
we will use the name of the file in two locations. Whenever that happens you
*really should* assign the file name to a variable first and then use that in
both locations. If you don't you will have to update the code in more than one
place everytime you want to modify the code to look at another file; and you
**will** eventually write the name differently resulting in code that doesn't
work.

```R
fname <- "programming_practical.md"
if(file.exists(fname)){
    source <- readLines( fname )
    print("The first line of the source:")
    print(source[1])
    print("The last line of the source:")
    print(source[ length(source) ])
}else{
    print("The file does not exist")
}
```

It should be fairly clear what is going on here. `file.exists` is an `R`
function that checks whether a file exists. If it does, `TRUE` is returned, if
not `FALSE` is returned. The contents of the brackets in the `if( ... )`
statement can be anything that evaluates to `TRUE` or `FALSE`. You can even
write:

```R
if(TRUE){
    ## do interesting things
    ## The code here will always be executed.
}

if(FALSE){
    ## the code within this block
    ## will never be executed
}
```

The above code might seem pointless since the code in the first statement will
always be executed and the code in the second one never. However, sometimes
when developing code it can be useful to exclude the execution of a block
temporarily and this is one of the ways in which this can be done.

The set of conditional statements available in different languages varies a
bit; there are also `if-elsif-else` and `switch` statements that allow the
program to choose between more than two alternative code paths. There is also
the rather useful ternary operator in many `C`-like languages:

```c
a = a < 0 ? b : c
```

The part to the left of the `?` should be an expression that evaluates to
`TRUE` or `FALSE`. If the expression evaluates to `TRUE` then the value to the
left of the `:` is returned; otherwise the value to the right of the `:` is
returned. This statement can be used to assign values to variables since one
of the two alternatives is returned. The statement above is just a short way
of writing:

```c
if(a < 0){
   a = b
}else{
   a = c
}
```

`R` does not have the ternary operator, but it does have the `ifelse`
statement that does the same thing:

```R
c <- ifelse( <cond>, a, b )
```

Here, `<cond>` should be replaced with an expression that evaluates to `TRUE`
or `FALSE`. If the statement is `TRUE` the second argument (`a`) will be
returned; otherwise the third argument (`b`) is returned. The `ifelse`
statement is actually more useful than the simple ternary operator in `C`,
because it works on vectors; all of the arguments can be vectors and the
statement can return a vector as well. This is often useful, but be aware that
`ifelse` can be a bit finnicky about the relative lengths of the three
arguments (you can sometimes get unexpected results, so you should always
confirm that the code does what you think it does).

### Loops

A loop is a block of code that is repeatedly executed as long as some
condition is `TRUE`. This isn't alway obvious in higher level languages where
the loop expression is abstracted away from the implementation. The two most
common types of loops are the `for` and `while` loops and I will describe
how to use them in `R`.

#### `while` loops in `R`

Execute a block of code while some condition is true. We can modify the little
script above to print out every line of the source file instead of just the
first and the last one. To make it moderately useful we will prefix each line
with its line number:

```R
fname <- "programming_practical.md"
if(file.exists(fname)){
    source <- readLines( fname )
    i <- 1
    while(i <= length(source)){
        ## paste allows us to make a string from several components
        ## including non string elements like integers
        print( paste(i, source[i]) )
        i <- i + 1
    }
}else{
    print("The file does not exist")
}
```

Here I simply set the value of an index variable (`i`) to 1, and then execute
the body of the loop as long as the value of `i` is smaller than or equal to
the length of the `source` character vector. The last statement of the loop
body increments the value of `i` by 1. If it did not then the loop would never
terminate and it would simply print out the first line of the input file until
killed.

#### `for` loops in `R`

The `for` loop executes once for every value or instance in a range. It is
short for 'for each'. In `R` it is used as:

```R
for(e in v){
    ## do something with e (an element of vector v)
}
```

You can read the above as, 'for each element e in v do something'. Here, `e`
will be assigned first to the first element of `v` and the code executed, then
`e` is assigned to the second element the code run, and then the third and so
on.

There is a variant of this that uses an index to subset `v` with instead of
assigning the elements of `v` directly to `e`:

```R
for(i in 1:length(v)){
    ## do something with v[i]
}
```

Although the latter form might seem pointlessly more verbose it can be used to
modify the values of `v`. Try the following:

```R
v <- 1:10
for(e in v){
    e <- e * 2
    print(e)
}

## did anything happen to v?
v

## but if we do:
v <- 1:10
for(i in 1:length(v)){
    v[i] <- v[i] * 2
    print(v[i])
}

## and note the difference
v
```

#### A slighlty less pointless example of using a loop in R

Instead of simply printing out each line of the source file we can try to ask
if there are any interesting[^interesting] patterns in it. One of the simplest
things we can ask is if there are any patterns in the number of characters per
line. To do this we can use the `nchar` function that tells us the number of
characters in a string.

[^interesting]: Well maybe not that interesting, but it's a start.

```R
fname <- "programming_practical.md"
if(file.exists(fname)){
    source <- readLines( fname )
    char.count <- vector(mode='numeric', length=length(source))
    for(i in 1:length(source)){
        char.count[i] <- nchar( source[i] )
    }
    plot( char.count )
}else{
    print("The file does not exist")
}
```

Not terribly interesting, though there are actually some patterns there that
indicate something about the structure of the file.

Of course `R` is a very high-level language (has lots of convenient functions
built in) and so you wouldn't actually do the above. Instead you would do
something like:

```R
## this excludes the reading from the file and assumes that the lines of the file
## are in the vector source

## general version..
char.count <- sapply( source, nchar )
plot( char.count )
## looks the same

## but since nchar is vectorised we can actually do:
char.count <- nchar(source)
plot(char.count)
## looks the same as well... and is probably much faster
```

Generally in `R` one doesn't use loops that much; instead you are more likely
to use the `apply` family of functions to apply a function to all elements (or
slices) of a data structure. However, to understand how to use `apply` it is
necessary to understand how to create functions in `R` and how they work. I
have introduced loops here for two reasons:

1. Loops are kind of universal; you find them in all programming languages
   that I know of.
2. There are things that are difficult to do with the `apply` statements; for
   example it is difficult to modify the contents of the data that is applied
   to the functions. This is because of *scope*. This will be explained in the
   function section.
   
### Functions

The conditional statements and the loops described above allow you to escape
the strict linearity of execution of the program code. However, all they
really allow is for you to either skip certain parts of the code, or to repeat
execution of a given block of code. Functions allow you to execute code that
is defined at locations distant from the code that calls the function (this
can be somewhere else in your code base, or code defined in a package that you
have loaded, or even part of the core language).

However, after the code in the function has finished executing, the control
will return to the statement immediately after the position where the function
was called. This maintains the linearity of the program and is in contrast to
the `goto` statement which simply allows you to move to another part of your
program. `R` doesn't have `goto`, and in general use of `goto` is discouraged
in most languages that has the `goto` function[^goto]. This means that you can
sort of think of the code in the function body as being more or less copied
and pasted into your code when you call it[^copied].

[^goto]: One of the most famous of all articles in computer science (by Edsger
    Dijkstra) was given the title, 'Goto considered harmful'.
	
[^copied]: But there are important differences; in particular with respect to
    scope.

In `R`, functions are treated as normal data objects that can be called by
appending brackets to the end of their name. The brackets can also contain
arguments to the function; these are variables whose values will be copied to
the function scope and can be used within the body of the function. To create
a function we call `function` and assign the return value to a variable
bearing the name of our function. Consider:

```R
## a function called print.hello
## print.hello takes no arguments but simply prints out a message
print.hello <- function(){
    print("Hello World!")
}

## the above code created a function called print.hello
## to call it, simply do:
print.hello()

## you can repeat that as many times as you like:
print.hello()
print.hello()

## For this function, the effect is exactly the same as if the body of
## the function had been included into the main body of the text

print("Hello World!")
## so this function is largely useless.
```

A slightly more complicated example might read in the contents of a file and
print it out line by line, like we did above:

```R
print.lines <- function(){
    fname <- "programming_practical.md"
    if(file.exists(fname)){
		source <- readLines(fname)
        char.count <- vector(mode='numeric', length=length(source))
        for(i in 1:length(source)){
            print(paste(i, source[i]))
        }
    }else{
        print("The file does not exist")
    }
}

## when I call the function, it is almost as if I had included the function body
## where I make the call:
print.lines()

## this is convenient if I want to call the function many times
print.lines()
print.lines()

### IMPORTANT NOTE ABOUT THE THING WE CALL SCOPE ###
## in the function we set the value of the variable fname
## This will not change the value of a variable of the same name
## that I have defined here..

fname <- "something_else"
print.lines()

## the value of fname has not changed here!
## THIS IS ONE OF THE THINGS YOU MUST UNDERSTAND
fname
```

In the example above I modified some of the code I used to demonstrate the use
of loops and put that code into a function. That function can now be called
many times. Since the function will do exactly the same thing every time, that
is not that useful. But you should notice that although there are several
variables that are defined and set in the function body, this will not change
or create any new variables outside of the function body. This is what we
refer to as *scope*; that is to say that we can have variables that are
available within different parts of the code, but not others. In this case the
variables defined in the function are *local* to the function. Variables that
are defined in the *main* part of the code (the stuff that is run directly)
have *global* scope in `R`. This means that they can be accessed within
functions; however, in `R` you cannot modify *global* variables within functions[^scope2].
The rules of scope varies between programming languages, and it is important
that you understand how they work in order to be able to effectively program.

[^scope2]: You will not be surprised to read that this is a bit of an
    oversimplification and not strictly true. It is actually possible to
    modify the variables of the calling scope, but it is frowned upon and
    unusual in `R`. So we will ignore that for the time being.
	
As mentioned, variables in the *global* (or *parental*) scope are accessible
within functions. This means that we *could* rewrite our previous function to
make use of a *global* `fname` variable:

```R
### DON'T DO AS SHOWN HERE. THIS IS AN EXAMPLE OF BAD CODE ### 
print.lines <- function(){
    if(file.exists(fname)){
        source <- readLines( fname )
        for(i in 1:length(source)){
            print(paste(i, source[i]))
        }
    }else{
        print("The file does not exist")
    }
}

fname <- "programming_practical.md"
## this will still work
print.lines()

## we can now modify the value of fname:
## to print out the lines of another file
fname <- "another_file_somewhere"
print.lines()
```

In the above example I have defined a variable by the name of `fname` in the
main body of the code. This can be read in functions and I have made use of
this to make the function read the contents of a different file name. So this
is now somewhat more flexible than before as the function can now do different
things. However, it is A BAD IDEA to write functions like these. This is
because this function now relies on the presence of a specific variable in the
global scope; this means that the function is not portable. I cannot make this
function available to other code bases, but will have to rewrite it every
time. It is also fragile as I might accidentally overwrite the value of
`fname` in the my main code.

A better alternative to make the function useful is to specify that the file
name should be specified as an argument to the function. We do this by
modifying the call to the `function` function:

```R
### THIS IS OK THOUGH !!
print.lines <- function(fname){
    ## the caller specifies the argument; this will now be accessible within the body of the
    ## function using the name given as the argument to function( )
    if(file.exists(fname)){
        source <- readLines( fname )
        for(i in 1:length(source)){
            print(paste(i, source[i]))
        }
    }else{
        print("The file does not exist")
    }
}

## we can no do:
print.lines("programming_practical.md")
print.lines("another_file_somewhere")
```

We now have a flexible function that will print out the lines of any file that
we specify in the argument to it. But maybe we don't want to print out lots
and lots of lines if the file is large. We can provide another argument that
specifies a maximum number of lines to print out:

```R
### THIS IS OK THOUGH !!
print.lines <- function(fname, max.lines){
    ## the caller specifies the argument; this will now be accessible within the body of the
    ## function using the name given as the argument to function( )
    if(file.exists(fname)){
        source <- readLines( fname )
        ## sorry about the next line it is a bit of a mouthful
        for(i in 1:min(max.lines, length(source))){
            print(paste(i, source[i]))
        }
    }else{
        print("The file does not exist")
    }
}

## we can no do:
print.lines("programming_practical.md", 10)
print.lines("another_file_somewhere", 15)
```

Here I have added another argument to the function definition (`max.lines`)
that we use to limit the number of lines that are printed out. To use it have
simply modified the expression `1:length(source)` to `1:min(max.lines,
length(source))`; i.e. the end point of the range should the whichever is
smaller of the number of lines in the file (`length(source)`) or the value of
`max.lines`.

In general when calling functions you have to be careful of the order of the
arguments; `print.lines("filename", 10)` is not the same as `print.lines(10,
"filename")`. This can be difficult when there are many arguments. To help
with this `R` allows you to use named arguments. So you can write:

```R
print.lines(fname="some_file", max.lines=10)
## or
print.lines(max.lines=10, fname="some_file")
```

and both should work. You can also use combinations, where you specify the
name of some of the arguments, but not others, but be careful with this until
you are comfortable with how it works.

In `R` you can also specify default arguments. That is, the values that
variables will take if they are not specified by the caller. Consider the
following set of functions that add two numbers together:

```R
sum.1 <- function(a, b){
    print(paste(a, "+", b, "=",  a + b ))
}

sum.2 <- function(a, b=1){
    print(paste(a, "+", b, "=",  a + b ))
}

## this is a bit unusual
sum.3 <- function(a=10, b){
    print(paste(a, "+", b, "=",  a + b ))
}

sum.4 <- function(a=10, b=1){
    print(paste(a, "+", b, "=",  a + b ))
}

sum.1(10, 20)
sum.2(10)
sum.2(a=10)

sum.2(b=10) ## oops

sum.3(20)   ## oops
sum.3(a=20) ## oops
sum.3(b=20) ## better

sum.4()  ## pointless...
```

If you run the code the meaning should be sort of self-explanatory.

#### Return values

Given that functions in `R` cannot modify the values of variables in the
calling scope they cannot directly do useful things. For example I can read a
file into some data structure in a function, but when the function is finished
executing and control returns to where the function was called from, all those
variables will be lost[^tears]. When this happens, control of execution is said to
have *returned* to the calling scope; and when a function reaches the end of
execution we say that it has *returned*. When a function *returns* it can pass
the value of a variable back to the calling function; this is called the
*return* value.[^return]

[^tears]: ... like tears in rain.

[^return]: Yeah, I know, there are too many 'returns' in that paragraph. But
    it's as good as I can make it.

In `R`, the values can be returned from functions either by explicitly calling
`return()` in the function, or through the evaluation of the last statement in
the function. Again, we can modify the function we used for printing lines and
get it to return the lines of the file.

```R
### Note the use of the default value of max.lines!!!
print.lines <- function(fname, max.lines=10){
    ## the caller specifies the argument; this will now be accessible within the body of the
    ## function using the name given as the argument to function( )
    if(file.exists(fname)){
        source <- readLines( fname )
        ## sorry about the next line it is a bit of a mouthful
        for(i in 1:min(max.lines, length(source))){
            print(paste(i, source[i]))
        }
        return(source)
    }else{
        print("The file does not exist")
    }
    ## we don't have to put anything here, but we can explicitly return
    ## a NULL value that can be checked afterwards
    NULL
}

## now when we call the function we can assign the content of the file to a
## local variable

lines.1 <- print.lines( "programming_practical.md")
lines.2 <- print.lines( "some_other_file", max.lines=1) ## print only one line

## but note that this will give you an error. See if you can work out why
lines.2 <- print.lines( "some_other_file", max.lines=0) ## try not to print out..
```

To some extent you should be familiar with the concept of *return* values,
since you use them almost everytime you call a function (`vector` returns a
vector, `matrix` returns a matrix, `c` concatenates vectors and returns a
vector, and so on). So the above is mainly about what we call the action and
how you can implement it in your own functions.

#### The ellipsis `...`

If you have looked up `R` manual pages you may have noticed '`...`' as one of
the arguments. This is called the ellipsis and represents arguments not
handled by the function but which will be passed to another function within the
function body. Consider the rather pointless `draw.text` function defined below:

```R
## a function that simply wraps the text function
draw.text <- function(txt, x, y, clear=TRUE, ...){
    if(clear)
        plot.new()
    text(x, y, txt, ...)
}

## The draw.text function has three named arguments; any other arguments
## passed to draw.text will simply be passed on to text.

## To use the function we have to set up somewhere to do the drawing
plot.new()
plot.window(xlim=c(0,1), ylim=c(0,1))

draw.text("Hello World!", 0.5, 0.5)
## in the next call to draw.text we specify the argument adj;
## this will be passed to the text function through the
## ellipsis (...)
draw.text("Goodbye World!", 1, 0, FALSE, adj=c(1,1))
    
## We can also pass other arguments to text
draw.text("Hello World!", 0, 1, adj=c(0,1), cex=3)
draw.text("Goodbye World!", 1, 0, FALSE, adj=c(1,0), cex=3)

## We don't have to name the arguments in ..., but it is really
## advisable to do so.
## The following works, because adj is the
## first argument to text() after x,y, and labels.
draw.text("Hello World!", 0, 1, TRUE, c(0,1), cex=3)

## but the following will do the wrong thing:
## since c(0,1), the argument that we want to pass on to
## text will instead be assigned to the clear argument (which
## should take a boolean value).
draw.text("Hello World!", 0, 1, c(0,1), cex=3)

## so _do_ name the arguments that are passed to ...
```

When a function takes `...` as part of its arguments the documentation should
state explicitly the function to which the the ellipsis arguments will be
passed.

#### Object orientation in `R`

`R` is kind of object orientated[^oop]; for `R` this means that
there can be several different functions that appear to have the same name. Such
functions are referred to as *generic* functions and `R` will look at the
arguments passed in the function call and determine what function will be
used. This is called *function* dispatch, and is too complicated a topic to
cover here. If you are feeling particularly ambitious you can find more details
at <https://www.stat.umn.edu/geyer/3701/notes/generic.html> .


[^oop]: There isn't much of an agreement of exactly what is needed in order to
    consider a computing language object. This doesn't really matter, but what
    does matter is how `R` determines what function to run.

As an example, consider the function `plot`. Remember that when you simply
type the name of a function you will usually see the source code for the
function. If you do this for `plot` you will get the following:

```R
function (x, y, ...) 
UseMethod("plot")
<bytecode: 0x5600db691618>
<environment: namespace:graphics>
```

The `UseMethod("plot")` bit tells us that `R` will do something magical to
choose a plotting function that depends on the class of the arguments passed
to the plot function. You can use the `methods` function to find what the
currently available[^avail] options are:

```R
methods(plot)
 [1] plot.acf*           plot.data.frame*    plot.decomposed.ts*
 [4] plot.default        plot.dendrogram*    plot.density*      
 [7] plot.ecdf           plot.factor*        plot.formula*      
[10] plot.function       plot.hclust*        plot.histogram*    
[13] plot.HoltWinters*   plot.isoreg*        plot.lm*           
[16] plot.medpolish*     plot.mlm*           plot.ppr*          
[19] plot.prcomp*        plot.princomp*      plot.profile.nls*  
[22] plot.raster*        plot.spec*          plot.stepfun       
[25] plot.stl*           plot.table*         plot.ts            
[28] plot.tskernel*      plot.TukeyHSD*     
see '?methods' for accessing help and source code
 ```

For `S3` class objects `R` will look for a method that has a `.classname`
appended to it. Eg, if you do:

```R
h1 <- hist( rnorm( 200 ) )
class(h1)  ## histogram

plot(h1)
## this should be equivalent to

plot.histogram(h1)
## except for the fact that the "*" at the end of the name in the
## output of methods(plot) means that the function is not exported
## explicitly, and that we can't call it directly.
##
## you can get the code to the function with:

getS3method("plot", "histogram")
```

If no suitable method is found then `R` will use the `.default` method.

Although the source of the function may not be exported (as for the class
`histogram`) you should still be able to get the documentation for the
specific functions:

```R
?plot.histogram

plot.histogram            package:graphics             R Documentation

Plot Histograms

Description:

     These are methods for objects of class ‘"histogram"’, typically
     produced by ‘hist’.

Usage:

     ## S3 method for class 'histogram'
     plot(x, freq = equidist, density = NULL, angle = 45,
                    col = NULL, border = par("fg"), lty = NULL,
                    main = paste("Histogram of",
                                 paste(x$xname, collapse = "\n")),
                    sub = NULL, xlab = x$xname, ylab,
                    xlim = range(x$breaks), ylim = NULL,
                    axes = TRUE, labels = FALSE, add = FALSE,
                    ann = TRUE, ...)
     
     ## S3 method for class 'histogram'
     lines(x, ...)
...
```

[^avail]: The set of options changes when you load packages and is not constant.

### The `apply` family of functions

As mentioned many times, you are not usually expected to make use of loops in
`R`. It is better to use vectorised statements when possible, and when not you
can use one of the `apply` like statements. In general the apply statements
take two arguments: a collection of objects (like a vector, a list, or a
matrix) and a function that should be applied to either all elements of the
collection or groups of the collection. The function can either be a named
function that exists, or a function defined in the call to the apply function
itself.

There are a number of `apply` statements:

#### `sapply` and `lapply`
`sapply` and `lapply` are used with vectors or lists to run a function on
every element of the list. `lapply` will always return a `list` whereas
`sapply` will try to return a `vector` or a `matrix` when possible (possible
when each call returns an object of the same length and type). The first
argument is the vector of data to which the function should be applied, and
the function is applied as the second argument. Consider:

```R
## a function that returns a variable length vector
f1 <- function(x){
    fct <- c(x)
    while( (x / 2) == as.integer(x / 2) ){
        x <- x / 2
        fct <- c(fct, x)
    }
    fct
}

## a function that always returns exactly two valuees
f2 <- function(x){
    c(x, x/2)
}

## try the functions on some numbers:
f1(54)
f1(128) ## maybe useful..
f2(128)

## define a vector we will apply these functions to:
v1 <- 1:100

## try first with lapply
v1.f1 <- lapply( v1, f1 )

### AND LOOK AT THE RESULT.
v1.f2 <- lapply( v1, f2 )

## sapply will simplify the ouput of the second of these functions

v1.sf1 <- sapply( v1, f1 )
v1.sf2 <- sapply( v1, f2 )

## look at the resulting objects...
v1.sf1
v1.sf2

class(v1.sf1)
class(v1.sf2)

dim(v1.sf1) ## NULL
dim(v1.sf2) ## something better
```

We often use apply like statements for small functions that we do not wish to
define in the main namespace (easy to get too many variables). So we can
define the function in the call to the apply function.

```R
v1 <- 1:100

v1.f1 <- lapply( v1, function(x){
    fct <- c(x)
    while( (x / 2) == as.integer(x / 2) ){
        x <- x / 2
        fct <- c(fct, x)
    }
    fct
})

## and that should do the same thing as above.. 
```

The argument to the function will be assigned from the elements of the
vector. The first time the function is executed it will take the first value
of the vector, then the second and so on. Whatever is returned by the function
will be combined into a list, vector or matrix depending on whether `lapply`
or `sapply` was called and whether the data can be simplified.

#### `apply`

`apply` is used to run functions on either the rows or columns of a
matrix. The first argument is the matrix, the second (MARGIN) defines whether
the function will act on the columns or the rows of the matrix, and the third
is the function that should be run. Like other `apply` like statements you can
also specify optional arguments to the function (i.e. ones that have default
values).

```R
## define a matrix
m <- matrix(1:12, nrow=4)
## look at it!!!
m

## calculate sums of rows
apply( m, 1, sum )

## and columns
apply( m, 2, sum )

## can also get the range for each colum and row
## range returns two values (min and max) so we will get a matrix back
## when we do this:
apply(m, 1, range)
apply(m, 2, range)

## do look to make sure that the numbers make sense!!!

## in R there is a special value NA for missing or nonsensical numbers
## these cause problems with arithmetic:

m[3,3] <- NA
m

apply(m, 1, sum)
apply(m, 2, sum)

## we can fix that (sort of) by doing
apply(m, 1, sum, na.rm=TRUE)
apply(m, 2, sum, na.rm=TRUE)

## na.rm is an optional argument to the sum function. You should be able to work out
## what it does, from the name and the above. If not, do some experimentation.
```

#### `tapply`

`tapply` is one of the more complicated apply statements. It operates on two
vectors and calls the specified function on values of the first argument
grouped by the second argument. Perhaps most easily explained by an example:

```R
## define a vector
v1 <- 1:100

## define a vector that can be used to group v1 into
## groups. We can call this an index

g1 <- sample(LETTERS[1:5], size=length(v1), replace=TRUE)

## too look at both at the same time you can do:
## note that in this v1 will be coerced to a character vector
## but that this does not affect v1 here.. 
head( cbind( v1, g1 ), n=20 )

## to calculate the sum of all elements with an A, B, C, D, or E in the same position we
## use tapply

tapply( v1, g1, sum )

## this is the same as if I had done..

sum( v1[ g1 == 'A' ] )
sum( v1[ g1 == 'B' ] )
## and so on.
```

#### `mapply`

All the previous `apply` functions apply a function to each element of a
single data structure. `mapply` allows you to pass elements from several lists
to a single function in parallel. For `mapply` you must first specify the
function and then the arguments that will be passed to it. I will just give a
very minimal example. If you wish to use `mapply`, then please read the
manual!

```R
v1 <- 1:10
v2 <- 20:11

mapply(sum, v1, v2)
```

Again, you should by now be able to work out what is going on here.

Of course, you can also give a function defintion as an argument to mapply.

## Drawing stuff in `R`

`R` contains many functions for standard plot types. These include (among others):

1. `plot.default`: x-y scatter plot
2. `barplot`: plot quantities as bars
3. `boxplot`: box-and-whisker plot(s) of the ranges of grouped values
4. `strpchar`: one dimensional scatter charts
5. `hist` : calculates and plots distributions as barplots
6. `dotchart`: plot quantites as points

These allow you to visualise numbers in fairly basic, not-so-pretty
ways. There are also many packages that provide much more complex and
appealing data visualisation (eg. see `ggplot2`). For a range of examples,
have a look at: <https://www.r-graph-gallery.com/all-graphs.html>.

Remember that `plot` is a generic function, and that calling the outcome of
`plot(X)` depends on the class of `X`. It can be almost anything, and this is
one of the reasons I will cover only the lower level drawing functions here.

There are countless examples on the internet detailing how to use these
different plot types and packages and I will not cover them here. Instead I
want to cover the low level graphics functions that are built in to `R` and
which allow you to build your own data visualisation. Using these is more akin
to *drawing with numbers* than to simply plotting values in more and more
pretty ways.


To draw in `R` you need to understand a rather small set of functions that
are used to either set up the plotting surface and it's coordinates, draw
things on that surface and to modify the way in which things are drawn:

1. Initialisation
   1. `plot.new` clear any previous plots
   2. `plot.window`  set up the coordinate system
2. Drawing functions
   1. `points` draw points
   1. `lines` draw lines
   2. `segments` draw lines specified in a different way
   3. `rect` draw rectangles
   4. `polygon` draw polygons
   5. `abline` draw lines specified by equations
   6. `text`, `mtext` draw text objects
   7. `legend` add a legend to a function
3. Modifying functions
   1. `par` changes graphics parameters in a large number of ways
   2. `rgb` specify colours using red-green-blue
   3. `hsv` specify colours using hue-saturation-value
   4. `layout` draw several plots on one page
  
When combined with fairly simple mathematical constructs you can, reasonably
efficiently draw more or less anything you like with this set of
functions. Note that many of the plot commands will pass arguments to `par`;
this means that you often have to look at the help page for the `par` function
in order to work out how to use the various drawing functions.

### Initialising a plotting area

To start a new plot or drawing simply call `plot.new()`. This will clear the
current graphics device (which may be a window on your screen, or a device
that draws the commands to a file). If you are using `RStudio` I would
probably recommend calling the `x11()` function first before calling
`plot.new()`. This will create a graphics device associated with a free
floating window that you can resize to your heart's content. The little
plotting tab that you get in `RStudio` is frequently too small for more
complicated graphs.

After calling `plot.new()` you can can call `plot.window()` to set up the
coordinate system. All that does is to define what the minimum and maximal
coordinates within your plotting system will be. Try the following:

```R
x11()
## but call x11 only once unless you know what you are doing, or if
## you inadvertently close the plotting window.

plot.new()
plot.window(xlim=c(0,100), ylim=c(0,100))

## let use draw some points at the indicated x and y positions:
points(x=c(0, 0, 100, 100), y=c(0, 100, 0, 100))

## the graphical parameter 'usr' defines the limits of the current plot
## window. Try:
par('usr')

## try the following:
plot.new()
plot.window(xlim=c(0,100), ylim=c(0,100), xaxs='i', yaxs='i')
points(x=c(0, 0, 100, 100), y=c(0, 100, 0, 100))
par('usr')
## do you see the difference?

## plot.window has an additional option asp:
## try
plot.new()
plot.window(xlim=c(0,100), ylim=c(0,100), asp=1)
points(x=c(0, 0, 100, 100), y=c(0, 100, 0, 100))

## then resize the window and see if you can work out what
## is different.
```

### Drawing points

To draw any number of points on a plot, simply specify their x and y
coordinates in two vectors:

```R
plot.new()
plot.window(xlim=c(0,100), ylim=c(0,100))

x <- sort(sample(1:100, 30))
y <- sort(sample(1:100, 30))
## if you don't understand the two lines above, then
## have a look at x and y. If you still don't understand
## then look at the help page for sample and sort and play
## around with them until things make sense

points(x, y)

## for points we can change the size of the points using the cex parameter:
points(x, y, cex=3) ## note how things are just added to the plot.

## we can specify different point sizes for every point:
points(x, y, cex=1:10/2)  ## the values are recycled..

## there are also different plotting characters (pch)
points(x, y, cex=1:10/2, pch=1:30)  ## the values are recycled..

## oops note the warnings.. seems only 25 different plotting characters
points(x, y, cex=1:10/2, pch=1:25)  ## the values are recycled..

## we can also specify the colours. There are many packages and functions
## that generate colour gradients; but we can also use the base hsv and rgb
## functions. Here the hsv one seems good..

points(x, y, cex=1:10/2, pch=1:25, col=hsv( 1:30/30, 0.8, 0.8 ))  ## the values are recycled..
## see if you can work out what is going on there:
## with hsv we specify hue (colour), saturation and the value (brightness) of each colour
## Each value should be between 0 and 1. hsv( x, 0, 1 ) is white for any value
## of x, whereas hsv(x, y, 0) is black for all values of x and y.

## plot characters are composed of lines. We can set the line width of these using
## lwd:
points(x, y, cex=1:10/2, pch=1:25, col=hsv( 1:30/30, 0.8, 0.8 ), lwd=3)

## all of the above is quite similar to simply doing:

plot(x, y, cex=1:10/2, pch=1:25, col=hsv( 1:30/30, 0.8, 0.8 ), lwd=3)
## except you then get a box and some axes and other things.

## note you could also do:
plot(x, y, cex=1:10/2, pch=1:25, col=hsv( 1:30/30, 0.8, 0.8 ), lwd=3, type='l')
points(x, y, cex=1:10/2, pch=1:25, col=hsv( 1:30/30, 0.8, 0.8 ), lwd=3)
## to gradually build up a plot in any way you like.

```

### Drawing lines with `lines` and `segments`

`lines` draws a line through a number of points specified by two vectors
specifiying the x and y coordinates. `segments` draws individual lines
specified by pairs of coordinates; it is useful when you want disconnected
lines specified by some set of numbers. Using the x and y coordinates we
created in the previous section (they will still be available as long as you
have not deleted them).

```R
plot.new()
plot.window(xlim=c(0,100), ylim=c(0,100))

## draw a line through x and y:
lines(x, y)

## lets make it a nicer colour and thickness
lines(x, y, col=hsv( 1:30/30, 0.8, 0.8 ), lwd=3)
## hmm that only gives us a single colour..

## lets put some points on the line
points(x, y, col=hsv( 1:30/30, 0.8, 0.8 ), lwd=3)

## to draw standard errors or deviations associated
## with these numbers:
## we can use segments..
p.se <- rnorm( length(x), sd=3 )

## speficy the x0, y0, x1, y1
## with lines drawn from (x0[i],y0[i]) to (x1[i],y1[i])
segments( x0=x, y0=y-p.se, x1=x, y1=y+p.se, lwd=3,
         col=hsv( 1:30/30, 0.8, 0.8 ) )
## hey we get individual colours again..

## actually if you really want to draw error bars you can use
## arrows.. it behaves very similarly

arrows( x0=x, y0=y-p.se, x1=x, y1=y+p.se, lwd=3,
         col=hsv( 1:30/30, 0.8, 0.8 ), angle=90, length=0.05, code=3 )
```

### Drawing rectangles with `rect`

`rect` draws rectangles specified in manner very similar to how lines are
specified in `segments`; the arguments are `left`, `bottom`, `right` and
`top`. We can use `rect` to make a sort of fancy barplot with the same data.

```R
plot.new()
plot.window(xlim=c(0,100), ylim=c(0,100))

## lets add a thick line to start with..
lines(x, y, lwd=3)
points(x, y, pch=3, lwd=3, cex=2 )

## we need to specify left and right coordinates. This is most easily done
## by using x[1:(n-1)] as the left coordinates and x[2:n] as the right coordinates
## we can make the tops of the rectangles the mean of neighbouring points
## (note that 1:(n-1) is the same as 2:n - 1
n <- length(x)
top <- (y[2:n-1] + y[2:n]) / 2

rect( x[1:(n-1)], 0, x[2:n], top )
## we can specify a transparent colour (the last 0.5 in the hsv command)
rect( x[1:(n-1)], 0, x[2:n], top, col=hsv( 1:30/30, 0.8, 0.8, 0.5 ) )

## note what happens when you repeat the rect call:
rect( x[1:(n-1)], 0, x[2:n], top, col=hsv( 1:30/30, 0.8, 0.8, 0.5 ) )
## the colours get more distinct.
```

### Drawing with `polygon`

We can do the same plot using the polygon function. This lets us draw a filled
polygon represening the area under the curve.

```R
plot.new()
plot.window(xlim=c(0,100), ylim=c(0,100))

## lets add a thick line to start with..
lines(x, y, lwd=3)

## for this we want to define new start and end points with 0 at the end.
n <- length(x)
polygon( c(x[1], x, x[n]), c(0, y, 0), col=hsv( 0.3, 0.8, 0.8, 0.5 ) )

## we can do something more fancy using hist and some more numbers.
h1 <- hist( rnorm( 10000, mean=0.5, sd=1 ) )
h2 <- hist( rnorm( 10000, mean=-0.5, sd=2 ))

## h1 and h2 are named lists. Check ?hist for more details. I will use some
## of the members of h1 and h2 here

plot.new()
plot.window( xlim=range( c(h1$mids, h2$mids)), ylim=c(0, max(c(h1$counts, h2$counts))))
n1 <- length( h1$mids )
n2 <- length( h2$mids )

polygon( c(h2$mids[1], h2$mids, h2$mids[n2]), c(0, h2$counts, 0), col=hsv(0.2, 0.8, 0.8, 0.5) )
polygon( c(h1$mids[1], h1$mids, h1$mids[n1]), c(0, h1$counts, 0), col=hsv(0.5, 0.8, 0.8, 0.5) )

## we can add some axes using axis:
axis(1)
axis(2)
axis(3)
axis(4)

## we can use legend to add a legend..
legend('topright', legend=c("data 1", "data 2"), fill=hsv( c(0.5, 0.2), 0.8, 0.8, 0.5 ))
```

### `abline` to draw equations

Use `abline` for drawing lines by equations:

```R
## 
plot.new()
plot.window(xlim=c(0,100), ylim=c(0,100))

## lets add some points..
points(x, y, lwd=3)

## use lm (linear model) to make a linear regression
lm1 <- lm(y ~ x)
abline(lm1)
## but the more general form is to specify
## a and b in
## y = a + bx
a <- lm1$coefficients[1]
b <- lm1$coefficients[2]
abline(a, b, col='red', lwd=2)

## you can also specify vertical and horizontal lines:
## and as for any line we can also specify the line type
## with lty
abline(h=seq(0, 100, 10), lty=2)
abline(v=seq(0, 100, 10), lty=2)

## there is also the grid command if that is what you
## want:
grid( col='blue', lwd=3)
```

#### Drawing text with `text` and `mtext`

As usual with R the most difficult thing to do well is to handle textual
data. The basic text command is `text`, and that takes a set of x and y
coordinates as well as a set of labels to draw and information about how to
draw the text.

```R
## again we will use our previously generated x and y coordinates:
plot.new()
plot.window(xlim=c(0,100), ylim=c(0,100))

## instead of drawing points, lets use letters..
n <- length(x)
text(x, y, letters[1:n] )
points(x, y) ## hmm. we don't have enough letters

## lets make some other labels;
lab <- paste("sample", 1:n, sep="_")

plot.new()
plot.window(xlim=c(0,100), ylim=c(0,100))
text(x, y, lab)

## that has a lot overlap. We can try to rotate the text with
## the parameter 'srt'
plot.new()
plot.window(xlim=c(0,100), ylim=c(0,100))
text(x, y, lab, srt=90)

## see what happens when you resize the window..
## the text size stays the same; R always tries to
## set the physical size of the text, whether that
## be on the screen or in a pdf.

## to change the size, use cex as usual:
plot.new()
plot.window(xlim=c(0,100), ylim=c(0,100))
text(x, y, lab, srt=90, cex=2)

## you can specify colour and others as usual.
## Use font to specify italic or bold, or both:
## bold
plot.new()
plot.window(xlim=c(0,100), ylim=c(0,100))
text(x, y, lab, srt=90, cex=2, font=2)

## italic
plot.new()
plot.window(xlim=c(0,100), ylim=c(0,100))
text(x, y, lab, srt=90, cex=2, font=3)

## bold italic
plot.new()
plot.window(xlim=c(0,100), ylim=c(0,100))
text(x, y, lab, srt=90, cex=2, font=4)

## if the font is available on your system you may be able to
## also specify a specific font using the family option:

plot.new()
plot.window(xlim=c(0,100), ylim=c(0,100))
text(x, y, lab, srt=90, cex=2, font=1, family="Times")

plot.new()
plot.window(xlim=c(0,100), ylim=c(0,100))
text(x, y, lab, srt=90, cex=2, font=1, family="Arial")
```

### Several plots on one page

There are several ways in which you can combine several plots on one page. The
simplest of these is to use `par(mfrow=c(nrow, ncol))` to split the page into
equally sized rows and columns. The `layout` function can do more complex
layouts where you can specify widths and heights of the rows as well as
merging the columns.

```R
## use par('mfrow'=c(..)) to set up a plotting surface:
## note that this will be set for your currently active graphics
## device.

par(mfrow=c(2,2))
hist(x, main="histogram of x", xlab="X-values" )
hist(y, main="histogram of y", xlab="Y-values" )
plot( x, y, type='l', lwd=3 )
plot( x, y, type='b', lwd=4, col=hsv(1:30/30, 0.5, 0.5, 0.8), pch=19, cex=3 )
```

To use `layout` you must specify a matrix as the first argument; that matrix
should contain values which specify the order in which the plots will be
plotted. Values of 0 will cause an empty section. 

```R
## our layout matrix:
m1 <- matrix(1:4, nrow=2, byrow=FALSE)
## look at m1!

layout(m1, widths=c(10, 5), heights=c(5, 10))

hist(x, main="histogram of x", xlab="X-values" )
hist(y, main="histogram of y", xlab="Y-values" )
plot( x, y, type='l', lwd=3 )
plot( x, y, type='b', lwd=4, col=hsv(1:30/30, 0.5, 0.5, 0.8), pch=19, cex=3 )

## you can also merge two or more panels:
m2 <- rbind(c(1,1,1), 2:4)
## look at m2!

## not terribly pretty, but hopefully you understand what is going on here:
layout(m2)
hist(x, main="histogram of x", xlab="X-values" )
hist(y, main="histogram of y", xlab="Y-values" )
plot( x, y, type='l', lwd=3 )
plot( x, y, type='b', lwd=4, col=hsv(1:30/30, 0.5, 0.5, 0.8), pch=19, cex=3 )
```

### A bit of trigonometry to specify positions

In all of the above we have only drawn straight lines and edges. We can use a
bit of simple trigonometry to draw circles, arcs and spirals.

```R
## to draw a circle we need to specify the x and y coordinates of a
## circle. That is sin(a) and cos(a) where a is an angle
## (I never remember which should be cos and sin, but it doesn't really
## matter here).

## angles in radians..
## a full circle
a1 <- seq(0, 2*pi, length.out=1000)

par(mfrow=c(1,1))
plot.new()
plot.window(xlim=c(0, 100), ylim=c(0,100), asp=1)

## specify some different radius values
r <- seq(10, 50, 10)
ori <- c(50,50)

## draws a circle
lines( ori[1] + r[1] * cos(a1), ori[2] + r[1] * sin(a1) )

## draws a bigger circle
lines( ori[1] + r[2] * cos(a1), ori[2] + r[2] * sin(a1) )

polygon( ori[1] + r[2] * cos(a1), ori[2] + r[2] * sin(a1), col=hsv(0.3, 0.7, 0.8, 0.3) )
polygon( ori[1] + r[1] * cos(a1), ori[2] + r[1] * sin(a1), col=hsv(0.3, 0.7, 1, 0.3) )

## in a different location:
polygon( 20 + r[1] * cos(a1), 80 + r[1] * sin(a1), col=hsv(0.3, 0.7, 1, 0.3) )
## draw a slice
a2 <- seq(0, 0.73*pi, length.out=1000)
polygon( 20 + r[1] * cos(a2), 80 + r[1] * sin(a2), col=hsv(0.3, 0.7, 1, 0.3) )

## to make an arc, simply connect to the origin..
polygon( c(20, 20 + r[1] * cos(a2), 20), c(80, 80 + r[1] * sin(a2), 80), col=hsv(0.8, 0.7, 1, 1) )

## to make a spiral simply specify more angles and an increasing radius
a3 <- seq(0, 10*pi, length.out=1000)
r2 <- seq(1, 50, length.out=1000)

lines( 50 + r2 * cos(a3), 50 + r2 * sin(a3), lwd=3 )
```

### Drawing summary

For most situations you should probably consider to make use of built in
plotting functions. However, if you know how to use `R`'s drawing primitives
it is possible for you to consider what *you want to do* rather than what
*functions are available*. And there are lots of situations where this can be
very useful. Although it can take some time to create your own plotting
functions you should not discount the time it takes to read and understand
package documentation.

## Additional stuff

`R` is a big complicated beast of a language. The above gives you some
material for getting started with it. There are many, many more useful things
within `R`, that we don't have time to go into. However, here are a few useful
functions or operators that you may want to look up in the future:

1. `%in%` Used as an operator, the expression `a %in% b` returns `TRUE` for
   every element of `a` which can be found in `b`.
2. `match(a, b)` Returns an index for values in `b` that are identical to the
   values in `a`. This function is useful, but confusing.
3. `sort` Returns a sorted vector (small to big by default).
4. `order` Returns a permutation (an index) that can be used to sort a
   vector. This is very useful if you want to sort a matrix or dataframe by
   one of its columns.
5. `%%` The modulus operator. Used as an operator `a %% b`, gives the integral
   remainder of the division of `a` by `b`.

