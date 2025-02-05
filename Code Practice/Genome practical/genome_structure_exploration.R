##
## Download genome information from Ensembl databases
## Visualise genome annotation at different genome loci
## and species
##
## Compile annotation statistics from different tables.

## probably you will need to:
install.packages("RMariaDB")
install.packages("DBI")

require(RMariaDB)
require(DBI)


## A function that connects to a public database:
## cred should be a named list including the following elements
## host, user, port, password.
ensembl.connect <- function(cred, db.name=NULL){
    if(!is.null(db.name))
        return( dbConnect(MariaDB(), user=cred$user, dbname=db.name, port=as.integer(cred$port),
                          host=cred$host, password='') )
    return( dbConnect(MariaDB(), user=cred$user, port=as.integer(cred$port),
                      host=cred$host, password='') )
    ## dbConnect(MySQL(), user=cred$user, host=cred$host,
    ##           password=cred$password, port=as.integer(cred$port), dbname=db.name )
}

## A function that connects to a database and performs a query
## cred:  a list as defined for ensembl connect
## query: a character vector specifying an SQL query
## db.name: the name of the database to query from (optional)
ensembl.query <- function(cred, query, db.name=NULL){
    db <- ensembl.connect(cred, db.name)
    data <- dbGetQuery(db, query)
    dbDisconnect(db)
    data
}

## From:
## https://metazoa.ensembl.org/info/data/mysql.html
##
##
## Mysql access from:
##
## Ensembl Genomes, all databases 	mysql-eg-publicsql.ebi.ac.uk 	4157
## Ensembl 	ensembldb.ensembl.org 	5306
## Ensembl Mart 	martdb.ensembl.org 	5316

## Define a list containing the host details for the all databases?
cred.all <- list(host='mysql-eg-publicsql.ebi.ac.uk', port=4157, user='anonymous', password="")

## to get a full list of databases do:
## This simply peforms the "show databases;" query at the
## database server specified by cred.all
db.all <- ensembl.query(cred.all, "show databases;")

## Have a look at the output of db.all!!

## don't execute the commented line. It was used by
## me to create a file that can be loaded by R.
## for you in case the db lookup didn't work.
## saveRDS( db.all, "db_all_1.rds")

## but do have a look at the databases. Try to work out what
## you get.
summary(db.all)
length(db.all) ## 1 ?
class(db.all) ## data.frame
dim(db.all)  ## [1] 7382    1
head(db.all)

## we have access to a total of 7382 different databases.
## That's a lot, but you may have noted that there is more than
## one version for each database.

## Genome annotation is normally present in the _core_ dbs.
## we can use grep to select:

## note that db.all is a data.frame, so we have to extract the first
## (and only) column from it in order to get a suitable data for grep
db.core <- grep('_core_', db.all[,1], value=TRUE)

summary(db.core)
## Length     Class      Mode 
##   5824 character character 

## we still have 5824 databses. To get only the latest version of the
## databases we can make use of tapply
## tapply can be used to split a vector into groups of elements based on an index
## Here we want to group all databases from the same species
## and select the one with the highest version number
## we can use the 'sub' function to create an index vector

## (i for index)
db.core.i <- sub("(.+?)_core_.+$", "\\1", db.core )

## look at the ouput!!
db.core.i
## In the sub statement above we match every element in db.core to:
## (.+?)_core_.+$
## .   matches any character
## +   means one or more of the preceding character class
## ?   don't be greedy (.+ could match the entire string, we don't want that)
## _core_  matches the exact string
## .+  matches one or more of anything.
##
## The brackets are used to capture the part of the string that matches;
## in this case we want to match everything before the '_core_' as that contains
## the species name.
##
## The second argument to 'sub' is what we want to substitute the string with;
## here we use, "\\1". This evaluates to whatever was captured in the first set
## of brackets. In this case the species name.

## We can use the resulting vector of species as an index to
## split the character vector by species.

## The function that we specify in the tapply expression
## will act on every group of elements
## that have a common index value (in this case, having the same species
## code). Hence the function needs to take as an argument a single
## group specified by the index. Here the function only returns the
## group itself. This is for practice only..

tmp <- tapply( db.core, db.core.i, function(x){ x })

## execute these one at a time!!!
summary(tmp)  ## prints too much stuff out
class(tmp) ## array ?
is.list(tmp) ## TRUE (so it's a list as well)
length(tmp)  ## 762
head(names(tmp), n=20) ## looks like species names

head(tmp)
## databases grouped by species ?

## instead of simply returning the argument in the function given to
## tapply, we can use it in any way we want...
tmp <- tapply( db.core, db.core.i, function(x){ length(x) })
head(tmp)  ## just a single number. The number of databases per species

## what we want is to get the database with the highest version number
## the easiest way to do this is to sort the elements in decreasing order
## and then picking the first element:

db.core.current <- tapply( db.core, db.core.i, function(x){
    sort(x, decreasing=TRUE)[1]
})

## but there is a problem with the above:
## xx_core_101_1 < xx_core_95_1
## so we need to modify the above. We need to extract the version numbers
## and treat them as numbers, not strings:

db.core.current <- tapply(db.core, db.core.i, function(x){
    ## x is now a character vector. We want to get the version numbers
    ## from it. We can use a combination of strsplit and sapply:
    v.no <- sapply( strsplit(x, "_"), function(y){
        l <- length(y)
        ## this will return the last three elements
        as.numeric( y[(l-2):l] )
    })
    v.no
    o <- order( v.no[1,], v.no[2,], v.no[3,], decreasing=TRUE )
    x[o[1]]
})

## the code above is probably a bit difficult to understand. If you have problems
## consider the following

x <- db.core[ db.core.i == db.core.i[1] ]
## this mimics what the tapply statement does. It will apply the specified
## function to sets of db.core grouped by the value of db.core.i. Here
## I'm making such a set; Have a look at the value of x to
## work out what is happening. If confused, look at the value of db.core.i[1]
## and so on.

## then consider what the following expressions do:
strsplit(x, "_")

## note the class of the object returned by strsplit
## (it's a list)
class( strsplit(x, "_"))

## You should by now understand what sapply does, and you should be able
## to guess what as.numeric does. But look at the output of:
sapply( strsplit(x, "_"), function(y){
        l <- length(y)
        ## this will return the last three elements
        as.numeric( y[(l-2):l] )
})

## compare to the value of x
## and see what the expression has done. The rest you should be able to guess
## and if not, you can always try to break up the individual pieces as above.

## if you made an object x, you may want to delete it at this point
## eg., rm(x)

## compare to db.core
head(db.core, n=30)
head(db.core.current)

length(db.core.current)
## 762 different genomes. Not bad.

## but we don't know what species we have.
## you can try googling the names given to the databases
## by the above code:

names(db.core.current)

## we can also look for specific species / genus names:

grep("homo", db.core.current, value=TRUE)
## nothing
grep("mus", db.core.current, value=TRUE)
## not the mouse I was looking for
grep("sacc", db.core.current, value=TRUE)
## at least we have some yeast
grep("danio", db.core.current, value=TRUE)

## but it doesn't look like we have any vertebrates.
## somehow the ensembl genomes doesn't actually have all the
## databases. We probably have to get these from the original
## ensembl site.
##
## to get the vertebrate databases we want to repeat all of the
## steps above.
##
## The whole point of programs are to avoid repeating yourself.
## Here we have identified a pattern of code that we may want
## to run more than once. Such patterns are better off embedded
## in a function. Here the function should simply take the
## database credentials needed to execute the query.
## We will return a list of three database lists
## all, core, core.current:

get.core.dbs <- function(cred){
    ## here I have simply copied and pasted the exact lines that we used above
    ## except for the variable name in the call to the ensembl.query function
    ## which is now 'cred' instead of 'cred.all'
    ## Since these are executed within a function they will not
    ## modify the global variables defined above (db.all, db.core, db.core.i,
    ## db.core.current)
    ## This is related to the 'scope' of the variables.
    db.all <- ensembl.query(cred, "show databases;")
    db.core <- grep('_core_', db.all[,1], value=TRUE)
    db.core.i <- sub("(.+?)_core_.+$", "\\1", db.core )
    db.core.current <- tapply(db.core, db.core.i, function(x){
        ## x is now a character vector. We want to get the version numbers
        ## from it. We can use a combination of strsplit and sapply:
        v.no <- sapply( strsplit(x, "_"), function(y){
            l <- length(y)
            ## this will return the last three elements
            as.numeric( y[(l-2):l] )
        })
        v.no
        o <- order( v.no[1,], v.no[2,], v.no[3,], decreasing=TRUE )
        x[o[1]]
    })
    list('all'=db.all[,1], 'core'=db.core, 'core.current'=db.core.current)
}

## we need then to construct a different cred list
## all that is different is the host and port number
cred.vert <- list(host='ensembldb.ensembl.org', port=5306, user='anonymous', password="")

db.vert <- get.core.dbs( cred.vert )
## Note the warning messages!
## see if you can work out:
## 1. Why we get the warnings
## 2. If the warnings matter
## to do this, look at the values of db.vert$core and db.vert$core.current
## and determine whether we are getting the latest version or not.
## I leave it to you to work out how to do that.

## to be consistent we can also create a data structure like db.vert for the
## db.all

db.all <- get.core.dbs( cred.all )
## and without warning messages.

## oops., that just overwrote our old db.all list. So now we have a bit of a mess
## with db.core, db.core.current, and db.all$all, db.all$core, db.all$core.current,
## db.vert$all, db.vert$core, db.vert$core.current. Quite a mess. We could consider
## tidying things up to clean up our workspace. Or better yet; tidying up the code that made
## the previous versions and then rerunning everything after restarting R (without saving
## the .RData file
## In a real project that would be good, because you always want to keep your code as simple
## as possible. Here we will ignore it for the time being.

## we can use the $ notation to subset named lists:
length(db.vert$core.current)  ## 336
db.vert$core.current


##########################################################################################
### We have some database names; we can now start to explore some of the
### the genomes
###
### Information about genes, transcripts and exons are found in the tables
###
### gene
### transcript
### exon
###
### in relational databases we use the 'SELECT' command to select from tables and
### combinations of tables. We can either select specific columns from the tables or
### we can select all columns using the glob character (*). We can limit the selections
### in various ways; the most simple is the 'limit' command.
### We will first have a look at some of the tables and see how we can do something useful:

## first select a database that you would like to look at
## to start with it is probably easier to look at a well studied vertebrate:
## you can specify a couple of core databases in a couple of variables:

db.dr <- db.vert$core.current['danio_rerio']
db.hs <- db.vert$core.current['homo_sapiens']
db.mm <- db.vert$core.current['mus_musculus']

### LOOK at what these variable contain!


### lets have a look at the first 10 rows of the gene table in mouse
tmp.gene <- ensembl.query(cred.vert, "select * from gene limit 10;", db.mm)
### you may get lots of warnings.. you can probably ignore.

## only 10 rows, we can look at the full table:
tmp.gene

## There are some columns that make sense:
## description, seq_region_start, seq_region_end, seq_region_strand
## these all give us text or numbers that make sense.
## but what about,
## analysis_id, seq_region_id, display_xref_id?

## We have coordinates for the begin and end of the gene, but
## we don't know on which chromosome:
##
## If you go go www.ensembl.org
## and search using one of the entries in the 'stable_id' column you can find
## out more...
##
## On the web page you will get a load of additional information;
## there is a symbol, synonyms and a whole lot more.
## This information (or almost all of it is inside the database) and we will
## obtain it later on.

## One of the columns gives a biotype. Clearly there will be more than one
## We can look at how many types by selecting unique elements from that columns.

tmp.biotype <- ensembl.query(cred.vert, "select distinct biotype from gene", db.mm)

## look at what you got...
tmp.biotype
## clearly a log of different ones. We can now have a look at a specific biotype.

tmp.gene <- ensembl.query(cred.vert, "select * from gene where biotype='protein_coding' limit 10", db.mm)

## and again
## look at it:
tmp.gene

## we can get the same from a different species by simply specifying a different database:
tmp.gene <- ensembl.query(cred.vert, "select * from gene where biotype='protein_coding' limit 10", db.dr)

tmp.gene

##### Gene location #####
### In the gene table you have the columns:
##
## seq_region_start
## seq_region_end
##
## these tell you the beginning and end of the gene
## but they do not tell you which chromosome; for that we just get a number:
##
## seq_region_id
##
## That number is the id for the seq_region in a different table; The reason for having that
## information in a different table is because we wish to store more information about the seq_region
## than just its name; and we do not wish to look up a seq region by a name as that is less efficient than
## using a number.
##
## You might be able to guess that the name of the table containing this information is 'seq_region' and that
## it has a column called seq_region_id.
##
## To get a table containing the names of the seq_regions and the genes we have to do a 'join' of these
## two tables.
##
## first we can have a look at the seq_region table:

ensembl.query(cred.vert, "describe seq_region;", db.dr)
## this gives us:
##
##             Field             Type Null Key Default          Extra
## 1   seq_region_id int(10) unsigned   NO PRI    <NA> auto_increment
## 2            name     varchar(255)   NO MUL    <NA>               
## 3 coord_system_id int(10) unsigned   NO MUL    <NA>               
## 4          length int(10) unsigned   NO        <NA>               

## from which we can see that we have a column name and length. These seem useful
## we also have a column called coord_system_id; which suggests that there is an
## additional table called coord_system that has more information about the seq_regions

### To do a join we need to specify some criterion that should be the same for
### every row of the table. If we don't do that, we will end up selecting a combination
### of every row of the two tables. In the absence of any safeguards to prevent this
### the database engine may very well crash, or the server could run out of memory.
###
### Having direct access to SQL databases is thus something that is a bit dangerous
### hopefully Ensembl have some safeguards on their databases...

## To join the gene and seq_region tables we will make use of aliases for the
## two tables; I will use a for gene and b for seq_region. I will also specify
## the columns that I am interested in.
##
## I will also use paste to build the query in several steps as I think this
## makes it easier to use:

query <- paste("select a.biotype, b.name, a.seq_region_start, a.seq_region_end, a.seq_region_strand from",
               "gene a",
               "inner join seq_region b on a.seq_region_id=b.seq_region_id",
               "limit 20;")
## again you should look at what query is at this point..

tmp.gene <- ensembl.query(cred.vert, query, db.dr)

## here you can try:
sapply( tmp.gene, class )
## this works, because ensembl.query returns a dataframe, and that is
## treated as a list of lists, with one element for each entry.
## Note that class of the "name" column.

## Try to modify the above query to include the description column.

## We can also get a gene symbol if we make use of the display_xref_id column
## this is not very obvious.. but we can expand the above query to include a symbol:

query <- paste("select c.display_label, a.biotype, b.name, a.seq_region_start, a.seq_region_end, a.seq_region_strand from",
               "gene a",
               "inner join seq_region b on a.seq_region_id=b.seq_region_id",
               "inner join xref c on a.display_xref_id=c.xref_id",
               "limit 20;")

tmp.gene <- ensembl.query(cred.vert, query, db.dr)

tmp.gene

## things are starting to make sense. This gives us enough to download a full gene table.
## we can do this from a couple of species:
##
## lets expand the query to also include the stable id and the description of the gene:
## and we can also change the column headers by using the as keyword:

query <- paste("select a.stable_id as id, c.display_label as gene, a.biotype, b.name as chr, a.seq_region_start as start,",
               "a.seq_region_end as end, a.seq_region_strand, a.description as strand from",
               "gene a",
               "inner join seq_region b on a.seq_region_id=b.seq_region_id",
               "inner join xref c on a.display_xref_id=c.xref_id",
               "limit 20;")

tmp.gene <- ensembl.query(cred.vert, query, db.dr)

## again: look at what tmp.gene, and see how it differs from the previous version

## that looks OK. Lets remove the limit statement from the query:
query <- paste("select a.gene_id, a.stable_id as id, c.display_label as gene, a.biotype, b.seq_region_id, b.name as chr, a.seq_region_start as start,",
               "a.seq_region_end as end, a.seq_region_strand, a.description as strand from",
               "gene a",
               "inner join seq_region b on a.seq_region_id=b.seq_region_id",
               "inner join xref c on a.display_xref_id=c.xref_id;")

## dr.gene, because we are selecting from the db.dr database (i.e. danio rerio)
## note: this may take a little time
dr.gene <- ensembl.query(cred.vert, query, db.dr)

dim(dr.gene)
##
## [1] 24329    10
## Interestingly, last year, we got:
## [1] 37140     10
##
## so, the number of genes has gone down!

head(dr.gene)

## and the beauty of having a database means that I can also get the human and mouse
## gene information. I will call these hs.gene and mm.gene (for the obvious reasons).

hs.gene <- ensembl.query(cred.vert, query, db.hs)
mm.gene <- ensembl.query(cred.vert, query, db.mm)

dim(hs.gene)
## [1] 49635    10
## Last year the numbers were:
## [1] 49100     10 !!!
## so the human gene count has gone up?

dim(mm.gene)
## [1] 56651    10
## Last year's result:
## [1] 56631     10
## The mouse count has also increased.
##
## You should be sceptical of these counts!

## So many more genes from mouse and human than danio rerio (zebra fish).
## what's going on. Well first we can have a look at the counts of different
## biotypes.
## R has a wonderful function called table that does this:
## 

## do:
table(dr.gene$biotype)

## wrapping the call to table in a call to sort( )
## can help us to see which classes are large:

sort(table(dr.gene$biotype))
sort(table(hs.gene$biotype))
sort(table(mm.gene$biotype))

## Note that the number of processed pseudogenes is much, much lower
## in D. rerio. This might an annotation artefact; but it could also
## be related to teleost genome biology. I _don't_ know, but I would
## suspect that someone will have published a story about this.
##
## Because there is a good reason why this might be the case.

## bit difficult to just look at these. We could try something like:
barplot( table(dr.gene$biotype), horiz=TRUE, las=2 )
## hmm, we need more space for the labels. Lets try again..

## to see the current margings (in numbers of lines)
par('mar')

## to set:
## (note that you may need to make the plot window bigger;
## you can try x11()
## or clicking some magnifying glass somewhere..
par('mar'=c(5.1, 20.1, 4.1, 2.1))
## Do one at a time!
barplot( table(dr.gene$biotype), horiz=TRUE, las=2, main="Danio rerio" )
barplot( table(mm.gene$biotype), horiz=TRUE, las=2, main="Mus musculus" )
barplot( table(hs.gene$biotype), horiz=TRUE, las=2, main="Homo sapiens" )

## this gives us a better idea as to why we see so many genes. 
## Last year, we obtained much, much larger numbers of protein
## coding genes. 
## The reason was related to the presence of alternate chromosomes.
## These are still there, but and you can see it from:

sort( table( dr.gene$chr ) )
sort( table( mm.gene$chr ) )
sort( table( hs.gene$chr ) )

## There appear to lots of weird chromosome names; in the case of Danio rerio
## these have the name CHR_ALT_... for human they seem to contain, HG or CHR
## Well for mouse we don't have any strangeness.
## These (at least in the case of Danio rerio) represent alternative chromosomes
## containing mutations of sorts. Well, at least that is my guess.
## we can remove these from the counts by doing..

## compare the output of:
barplot( table(dr.gene$biotype), horiz=TRUE, las=2, main="Danio rerio" )
## with
barplot( table(dr.gene$biotype[ grep("CHR", dr.gene$chr, invert=TRUE) ]), horiz=TRUE, las=2, main="Danio rerio" )

## Last year that still gave us more than 25,000 protein coding genes for Danio rerio. More than I told you on
## in the lecture, but actually consistent with what we find on:
## http://www.ensembl.org/Danio_rerio/Info/Annotation
##
## So as I said, gene numbers can go down as well as up
## here they have gone up and then down since when I made the slide.
##

## Lets do the same for human (no point for mouse)
## note that the | symbol in the regular expression means OR
barplot( table(hs.gene$biotype[ grep("CHR|HG", hs.gene$chr, invert=TRUE) ]), horiz=TRUE, las=2, main="Homo sapiens" )

## so actually of these three human would appear to have the fewest protein coding
## genes. Remember that these numbers can change.

## It's quite a pain to keep writing grep("CHR", mm.gene$chr, invert=TRUE)
## so instead we can define some logical vectors
##
## grep returns the indices of the rows we wish to keep.
## grepl returns a logical vector of TRUE / FALSE values instead
##
## this is more useful since we can easily combine several logical vectors using
## boolean logic.

## then we have to think of a name for these vectors..
## nchr for normal chromosome ?
dr.nchr <- !grepl("CHR", dr.gene$chr)
mm.nchr <- !grepl("CHR", mm.gene$chr)
hs.nchr <- !grepl("CHR|HG", hs.gene$chr)

## the ! in front of the grepl statement means NOT
## it converts TRUE values to FALSE and vice versa.


## to see what proportion are TRUE for the above:
## (when you call sum() on a logical vector R will treat TRUE
## values as 1, and FALSE values as 0)
sum(dr.nchr) / length(dr.nchr) ## and similar for the others


## check that that worked..
sort( table( dr.gene$chr[dr.nchr] ))
sort( table( mm.gene$chr[mm.nchr] ))
sort( table( hs.gene$chr[hs.nchr] ))

## we still have some weirdness
## these are probably small contigs or scaffolds that could
## not be assembled into the main chromosome sequences.
## We can ignore them for the time being as the numbers of rows
## (and hence genes) are small.

### having these tables we can consider to visualise the location of
### the genes on the genome. To do this we want to know the lengths of the
### chromosomes. This we can get from seq_region. In doing this I will pull
### in a little bit of extra information from the coord_system table
### (hopefully you remember that coord_system_id column in the seq_region table).

query <- paste("select a.seq_region_id, a.name, a.length, b.name as type, b.rank",
               "from seq_region a",
               "inner join coord_system b on a.coord_system_id=b.coord_system_id;")

dr.chr <- ensembl.query(cred.vert, query, db.dr)
mm.chr <- ensembl.query(cred.vert, query, db.mm)
hs.chr <- ensembl.query(cred.vert, query, db.hs)

### LOOK at these tables to get a feeling for the information that they hold!

## we may want to access the chromosmes by their id at some
## stage in the future. To make this simple we can set the
## rownames attribute:

rownames( dr.chr ) <- dr.chr$seq_region_id
rownames( mm.chr ) <- mm.chr$seq_region_id
rownames( hs.chr ) <- hs.chr$seq_region_id

#### WARNING:
#### if you do
head( rownames(dr.chr) )

### You will see numbers in quotes. The quotes means that the numbers
### are actually strings. This is unfortunate, because if you try to
### use the seq_region_id from the gene tables R may treat those
### as numbers and select rows by their numeric value. That might work
### in some of the cases here, but it is not guaranteed to work. To overcome
### this we need to make sure that R treats thes as 'names' and not 'numbers'
### to do this we will later use as.character( ) to coerce numbers into strings

## have a look at these tables as you would usually do;
## then:

table( dr.chr$type )
table( mm.chr$type )
table( hs.chr$type )

## what a lot of chromosomes..
## you can also try:

plot( sort( dr.chr$length ))
plot( sort( mm.chr$length ))
plot( sort( hs.chr$length ))

## you can try to colour by type:
## but then we have to reorder the type values as well
## to do this we can define an order using the order()
## function:
o <- order(dr.chr$length)
plot( dr.chr$length[o], col=as.factor(dr.chr$type[o]) )

## you can do the same with the others as well.

## what you are seeing here are the issues related to trying to deal with assemblies and
## variant assemblies. These things are a bit of a pain.
## clearly if we do the analysis we will need to not select some sort of chromosomes
## In a similar way to the use of the nchr used above.

### A function to draw a genome...
## first you have to consider what the function needs to know
## here the gene locations, the chromosome lengths and some way to select the chromosomes
## that we are interested in should suffice.
## we will also only use rank == 1
## the arguments:
## genes:   a dataframe with gene information as obtained above
## chr:     a dataframe containing seq_region information
## min.length: the minimum length that a region needs to be included in the drawing
## chr.type:   the type of region
## exclude:    a pattern to exclude regions that we don't want
## main:       the title for the plot
## biotype:    specific types of genes to include
## border:     draw the border around gene representations; specify the colour of the border or NA
draw.genome <- function(genes, chr, min.length, chr.type="chromosome", exclude="CHR", main="Genome",
                        biotype=NULL, border=NA){
    ## a boolean for which chromosomes we will include
    ## & means AND
    ## I tend to use b for booelan; chromosome boolean
    chr.b <- chr$type == chr.type & chr$length >= min.length & !grepl(exclude, chr$name) & chr$rank == 1
    chr <- chr[chr.b,]
    ## the expression: 'a %in% b'
    ## allows us to check which members of set a are present in set b
    genes <- genes[ genes$chr %in% chr$name, ]
    if(!is.null(biotype))
        genes <- genes[ genes$biotype %in% biotype, ]
    ## we then need to work out the y coordinates for the chromosomes
    ## We can sort the chromosomes by length?
    chr <- chr[ order( chr$length ), ]
    ## and then make a named vector that can be used to look up positions
    y <- 1:nrow(chr)
    names(y) <- chr$name
    ## we now know how much space we will need. lets set up a plotting surface:
    plot.new()   ## just makes a new plot
    ## and then the coordinates:
    plot.window( xlim=c(0, max(chr$length)), ylim=c(-2, max(y)+2), xaxs='i', yaxs='i' )
    ## draw the chromosomes as lines using the segments function
    ## segments( x0, y0, x1, y1)
    segments( 1, y, chr$length, y )
    ## and then draw the genes as rectangles..
    ## here we have to use the name of the chromosome to look up the y positions
    ## we can do this once only..
    genes.y <- y[ genes$chr ]
    rect( genes$start, genes.y - 0.2, genes$end, genes.y + 0.2, col='blue', border=border )
    mtext(chr$name, side=2, at=y, las=2, line=1)
    mtext(main, side=3, line=0, cex=2)
}

draw.genome(dr.gene, dr.chr, min.length=1e6 )
draw.genome(mm.gene, mm.chr, min.length=1e6 )
draw.genome(hs.gene, hs.chr, min.length=1e6 )

## Try to change the text at the top of the plot to something more
## reasonable.

## also consider the output of:
par("mar")
## Would it make sense to change the value of "mar"?

## There are a few interesting things that you should be able to work out from these plots
## discuss amongst each other.. 

## we can do plots for specific biotypes:
draw.genome(dr.gene, dr.chr, min.length=1e6, biotype='protein_coding')
draw.genome(dr.gene, dr.chr, min.length=1e6, biotype='ribozyme') ## nothing?
draw.genome(dr.gene, dr.chr, min.length=1e6, biotype='snRNA') ## so few?
draw.genome(dr.gene, dr.chr, min.length=1e6, biotype='rRNA') ## so few?
draw.genome(dr.gene, dr.chr, min.length=1e6, biotype='antisense') 
draw.genome(dr.gene, dr.chr, min.length=1e6, biotype='miRNA') ## so few locations?

## to work out why so few locations for some of these you can try:
debug(draw.genome)
draw.genome(dr.gene, dr.chr, min.length=1e6, biotype='miRNA')
undebug(draw.genome)

## and maybe try with
draw.genome(dr.gene, dr.chr, min.length=1e6, biotype='miRNA', border='black')
## ahh, too short to be seen without a border.

draw.genome(mm.gene, mm.chr, min.length=1e6, biotype='miRNA', border='black')
draw.genome(hs.gene, hs.chr, min.length=1e6, biotype='miRNA', border='black')

### so we can see some interesting things there.

## we can have a look at the distribution of gene lengths in the different species:
## we can use the with()
## function here. It imports the namespace of named lists and dataframes
## this means that instead of writing..

hist( dr.gene$end - dr.gene$start, breaks=50 )

## I can write:
with(dr.gene, hist(end - start, breaks=50))
with(mm.gene, hist(end - start, breaks=50))
with(hs.gene, hist(end - start, breaks=50))

## probably better to do:
with(dr.gene, hist(log10(end - start), breaks=50))
with(mm.gene, hist(log10(end - start), breaks=50))
with(hs.gene, hist(log10(end - start), breaks=50))

## and there you should see something quite interesting.

## to plot these on the same plot we can do:
par(mfrow=c(2,2))
with(dr.gene, hist(log10(end - start), breaks=50, main="Danio rerio"))
with(mm.gene, hist(log10(end - start), breaks=50, main="Mus musculus"))
with(hs.gene, hist(log10(end - start), breaks=50, main="Homo sapiens"))

## These do not actually look as one would expect. What might
## be the reason for this?
## What do you think is going on here. What about biotypes?
## Or maybe we should restrict ourselves to specific chromosomes like
## we did in the drawing function. There are lots of options.

### We can also consider the proportion of the genome that is genic (containing
### genes) and that which is intergenic. Since there is some uncertainty
### about which chromosomes that we want to use we can instead
### calculate the values for all individual chromosomes. Afterwards we can
### select those that we are interested in. This means that we can also see
### if there are chromosomes that are different from each other:

## To perform operations on all genes for a specific chromosome we can
## use tapply...
## note that we use dr.gene$seq_region_id for the index rather than the
## chr name (in the chr column). This is because the name is unfortunately
## not unique.. i.e. it appears more than once..
## using an index for the rows:

dr.gene.chr.g <- tapply( 1:nrow(dr.gene), dr.gene$seq_region_id, function(i){
    with(dr.gene, sum(end[i] - start[i])) })
## that gives us the total genic region per chromosome

## have a look
head(dr.gene.chr.g)
length(dr.gene.chr.g)
head(names(dr.gene.chr.g))
## note the quotes around the numbers. This means that R should treat these
## as character vectors rather than numbers. But R might decide to convert
## them to numbers. If it does we would hope to get some errors.
## Maybe we would get errors because:

nrow(dr.chr)
range(dr.chr$seq_region_id)
## and
range(dr.gene$seq_region_id)
## but since 36565 < 36615
## R might still do some stupid conversion.
## in order to check we should check some of the numbers coming out.

head(dr.gene.chr.g)
head( dr.chr[ names(dr.gene.chr.g), ] )
which( dr.chr$seq_region_id == 31640 ) ## 32767
### these suggest that R has not done anything strange

## to get the ratios:
hist( dr.gene.chr.g / dr.chr[ names(dr.gene.chr.g), 'length' ] )

## but it's probably more informative to look at log ratios:
hist( log10(dr.gene.chr.g / dr.chr[ names(dr.gene.chr.g), 'length' ]), breaks=50 )

## we can also plot against length. Again, probably better to log transform:
plot( log10( dr.chr[ names(dr.gene.chr.g), 'length' ] ),
     log10(dr.gene.chr.g / dr.chr[ names(dr.gene.chr.g), 'length' ]) )

## These results are a bit weird. Try to think what might be causing
## the observation and how to deal with it.

## We can also make a subset of the rows that we are interested in.
## Generally we are more interested in seq_regions that have a rank of 1. We can get these by
## defining a logical vector:

dr.r1.b <- dr.chr[ names(dr.gene.chr.g), 'rank' ] == 1
sum(dr.r1.b)  ## 816
## That's still a tol of contigs. Lets select that the type should be chromosome as well

dr.chr.b <- dr.chr[ names(dr.gene.chr.g), 'type' ] == 'chromosome'
sum(dr.chr.b)  ## 159

### maybe we can remove the ALT_CHR instead
dr.chr.b <- !grepl("CHR", dr.chr[ names(dr.gene.chr.g), 'name' ])

## check the intersection of rank 1 and dr.chr.b:
sum(dr.chr.b) ## 159
sum(dr.chr.b & dr.r1.b) ## 26 .. I like this one..

## for convenience define a short term variable:
b <- dr.chr.b & dr.r1.b
## and use to select specific names.. 
plot( log10( dr.chr[ names(dr.gene.chr.g[b]), 'length' ] ), log10(dr.gene.chr.g[b] / dr.chr[ names(dr.gene.chr.g[b]), 'length' ]) )

## here a linear scale may actually work better... 
plot(  dr.chr[ names(dr.gene.chr.g[b]), 'length' ], dr.gene.chr.g[b] / dr.chr[ names(dr.gene.chr.g[b]), 'length' ] )

## you can add labels to your graph using the text function:
## text(x, y, labels)
text( dr.chr[ names(dr.gene.chr.g[b]), 'length' ], dr.gene.chr.g[b] / dr.chr[ names(dr.gene.chr.g[b]), 'length' ],
     dr.chr[ names(dr.gene.chr.g[b]), 'name'] )

## in which case it can be nicer to do:
plot(  dr.chr[ names(dr.gene.chr.g[b]), 'length' ],
     dr.gene.chr.g[b] / dr.chr[ names(dr.gene.chr.g[b]), 'length' ], type='n' )
## which sets up the plotting surface but doesn't plot anything
text( dr.chr[ names(dr.gene.chr.g[b]), 'length' ], dr.gene.chr.g[b] / dr.chr[ names(dr.gene.chr.g[b]), 'length' ],
     dr.chr[ names(dr.gene.chr.g[b]), 'name'] )

## To make the plot easier to look at, try setting the x and y labels to something
## more readable. Use the "xlab" and "ylab" arguments to plot()

## the plot should be reasonable meaningful to you.
## We can of course do the same thing for the human and mouse tables:
mm.gene.chr.g <- tapply( 1:nrow(mm.gene), mm.gene$seq_region_id, function(i){
    with(mm.gene, sum(end[i] - start[i])) })
## that gives us the total genic region per chromosome

hs.gene.chr.g <- tapply( 1:nrow(hs.gene), hs.gene$seq_region_id, function(i){
    with(hs.gene, sum(end[i] - start[i])) })
## that gives us the total genic region per chromosome

## at this point you may notice that we are kind of repeating ourselves here.
## This should make you think: 'maybe I should create some functions that create
## these data structures'.
##
## but for now we will plough on making a bigger and bigger mess...

## lets look at the density distributions:
par(mfrow=c(2,2))
hist( log10(dr.gene.chr.g / dr.chr[ names(dr.gene.chr.g), 'length' ]), breaks=50, main='Danio rerio', xlab='gene density' )
hist( log10(mm.gene.chr.g / mm.chr[ names(mm.gene.chr.g), 'length' ]), breaks=50, main='Mus musculus', xlab='gene density' )
hist( log10(hs.gene.chr.g / hs.chr[ names(hs.gene.chr.g), 'length' ]), breaks=50, main='Homo sapiens', xlab='gene density' )

## lets define some subsets as before..
mm.r1.b <- mm.chr[ names(mm.gene.chr.g), 'rank' ] == 1
mm.chr.b <- !grepl("CHR", mm.chr[ names(mm.gene.chr.g), 'name' ])

hs.r1.b <- hs.chr[ names(hs.gene.chr.g), 'rank' ] == 1
hs.chr.b <- !grepl("CHR|HG", hs.chr[ names(hs.gene.chr.g), 'name' ])

dr.b <- dr.r1.b & dr.chr.b
mm.b <- mm.r1.b & mm.chr.b
hs.b <- hs.r1.b & hs.chr.b

## to get the number of regions passing the filters:
sum(dr.b) 
sum(mm.b)
sum(hs.b)

## in which case it can be nicer to do:
par(mfrow=c(2,2))

plot(  dr.chr[ names(dr.gene.chr.g[dr.b]), 'length' ], dr.gene.chr.g[dr.b] / dr.chr[ names(dr.gene.chr.g[dr.b]), 'length' ],
     type='n', main='Danio rerio', xlab='chr length', ylab='gene density' )
text( dr.chr[ names(dr.gene.chr.g[dr.b]), 'length' ], dr.gene.chr.g[dr.b] / dr.chr[ names(dr.gene.chr.g[dr.b]), 'length' ],
     dr.chr[ names(dr.gene.chr.g[dr.b]), 'name'] )

plot(  mm.chr[ names(mm.gene.chr.g[mm.b]), 'length' ], mm.gene.chr.g[mm.b] / mm.chr[ names(mm.gene.chr.g[mm.b]), 'length' ],
     type='n', , main='Mus musculus', xlab='chr length', ylab='gene density' )
text( mm.chr[ names(mm.gene.chr.g[mm.b]), 'length' ], mm.gene.chr.g[mm.b] / mm.chr[ names(mm.gene.chr.g[mm.b]), 'length' ],
     mm.chr[ names(mm.gene.chr.g[mm.b]), 'name'] )

plot(  hs.chr[ names(hs.gene.chr.g[hs.b]), 'length' ], hs.gene.chr.g[hs.b] / hs.chr[ names(hs.gene.chr.g[hs.b]), 'length' ],
     type='n', main='Homo sapiens', xlab='chr length', ylab='gene density' )
text( hs.chr[ names(hs.gene.chr.g[hs.b]), 'length' ], hs.gene.chr.g[hs.b] / hs.chr[ names(hs.gene.chr.g[hs.b]), 'length' ],
     hs.chr[ names(hs.gene.chr.g[hs.b]), 'name'] )

## Do the observations fit our expectations? Consider the number of genes and the
## size of the respective genomes.

##################################################################################
#### A comment about the code used here  #########################################
###
### By this time you should start to get a bit confused about the large
### numbers of variables (data objects) that we have defined. Which one
### does what, and am I about to clobber an old variable when I define a new
### one.
###
### you should also start to feel that the code feels very repetitive. Many of the
### statements look very similar with only small differences (eg. mm instead of dr).
###
### This is the sign of bad coding discipline; such statements are _very_ easy to screw
### up. In general writing the code as I have done here is bad practice. However, it is
### fairly common when exploring data; usually you don't know enough about the data beforehand
### in order to know exactly what you want to do. This makes it rather common to
### split analyses into an exploratory phase and a final analytical phase done
### completely separately and with much better coding discipline.
###
### So how to improve the code. In R, we have two primary ways to make the code better:
###
### 1. Create named functions for actions that you may use more than once. Sometimes it is
###    even nice to create a function for a single use purpose, since this allows you to
###    use a new namespace where you can create lots of new variables without worrying that
###    you will be overwriting variables in the parent name space.
###
### 2. Store data as part of collections; in practice this usually means that instead of having
###    a separate variable for every data object, I instead have a variable that contains a
###    list of several data objects. For example, instead of having:
###
###    mm.gene, hs.gene, dr.gene
###
###    I make a list containing those three elements. At this point I could do that by:
###    gene <- list('mm'=mm.gene, 'hs'=hs.gene, 'dr'=dr.gene)
###    but it would be better to make the list from the beginning.
###
###    Once you have a list, you can either, loop through the list calling code on each element of the list
###    eg.:
###    for(i in 1:length(gene)){
###         some_function( gene[[i]] )
###    }
###
###    Or use lapply() to apply some code to every element of the list:
###    lapply( gene, function(x){  <do something to x> }
###
### for example, to get the mm.gene, and so on as a list, I could do the following:
###

db.sel <- db.vert$core.current[ c('danio_rerio', 'homo_sapiens', 'mus_musculus') ]

## but what is the value of query at this time?
## we have to reset the value of query; this is ugly, there should be a better way?
query <- paste("select a.gene_id, a.stable_id as id, c.display_label as gene, a.biotype, b.seq_region_id, b.name as chr, a.seq_region_start as start,",
               "a.seq_region_end as end, a.seq_region_strand, a.description as strand from",
               "gene a",
               "inner join seq_region b on a.seq_region_id=b.seq_region_id",
               "inner join xref c on a.display_xref_id=c.xref_id;")

gene <- lapply( db.sel, function(x){ ensembl.query(cred.vert, query, x ) })

## check that it worked
sapply( gene, dim )

### For the plotting you can then simply do, lapply(gene, function(x){ some_plottting_function(x) })
### instead of rewriting stuff all the time.
###
### We may cover some of this at a later point. For now we will act in exploratory mode and
### make a fantastic mess.

##############################################################################################
#### Genes --> transcripts --> exons ##############
####
#### Every gene can have mutiple transcripts; each transcript contains a selection of
#### the exons associated with that gene. That is a transcript contains many exons, but
#### exons can also be found in many transcripts. This is what is known as a many to many
#### relationship. These can be a bit of a pain to represent efficiently in code. Here
#### we will do the simplest thing we can; which is also enormously inefficent, and
#### really horrible when you think about it.
####
#### We will make two tables:
#### transcript, containing gene id, transcript id and the start and end points
#### exon containing the exons and their start and end points for every
###  transcript
###
### this is a bit monstrous, because it means that if two transcripts share
### a number of exons, then those exons will be repeated in the table.

## here we don't care about the chromosome name. We will just use the seq_region_id
tr.query <- paste("select a.gene_id as g_id, b.transcript_id as tr_id, a.stable_id as gene, b.stable_id as transcript,",
                  "b.biotype,",
                  "b.seq_region_id, b.seq_region_start, b.seq_region_end, b.seq_region_strand",
                  "from gene a",
                  "inner join transcript b on a.gene_id=b.gene_id;")

## now I can use my db.sel 

transcript <- lapply( db.sel, function(x){ ensembl.query(cred.vert, tr.query, x) })

sapply( transcript, dim )
## looks good. And it's quite reasonable..
## use head as well; but note that head(transcript) probably
## isn't a good idea. Try head(transcript[[1]]) instead.
## What do YOU think?

## how about:
names(transcript)

## to get the exons is much more complicated.
ex.query <- paste("select a.transcript_id, a.stable_id as transcript, b.rank,",
                  "c.exon_id, c.stable_id as exon, c.seq_region_id, c.seq_region_start, c.seq_region_end, c.seq_region_strand",
                  "from transcript a",
                  "inner join exon_transcript b on a.transcript_id=b.transcript_id",
                  "inner join exon c on b.exon_id=c.exon_id",
                  "order by a.transcript_id, b.rank;")

## again I use my db.sel list
## Note: this may take some time to run.. 
exon <- lapply( db.sel, function(x){ ensembl.query(cred.vert, ex.query, x) })

## and again:
sapply( exon, dim )

## and have a look using head, tail, whatever you can think of.

## We can quite easily ask how many exons are associated with individual transcripts
## using the exon table. The number of exons derived from a given gene
## will be equal to the max rank in the table. So we can plot the distribution 
## we can use a tapply wrapped in an lapply to do this for the full set of species:

par(mfrow=c(2,2))
lapply( exon, function(x){
    ## note that we don't need to define a function for
    ## the tapply statement as we can simply use max:
    hist(tapply( x$rank, x$transcript_id, max ), breaks=50)
})

### there is some weirdness there. Maybe we should modify to look for specific
### biotypes? To consider only protein coding we need to consider both the
### transcript and exon tables. We can use mapply for this. In mapply, the function
### comes first.. followed by the vector whose elements are used in the function:

par(mfrow=c(2,2))
mapply( function(tr, ex){
    pc.b <- tr$biotype == 'protein_coding'
    ex.b <- ex$transcript_id %in% tr$tr_id[ pc.b ]
    hist( tapply( ex$rank[ex.b], ex$transcript_id[ex.b], max ), breaks=50)},
    transcript, exon)

### sometimes a quantile plot is clearer:
par(mfrow=c(2,2))
mapply( function(tr, ex){
    pc.b <- tr$biotype == 'protein_coding'
    ex.b <- ex$transcript_id %in% tr$tr_id[ pc.b ]
    plot( log2( sort(tapply( ex$rank[ex.b], ex$transcript_id[ex.b], max )) ))},
    transcript, exon)

### having obtained these datatables we can now consider how to visualise the
### transcripts associated with a given gene.
### for this you need to have some way of finding the name of a gene that you wish to specify.

### you can look at:
head(dr.gene$gene, n=40)

## To find some example gene symbols you can try:
## 
grep("hox", gene$danio_rerio$gene, ignore.case=TRUE, value=TRUE)
grep("hox", gene$mus_musculus$gene, ignore.case=TRUE, value=TRUE)
grep("hox", gene$homo_sapiens$gene, ignore.case=TRUE, value=TRUE)
grep("znf", gene$homo_sapiens$gene, ignore.case=TRUE, value=TRUE)
grep("rnf", gene$homo_sapiens$gene, ignore.case=TRUE, value=TRUE)
grep("col", gene$homo_sapiens$gene, ignore.case=TRUE, value=TRUE)
grep("mcm", gene$homo_sapiens$gene, ignore.case=TRUE, value=TRUE)
grep("plx", gene$homo_sapiens$gene, ignore.case=TRUE, value=TRUE)
grep("sem", gene$homo_sapiens$gene, ignore.case=TRUE, value=TRUE)
grep("sox", gene$homo_sapiens$gene, ignore.case=TRUE, value=TRUE)
grep("tnf", gene$homo_sapiens$gene, ignore.case=TRUE, value=TRUE)
grep("pou", gene$homo_sapiens$gene, ignore.case=TRUE, value=TRUE)
grep("mad", gene$homo_sapiens$gene, ignore.case=TRUE, value=TRUE)
grep("tbx", gene$homo_sapiens$gene, ignore.case=TRUE, value=TRUE)

## Define a function that takes a gene symbol, a table of genes (genes),
## a table of transcripts and a table of exons as arguments:
draw.gene <- function(g.symbol, genes, transcripts, exons, sp=""){
    ## first identify a gene associated with the exact name
    g.i <- which( genes$gene == g.symbol )
    if(length(g.i) == 0){
        stop(paste("No gene with symbol", g.symbol))
    }
    if(length(g.i) > 1)
        warning(paste("Found", length(g.i), "genes with symbol: ", g.symbol, "\nWill use the first one"))
    ## then get the transcripts associated with the gene;
    ## note that we have ended up using the gene id here. Would have been 
    tr.i <- which( transcripts$g_id== genes[g.i[1], 'gene_id'] )
    ## now we will do this the slow simple way.
    ## first set up a plotting surface. We will plot transcript on a given y level:
    y <- 1:length(tr.i)
    plot.new()
    plot.window(xlim=range(genes[g.i[1], c('start', 'end')]), ylim=c(0, max(y)+1))
    axis(1)
    mtext(paste(sp, g.symbol), cex=2)
    ## draw a line for each transcript
    with(transcripts, segments( seq_region_start[tr.i], y, seq_region_end[tr.i], y) )
    ## lets loop through the different transcripts and find the exons associated with
    ## each transcript...
    for(i in 1:length(tr.i)){
        ex.i <- which(exons$transcript_id == transcripts$tr_id[ tr.i[i] ])
        with(exons, rect( seq_region_start[ex.i], y[i]-0.25, seq_region_end[ex.i], y[i]+0.25, col=rgb(0.2, 0.2, 0.8, 0.6)))
    }
    if( genes$seq_region_strand[g.i[1]] == 1 ){
        with(genes, arrows( start[g.i[1]], 0, end[g.i[1]], 0 ))
    }else{
        with(genes, arrows( end[g.i[1]], 0, start[g.i[1]], 0 ))
    }
}



dr.col <- sort(grep("col", gene$danio_rerio$gene, ignore.case=TRUE, value=TRUE))
hs.col <- sort(grep("col", gene$homo_sapiens$gene, ignore.case=TRUE, value=TRUE))
mm.col <- sort(grep("col", gene$homo_sapiens$gene, ignore.case=TRUE, value=TRUE))
## look at what these are:

## and then.. 
par(mfrow=c(2,2))
draw.gene( "col2a1a", gene$danio_rerio, transcript[[1]], exon[[1]], sp="D. rerio")
draw.gene( "Col2a1", gene$mus_musculus, transcript[[3]], exon[[3]], sp="M. musculus")
draw.gene( "COL2A1", gene$homo_sapiens, transcript[[2]], exon[[2]], sp="H. sapiens")

### Play around a bit with different genes and
### Consider how you can improve these visualisations. 

### You can save R data structures;
### these can then be reloaded in other R sessions if you
### should wish to do so.

saveRDS( db.vert, "db_vert.rds")
saveRDS( db.all, "db_all_2.rds")
saveRDS( dr.gene, "dr_gene.rds")
saveRDS( mm.gene, "mm_gene.rds")
saveRDS( hs.gene, "hs_gene.rds")
saveRDS( dr.chr, "dr_chr.rds")
saveRDS( mm.chr, "mm_chr.rds")
saveRDS( hs.chr, "hs_chr.rds")
saveRDS( gene, "gene.rds")
saveRDS( transcript, "transcript.rds")
saveRDS( exon, "exon.rds" )
