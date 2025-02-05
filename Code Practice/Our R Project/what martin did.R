> debug(find_CpG_mp)
> # Call all CpG islands and save in a list
  > found_CpG_mp <- find_CpG_mp(dna_seq, matches)
  debugging in: find_CpG_mp(dna_seq, matches)
  debug at CpG_project_functions.R#175: {
  dna_seq <- unlist(strsplit(as.character(dna_seq), split = ""))
  calculate_cg_density_fast <- function(cg_positions, window_size = 200, 
                                        step_size = 1) {
    cg_positions <- sort(cg_positions)
    results <- list()
    seq_length <- length(dna_seq)
    print(paste("Sequence length:", seq_length))
    if (seq_length < window_size) {
      print("Sequence is too short for the window size")
      stop("The sequence is shorter than the window size.")
    }
    for (i in start(cg_positions)) {
      print(paste("Processing window starting at position:", 
                  i))
      if ((i + window_size - 1) > seq_length) {
        print("Window exceeds sequence length, skipping")
        next
      }
      window_seq <- dna_seq[i:(i + window_size - 1)]
      gc_content <- sum(window_seq == "G" | window_seq == 
                          "C")/window_size
      observed_cpg <- sum(window_seq[1:(window_size - 1)] == 
                            "C" & window_seq[2:window_size] == "G")
      num_c <- sum(window_seq == "C")
      num_g <- sum(window_seq == "G")
      expected_cpg <- (num_c * num_g)/window_size
      cpg_ratio <- ifelse(expected_cpg > 0, observed_cpg/expected_cpg, 
                          NA)
      print(paste("Window:", i, "GC content:", gc_content, 
                  "CpG ratio:", cpg_ratio, "Observed CpG:", observed_cpg))
      if (!is.na(gc_content) && !is.na(cpg_ratio)) {
        if (gc_content >= gc_threshold && cpg_ratio >= 
            cpg_ratio_threshold) {
          results <- c(results, list(list(start = i, 
                                          end = i + window_size - 1, gc_content = gc_content, 
                                          cpg_ratio = cpg_ratio)))
          print(paste("CpG island found at window:", 
                      i))
        }
      }
    }
    return(results)
  }
  }
Browse[1]> n
debug at CpG_project_functions.R#178: dna_seq <- unlist(strsplit(as.character(dna_seq), split = ""))
Browse[1]> n
debug at CpG_project_functions.R#180: calculate_cg_density_fast <- function(cg_positions, window_size = 200, 
step_size = 1) {
  cg_positions <- sort(cg_positions)
  results <- list()
  seq_length <- length(dna_seq)
  print(paste("Sequence length:", seq_length))
  if (seq_length < window_size) {
    print("Sequence is too short for the window size")
    stop("The sequence is shorter than the window size.")
  }
  for (i in start(cg_positions)) {
    print(paste("Processing window starting at position:", 
                i))
    if ((i + window_size - 1) > seq_length) {
      print("Window exceeds sequence length, skipping")
      next
    }
    window_seq <- dna_seq[i:(i + window_size - 1)]
    gc_content <- sum(window_seq == "G" | window_seq == "C")/window_size
    observed_cpg <- sum(window_seq[1:(window_size - 1)] == 
                          "C" & window_seq[2:window_size] == "G")
    num_c <- sum(window_seq == "C")
    num_g <- sum(window_seq == "G")
    expected_cpg <- (num_c * num_g)/window_size
    cpg_ratio <- ifelse(expected_cpg > 0, observed_cpg/expected_cpg, 
                        NA)
    print(paste("Window:", i, "GC content:", gc_content, 
                "CpG ratio:", cpg_ratio, "Observed CpG:", observed_cpg))
    if (!is.na(gc_content) && !is.na(cpg_ratio)) {
      if (gc_content >= gc_threshold && cpg_ratio >= cpg_ratio_threshold) {
        results <- c(results, list(list(start = i, end = i + 
                                          window_size - 1, gc_content = gc_content, cpg_ratio = cpg_ratio)))
        print(paste("CpG island found at window:", i))
      }
    }
  }
  return(results)
}
Browse[1]> n
exiting from: find_CpG_mp(dna_seq, matches)
> class(matches)
[1] "XStringViews"
attr(,"package")
[1] "Biostrings"
> length(matches)
[1] 1229881
> matches[1]
Views on a 77842275-letter DNAString subject
subject: TGTTATTCGTGATGACTGGGCGTAGGGTACGTAGGCGGATAATC...GGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGG
views:
  start end width
[1]     8   9     2 [CG]
> tmp <- as.matrix(matches)
Warning message:
  In .local(x, ...) :
  as.matrix() on an XStringViews object 'x' has changed behavior: now the
views in 'x' must be of equal width and each view is converted into a row of
single characters. To achieve the old behavior, do 'as.matrix(ranges(x))'.
To supress this warning, do 'suppressWarnings(as.matrix(x))'.
This warning will be removed in BioC 2.12.
> tmp <- as.matrix(ranges(matches))
> head(tmp)
[,1] [,2]
[1,]    8    2
[2,]   21    2
[3,]   30    2
[4,]   36    2
[5,]  114    2
[6,]  163    2
> dim(tmp)
[1] 1229881       2
> nrow(tmp)/1e6
[1] 1.229881
> head(tmp)
[,1] [,2]
[1,]    8    2
[2,]   21    2
[3,]   30    2
[4,]   36    2
[5,]  114    2
[6,]  163    2
> hist(diff(tmp[,1]))
> range(tmp[,1])
[1]        8 77837386
> hist(log10(diff(tmp[,1])))
> mtch <- lapply( dna_seq, function(x){ matchPattern("CG", x)})
> tmp <- matchPatterns("CG", dna_seq[1], max.mismatch=0, with.indels=FALSE)
Error in matchPatterns("CG", dna_seq[1], max.mismatch = 0, with.indels = FALSE) : 
  could not find function "matchPatterns"
> tmp <- matchPattern("CG", dna_seq[1], max.mismatch=0, with.indels=FALSE)
> tmp <- matchPattern("CG", dna_seq[1], max.mismatch=0)
> mtch <- lapply( dna_seq, function(x){ matchPattern("CG", x, max.mismatch=0, with.indels=FALSE)})
> length(mtch)
Error: object 'mtch' not found
> length(dna_seq)
[1] 77842275
> class(dna_seq)
[1] "DNAString"
attr(,"package")
[1] "Biostrings"
> width(dna_seq)
Error: unable to find an inherited method for function ‘width’ for signature ‘x = "DNAString"’
> # Read the DNA content of the slime mold
  > fasta_data <- readDNAStringSet("assembly.fasta")
> plot(sort(width(dna_seq)))
Error in h(simpleError(msg, call)) : 
  error in evaluating the argument 'x' in selecting a method for function 'plot': error in evaluating the argument 'x' in selecting a method for function 'sort': unable to find an inherited method for function ‘width’ for signature ‘x = "DNAString"’
> plot(sort(width(fasta_data)))
> plot(cumsum(sort(width(fasta_data))))
> x11()
> plot(cumsum(sort(width(fasta_data))))
> plot(cumsum(sort(width(fasta_data))))
> s <- sort(width(fasta_data))
> plot(s, cumsum(s))