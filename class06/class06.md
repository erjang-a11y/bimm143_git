# Class 6: R functions
Erin (PID: A17713992)

- [Background](#background)
- [A first function](#a-first-function)
- [Q2. A generate_dna() function](#q2-a-generate_dna-function)
- [A generate_protein() function](#a-generate_protein-function)
- [Are our peptides “unique in
  nature”?](#are-our-peptides-unique-in-nature)
- [Connecting our findings to
  immunology](#connecting-our-findings-to-immunology)

## Background

All functions in R have at least 3 things:

- a *name* (we pick that and use it to call the function)
- input *arguments* (one or more comma separated inputs that go inside
  the brackets when we call the function)
- the *body* (the lines of R code that do the work of the function)

## A first function

Here we will create a function to add some numbers. Let’s call it
`add()`

Input arguments can be either **“required”** or **“optional”**. The
latter have fall-back default values that will be used if the user does
not specify them.

> Q1a. Your first version of the function should add two input numbers
> together. For example, add(4, 7) should return 11.

``` r
add <- function(x, y=0) {
  x + y
}
```

Can we use our new function:

``` r
add(10, 100)
```

    [1] 110

``` r
add(10)
```

    [1] 10

> Q1b. For you second version, adapt your first function so it can take
> a single input vector or two inputs as before. For example, add(4, 7)
> and add( c(4,7,10) ).

``` r
add <- function(x, y=0) {
  sum(x, y)
}
```

``` r
add(c(4, 7, 10))
```

    [1] 21

> Q1c. Finally, on your own (outside of class) create a third version of
> your function that can add any number of inputs that the user
> provides. For example, add(1, 2, 3, -4) should return 2. Hint: explore
> the … (dots) argument or the base R function sum().

``` r
add <- function(x, y, ...) {
  sum(x, y, ...)
}
```

``` r
add(1, 2, 3, 4)
```

    [1] 10

We can explicitly set a **return** value output from a funcion (rather
than the default last line) by using the `return()` function call.

``` r
add <- function(x, y=0, z=0) {
  sum(x, y, z)
  cat("Is it break time yet?\n")
}

add(10, 100)
```

    Is it break time yet?

``` r
add <- function(x, y=0, z=0) {
  return(sum(x, y, z))
  cat("Is it break time yet?\n")
}

add(10, 100)
```

    [1] 110

## Q2. A generate_dna() function

A useful function here is the “base R” `sample()` function:

``` r
sample(1:5, size=60, replace=TRUE)
```

     [1] 3 4 4 4 5 1 1 4 5 2 5 4 2 1 2 4 1 5 2 2 2 4 5 4 1 3 1 1 1 3 4 5 2 2 1 3 4 4
    [39] 2 5 1 2 5 1 3 5 5 1 5 1 4 5 1 5 5 3 3 3 2 3

We can use this to make a random nucleotide sequence if we draw from
“A”, “C”, “G”, and “T”…

``` r
sample(x=c("A", "C", "G", "T"), size = 10, replace = TRUE)
```

     [1] "C" "G" "A" "G" "C" "C" "A" "T" "G" "A"

> **Q2a**. Write a function generate_dna() that returns a random DNA
> sequence of a length specified by the user. Your first version should
> return a multi-element vector of single character nucleotides. For
> example generate_dna(6) might return “A”, “T”, “T”, “G”, “A”, “C”.

``` r
generate_dna <- function(len) {
  sample(x=c("A", "C", "G", "T"), size = len, replace = TRUE)
}
```

``` r
generate_dna(len = 6)
```

    [1] "T" "G" "C" "A" "G" "G"

> **Q2b**. Your second version should *optionally* be able to return
> either a multi-element vector of single character nucleotides (as
> before) or a **single character string** (not a vector of single
> letters but a singe vector of multiple letters). For example
> “AAGGCTG”.

``` r
generate_dna <- function(len, single.element = TRUE) {
  ans <- sample(x=c("A", "C", "G", "T"), size = len, replace = TRUE)
  if(single.element) {
    ans <- paste(ans, collapse = "")
  }
  return(ans)
}
```

``` r
generate_dna(len = 6)
```

    [1] "CACGTA"

Functions that could be useful here are `paste()`, `if()`, `cat()`, and
`return()`.

``` r
paste( c("A", "C", "G"), collapse = "---")
```

    [1] "A---C---G"

> **Q2c**. Finally, create a final version of your function that prints
> out a FASTA format sequence with an id line indicating the sequence
> length.

    >len9
    CGAAGGCTG

``` r
generate_dna <- function(len, single.element = TRUE) {
  ans <- sample(x=c("A", "C", "G", "T"), size = len, replace = TRUE)
  if(single.element) {
    ans <- paste(ans, collapse = "")
  }
  
  ## Format as FASTA with an ID line
  cat( paste(">len", len, "\n", sep = "") )
  cat(ans)
  cat("\n")
  ##
  
  return(ans)
}
```

``` r
x <- generate_dna(4)
```

    >len4
    TCAT

## A generate_protein() function

> **Q3**. Write a function `generate_protein()` that returns a random
> peptide/protein sequence of a length specified by the user. For
> example `generate_protein(6)` might return `"WQRTAG"`.

``` r
generate_protein <- function(len) {
  ## choose amino acids at random according to specified length
  ans <- sample(x=c("A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S",
"T", "W", "Y", "V"), size = len, replace = TRUE)
  
  ## concatenate the sequence
  ans <- paste(ans, collapse = "")
  return(ans)
}
```

``` r
generate_protein(len = 6)
```

    [1] "MQMVAQ"

> **Q4**. Generate random protein sequences of length 6 to 13. Adapt and
> use your `generate_protein()` function to generate a series of random
> protein sequences ranging from 6 to 13 amino acids in length (one
> sequence per length). Take advantage of the base R function `for()` or
> `sapply()` so that you do not have to call `generate_protein()` eight
> times by hand.

``` r
for(l in 6:13) {
  cat(">", l, "\n", sep = "")
  cat( generate_protein(l), "\n" )
}
```

    >6
    VDWKDC 
    >7
    VEPYTHL 
    >8
    GGCNHEMA 
    >9
    VFWFNPCAK 
    >10
    WWVHRIMQSL 
    >11
    LFETRHYEISI 
    >12
    QGDYWPHHKATR 
    >13
    IATFFCVMWLGTE 

## Are our peptides “unique in nature”?

> **Q5**. Take your FASTA-formatted peptides from Q4 and run them as a
> single BLASTp search against the Non-redundant protein sequences (nr)
> database at https://blast.ncbi.nlm.nih.gov/. For this question do not
> restrict the organism (leave the Organism field blank so that the full
> nr database is searched).

| Length (aa) | Best hit % identity | Best hit % coverage | Unique |
|-------------|---------------------|---------------------|--------|
| 6           | 100                 | 100                 | N      |
| 7           | 100                 | 100                 | N      |
| 8           | 88                  | 100                 | Y      |
| 9           | 89                  | 100                 | Y      |
| 10          | 90                  | 100                 | Y      |
| 11          | 100                 | 50                  | Y      |
| 12          | 92                  | 75                  | Y      |
| 13          | 77                  | 90                  | Y      |

> **Q5a**. At which sequence length do your randomly generated peptides
> start to look “unique in nature” (i.e. no 100% coverage + 100%
> identity hit)?

At sequence length 8, my randomly generated peptides start to look
“unique in nature”.

> **Q5b**. Speculate why very short random peptides are almost always
> found in nr, while longer ones typically are not. Your answer should
> refer both to the size of the sequence space (20𝐿 for a peptide of
> length 𝐿) and to the size of the known protein universe.

The total possible combinations of amino acids for shorter sequences are
much smaller than for longer sequences, as supported by the exponential
growth formula 20^L. As a result, there is statistically a higher chance
of short peptides appearing in nature by chance, whereas longer peptides
would be much harder to find in nature. In addition, the size of the
known protein universe is approximately 10<sup>8-10</sup>9, which
displays the great number of possible peptide sequences that can exist.

## Connecting our findings to immunology

> **Q6**. In 3–6 sentences total and using your Q5 data and the
> reasoning from Q5b, what do you think this minimum length is and why
> might it be a bad design choice for the immune system to present very
> short peptides?

Based on the reasoning from the previous question, the minimum length
most likely is around at least 10 amino acids. Presenting shorter
peptides would be a bad design choice for the immune system because
these lengths would not be unique in nature. This would ultimately
result in a high binding affinity for self immune cells, causing the
body to attack its own cells. Having longer amino acid lengths ensures
stability and the ability to bind to unique cells.
