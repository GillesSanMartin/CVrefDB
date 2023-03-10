% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Misc.R
\name{blastn}
\alias{blastn}
\title{Simple wraper function using system calls to \code{blastn} and \code{makeblastdb}}
\usage{
blastn(
  db,
  query,
  out = "Blast_output.tsv",
  outfmt = "6 qacc saccver evalue bitscore length pident",
  retrieve_nonmatching = TRUE,
  num_threads = parallel::detectCores(),
  max_target_seqs = 20,
  evalue = 10,
  verbose = 0,
  clean = FALSE
)
}
\arguments{
\item{db}{A reference database in fasta format}

\item{query}{Sequences to blast in fasta format or in blast format}

\item{out}{Name of the output file where the blast results will be saved}

\item{outfmt}{Output format string to be passed to the \code{blastn} function.
Defaults to "6 qacc saccver evalue bitscore length pident".
Type \code{blastn -help} in your system shell for more options.}

\item{retrieve_nonmatching}{If TRUE (default) : will keep sequence IDs with no blast match}

\item{num_threads}{Number of cores to be used (default = all cores)}

\item{max_target_seqs}{Maximum number of best blast hits to return (default = 20)}

\item{evalue}{Minimum E value to be considered by \code{blastn} (default = 10)}

\item{verbose}{0 for no messages (default), 1 for time used, 2 for blast DB creation messages}

\item{clean}{If TRUE, removes the intermediate files created (i.e. blast db based on a fasta file). Default to FALSE}
}
\value{
A file saved on the disk containing for each query sequence the
\code{max_target_seqs} top hits found by blast (or an empty line if blast didn't
find any hit with E value >  \code{evalue}).
With the default value for \code{outfmt} this file will be a 6 column file (without
column names) containing in order : the Query ID, the reference database ID,
E-value, bit score, length of the alignment and \% of identity.
}
\description{
\code{blastn} will use \code{system()} calls to a locally installed version of blastn
to blast query sequences against a reference database. If the reference
database is in fasta format, the function will call \code{makeblastdb} to create
a database in blast format. Unlike the standard blast results, this function
will add in the output the IDs of the sequences for which blast
didn't find any hit (i.e. hits with Expect value < \code{evalue}).
}
\section{Warning }{

You need to have \code{blastn} installed on your machine and available from the command line.
You can test if this is the case by typing in the R console : \code{system("blastn -version")}
which should print the version of \code{blast} available on your system.
}

\examples{

# Create a temporary file for a toy reference database
temp_refDB <- tempfile("temp_refDB", fileext = ".fasta")

# Add a few sequences to this temporary file
writeLines(
    ">AB020426.1
TGCAGAATCCCGTGAACCATCGAGTCTTTGAACGCAAGTTGCGCCCGAGGCCTTTTGGTTGAGGGCACGCCTGCCTGGG
CGTCACGCCTTGTTTCGCTCTGTGCCCATGCTCTTTCGGGGGCGGTCATGGATGCGGAGATTGGCCCTCCGTGCCTAGT
>AB020453.1
TGCAGAATCCCGTGAACCATCGAGTCTTTGAACGCAAGTTACGCCCGAGTCCTTTCGGTTGAGGGCACGCCTGCCTGGG
CGTCACGCCTTGTTTCGCTCTATGCCTATGCTCTTTCGGGGGCGGTCATGGATGCGGAGATTGGCTCTCTGTGCCTCGT
>AB020454.1
TGCAGAATCCCGTGAACCATCGAGTCTTTGAACGCAAGTTGCGCCCGAGGCCTTTCGGTTGAGGGCACGCCTGCCTGGG
CGTCACGCCTTGTTTCGCTCTGTGCCCATGCTCTTTCCGGGGCGGTCATGGATGCGGAGATTGGCCCTCCGTGCCTCGT
>AB020456.1
TGCAGAATCCCGTGAACCATCGAGTCTTTGAACGCAAGTTGCGCCCGAGGCCTTTCGGTTGAGGGCACGCCTGCCTGGG
CGTCACGCCTTGTTTCGCTCTGTGCCCATGCTCTTTCGGGGGCGGTCATGGATGCGGAGATTGGCCCTCCGTGCCTCGT",
    temp_refDB)


# Create a temporary file with a few target sequences
temp_query_sequences <- tempfile("temp_query_sequences", fileext = ".fasta")

# Add a few sequences to this temporary file
writeLines(
    ">ID001
TGCAGAATCCCGTGAACCATCGAGTCTTTGAACGCAAGTTGCGCCCGAGGCCTTTTGGTTGAGGGCACGCCTGCCTGGGC
GTCACGCCTTGTTTCGCTCTGTGCCCATGCTCTTTCGGGGGCGGTCATGGATGCGGAGATTGGCCCTCCGTGCCTAGTGT
>ID002
ATCGACCTTTCGGTTGAGGGCACGCCTGCCTGGGCGTCACGCCTTGTTTCGCTCTATGCCTGTGCAGAATCCCGTGAACC
GGAGATTGGCTCTCTGTGCCTCGTGTGCGGCGGGCTTAAGCGCGGGCTGTCGGCGTCGGGATGGGCACGGTCATGGATGC
>ID003
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
    temp_query_sequences)

# These fasta files are temporarily stored here :
# list.files(tempdir())


# Blast the query sequences against the toy reference database

blastn(
    db = temp_refDB, # in typical use, this would be a fasta file on the disk
    query = temp_query_sequences, # in typical use, this would be a fasta file on the disk
    out = paste0(tempdir(), "/Blast_output.tsv"),
    retrieve_nonmatching = TRUE,
    num_threads = parallel::detectCores(),
    max_target_seqs = 20,
    evalue = 10,
    verbose = 1, # to print an estimate of the computing time
    clean = TRUE # remove the files created (blast db based on a fasta file)
)


# The results are saved on the disk in a tab separated values file.
# Note that blast did not find any match for the third query sequence (ID003)
# But we keep a line for this query sequence with NAs and an empty string for
# the ID of the reference database...

blast_results <- read.table(paste0(tempdir(), "/Blast_output.tsv"), sep = "\t")
colnames(blast_results) <- c("TaxID_query", "TaxID_blast", "E_value",
                             "Bit_score", "Length", "Identity")
blast_results


# Delete the temporary files

unlink(temp_refDB)
unlink(temp_query_sequences)


}
