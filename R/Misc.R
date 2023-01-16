




# /*
# ----------------------------- split_taxonomy function ---------------------------
# */



#'
#' Split the taxonomic lineages in classical NCBI format into separate columns for each taxonomic level.
#'
#' `split_taxonomy()` takes in a 7 levels NCBI taxonomic lineage string and returns a 8
#' columns data.frame, 1 column for each taxonomic level (eg Kingdom, Phylum,
#' Class, Order, Family, Genus and species). The column "species" contains the
#' binomial specific name (i.e. genus + species) while and 8th column called
#' "species_orig" contains only the specific epithet (as in the original NCBI string).
#' If the argument clean = TRUE, the "species" column will be cleaned (see details)
#' while the column "species_orig" remains unchanged.
#'
#' @details
#' If the argument clean = TRUE, the function tries as best as possible to
#' "clean" the "species" column to remove any subspecific taxonomic information
#' (eg subspecies, variety, form) and to transform the uncertain specific ID (aff., cf. x)
#' by the generic "sp." (see examples). The column "species_orig" remains always
#' unchanged so that you can still retrieve the untransformed information.
#'
#' @param taxonomy A character vector with taxonomic lineages in NCBI format (see examples)
#' @param clean A logical. If TRUE (default) it will try to "clean" the species
#' level taxonomy e.g. by removing any subspecific information (see details).
#'
#' @return A data.frame with 8 columns for each taxonomic level between Phylum
#' and species + the binomial specific name (cleaned or not).
#'
#' @examples
#'
#' # example of NCBI taxonomic lineages
#'
#' ex <- c(
#' "k__Viridiplantae; p__Streptophyta; c__Magnoliopsida; o__Sapindales; f__Sapindaceae; g__Acer; s__tataricum subsp. aidzuense",
#' "k__Viridiplantae; p__Streptophyta; c__Magnoliopsida; o__Asterales; f__Asteraceae; g__Amphipappus; s__fremontii var. spinosus",
#' "k__Viridiplantae; p__Streptophyta; c__Magnoliopsida; o__Poales; f__Poaceae; g__Helictochloa; s__pratensis subsp. aff. pratensis GW-2014",
#' "k__Viridiplantae; p__Streptophyta; c__Magnoliopsida; o__Rosales; f__Moraceae; g__Artocarpus; s__nitidus cf. subsp. humilis EMG-2016",
#' "k__Viridiplantae; p__Streptophyta; c__Magnoliopsida; o__Rosales; f__Moraceae; g__Artocarpus; s__cf. nitidus",
#' "k__Viridiplantae; p__Streptophyta; c__Magnoliopsida; o__Rosales; f__Moraceae; g__Artocarpus; s__aff.nitidus",
#' "k__Viridiplantae; p__Streptophyta; c__Magnoliopsida; o__Rosales; f__Moraceae; g__Artocarpus; s__aff.",
#' "k__Viridiplantae; p__Streptophyta; c__Magnoliopsida; o__Rosales; f__Moraceae; g__Artocarpus; s__x",
#' "k__Viridiplantae; p__Streptophyta; c__Magnoliopsida; o__; f__Moraceae; g__x; s__",
#' "k__Viridiplantae; p__Streptophyta; c__Magnoliopsida; o__Malvales; f__Tiliaceae; g__Tilia; s__x europaea"
#' )
#'
#' (my_taxonomy <- split_taxonomy((ex)))
#'
#' # It might still be useful to display the species with more than 1 space
#' # character in their name after cleaning in order to detect hybrids + other
#' # potentially problematic names that could need further cleaning :
#'
#' unique(my_taxonomy[nchar(gsub("[^ ]*", "", my_taxonomy$Species))>1,
#'                    c("Species", "Species_orig")])
#'
#'
#'
#' # Read a true NCBI taxonomy file from a file on the disk
#' taxonomy_path <- system.file("extdata/ITS2_Rosales_Restricted.tsv",
#'                              package = "CVrefDB")
#' taxo <- read.table(taxonomy_path, header = TRUE, sep =  "\t")
#'
#' # Split and clean the taxonomy and append it to the original data
#' taxo <- cbind(taxo, split_taxonomy(taxo$Taxonomy))
#' head(taxo, n = 15)
#'
#'
#' @export
#'
split_taxonomy <- function(taxonomy, clean = TRUE){

    # NB : it is important to split first and gsub after to keep the right
    # number of columns even when some taxonomic levels are empty

    # To handle cases where taxonomy is read as a factor and not a character
    if(!is.character(taxonomy)){taxonomy <- as.character(taxonomy)}

    tmp <- do.call(rbind, strsplit(taxonomy, split = ";"))
    tmp <- gsub(" ?[a-z]{1,1}__", "", tmp)

    tmp <- as.data.frame(tmp)

    colnames(tmp) <- c("Kingdom", "Phylum", "Class", "Order",
                       "Family", "Genus", "Species")

    # keep the full unmodified name
    tmp$Species_orig <- tmp$Species
    tmp$Species <- paste(tmp$Genus, tmp$Species, sep = " ")

    if(clean == TRUE) {

        # cleaning of the species name

        # Remove everything that is after sp.
        tmp$Species <- gsub("(.* sp\\.) .*", "\\1", tmp$Species)

        # Remove the subspecies, variety and forms information --> keep the species level
        tmp$Species <- gsub("(.*) cf\\. var\\. .*", "\\1", tmp$Species)
        tmp$Species <- gsub("(.*) cf\\. subsp\\. .*", "\\1", tmp$Species)
        tmp$Species <- gsub("(.*) aff\\. var\\. .*", "\\1", tmp$Species)
        tmp$Species <- gsub("(.*) aff\\. subsp\\. .*", "\\1", tmp$Species)
        tmp$Species <- gsub("(.*) var\\. .*", "\\1", tmp$Species)
        tmp$Species <- gsub("(.*) subsp\\. .*", "\\1", tmp$Species)
        tmp$Species <- gsub("(.*) f\\. .*", "\\1", tmp$Species)

        # aff. and cf. are uncertain IDs --> converted to sp.
        # tmp$Species <- gsub("(.*) aff\\. .*", "\\1 sp.", tmp$Species)
        # tmp$Species <- gsub("(.*) cf\\. ?.*", "\\1 sp.", tmp$Species)
        tmp$Species <- gsub("([^ ]*) ?aff\\..*", "\\1 sp.", tmp$Species)
        tmp$Species <- gsub("([^ ]*) ?cf\\..*", "\\1 sp.", tmp$Species)
        # tmp$Species <- gsub("([^ ]*) ?x$", "\\1 sp.", tmp$Species)
        tmp$Species <- gsub("([^ ]*) +x$", "\\1 sp.", tmp$Species)




    }
    return(tmp)
}




# /*
# ----------------------------- blastn R function ---------------------------
# */


# arguments for testing

# db = "Fold01_DB.fasta"
# query = "Fold01_Query.fasta"
# num_threads = detectCores()
# max_target_seqs = 20
# evalue = 10^300
# outfmt = "6 qacc saccver evalue bitscore length pident"
# out = "Blast_output.tsv"
# verbose = 1
# retrieve_nonmatching = TRUE





#' Simple wraper function using system calls to `blastn` and `makeblastdb`
#'
#' `blastn` will use `system()` calls to a locally installed version of blastn
#' to blast query sequences against a reference database. If the reference
#' database is in fasta format, the function will call `makeblastdb` to create
#' a database in blast format. Unlike the standard blast results, this function
#' will add in the output the IDs of the sequences for which blast
#' didn't find any hit (i.e. hits with Expect value < `evalue`).
#'
#' @param db A reference database in fasta format
#' @param query Sequences to blast in fasta format or in blast format
#' @param out Name of the output file where the blast results will be saved
#' @param outfmt Output format string to be passed to the `blastn` function.
#'   Defaults to "6 qacc saccver evalue bitscore length pident".
#'   Type `blastn -help` in your system shell for more options.
#' @param retrieve_nonmatching If TRUE (default) : will keep sequence IDs with no blast match
#' @param num_threads Number of cores to be used (default = all cores)
#' @param max_target_seqs Maximum number of best blast hits to return (default = 20)
#' @param evalue Minimum E value to be considered by `blastn` (default = 10)
#' @param verbose 0 for no messages (default), 1 for time used, 2 for blast DB creation messages
#' @param clean If TRUE, removes the intermediate files created (i.e. blast db based on a fasta file). Default to FALSE
#'
#' @return
#'   A file saved on the disk containing for each query sequence the
#'   `max_target_seqs` top hits found by blast (or an empty line if blast didn't
#'   find any hit with E value >  `evalue`).
#'   With the default value for `outfmt` this file will be a 6 column file (without
#'   column names) containing in order : the Query ID, the reference database ID,
#'   E-value, bit score, length of the alignment and % of identity.
#'
#' @section Warning :
#'  You need to have `blastn` installed on your machine and available from the command line.
#'  You can test if this is the case by typing in the R console : `system("blastn -version")`
#'  which should print the version of `blast` available on your system.
#'
#' @importFrom parallel detectCores
#'
#' @export
#'
#' @examples
#'
#' # Create a temporary file for a toy reference database
#' temp_refDB <- tempfile("temp_refDB", fileext = ".fasta")
#'
#' # Add a few sequences to this temporary file
#' writeLines(
#'     ">AB020426.1
#' TGCAGAATCCCGTGAACCATCGAGTCTTTGAACGCAAGTTGCGCCCGAGGCCTTTTGGTTGAGGGCACGCCTGCCTGGG
#' CGTCACGCCTTGTTTCGCTCTGTGCCCATGCTCTTTCGGGGGCGGTCATGGATGCGGAGATTGGCCCTCCGTGCCTAGT
#' >AB020453.1
#' TGCAGAATCCCGTGAACCATCGAGTCTTTGAACGCAAGTTACGCCCGAGTCCTTTCGGTTGAGGGCACGCCTGCCTGGG
#' CGTCACGCCTTGTTTCGCTCTATGCCTATGCTCTTTCGGGGGCGGTCATGGATGCGGAGATTGGCTCTCTGTGCCTCGT
#' >AB020454.1
#' TGCAGAATCCCGTGAACCATCGAGTCTTTGAACGCAAGTTGCGCCCGAGGCCTTTCGGTTGAGGGCACGCCTGCCTGGG
#' CGTCACGCCTTGTTTCGCTCTGTGCCCATGCTCTTTCCGGGGCGGTCATGGATGCGGAGATTGGCCCTCCGTGCCTCGT
#' >AB020456.1
#' TGCAGAATCCCGTGAACCATCGAGTCTTTGAACGCAAGTTGCGCCCGAGGCCTTTCGGTTGAGGGCACGCCTGCCTGGG
#' CGTCACGCCTTGTTTCGCTCTGTGCCCATGCTCTTTCGGGGGCGGTCATGGATGCGGAGATTGGCCCTCCGTGCCTCGT",
#'     temp_refDB)
#'
#'
#' # Create a temporary file with a few target sequences
#' temp_query_sequences <- tempfile("temp_query_sequences", fileext = ".fasta")
#'
#' # Add a few sequences to this temporary file
#' writeLines(
#'     ">ID001
#' TGCAGAATCCCGTGAACCATCGAGTCTTTGAACGCAAGTTGCGCCCGAGGCCTTTTGGTTGAGGGCACGCCTGCCTGGGC
#' GTCACGCCTTGTTTCGCTCTGTGCCCATGCTCTTTCGGGGGCGGTCATGGATGCGGAGATTGGCCCTCCGTGCCTAGTGT
#' >ID002
#' ATCGACCTTTCGGTTGAGGGCACGCCTGCCTGGGCGTCACGCCTTGTTTCGCTCTATGCCTGTGCAGAATCCCGTGAACC
#' GGAGATTGGCTCTCTGTGCCTCGTGTGCGGCGGGCTTAAGCGCGGGCTGTCGGCGTCGGGATGGGCACGGTCATGGATGC
#' >ID003
#' AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
#' AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
#'     temp_query_sequences)
#'
#' # These fasta files are temporarily stored here :
#' # list.files(tempdir())
#'
#'
#' # Blast the query sequences against the toy reference database
#'
#' blastn(
#'     db = temp_refDB, # in typical use, this would be a fasta file on the disk
#'     query = temp_query_sequences, # in typical use, this would be a fasta file on the disk
#'     out = paste0(tempdir(), "/Blast_output.tsv"),
#'     retrieve_nonmatching = TRUE,
#'     num_threads = parallel::detectCores(),
#'     max_target_seqs = 20,
#'     evalue = 10,
#'     verbose = 1, # to print an estimate of the computing time
#'     clean = TRUE # remove the files created (blast db based on a fasta file)
#' )
#'
#'
#' # The results are saved on the disk in a tab separated values file.
#' # Note that blast did not find any match for the third query sequence (ID003)
#' # But we keep a line for this query sequence with NAs and an empty string for
#' # the ID of the reference database...
#'
#' blast_results <- read.table(paste0(tempdir(), "/Blast_output.tsv"), sep = "\t")
#' colnames(blast_results) <- c("TaxID_query", "TaxID_blast", "E_value",
#'                              "Bit_score", "Length", "Identity")
#' blast_results
#'
#'
#' # Delete the temporary files
#'
#' unlink(temp_refDB)
#' unlink(temp_query_sequences)


#'
#'
blastn <- function(db, # reference database in fasta format
                   query, # sequences to blast in fasta format
                   out = "Blast_output.tsv",
                   outfmt = "6 qacc saccver evalue bitscore length pident",
                   retrieve_nonmatching = TRUE, #
                   num_threads = parallel::detectCores(),
                   max_target_seqs = 20,
                   evalue = 10,
                   verbose = 0,  # 0 for no messages, 1 for time used, 2 for blast DB creation messages
                   clean = FALSE # remove the files created (blast db based on a fasta file)
){

    t0 <- Sys.time()

    blast_db <- db


    # if db is a fasta file : create a blast database

    if(grepl("\\.fasta$", db)){

        blast_db <- gsub("\\.fasta$", "", db)

        # build the blast DB based on a fasta file
        system(paste0('makeblastdb -in ', db,
                      ' -parse_seqids -blastdb_version 5 -dbtype nucl',
                      ' -title "', blast_db, '"',
                      ' -out "',blast_db,'"'),
               ignore.stdout = (verbose < 2) )
    }


    system(paste0(
        'blastn',
        ' -db ', blast_db,
        ' -query ', query,
        ' -num_threads ', num_threads,
        ' -max_target_seqs ', max_target_seqs,
        ' -evalue ', evalue,
        ' -outfmt "', outfmt, '"',
        ' -out ', out
    ))


    if(retrieve_nonmatching){

        # retrieve IDs from the query file
        # Use Biostrings if available.
        if(requireNamespace("Biostrings", quietly = TRUE)){
            ID <- data.frame(ID = names(Biostrings::readDNAStringSet(query)))
            # keep only the identifier part of the name (remove comments after withe space)
            ID$ID <- gsub("^([^ ]*).*", "\\1", ID$ID)
        } else {
            ID <- scan(query, what = "character")
            ID <- ID[grep("^>", ID)]
            ID <- data.frame(ID = gsub("^>([^ ]*).*", "\\1", ID))
        }

        # merge with blast outputs
        # use data.table if available for faster read/write/merge.
        # If not available use base R equivalents

        if(requireNamespace("data.table", quietly = TRUE)){

            blast_output <- data.table::fread(out, sep = "\t", header = F)
            blast_output <- merge(ID, blast_output, all = T,
                                  by.x = "ID", by.y = "V1", sort = FALSE)
            data.table::fwrite(blast_output, out, sep = "\t", row.names = F,
                               col.names = F)

        } else {

            blast_output <- read.csv(out, sep = "\t", header = F)
            blast_output <- merge(ID, blast_output, all = T,
                                  by.x = "ID", by.y = "V1", sort = FALSE)
            write.csv(blast_output, out, sep = "\t", row.names = F, header = F)

        }
    }

    if(clean & grepl("\\.fasta$", db)){
        file.remove(paste0(blast_db, c(".ndb", ".nhr", ".nin", ".nog", ".nos",
                                       ".not", ".nsq", ".ntf", ".nto")))
    }

    t1 <- Sys.time()
    time_count = round(difftime(t1, t0, units = "mins"),2)
    if(time_count == 0) {time_count = round(difftime(t1, t0, units = "secs"),2)}

    if(verbose >= 1){
        # cat("\nFold", fold, "\n")
        cat("Blast finished in", time_count,  attr(time_count, "units"))
    }


}



# /*
# ----------------------------- F-score function ---------------------------
# */


#' Function to compute Recall (R), Precision (P) and F-score (F) for each taxon
#'
#'
#' When we don't take all sequences into account, there are several ways to compute the
#' % of correctly identified sequences :
#' * P = Precision = % of taxon_predicted which is correct `P = TP/(TP+FP)`
#' * R = Recall = % of taxon_true for which the prediction is correct `R = TP/(TP+FN)`
#' * F = `F-score` = harmonic mean of F and P `F = 2*P*R/(P+R)`
#'
#' The problem is that there are in fact several ways to combine these values
#' to obtain a global estimates.
#'
#' For example we can compute F, R and P for each species and genus and then
#' compute their average for each family to obtain a global estimate of the
#' the prediction quality at at the genus and species level for each family.
#' The difficulty with this approach is that there are lots of NAs and 0
#' (for example is a genus is never assigned, its Precision will be NA) and it
#' is difficult to choose how to combine them. This typically occurs in species
#' poorly represented in the database and can have a huge influence on the results.
#'
#' Another possibility (not available with this function)
#' is to compute R and F across all sequences for
#' each family. Eg for each family we compute the % of sequences corresponding to
#' the genus of each family to be correctly predicted (Recall) or the % of sequences
#' that are correctly predicted to be a genus of each family (Precision).
#' We can then compute their harmonic mean to obtain the F score.
#' NB : with this approach all metrics (Recall, Precision and F score) are almost
#' identical and can therefore be interpreted as a % of correct predictions.
#' This is also a lot more complex to compute.
#'
#' The function[score_per_taxon()] provides yet another (more simple) approach.
#' Seee in particular the examples of this function's help for a comparion of the
#' outputs provided by F_score() and score_per_taxon().
#
#'
#' @param taxon_true A character vector with the true taxon
#' @param taxon_predicted A character vector with the corresponding assigned/predicted taxon
#' @param predicted_NA_wrong A logical. If TRUE : the unassigned taxa
#'   (ie with missing values) will be considered as not correctly identified.
#'   Imagine you have 3 sequences, 1 assignment is correct, 1 is incorrect and
#'   1 is NA (no assignment). If predicted_NA_wrong = FALSE, we compute the
#'   proportion of correct prediction as 1/2 (the NA is ignored).
#'   If predicted_NA_wrong = TRUE, we compute the
#'   proportion of correct prediction as 1/3 (the NA is counted in the total.
#'   This is not an incorrect assignment but the assignment is still unsuccessful).
#'
#' @return a data.frame with for each taxon the F-score (F), Recall (R) and Precision (P)
#' and the number of sequences used to compute R (Rnb) and P (Pnb)
#' @export
#'
#' @examples
#'
#' # retrieve path to two example datasets
#' fasta_path <- system.file("extdata/ITS2_Rosales_Restricted.fasta",
#'                           package = "CVrefDB")
#' taxonomy_path <- system.file("extdata/ITS2_Rosales_Restricted.tsv",
#'                              package = "CVrefDB")
#'
#' # Read the fasta and tsv files
#' fasta <- Biostrings::readDNAStringSet(fasta_path)
#' reftaxo <- read.table(taxonomy_path, sep = "\t", header = TRUE)
#'
#' # Cross validate and imput taxonomies with 4 methods
#' output_leakedCV <- CV_blastn(fasta_db = fasta, taxo = reftaxo,
#'                              k = 10, seed = 12, verbose = TRUE)
#' assigned_long <- assign_taxonomy(output_leakedCV, taxo = reftaxo, Order = NA)
#'
#' # Recall, Precision and F-score for TopHitPlus method and family level
#' tmp <- assigned_long[assigned_long$Method == "TopHitPlus" &
#'                          assigned_long$Tax_level == "Family",]
#' F_score(taxon_predicted = tmp$Taxon, taxon_true = tmp$Taxon_true)
#'
#' # Compute these statistics for the 4 methods and all taxonomic levels
#' # Note that in many cases we have NAs or evn 0 values that are poorly exploitable
#' # due eg to the limited number of sequences
#' library(tidyr)
#' library(dplyr)
#' metrics <-
#'     assigned_long %>%
#'     mutate(Taxon = replace(Taxon, Taxon == "", NA)) %>%
#'     group_by(Method, Tax_level) %>%
#'     nest() %>%
#'     mutate(F_score = lapply(data, function(x) {
#'         F_score(taxon_true = x$Taxon_true,
#'                 taxon_predicted = x$Taxon,
#'                 predicted_NA_wrong = TRUE)
#'     })) %>%
#'     select(-data) %>%
#'     unnest(F_score) %>%
#'     ungroup()
#'
#' metrics %>%
#'     slice(1:15)


F_score <- function(taxon_true, taxon_predicted,
                    # Should NA in the predicted taxon be considered as wrong predictions ?
                    predicted_NA_wrong = FALSE){

    df = data.frame(taxon_true, taxon_predicted) %>%
        dplyr::mutate(CorrectID = (taxon_true == taxon_predicted))

    if(predicted_NA_wrong){
        df[is.na(df$CorrectID), "CorrectID"] <- FALSE
    }

    # accuracy
    # A <- mean(df$CorrectID, na.rm = T)

    # Recall = % of correct taxon_true
    R <-
        df %>%
        dplyr::group_by(taxon_true) %>%
        dplyr::summarize( R = round(mean(CorrectID, na.rm = T)*100, 2),
                          Rnb = dplyr::n()) %>%
        dplyr::rename(Taxon = taxon_true)

    # Precision = % of correct taxon_predicted
    P <-
        df %>%
        dplyr::group_by(taxon_predicted) %>%
        dplyr::summarize( P = round(mean(CorrectID, na.rm = T)*100, 2),
                          Pnb = dplyr::n()) %>%
        dplyr::rename(Taxon = taxon_predicted)

    # Merge + F-score =  Harmonic mean of R and P
    res <- dplyr::full_join(R, P, by = "Taxon") %>%
        dplyr::mutate(F = round(2*(P*R)/(P+R), 2)) %>%
        dplyr::filter(!is.na(Taxon)) %>%
        dplyr::arrange(Taxon) %>%
        dplyr::select(Taxon, R, P, F, Rnb, Pnb)

    return(res)
}



# /*
# ----------------------------- score_per_taxon function ---------------------------
# */




#' Calculate a score (% of correctly identified sequences) per taxon
#'
#' This function calculates a score for each
#' taxon at a specified taxonomic level (e.g. family), based on the
#' percentage of sequences that were correctly identified. For example,
#' if the grouping taxonomic level is "Family", the function will compute
#' for the sequences of each true family the % of sequences correctly identified
#' at the Family, Genus, and Species level.
#'
#' NB : there are in fact many different ways to compute similar statistics which
#' will provide results with more or less different results and interpretation.
#' We chose here a rather simple and straightforward way to proceed. The code
#' of the function is pretty simple and would be easy to adapt for other use cases.
#' Computing these values by yourself is probably the best way to be certain
#' to understand exactly the meaning of the computations...
#'
#' The question we are asking here with the default values is :
#' if I know that a group of sequences are
#' in a given Family (true family), which proportion of these sequences are correctly
#' identified at the Family, Genus or Species level for each family ?
#' See the examples for a comparison of the scores computed here and the usual
#' Recall and Precision scores which are sometimes identical to the output tof
#' current function and sometimes not...
#'
#' @param assigned a data.frame, output form function assign_taxonomy()
#' @param grouping_tax_level a character string indicating the taxonomic level to use for grouping
#'   the results. Possible values are :"Kingdom", "Phylum", Class", "Order,
#'   "Family", "Genus", or "Species". Default is "Family". NB : the chosen value
#'   should be present in the `Tax_level` column of the `assigned` table.
#'   NB : the case is important)
#' @param grouping_taxon a character string with only one of two possible values :
#'   `Taxon_true` (default) : we will use as grouping taxon the true taxon corresponding to
#'   each sequence.
#'   `Taxon` : we will use the predicted/assigned taxon as grouping taxon
#'   (e.g. all sequences predicted to be in a given family will be grouped together).
#'   The results of both options will usually be very similar unless you chose
#'   a `grouping_tax_level` at which the error rate is high.
#' @param predicted_NA_wrong a logical value indicating whether sequences that were predicted to
#'   be "NA" should be considered as incorrect. Default is FALSE. Empty strings are considered as NAs.
#'   Imagine you have 3 sequences, 1 assignment is correct, 1 is incorrect and
#'   1 is NA (no assignment). If predicted_NA_wrong = FALSE, we compute the
#'   proportion of correct prediction as 1/2 (the NA is ignored).
#'   If predicted_NA_wrong = TRUE, we compute the
#'   proportion of correct prediction as 1/3 (the NA is counted in the total.
#'   This is not an incorrect assignment but the assignment is still unsuccessful)
#' @return a data.frame with 6 columns  :
#'   * `Method` : The assignment method always used as a grouping factor
#'   * `Tax_level` : Taxonomic level (eg if Tax_level = "Genus", Pct represnet the %
#'   of sequences correctly predicted at genus level).
#'   * `Level` : Taxonomic level code (lower case first letter of `Tax_level`)
#'   * `Grouping_taxon` : grouping taxon
#'   * `Pct` : % of of sequences correctly identified for a given grouping taxon at a
#'   given taxonomic level
#'   * `N` : Number of sequences used to compute `Pct`.
#'   Excluding NAs if `predicted_NA_wrong` = FALSE
#'   and including NAs in the total if `predicted_NA_wrong` = TRUE
#'
#' @export
#' @examples
#'
#' # retrieve path to two example datasets
#' fasta_path <- system.file("extdata/ITS2_Rosales_Restricted.fasta",
#'                           package = "CVrefDB")
#' taxonomy_path <- system.file("extdata/ITS2_Rosales_Restricted.tsv",
#'                              package = "CVrefDB")
#'
#' # Read the fasta and tsv files
#' fasta <- Biostrings::readDNAStringSet(fasta_path)
#' reftaxo <- read.table(taxonomy_path, sep = "\t", header = TRUE)
#'
#' # Cross validate and input taxonomies with 4 methods
#' output_leakedCV <- CV_blastn(fasta_db = fasta, taxo = reftaxo,
#'                              k = 10, seed = 12, verbose = TRUE)
#' assigned_long <- assign_taxonomy(output_leakedCV, taxo = reftaxo,
#'                                  Order = NA, Assignation_method = "TopHitPlus")
#'
#' # Compute the scores per family.
#'
#' # Interpretation : For the sequences that are known to be in the family Rosaceae (column
#' # (Grouping_taxon), the % of sequences correctly predicted (Pct) is 100%, 93.9%
#' # and 7.61% at Family, genus and species level respectively (Tax_level).
#' # The % are based on 460 sequences (column "N")
#'
#' score_per_taxon(assigned_long, grouping_tax_level = "Family")
#'
#' # compute the scores per genus and look at the % of sequences correctly predicted
#' # at the genus level (for genus with more than 5 sequences available)
#'
#' per_genus <- score_per_taxon(assigned_long, grouping_tax_level = "Genus")
#'
#' library(dplyr)
#' per_genus %>%
#'     filter(Tax_level == "Genus" & N > 5) %>%
#'     arrange(desc(Pct))
#'
#' # NB : in this particular case (at the genus level), Pct corresponds to the
#' # statistic usually called "Recall" (R) which can be computed with another function
#' # of this package : F_score()
#'
#' tmp <- assigned_long[assigned_long$Tax_level == "Genus", ]
#' F_score(taxon_true = tmp$Taxon_true,
#'         taxon_predicted = tmp$Taxon) %>%
#'     filter(Rnb>5) %>%
#'     arrange(desc(R))
#'
#' # You could compute the equivalent of Precision (P) with the argument grouping_taxon
#' # set to "Taxon" instead of of the default "Taxon_true"
#'
#' per_genus_precision <-
#'     score_per_taxon(assigned_long,
#'                     grouping_tax_level = "Genus",
#'                     grouping_taxon = "Taxon")
#'
#' # show a subset
#' per_genus_precision %>%
#'     filter(Tax_level == "Genus" & N > 5) %>%
#'     arrange(desc(Pct))
#'
#' # Compare with column P of F_score()
#' F_score(taxon_true = tmp$Taxon_true,
#'         taxon_predicted = tmp$Taxon) %>%
#'     filter(Pnb>5) %>%
#'     arrange(desc(P))
#'
#' # The whole point of the score_per_taxon() function is that for each of these genus
#' # it can also provide the % of correct prediction at the species level and this score
#' # may differ from what you could compute with F_score() at the species level...
#' # score_per_taxon() is hierarchical while F_score() is not. F_score() will
#' # compute the scores for each species independently of the genus in which this species is.
#'
#' # With score_per_taxon() we can for example compute for each genus the % of
#' # sequences correctly identified at the species level :
#' per_genus %>%
#'     filter(Tax_level == "Species" & N > 5) %>%
#'     arrange(desc(Pct))
#'
#' # With F_score() you may compute F, R and P for each species and then try to
#' # aggregate these values at the genus level but it is not straightforward and
#' # it will often cause problems because of the presence of many NA, NaN or 0
#' # and it is difficult to decide how to handle these values properly...
#' tmp <- assigned_long[assigned_long$Tax_level == "Species", ]
#' F_score(taxon_true = tmp$Taxon_true,
#'         taxon_predicted = tmp$Taxon) %>%
#'     filter(Rnb>5) %>%
#'     arrange(desc(R))



score_per_taxon <- function(
        assigned , # a data.frame, output form function assign_taxonomy()
        grouping_tax_level = "Family",
        grouping_taxon = "Taxon_true", # either "Taxon_true" or "Taxon"
        predicted_NA_wrong = FALSE
){
    if(predicted_NA_wrong) {
        assigned <-
            assigned %>%
            dplyr::mutate(Taxon = replace(Taxon, Taxon == "", NA)) %>%
            dplyr::mutate(CorrectID = Taxon_true == Taxon & !is.na(Taxon)) # NA are considered as incorrect
    } else{
        assigned <-
            assigned %>%
            dplyr::mutate(Taxon = replace(Taxon, Taxon == "", NA)) %>%
            dplyr::mutate(CorrectID = Taxon_true == Taxon)} # NAs remain as NA and will be ignored

    # retrieve the grouping taxonomy corresponding top each ID
    ID_taxonomy <-
        assigned %>%
        dplyr::filter(Tax_level == grouping_tax_level) %>%
        dplyr::select(ID, Grouping_taxon = !!grouping_taxon) %>%
        dplyr::distinct()

    assigned %>%
        dplyr::left_join(ID_taxonomy, by = "ID") %>%
        dplyr::group_by(Method, Tax_level, Level, Grouping_taxon) %>%
        dplyr::summarize(Pct = mean(CorrectID, na.rm = TRUE)*100, # % of correctly identified sequences
                         N = sum(!is.na(CorrectID))) # Number of sequences used for the computation of
}





# /*
# ----------------------------- declare unquoted variable names  ---------------------------
# */



# declare unquoted variable names to avoid Notes on R CMD check

utils::globalVariables(c("fasta_db", "taxo", "Label", "Fold", "TaxID_Query",
                         "Kingdom", "Phylum", "Class", "Order", "Family", "Genus",
                         "Species", "Species_orig", "TaxID_Blast", "E_Value",
                         "Bit_Score", "Length", "Identity", "CorrectID", "Taxon",
                         "Rnb", "Pnb", "Nb_tot_hits", "Tax_level", "Level",
                         "Cons", "Method", "Feature.ID", "Tax_level", "Method",
                         "ID","Grouping_taxon", "Taxon_true"))



