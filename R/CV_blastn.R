




# /*
# ----------------------------- function to perform the cross validation ---------------------------
# */




#' Cross-validation of taxonomic assignations by blast on a DNA reference database
#'
#' This function will repeatedly extract a random sample of sequences from a DNA barcode reference
#' database in fasta format, then blast these query sequences against the remainder
#' of the database (k-fold cross-validation) or the whole database (leaked cross-validation)
#' and return the output of the top blast hits with the true and predicted taxonomies and
#' blastn statistics (bit score, identity,...) in a data.frame or in a csv file on the disk.
#'
#' @param fasta_db Path to a reference database in fasta format or a DNAStringSet object (Biostring package)
#' @param taxo data.frame or path to a 2 columns tsv file (tab separated values) containing the
#'   taxonomic descriptions of the `fasta_db` : The first column must contain the
#'   Taxon ID corresponding to the fasta file name. The second column must contain
#'   the NCBI formatted taxonomies (see [split_taxonomy()] for examples).
#' @param output Name and path of the final output file to save on the disk.
#'   If NULL (default when fasta_db is a DNAStringSet), the function will return
#'   the results in a data.frame.
#'   Otherwise, the function returns nothing but saves the output on the disk. If fasta_db
#'   is a path, the default behavior is to save the output in the same directory
#'   in a file named "CV_blastn_output.csv".
#' @param label A string. The first column of the output file will contain label
#'   repeated on every row. Defaults to the output file name when it is not NULL
#'   or to an empty string otherwise.
#' @param folds_dir Where should the `__Folds__` directory be created  ? i.e.
#'   the temporary directory containing the query and reference DB for each fold + intermediary files.
#'   By default it will be the same directory as fasta_db when it is a path or
#'   the working directory when fasta_db is a DNAstringset. NB : any existing
#'   `__Folds__` direcory will be erased...
#' @param alternative_fasta_db If not NA or NULL (the default), `fasta_db` is
#'   blasted against this alternative db instead of
#'   blasted against itself. When k>1 the sequences from `fasta_db` in each fold
#'   are removed from `alternative_fasta_db` based on their IDs
#' @param alternative_taxo Taxonomic file corresponding to `alternative_fasta_db`
#' @param k Number of folds to create (default : k = 10). If `k==1` the function
#'   will use a leaked validation approach, i.e. the full database is blasted
#'   against itself without removing any sequence.
#'   The number of folds should not be >99 to avoid problems.
#' @param folds_subset Folds id to effectively blast against the database.
#'   e.g. :  if `folds_subset` = 3:5 , we will blast only the folds 3, 4 and 5.
#'   Useful mainly to test the code on a subset of the folds.
#'   Defaults to 1:k, i.e. all the folds.
#' @param sequences_subset A vector of integers : we will blast only these
#'   sequences for each fold. E.g. `sequences_subset` = 1:10 -->
#'   we will blast only the first 10 sequences for each fold.
#'   Useful mainly to test the code on a subset of the query sequences.
#'   Ignored for any non numeric value.
#' @param query_filter A vector of characters containing sequences IDs to filter out
#'    from the query sequences in each fold. Defaults to NULL (all sequences are used)
#' @param db_filter A vector of characters containing sequences IDs to filter out
#'    from the reference database used for each fold. Defaults to NULL
#'    (all sequences are used).
#' @param nb_cores An integer providing the number of cores to be used by `blastn`
#'    (defaults to all cores).
#' @param max_target_seqs An integer. Maximum number of hits to be returned by
#'    blast. Default = 15
#' @param verbose If TRUE (default) : prints messages during the process
#'    including computing time.
#' @param clean_species If TRUE (default), the species names will be cleaned as
#'    much as possible (remove subspecies, transform cf. into sp. ...)
#' @param folds_remove Logical. If TRUE (default), the directory `__Folds__`
#'    created during the cross-validation process will be removed.
#' @param seed An integer to chose if you want a perfectly reproducible result
#' (i.e. in the way the fold are created). Defaults to the system time, so
#' the results will change at each run.
#'
#' @return A csv file saved on the disk (with `;` separated columns and `,`
#' as decimal separators) or a data.frame containing the results of the
#' cross-validation process.
#' For each query sequence (`TaxID_query`) the table contains the known ("true") taxonomy
#' (i.e. columns `Kindgom_true`, `Phylum_true`,..., `Genus_true`, `Species_true`)
#' and the corresponding top 15 blast hits (by default) with their ID (`TaxID_blast`)
#' and corresponding taxonomic information
#' (i.e. columns `Kindgom`, `Phylum`,..., `Genus`, `Species`). For each hit
#' we also have 4 descriptive statistics provided by `blastn` : `E_value`,
#' `Bit_score`, `Length` of the alignment and `Identity` (% of identical bases).
#' For each line we also have an ID for the `Fold` (e.g. F01, F02,.... F10
#' for 10-fold Cross-Validation) and a `Label` which is identical for all
#' lines (by default, it is the output file name or an empty string).
#' The purpose of the label is for example to identify the strategy used when
#' CV_blast is run multiple times with different options/strategies and
#' the final results are merged (see examples).
#'
#' @importFrom parallel detectCores
#'
#' @export
#'
#' @examples
#'
#' # Example 1
#'
#' # Retrieve an example reference database (fasta file) and its corresponding
#' # taxonomic information from the package data examples.
#' # This example contains >7000 sequences for the ITS2 barcode of plants from the
#' # order "Rosales" restricted to a portion of the gene amplified by a given
#' # given primer used in metabarcoding.
#' #
#' fasta_path <- system.file("extdata/ITS2_Rosales_Restricted.fasta",
#'                           package = "CVrefDB")
#' taxonomy_path <- system.file("extdata/ITS2_Rosales_Restricted.tsv",
#'                              package = "CVrefDB")
#'
#' # Read the fasta and tsv files
#' fasta <- Biostrings::readDNAStringSet(fasta_path)
#' fasta
#' reftaxo <- read.table(taxonomy_path, sep = "\t", header = TRUE)
#' head(reftaxo)
#'
#' # Execute a 10 fold cross validation (k = 10 by default) on the reference database.
#' # NB to speed up the computation time we restrict ourselves to the first 3 folds
#' # (`folds_subset = 1:3`) and to the first 10 query sequences in each fold
#' # (`sequences_subset = 1:10`). With their default values, the process is executed
#' # on all sequences
#'
#' CV_results <-
#'     CV_blastn(fasta_db = fasta,
#'               taxo = reftaxo,
#'               label = "Test", # defaults to an empty string
#'               folds_subset = 1:3, # use only the 3 first folds (just for testing)
#'               sequences_subset = 1:10, # use onlythe first 10 sequences in each fold
#'               verbose = FALSE # silence the processing time logs
#'     )
#'
#' # The output contains 15 rows (by default max_target_seqs = 15) for each
#' # queried sequence (TaxID_query) with its "true" taxonomy.
#' # Each row contains one blast hit (TaxID_blast)
#' head(CV_results)
#'
#'
#'
#' # Example 2
#'
#' # Instead of reading the fasta file and its corresponding taxonomy, you can
#' # also simply pass their path to CV_blastn(). By default, the output will then
#' # be saved in a csv file on the disk (argument output) instead of returned
#' # as a data.frame by the function (unless output = NULL).
#'
#' # Retrieve the paths of 2 example files
#' #
#' # Remember that the taxonomy file MUST be a tab separated text file with 2 columns :
#' # col 1 = sequences IDs used in the fasta file, col 2 = full taxonomy in NCBI format
#'
#' fasta_path <- system.file("extdata/ITS2_Rosales_Restricted.fasta",
#'                           package = "CVrefDB")
#' taxonomy_path <- system.file("extdata/ITS2_Rosales_Restricted.tsv",
#'                              package = "CVrefDB")
#'
#'
#'
#' # First execute a 10 fold cross validation (k = 10) on the reference database.
#' # NB to speed up the computation time we restrict ourselves to the first 3 folds
#' # (`folds_subset = 1:3`) and to the first 10 query sequences in each fold
#' # (`sequences_subset = 1:10`).
#'
#' # File where the output will be saved. Here, we save it in the temporary
#' # "directory"
#'
#' output_path_10FoldCV <- paste0(tempdir(), "/ITS2_Rosales_Restricted_10FoldsCV.csv")
#'
#' CV_blastn(fasta_db = fasta_path,
#'           taxo = taxonomy_path,
#'           output = output_path_10FoldCV,
#'           k = 10, # 10-fold cross-validation
#'           label = "10FoldsCV",
#'           folds_subset = 1:3, # use only on the 3 first folds (just for testing)
#'           sequences_subset = 1:10, # use only the first 10  sequences in each fold
#'           seed = 123 # use any integer to have replicable results
#' )
#'
#'
#' # File where the output will be saved. Here, we save it in the temporary
#' # " directory"
#'
#' output_path_LeakedCV <- paste0(tempdir(), "/ITS2_Rosales_Restricted_LeakedCV.csv")
#'
#' CV_blastn(fasta_db = fasta_path,
#'           taxo = taxonomy_path,
#'           output = output_path_LeakedCV,
#'           k = 1, # Leaked cross validation
#'           label = "LeakedCV",
#'           sequences_subset = 1:30, # use only thefirst 30 query sequences
#'           seed = 123 # use any integer to have replicable results
#' )
#'
#'
#'
#' # Read the output files from the disk (NB with real large files, the use of
#' # faster approaches like `data.table::fread` is recommended).
#' # Note that the outputs are European csv files with semicolon as field separator
#' # and comma as decimal separator
#'
#' output_10FoldCV <- read.table(output_path_10FoldCV,
#'                               header = TRUE, sep =  ";", dec = ",")
#' output_leakedCV <- read.table(output_path_LeakedCV,
#'                               header = TRUE, sep =  ";", dec = ",")
#'
#' # Combine the 2 files. NB : the `Label` column allows you
#' # to keep track of the origin of the results once the outputs are combined
#'
#' full_output <- rbind(output_10FoldCV, output_leakedCV)
#' head(full_output)
#'
#' # Note that this output file does not assign a taxonomy. It just provides the
#' # n best blast hits. You can use the `assign_taxonomy()` function to do that.
#'
#' # A rough way to assign a taxonomy is to take the
#' # hit with highest Bit_score value and lowest Evalue (and in case of ties
#' # we simply take the first in the table...) :
#'
#' library(dplyr)
#' assigned <-
#'     full_output %>%
#'     group_by(Label, TaxID_query) %>%
#'     arrange(desc(Bit_score), E_value) %>%
#'     slice_head(n=1) %>%
#'     as.data.frame()
#' head(assigned)
#'
#'
#' # You can then for example compute the % of correct predictions
#' # at species, genus or family level for the leaked CV and 10 Folds CV
#' # approaches:
#'
#' assigned %>%
#'     group_by(Label) %>%
#'     summarize(
#'         Species = sum(Species_true == Species)*100/n(),
#'         Genus   = sum(Genus_true   == Genus)*100/n(),
#'         Family  = sum(Family_true  == Family)*100/n()
#'     )



CV_blastn <- function(

    fasta_db, # reference database in fasta format (path to a file or DNAStringSet object)
    taxo, # A 2 columns tsv file: Col 1 = Taxon ID corresponding to the fasta file name,  Col 2 = taxonomies
    output = if(!is.character(fasta_db)) {NULL} else {paste0(dirname(fasta_db), "/CV_blastn_output.csv")}, # name of the final output
    label = ifelse(is.null(output), "", output), # the first column of the output will contain label repeated

    # where should the '__Folds__' directory be created  ?
    # (containing the query and reference DB for each fold + intermediary files)
    folds_dir = ifelse(is.character(fasta_db), dirname(fasta_db), "."),

    # alternative fasta_db and the corresponding taxonomic file.
    # If not NA or NULL fasta_db is blasted against this alternative db instead of
    # blasted against itself. When k>1 the sequences from fasta_db in each fold
    # are removed from alternative_fasta_db based on their name
    alternative_fasta_db = NULL,
    alternative_taxo = NULL,

    # k = number of folds to create
    # if k==1 --> leaked validation the full database is blasted against itself
    # without removing any sequence
    k = 10,

    # fold = folds id to effectively blast against the database
    # eg :  if folds_subset = 3:5 , we will blast only the folds 3, 4 and 5
    folds_subset = 1:k,

    # sequences subset : ignored for any non numeric value
    # if a numeric vector is provided we will blast only these
    # sequences for each fold.
    # Eg sequences_subset = 1:10 --> we will blast only the first 10 sequences
    # for each fold
    sequences_subset = NULL,

    # filter out sequences from the query of the reference file based on their ID
    query_filter = NULL,
    db_filter = NULL,

    # Number of cores to be used by blastn
    nb_cores = parallel::detectCores(),

    # maximum number of sequences to be returned by blast
    max_target_seqs = 15,

    # id TRUE, prints messages with timing information
    verbose = TRUE,

    # should we try to clean the species names (remove subspecies, transform cf. into sp. ...)
    clean_species = TRUE,
    folds_remove = TRUE, # should the Folds directory be removed when finished ?
    seed = as.numeric(Sys.time()) # use any integer to have replicable results

) {

    if(verbose) {
        cat("\n--------------------------------------------------\n\n")
        cat(format(Sys.time(), "%Y-%m-%d %X")," : Start processing fasta file \n")
    }

    t00 <- Sys.time()


    # /*
    # ----------------------------- Read data ---------------------------
    # */


    # read fasta db
    if(is.character(fasta_db)){
        f <- Biostrings::readDNAStringSet(fasta_db)
    } else if(inherits(fasta_db, "DNAStringSet")) {
        f <- fasta_db
    } else {
        stop("fasta_db must be either a path to a fasta file or an object with class 'DNAStringSet' from package Biostrings")
    }

    # read alternative fasta bd if any
    if(!is.null(alternative_fasta_db)){

        if(is.character(fasta_db)){
            alt_f <- Biostrings::readDNAStringSet(alternative_fasta_db)
        } else if(inherits(fasta_db, "DNAStringSet")) {
            alt_f <- alternative_fasta_db
        } else {
            stop("alternative_fasta_db must be either NULL or a path to a fasta file or an object with class 'DNAStringSet' from package Biostrings")
        }
    }

    # read taxonomy reference database

    if(is.data.frame(taxo)){
        reftaxo <- taxo
    } else if(is.character(taxo)) {
        reftaxo <- unique(data.table::fread(taxo, sep = "\t", header = FALSE))
    } else {
        stop("taxo must be either a 2 columns data.frame or a character string providing a path to a tab separated text file")
    }
    colnames(reftaxo) <- c("Feature.ID", "Taxon")
    reftaxo <- cbind(reftaxo, split_taxonomy(reftaxo$Taxon,clean = clean_species))
    reftaxo$Taxon <- NULL


    if(!is.null(alternative_taxo)){
        if(is.data.frame(alternative_taxo)){
            alt_reftaxo <- unique(alternative_taxo)
        } else if(is.character(alternative_taxo)) {
            alt_reftaxo <- unique(data.table::fread(taxo, sep = "\t", header = FALSE))
        } else {
            stop("alternative_taxo must be either NULL or a 2 columns data.frame or a character string providing a path to a tab separated text file")
        }
        colnames(alt_reftaxo) <- c("Feature.ID", "Taxon")
        alt_reftaxo <- cbind(alt_reftaxo, split_taxonomy(alt_reftaxo$Taxon,
                                                         clean = clean_species))
        alt_reftaxo$Taxon <- NULL
    }


    # /*
    # ----------------------------- Create folds ---------------------------
    # */

    # create directory to store the query and reference DB for each fold
    # If this directory already exists delete it
    if(dir.exists(file.path(folds_dir, "__Folds__"))){
        unlink(file.path(folds_dir, "__Folds__"), recursive = T)
    }
    dir.create(file.path(folds_dir, "__Folds__"))


    # generate k folds indices: sequences selected at random among 100/k % of the database

    set.seed(seed)
    random_sample <- sample(1:length(f), replace = F)
    fold <- rep(paste0("Fold", sprintf("%02.0f", 1:k)),
                each = ceiling(length(f)/k))[1:length(f)]
    fold <- split(random_sample , fold)

    for (i in 1:length(fold)){

        # when k == 1, no real Cross validation but "leaked" cross-validation
        # : we blast all sequences against the full database without removing these
        # query sequences from the reference database

        if(is.null(alternative_fasta_db)){
            if(k==1){
                query <- f
                db <- f
            } else {
                query <- f[fold[[i]]] # subset of 100/k % random sequences
                db <- f[-fold[[i]]] # reference database of 100-k % of remaining sequences
            }
        } else {
            if(k==1){
                query <- f
                db <- alt_f
            } else {
                query <- f[fold[[i]]] # subset of 100/k % random sequences
                db <- alt_f[!names(alt_f) %in% names(query)] # keep only the names that are different from the query
            }
        }

        # filter out some of the query sequences based on their ID
        if(!any(is.na(query_filter)) & !any(is.null(query_filter))){
            query <- query[names(query) %in% query_filter]
        }

        # filter out some of the database sequences based on their ID
        if(!any(is.na(db_filter)) & !any(is.null(db_filter))){
            db <- db[names(db) %in% db_filter]
        }

        # filter out some of the query sequences based on their position in the fold
        if(is.numeric(sequences_subset) & length(sequences_subset) < length(query)){
            query <- query[sequences_subset]
        }

        if(length(query) == 0 | length(db) == 0){ next }

        Biostrings::writeXStringSet(query, paste0(folds_dir, "/__Folds__/", names(fold)[i], "_Query.fasta"))
        Biostrings::writeXStringSet(db, paste0(folds_dir,"/__Folds__/", names(fold)[i], "_DB.fasta"))
    }

    # t1 <- Sys.time()
    # time_count <- round(difftime(t1,t00, units = "mins"), 2)
    #
    # if(verbose) { cat("\n", k, "folds created in ", time_count, "mins \n\n") }



    # /*
    # ----------------------------- blast each fold ---------------------------
    # */

    # check that the folds_subset are present after filtering

    existing_folds <- list.files(paste0(folds_dir, "/__Folds__/"), "_Query.fasta")
    existing_folds <- gsub("_Query\\.fasta", "", existing_folds)

    folds_subset <- folds_subset[paste0("Fold", sprintf("%02.0f", folds_subset))
                                 %in% existing_folds]

    if(length(folds_subset)==0){
        stop("All the Folds selected with the option 'folds_subset' are empty after
         filtering with 'db_filter' and/or 'query_filter'")
    }

    for(i in folds_subset){

        fold <- sprintf("%02.0f", i)
        t0 <- Sys.time()

        if(verbose){
            cat("\n", format(Sys.time(), "%Y-%m-%d %X"), "- Blast Fold", fold)
        }

        fold_db <- paste0(folds_dir, "/__Folds__/Fold", fold, "_DB")
        fold_db_fasta <- paste0(folds_dir, "/__Folds__/Fold", fold, "_DB.fasta")
        fold_query <- paste0(folds_dir, "/__Folds__/Fold", fold, "_Query.fasta")
        fold_output <- paste0(folds_dir, "/__Folds__/Output_Fold", fold, ".tsv")


        blastn(
            db = fold_db_fasta,
            query = fold_query,
            out = fold_output,
            retrieve_nonmatching = TRUE,
            num_threads = nb_cores,
            max_target_seqs = max_target_seqs,
            evalue = 10,
            verbose = 0,  # 0 for no messages, 1 for time used, 2 for blast DB creation messages
            clean = FALSE # remove the files created (blast db based on a fasta file)
        )


        if(verbose){
            t1 <- Sys.time()
            time_count = round(difftime(t1, t0, units = "mins"),2)
            cat(" - Job finished in", time_count, "minutes")
        }
    }



    # /*
    # ---------------- Concatenate the results of the different folds  -------------------
    # */


    f <- list.files(paste0(folds_dir, "/__Folds__"), pattern = "Output_", full.names = T)

    d <- vector(length = length(f), mode = "list")
    for(i in 1:length(d)){
        tmp <- data.table::fread(f[i], sep = "\t", header = F, stringsAsFactors = F)
        colnames(tmp) <- c("TaxID_Query", "TaxID_Blast", "E_Value", "Bit_Score", "Length", "Identity")
        tmp <- data.frame(Label = label,
                          Fold = gsub("Output_Fold(\\d{2,4})\\.tsv", "F\\1", basename(f[i])),
                          tmp)
        d[[i]] <- tmp
    }

    d <- as.data.frame(do.call(rbind, d))




    # /*
    # -------------------------- Merge blast output with taxonomic information ---------------------------
    # */


    tmp <-
        d  %>%
        dplyr::left_join(reftaxo[, c("Feature.ID", "Kingdom", "Phylum", "Class",
                                     "Order", "Family", "Genus", "Species", "Species_orig")],
                         by = c("TaxID_Query" = "Feature.ID"))  %>%
        dplyr::select(Label, Fold, TaxID_Query, Kingdom_true = Kingdom,
                      Phylum_true = Phylum, Class_true = Class,
                      Order_true = Order, Family_true = Family,
                      Genus_true = Genus, Species_true = Species,
                      Species_orig_true = Species_orig,
                      TaxID_Blast, E_Value, Bit_Score, Length, Identity)

    if(is.null(alternative_fasta_db)){

        tmp <-
            tmp %>%
            dplyr::left_join(reftaxo[, c("Feature.ID", "Kingdom", "Phylum", "Class",
                                         "Order", "Family", "Genus", "Species", "Species_orig")],
                             by = c("TaxID_Blast" = "Feature.ID"))  %>%
            # uniformize the case n column names : first letter = uppercase, after_ = lower case
            dplyr::rename(TaxID_query = TaxID_Query,
                          TaxID_blast = TaxID_Blast,
                          E_value = E_Value,
                          Bit_score = Bit_Score)

    } else {

        tmp <-
            tmp   %>%
            dplyr::left_join(alt_reftaxo[, c("Feature.ID", "Kingdom", "Phylum", "Class",
                                             "Order", "Family", "Genus", "Species", "Species_orig")],
                             by = c("TaxID_Blast" = "Feature.ID")) %>%
            # uniformize the case n column names : first letter = uppercase, after_ = lower case
            dplyr::rename(TaxID_query = TaxID_Query,
                          TaxID_blast = TaxID_Blast,
                          E_value = E_Value,
                          Bit_score = Bit_Score)

    }



    if(!is.null(output)){
        data.table::fwrite(tmp, output,
                           row.names = F, sep = ";", dec = ",", quote = TRUE)
    }

    if(folds_remove){
        unlink(file.path(folds_dir, "__Folds__"), recursive = T)
    }


    t1 <- Sys.time()
    time_count <- round(difftime(t1,t00, units = "mins"), 2)

    if(verbose){
        cat("\n\nProcessing of fasta file finished in", time_count, "mins")
        cat("\n\n--------------------------------------------------\n")
    }

    if(is.null(output)){
        return(tmp)
    }


}





