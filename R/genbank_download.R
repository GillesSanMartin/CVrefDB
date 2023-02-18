

# /*
# ----------------------------- genbank_count ---------------------------
# */


#' Number of sequences available on GenBank for a given Entrez query
#'
#' Uses the NCBI API to retrieve the number of sequences available in GenBank
#' based on the provided Entrez query
#'
#' @param entrez_query A character string containing the Entrez query.
#' It can contain spaces and new lines (which
#' are replaced internally by `+` signs to make the query string valid).
#'
#' @return An integer
#' @export
#'
#' @examples
#'
#' genbank_count("Adalia[organism] AND (COI[Title] or CO1[Title])")
#'
genbank_count <- function(entrez_query){

    # replace spaces and new lines by "+"
    entrez_query <- gsub(" |\\\n", "\\+", entrez_query)

    # build url and add API key if available in .Renviron
    my_url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&rettype=count"
    # ENTREZ_KEY <- Sys.getenv("ENTREZ_KEY")
    # if(ENTREZ_KEY != "") my_url <- paste0(my_url, "&api_key=", ENTREZ_KEY)
    my_url <- paste0(my_url, "&term=", entrez_query)

    my_count <-
        httr::RETRY("GET", url = my_url, times = 20, pause_cap = 300) %>%
        xml2::read_xml() %>%
        xml2::xml_integer()
    return(my_count)
}
# genbank_count(entrez_query)



# /*
# ----------------------------- genbank_download ---------------------------
# */




#' Downloads FASTA sequences and taxonomy information from NCBI GenBank database
#' based on a given Entrez query.
#'
#' The sequences and their taxonomy are downloaded by batches of size = `retmax`
#' and saved on the disk in two separated files.
#' In case of server error, a new request is sent repeatedly with increasing
#' waiting time. Typical usage : build a local reference database for barcoding
#' or metabarcoding.
#'
#' @param entrez_query A character string containing the Entrez query used to
#' retrieve sequences from GenBank. It can contain spaces and new lines (which
#' are replaced internally by `+` signs to make the query string valid).
#' @param retmax An integer value indicating the maximum number of sequences
#' to retrieve per batch (default: 350). Higher values than 350 will
#' most likely cause systematically server errors. This values shoulp
#' not be changed in most of the cases.
#' @param seq_file_name A character string specifying the name of the output
#' FASTA file (default: "sequences.fasta").
#' @param taxo_file_name A character string specifying the name of the output
#' taxonomy file (default: "taxonomy.tsv").
#' @param append A logical value indicating whether to append the sequences
#' and taxonomy to existing files or create new ones (default: FALSE).
#' If FALSE, any existing file with the names indicated by `seq_file_name` and
#' `taxo_file_name` will be silently erased.
#' @param verbose A logical value indicating whether to display progress
#' messages and time elapsed (default: TRUE).
#' @param remove_tmp_files A logical value indicating whether to remove
#' 2 temporary files (default: TRUE) : 1) `__accver_taxid__.tsv` containing
#' the accession.version IDs and their corresponding taxid.
#' 2) `__taxonomy_long__.tsv`
#'
#' @return The function returns nothing but saves the results into 2 files
#' on the disk. The fasta file containing the sequences is progressively updated
#' at each batch. The taxonomy file is computed when all batches are completed.
#'
#' @examples
#'
#' \dontrun{
#' # COI sequences for the genus Pullex
#' # Note that we can have spaces and new lines inside the query
#'
#' my_query <-
#' "Pulex[Organism] AND
#' (CO1[Title] OR COI[Title] OR
#' cytochrome c oxidase subunit 1[Title] OR cytochrome c oxidase subunit I[Title] OR
#' cytochrome oxidase subunit 1[Title] OR cytochrome oxidase subunit I[Title]) AND
#' 100:10000000[Sequence Length]) NOT
#' (uncultured OR environmental sample OR incertae sedis OR unverified)"
#'
#' genbank_download(my_query)
#'
#' # We can download the sequences of another flea genus and append the results
#' # to the existing files :
#'
#' my_query <- "Xenopsylla[Organism] AND (CO1[Title] OR COI[Title])"
#' genbank_download(my_query, append = TRUE)
#'
#'
#' # Inspect the results
#'
#' fasta <- Biostrings::readDNAStringSet("sequences.fasta")
#' fasta
#'
#' tax_wide <- read.table("taxonomy.tsv",  sep= "\t", header = TRUE)
#' head(tax_wide)
#' tail(tax_wide)
#'
#' # NB : before the database can be used, there are typically a few
#' # cleaning steps to perform after the download (remove sequences with
#' # too much homopolymers or ambiguous sequences, dereplicate, ...).
#' # See vignette "Build_database" for guidance.
#' #
#'
#' }
#'
#' @export


genbank_download <- function(
        entrez_query,
        retmax = 350,
        seq_file_name = "sequences.fasta",
        taxo_file_name = "taxonomy.tsv",
        append = FALSE,
        verbose = TRUE,
        remove_tmp_files = TRUE
){

    # remove existing temporary files
    if(file.exists("__accver_taxid__.tsv")) {unlink("__accver_taxid__.tsv")}
    if(file.exists("__taxonomy_long__.tsv")) {unlink("__taxonomy_long__.tsv")}

    if(append == FALSE){
        if(file.exists(taxo_file_name)) {unlink(taxo_file_name)}
        if(file.exists(seq_file_name)) {unlink(seq_file_name)}
    }


    # replace spaces and new lines by "+"
    entrez_query <- gsub(" |\\\n", "\\+", entrez_query)
    # retrieve NCBI API key if available in .Renviron --> probably not useful
    ENTREZ_KEY <- Sys.getenv("ENTREZ_KEY")


    # Initialization of the first sequence to be retrieved
    # this value will be updated in each repetition of the while loop
    retstart = 0
    # initialization of the vector that will collect the accession numbers
    my_accession <- vector(length = retmax)

    if(verbose) {
        t0 <- Sys.time()
        my_count <- genbank_count(entrez_query)
    }

    while(length(my_accession) == retmax) {

        if(verbose) {cat("Downloading sequences", retstart,
                         "to", retstart + retmax, "from", my_count, "sequences")}

        # retrieve the accession numbers corresponding to this entrez query
        my_url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore"
        # if(ENTREZ_KEY != "") my_url <- paste0(my_url, "&api_key=", ENTREZ_KEY)

        my_url <- paste0(my_url,
                         "&rettype=uilist&idtype=acc&retmode=xml",
                         "&retmax=", retmax,
                         "&retstart=", retstart,
                         "&term=")

        my_url <- paste0(my_url, entrez_query)

        my_xml <- httr::RETRY("GET", url = my_url, times = 20, pause_cap = 300) %>%
            httr::content("text", encoding = "UTF-8") %>%
            xml2::read_xml()

        my_accession <-
            my_xml %>%
            xml2::xml_find_all("//IdList/Id") %>%
            xml2::xml_text()

        # retrieve the fasta and taxid corresponding to these accession numbers

        my_url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&retmode=xml"
        # if(ENTREZ_KEY != "") my_url <- paste0(my_url, "&api_key=", ENTREZ_KEY)

        my_url  <- paste0(my_url,
                          "&id=",
                          do.call(paste, as.list(c(my_accession, sep = ","))))

        my_xml <- httr::RETRY("GET", url = my_url, times = 20, pause_cap = 300) %>%
            httr::content("text", encoding = "UTF-8") %>%
            xml2::read_xml()

        my_sequences <-
            my_xml %>%
            xml2::xml_find_all("//TSeq_sequence") %>%
            xml2::xml_text()

        my_accver <-
            my_xml %>%
            xml2::xml_find_all("//TSeq_accver") %>%
            xml2::xml_text()

        my_taxid <-
            my_xml %>%
            xml2::xml_find_all("//TSeq_taxid") %>%
            xml2::xml_text()

        my_fasta <- paste0(">", my_accver, "\n", my_sequences)
        accver_taxid <- data.frame(accession = my_accver, taxid = my_taxid)


        # retrieve the taxonomy lineages corresponding to these taxid
        my_url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&rettype=xml"
        # if(ENTREZ_KEY != "") my_url <- paste0(my_url, "&api_key=", ENTREZ_KEY)

        my_url  <- paste0(my_url,
                          "&id=",
                          do.call(paste, as.list(c(unique(my_taxid), sep = ","))))

        my_xml <- httr::RETRY("GET", url = my_url, times = 20, pause_cap = 300) %>%
            httr::content("text", encoding = "UTF-8") %>%
            xml2::read_xml()

        # transform the xml taxonomy hierarchical structure into a 2 dimensions table
        taxonomy_long <-
            xml2::as_list(my_xml) %>%
            tibble::as_tibble() %>%
            tidyr::unnest_wider(TaxaSet) %>%
            dplyr::select(TaxId, ScientificName, LineageEx) %>%
            tidyr::unnest_longer(c(TaxId, ScientificName, LineageEx)) %>%
            dplyr::select(taxid = TaxId, name = ScientificName, lineage = LineageEx) %>%
            tidyr::unnest_wider(lineage) %>%
            tidyr::unnest_longer(c(ScientificName, Rank)) %>%
            dplyr::select(taxid, name, taxon = ScientificName, rank = Rank)

        # save the data of the current batch on the disk
        write(my_fasta, seq_file_name, append = TRUE)

        suppressWarnings(
            utils::write.table(accver_taxid, "__accver_taxid__.tsv",  append = TRUE,
                               sep= "\t", quote = TRUE, row.names = FALSE,
                               col.names = !file.exists("__accver_taxid__.tsv"))
        )
        suppressWarnings(
            utils::write.table(taxonomy_long, "__taxonomy_long__.tsv", append = TRUE,
                               sep= "\t", quote = TRUE, row.names = FALSE,
                               col.names = !file.exists("__taxonomy_long__.tsv"))
        )

        # increment
        retstart = retstart + retmax

        if(verbose) {
            t1 <- round(difftime(Sys.time(), t0, units = "mins"),2)
            cat("- Time elapsed :", t1, "minutes\n")
        }
    }



    # Merge accession numbers with the full taxonomy lineage in wide format

    accver_taxid <- utils::read.table("__accver_taxid__.tsv", sep= "\t",  header = TRUE)
    tax_long <- utils::read.table("__taxonomy_long__.tsv", sep= "\t",  header = TRUE)

    # we will keep only the main taxonomic ranks
    # --> you lose the ranks "Clade" and "no rank" and  some others
    taxonomic_ranks_of_interest <-
        c("superkingdom", "kingdom", "subkingdom",
          "superphylum", "phylum", "subphylum", "infraphylum",
          "superclass", "class", "subclass", "infraclass",
          # "cohort", "subcohort",
          "superorder", "order", "suborder", "infraorder", "parvorder",
          "superfamily", "family", "subfamily",
          "tribe", "subtribe",
          "genus", "subgenus",
          # "section", "series",
          "species", "subspecies")

    # We attach an empty taxon with all taxonomic levels.
    # The aim is to always obtain a final wide table with the
    # same columns in the same order...
    tax_long <- rbind(
        data.frame(taxid = NA, name = NA, taxon = NA,
                   rank = taxonomic_ranks_of_interest),
        tax_long
    )

    tmp <-
        tax_long %>%
        dplyr::filter(rank %in% taxonomic_ranks_of_interest) %>%
        dplyr::distinct() %>%
        tidyr::pivot_wider(names_from = rank, values_from = taxon ) %>%
        dplyr::slice(-1)

    tax_wide <-
        accver_taxid %>%
        dplyr::left_join(tmp, by = "taxid")

    suppressWarnings(
        utils::write.table(tax_wide, taxo_file_name, append = TRUE,
                           sep= "\t", quote = TRUE, row.names = FALSE,
                           col.names = !file.exists(taxo_file_name))
    )

    if(remove_tmp_files){
        unlink("__accver_taxid__.tsv")
        unlink("__taxonomy_long__.tsv")
    }

}



# For testing

# taxo_file_name = "taxonomy.tsv"
# seq_file_name = "sequences.fasta"
# append = FALSE
# verbose = TRUE
# remove_tmp_files = TRUE
# retmax = 350
#
# entrez_query <- "
# Pulex[Organism] AND
# (CO1[Title] OR COI[Title] OR
# cytochrome c oxidase subunit 1[Title] OR cytochrome c oxidase subunit I[Title] OR
# cytochrome oxidase subunit 1[Title] OR cytochrome oxidase subunit I[Title]) AND
# 100:10000000[Sequence Length]) NOT
# (uncultured OR environmental sample OR incertae sedis OR unverified)"
#









