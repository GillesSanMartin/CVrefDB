

# /*
# ----------------------------- assign_taxonomy() ---------------------------
# */



#' Assign taxonomy to sequences based on top n blast hits
#'
#' This function allows to assign a taxonomy to sequences based on their top
#' blastn hits. It can use different methods for the assignment based on the bit
#' score, the consensus score and it allows optionally to specify the minimum
#' identity, length and E-value to take into consideration for each taxonomic
#' level. The output contains a taxonomic assignment for each taxonomic level
#' (Species, genus, family,...) and identity and consensus scores that can be
#' used to decide to dismiss untrustworthy identifications.
#'
#'
#' @details
#' This function assigns a taxonomy for each query sequence and each requested.
#' taxonomic level. I.e. for a given sequence, you will get a proposed
#' identification for the family level, the genus level, the species level,...
#' For each assignment, the function also provides several statistics to help
#' you assess the reliability of the assignment, in particular, the Identity
#' score and the Consensus score (both expressed as %).
#'
#' The Identity score (%), directly provided by blast, gives the % of base pairs
#' (bp) that are identical between the target sequence and the sequence matched
#' in the reference database. Once a taxon is assigned, we attribute the
#' Identity score corresponding to the maximum bit score for this taxon. NB :
#' this is not necessarily the maximum identity score because the bit score
#' penalizes matches with low values for alignment length (i.e. it is better to
#' have 98% of Identity with a 400 pb alignment length than to have 100%
#' identity with a 50 pb alignment length. . . In such a case, 98% will be the
#' identity score assigned to the taxon.).
#'
#' The Consensus score is simply the % of the 10 top hits (by default) matching the same taxon.
#' If among the top 10 hits, blast matches 2 Brassica napus, 6 Brassica nigra
#' and 2 Sinapis alba, the consensus score is 60% for Brassica nigra,
#' 80% for the genus Brassica and 100% for the family Brassicaceae.
#'
#' Four different methods can be chosen :
#'
#' * `TopHit` : we simply assign the taxon with the best Bit score (default blast
#' output if we keep only 1 match).
#' * `TopHitPlus` : we also choose the taxon with the best bit score but we
#' break the potential ties by choosing the taxon with the highest consensus
#' score among the ties.
#' * `TopN` : we choose among the top 10 hits (by defaults `Top_n = 10`),
#' the taxon with the highest  consensus score at each taxonomic level.
#' With this approach the chosen species might be in a different genus
#' than the taxon chosen at the genus level for example.
#' * `TopNPlus` : we first eliminate the sequences which do not fulfill some
#' predefined requirements, then we compute the consensus score and choose the
#' taxon with the highest consensus score at each taxonomic level. With defaults
#' values, we eliminate the matched sequences with a length < 100 base pairs
#' (`Length_min = 100`) and with an E score > 10^-10 (`E_max = 10^-10`). The
#' identity score must also be at least 97% for the species level
#' (`Ident_min_s = 97`), 90% for the genus level (`Ident_min_g = 90`), 80% for the family
#' level (`Ident_min_f = 80`),... Also the difference of identity score between
#' the best hit and the matched sequence must be <5% for the Species level
#' (`Ident_max_diff_s = 5%`), 10% for the genus level (`Ident_max_diff_g = 10%`),
#' 20% for the family level (`Ident_max_diff_g = 20%`),...
#' (approach inspired e.g. by Milla et al.(2022) and
#' [https://github.com/Joseph7e/Assign-Taxonomy-with-BLAST](https://github.com/Joseph7e/Assign-Taxonomy-with-BLAST)
#'


#'
#' @param df a data frame with the top blastn hits for each query sequence ID
#'   (typically, the output of the function [CV_blastn()]. For each hit, it
#'   should contain 4 blastn statistics (Identity, Length, bit score, E-value)
#'   and a column for each taxonomic level of interest ( e.g. : `Species`,
#'   `Genus`, `Family`, `Order`, `Class`, `Phylum`, `Kingdom`).
#' @param taxo  Either a 2-column data.frame with
#'   sequence ID (col 1) and QIIME2 taxonomy (col2) or a a character string
#'   giving the path to a tsv file containing such information (tab separated
#'   text file). This file is optional and can
#'   be used to add the "true" taxon as a supplementary column in the output.
#'   Default is `NULL` and then, the output does not contain a "Taxon_true" column.
#' @param Assignation_method a character vector with at least one of the
#'   following values: "TopHit", "TopHitPlus", "TopN", "TopNPlus". It specifies
#'   the method used for the assignment (see details). By default, the results
#'   for the 4 methods are provided.
#' @param Top_n an integer giving the number of top hits to consider when using
#'   the "TopN" or "TopNPlus" methods and to compute the consensus score (for
#'   all methods). Default is 10.
#' @param ID a character string giving the name of the column with the query
#'   taxonomic ID. Default is "TaxID_query".
#'
#' @param Ident a character string giving the name of the column containing the
#'   identity. Default is "Identity".
#' @param Bit a character string giving the name of the column containing the
#'   bit score. Default is "Bit_score".
#' @param Length a character string giving the name of the column containing the
#'   length. Default is "Length". This value is only used by TopNPlus method but
#'   the Length is also displayed in the final output. Its value can also
#'   be set to `NA` if "TopNPlus" is not one of the methods (the column Length
#'   in the output is then a column of NAs)
#' @param E a character string giving the name of the column containing the
#'   E-value. This parameter is only used with `TopNPlus` method.
#'   Default is "E_value" if any of the Methods is `TopNPlus` and `NA` otherwise.
#'
#' @param Species a character string giving the name of the column containing
#'   the species level taxonomic identification. Default is "Species".
#'   Use `NA` instead to ignore this level.
#' @param Genus a character string giving the name of the column containing the
#'   genus level taxonomic identification. Default is "Genus".
#'   Use `NA` instead to ignore this level.
#' @param Family a character string giving the name of the column containing the
#'   family level taxonomic identification. Default is "Family".
#'   #'   Use `NA` instead to ignore this level.
#' @param Order a character string giving the name of the column containing the
#'   order level taxonomic identification. Default is "Order".
#'   Use `NA` instead to ignore this level.
#' @param Class a character string giving the name of the column containing the
#'   class level taxonomic identification. Default is `NA` and the class level
#'   is ignored.
#' @param Phylum a character string giving the name of the column containing
#'   the Phylum level taxonomic identification. Default is `NA` and the Phylum
#'   level is ignored.
#' @param Kingdom a character string giving the name of the column containing
#'   the kingdom level taxonomic identification. Default is `NA` and the kingdom
#'   level is ignored.
#'
#' @param E_max a numeric giving the maximum E-value to consider when using the
#'   "TopNPlus" method. Default is `10^-10`.
#' @param Length_min an integer giving the minimum length to consider when using
#'   the "TopNPlus" method. Default is 100.
#'
#' @param Ident_min_s a numeric giving the minimum identity for the species
#'   level to consider when using the "TopNPlus" method. Default is 97.
#' @param Ident_min_g a numeric giving the minimum identity for the genus level
#'   to consider when using the "TopNPlus" method. Default is 90.
#' @param Ident_min_f a numeric giving the minimum identity for the family level
#'   to consider when using the "TopNPlus" method. Default is 80.
#' @param Ident_min_o a numeric giving the minimum identity for the order level
#'   to consider when using the "TopNPlus" method. Default is 70.
#' @param Ident_min_c a numeric giving the minimum identity for the class level
#'   to consider when using the "TopNPlus" method. Default is 60.
#' @param Ident_min_p a numeric giving the minimum identity for the Phylum
#'   level to consider when using the "TopNPlus" method. Default is 50.
#' @param Ident_min_k a numeric giving the minimum identity for the kingdom
#'   level to consider when using the "TopNPlus" method. Default is 40.
#' @param Ident_max_diff_s  Minimum difference between the best identity and the
#'   identity of each hit when using the "TopNPlus" method. A different value
#'   can be used for each taxonomic level (cf next parameters). If Identity of
#'   best Bit score is 97% and Ident_max_diff_s = 5, we will consider only
#'   Identities > 92% for the species level identifications. Default is 5.
#'
#' @param Ident_max_diff_g same as `Ident_max_diff_s` for the genus level.
#'   Default is 10
#' @param Ident_max_diff_f same as `Ident_max_diff_s` for the family level.
#'   Default is 20
#' @param Ident_max_diff_o same as `Ident_max_diff_s` for the order level.
#'   Default is 30
#' @param Ident_max_diff_c same as `Ident_max_diff_s` for the class level.
#'   Default is 40
#' @param Ident_max_diff_p same as `Ident_max_diff_s` for the Phylum level.
#'   Default is 50
#' @param Ident_max_diff_k same as `Ident_max_diff_s` for the kingdom level.
#'   Default is 60
#'
#' @return A table with the following columns :
#'   * Method : assignment method
#'   used (TopHit, TopHitPlus, TopN, TopNPlus)
#'   * Tax_level : the taxonomic level of the prediction (eg Family, Genus,
#'   Species,...)
#'   * Level : an abreviated version using the forst letter of the taxonomic
#'   level (f = Family, g = genus, s = species,...)
#'   * ID : the ID of the query sequence
#'   * Taxon : the assigned (predicted) taxon for that sequence ID at that
#'   taxonomic level.
#'   * Bit : Bit score of the assigned taxon
#'   * Length : Length of the alignment corresponding to the assigned taxon
#'   * Ident : Identity score (%) of the assigned taxon
#'   * Cons : Consensus score of the assigned taxon (% of blastn hits agreeing
#'   on that taxon among the `Top_n` top hits)
#'   * Taxon_true : This column will be present only when `taxo` is not NULL.
#'   It contains the "true" taxon associated with this sequence ID according to
#'   the original reference database.
#' @export
#'
#' @examples
#'
#' library(dplyr)
#' library(tidyr)
#'
#' # 1. Example with a fake dataset
#'
#' # Fake top 15 blast hits for a sequence corresponding to Scirpoides holoschoenus
#' # and inspired by a real life example where the 4 assignment methods give
#' # different results
#'
#' example_top15_blast_hits <-
#'     data.frame(
#'         TaxID_query = "ID001",
#'         Species_true = "Scirpoides holoschoenus",
#'
#'         Bit_score = c(1075,704,682,651,630,630,549,549,529,525,525,520,520,518,512),
#'         Length = c(589,607,608,595,597,597,568,568,570,541,571,616,609,569,569),
#'         Identity = c(99.660,88.138,87.500,86.891,86.265,86.265,84.683,84.507,83.860,
#'                      84.843,83.713,82.630,82.594,83.480,83.304),
#'         E_value = 10^-14,
#'
#'         Phylum = "Streptophyta",
#'         Class = "Magnoliopsida",
#'         Order = "Poales",
#'         Family = "Cyperaceae",
#'         Genus = c("Scirpoides","Erioscirpus","Erioscirpus","Dracoscirpoides",
#'                   "Dracoscirpoides","Dracoscirpoides", "Cyperus","Cyperus","Cyperus",
#'                   "Cyperus","Cyperus","Cyperus","Cyperus","Cyperus","Cyperus"),
#'         Species = c("Scirpoides holoschoenus","Erioscirpus microstachyus",
#'                     "Erioscirpus comosus", "Dracoscirpoides ficinioides",
#'                     "Dracoscirpoides falsa", "Dracoscirpoides falsa",
#'                     "Cyperus pulchellus", "Cyperus leucocephalus"," Cyperus haspan",
#'                     "Cyperus spiralis", "Cyperus flaccidus", "Cyperus sp.",
#'                     "Cyperus diffusus", "Cyperus tenuispica", "Cyperus haspan")
#'     )
#'
#' # corresponding true taxonomy
#' example_taxonomy <-
#' data.frame(ID = "ID001",
#' Taxo = "k__Viridiplantae; p__Streptophyta; c__Magnoliopsida; o__Poales; f__Cyperaceae; g__Scirpoides; s__holoschoenus")
#'
#'
#' # Assign the taxonomy with mostly default options and in "long" format
#'
#' tax_pred <-
#'     assign_taxonomy(df = example_top15_blast_hits,
#'                     taxo = example_taxonomy,
#'                     Order = NA, Family = NA) # ignore the family and genus levels
#' tax_pred
#'
#'
#' # Present the same data in a different, wider, format
#'
#' tax_pred %>%
#'     select(-c(Bit, Length)) %>%
#'     pivot_assign_taxonomy()
#'
#'
#'
#'
#' # change some of the many available options
#' tax_pred <-
#'     assign_taxonomy(df = example_top15_blast_hits,
#'                     taxo = example_taxonomy,
#'                     # use only the TopHitPlus and TopNPlus methods
#'                     Assignation_method = c("TopHitPlus", "TopNPlus"),
#'                     # use the top 15 blast hits instead of the default 10
#'                     Top_n = 15,
#'                     # change the defaults for TopNPlus
#'                     Ident_min_s = 80, Ident_max_diff_s = 20,
#'                     # Consider also the Phylum and Class level by indicating the
#'                     # name of the columns containing this information
#'                     Class = "Class", Phylum = "Phylum"
#'     )
#' tax_pred
#'
#'
#'
#' # ----------------------------------------------
#'
#' # 2 Example on a real dataset
#'
#' # Retrieve an example reference database (fasta file) and its corresponding
#' # taxonomic information from the package data examples.
#' # This example contains >7000 sequences for the ITS2 barcode of plants from the
#' # order "Rosales" restricted to a portion of the gene amplified by a given
#' # given primer used in metabarcoding.
#' #
#' # Remember that the taxonomy file MUST be a tab separated text file with 2 columns :
#' # col 1 = sequences IDs used in the fasta file, col 2 = full taxonomy in QIIME2 format
#'
#' fasta_path <- system.file("extdata/ITS2_Rosales_Restricted.fasta",
#'                           package = "CVrefDB")
#' taxonomy_path <- system.file("extdata/ITS2_Rosales_Restricted.tsv",
#'                              package = "CVrefDB")
#'
#'
#'
#' # First, execute a 10 fold cross validation (k = 10) on the reference database.
#'
#' # File where the output will be saved. Here, we save it in the temporary
#' # " directory"
#'
#' output_path_10FoldCV <- paste0(tempdir(), "/ITS2_Rosales_Restricted_10FoldsCV.csv")
#'
#' CV_blastn(fasta_db = fasta_path,
#'           taxo = taxonomy_path,
#'           output = output_path_10FoldCV,
#'           k = 10, # 10-fold cross-validation
#'           label = "Restricted_10FoldsCV",
#'           folds_subset = 1:5, # use only on the 5 first folds (just for testing)
#'           sequences_subset = 1:25, # use only the first 25  sequences in each fold
#'           folds_remove = TRUE, # should the __Folds__ dir be removed when finished ?
#'           seed = 123 # use any integer to have replicable results
#' )
#'
#' # Read the output files from the disk (NB with real large files, the use of
#' # faster approaches like `data.table::fread` is recommended)
#'
#' output_10FoldCV <- read.table(output_path_10FoldCV,
#'                               header = TRUE, sep =  ";", dec = ",")
#'
#' # Assign the taxonomy
#' assigned <- assign_taxonomy(df = output_10FoldCV,
#'                             taxo = taxonomy_path, # to add the "true" taxon
#'                             Assignation_method = c("TopHit", "TopHitPlus",
#'                                                    "TopN", "TopNPlus"),
#'                             Top_n = 10, # number of top hits to consider
#'                             Order = NA # skip the Order level
#' )
#'
#' # The ouput contains the assigned/predicted taxon (`Taxon`) and real taxon
#' # (`Taxon_true`) for each query sequence ID, each Method and each taxonomic
#' # level along with their identity and consensus scores
#' head(assigned)
#'
#'
#'
#' # We can now explore this results eg by checking the % of correct identifications
#' # for each taxonomic level and each method (NB the missing assignments are
#' # not counted here as correct identification)
#' library(dplyr)
#' assigned %>%
#'     mutate(CorrectID = Taxon == Taxon_true & !is.na(Taxon)) %>%
#'     group_by(Method, Tax_level) %>%
#'     summarize(CorrectID = sum(CorrectID)*100/n())
#'
#'
#' # We can explore also the relationship between the identity or consensus score
#' # with the proportion of correct identification
#' library(ggplot2)
#'
#' assigned %>%
#'     mutate(CorrectID = Taxon == Taxon_true & !is.na(Taxon)) %>%
#'     filter(Method == "TopHit") %>%
#'     ggplot(aes(y = as.numeric(CorrectID), x = Ident, group = Method)) +
#'     facet_wrap(~Tax_level, ncol = 4) +
#'     geom_point(alpha = 0.2,
#'                position = position_jitter(width = 0, height = 0.05)) +
#'     stat_smooth(method = "glm", se = FALSE,  method.args = list(family = binomial))+
#'     ylab("Proportion of correct identification") +
#'     xlab("Identity score (%)") +
#'     theme_bw()
#'
#' assigned %>%
#'     mutate(CorrectID = Taxon == Taxon_true & !is.na(Taxon)) %>%
#'     filter(Method == "TopHit") %>%
#'     ggplot(aes(y = as.numeric(CorrectID), x = Cons, group = Method)) +
#'     facet_wrap(~Tax_level, ncol = 4) +
#'     geom_point(alpha = 0.2,
#'                position = position_jitter(width = 0, height = 0.05)) +
#'     stat_smooth(method = "glm", se = FALSE,  method.args = list(family = binomial))+
#'     ylab("Proportion of correct identification") +
#'     xlab("Consensus score (%)") +
#'     theme_bw()
#'




assign_taxonomy <- function(

    df, # data frame or data.table ?
    taxo = NULL, # An optional 2 columns file with sequence ID (col 1) and QIIME2 taxonomy (col2)

    # a character vector with at least one of these 4 values
    Assignation_method = c("TopHit", "TopHitPlus", "TopN", "TopNPlus"),
    Top_n = 10, # number of top hits to consider

    # Name of the columns with the query taxonomic ID
    ID = "TaxID_query",

    # Name of the columns containing 4 descriptive statistics
    # Length and E are only used for TopNPlus. But Length is displayed in the final results
    Ident = "Identity",
    Bit = "Bit_score",
    Length = "Length", # can also be "NA"
    E = ifelse(any(Assignation_method == "TopNPlus"), "E_value", NA),

    # Name of the columns containing the various taxonomic level IDs
    # If = NA --> this taxonomic level is ignored (not yet functional though...)
    Species = "Species",
    Genus = "Genus",
    Family = "Family",
    Order = "Order",
    Class = NA,
    Phylum = NA,
    Kingdom = NA,

    # Thresholds values used only with TopNPlus method

    E_max = 10^-10,
    Length_min = 100,

    # Minimum identity for each taxonomic level
    # s = species, g = genus, f = family,...
    Ident_min_s = 97,
    Ident_min_g = 90,
    Ident_min_f = 80,
    Ident_min_o = 70,
    Ident_min_c = 60,
    Ident_min_p = 50,
    Ident_min_k = 40,

    # Minimum difference between the best identity and the identity of each
    # hit. A different value can be used for each taxonomic level
    # If Identity of best Bit score is 97% and Ident_max_diff_s = 5, we will
    # consider only Identities > 92%
    Ident_max_diff_s = 5,
    Ident_max_diff_g = 10,
    Ident_max_diff_f = 20,
    Ident_max_diff_o = 30,
    Ident_max_diff_c = 40,
    Ident_max_diff_p = 50,
    Ident_max_diff_k = 60

){

    # Select and rename the columns of interest

    orig_names <- c(ID, Kingdom, Phylum, Class, Order, Family, Genus, Species,
                    E, Bit, Ident, Length)
    my_names <- c("ID", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species",
                  "E", "Bit", "Ident", "Length")

    my_names <- my_names[!is.na(orig_names)]
    orig_names <- orig_names[!is.na(orig_names)]

    df <- as.data.frame(df)[,orig_names]
    colnames(df) <- my_names

    TopHit <- NA
    TopHitPlus <- NA
    TopN <- NA
    TopNPlus <- NA

    taxonomic_levels <- c(Kingdom, Phylum, Class, Order, Family, Genus, Species)
    names(taxonomic_levels) <- c("Kingdom", "Phylum", "Class", "Order",
                                 "Family", "Genus", "Species")
    taxonomic_levels <- taxonomic_levels[!is.na(taxonomic_levels)]


    # full list of IDs to retrieve the ones that have been dropped
    IDs <- df %>% dplyr::select(ID) %>% dplyr::distinct() %>% dplyr::pull()

    IDs <- data.frame(Tax_level = names(taxonomic_levels),
                      Level = tolower(substring(names(taxonomic_levels),1,1)),
                      ID = rep(IDs, each = length(taxonomic_levels))
    )


    # fill missing values with 0 or 9999 for an easier handling

    IDs_with_NA <- df %>% dplyr::filter(is.na(Bit)) %>% dplyr::select(ID) %>% dplyr::pull()
    df[is.na(df$Bit), c("Bit", "Ident", "Length")] <- 0
    df[is.na(df$E), "E"] <- 9999

    # General data preparation

    tmp <-
        df %>%
        dplyr::group_by(ID) %>%
        dplyr::arrange(ID, -Bit, -Ident, -Length) %>%
        dplyr::slice_head(n = Top_n) %>%
        dplyr::mutate(Nb_tot_hits = dplyr::n())


    # prepare data for the 3 simple methods : TopHit, TopHitPlus, TopN

    if(any(Assignation_method != c("TopNPlus"))){

        # Compute Consensus scores and Identity for best Bit score for each taxonomic
        # level

        k <- p <- c <- o <- f <- g <- s <- NA

        if(!is.na(Species)){
        s <-
            tmp %>%
            dplyr::group_by(ID, Species) %>%
            dplyr::summarise(Bit = max(Bit),
                      Length = Length[which.max(Bit)],
                      Ident = Ident[which.max(Bit)],
                      Nb_tot_hits = mean(Nb_tot_hits),
                      Cons = round(dplyr::n() *100/ Nb_tot_hits, 2),
                      .groups = "drop_last") %>%
            dplyr::mutate(Tax_level = "Species", Level = "s") %>%
            dplyr::select(Tax_level, Level, ID, Taxon = Species, Bit, Length, Ident, Cons)
        }

        if(!is.na(Genus)){
        g <-
            tmp %>%
            dplyr::group_by(ID, Genus) %>%
            dplyr::summarise(Bit = max(Bit),
                      Length = Length[which.max(Bit)],
                      Ident = Ident[which.max(Bit)],
                      Nb_tot_hits = mean(Nb_tot_hits),
                      Cons = round(dplyr::n() *100/ Nb_tot_hits, 2),
                      .groups = "drop_last") %>%
            dplyr::mutate(Tax_level = "Genus", Level = "g") %>%
            dplyr::select(Tax_level, Level, ID, Taxon = Genus, Bit, Length, Ident, Cons)
        }

        if(!is.na(Family)){
        f <-
            tmp %>%
            dplyr::group_by(ID, Family) %>%
            dplyr::summarise(Bit = max(Bit),
                      Length = Length[which.max(Bit)],
                      Ident = Ident[which.max(Bit)],
                      Nb_tot_hits = mean(Nb_tot_hits),
                      Cons = round(dplyr::n() *100/ Nb_tot_hits, 2),
                      .groups = "drop_last") %>%
            dplyr::mutate(Tax_level = "Family", Level = "f") %>%
            dplyr::select(Tax_level, Level, ID, Taxon = Family, Bit, Length, Ident, Cons)
        }

        if(!is.na(Order)){
        o <-
            tmp %>%
            dplyr::group_by(ID, Order) %>%
            dplyr::summarise(Bit = max(Bit),
                      Length = Length[which.max(Bit)],
                      Ident = Ident[which.max(Bit)],
                      Nb_tot_hits = mean(Nb_tot_hits),
                      Cons = round(dplyr::n() *100/ Nb_tot_hits, 2),
                      .groups = "drop_last") %>%
            dplyr::mutate(Tax_level = "Order", Level = "o") %>%
            dplyr::select(Tax_level, Level, ID, Taxon = Order, Bit, Length, Ident, Cons)
        }

        if(!is.na(Class)){
            c <-
                tmp %>%
                dplyr::group_by(ID, Class) %>%
                dplyr::summarise(Bit = max(Bit),
                                 Length = Length[which.max(Bit)],
                                 Ident = Ident[which.max(Bit)],
                                 Nb_tot_hits = mean(Nb_tot_hits),
                                 Cons = round(dplyr::n() *100/ Nb_tot_hits, 2),
                                 .groups = "drop_last") %>%
                dplyr::mutate(Tax_level = "Class", Level = "c") %>%
                dplyr::select(Tax_level, Level, ID, Taxon = Class, Bit, Length, Ident, Cons)
        }

        if(!is.na(Phylum)){
            p <-
                tmp %>%
                dplyr::group_by(ID, Phylum) %>%
                dplyr::summarise(Bit = max(Bit),
                                 Length = Length[which.max(Bit)],
                                 Ident = Ident[which.max(Bit)],
                                 Nb_tot_hits = mean(Nb_tot_hits),
                                 Cons = round(dplyr::n() *100/ Nb_tot_hits, 2),
                                 .groups = "drop_last") %>%
                dplyr::mutate(Tax_level = "Phylum", Level = "p") %>%
                dplyr::select(Tax_level, Level, ID, Taxon = Phylum, Bit, Length, Ident, Cons)
        }

        if(!is.na(Kingdom)){
            k <-
                tmp %>%
                dplyr::group_by(ID, Kingdom) %>%
                dplyr::summarise(Bit = max(Bit),
                                 Length = Length[which.max(Bit)],
                                 Ident = Ident[which.max(Bit)],
                                 Nb_tot_hits = mean(Nb_tot_hits),
                                 Cons = round(dplyr::n() *100/ Nb_tot_hits, 2),
                                 .groups = "drop_last") %>%
                dplyr::mutate(Tax_level = "Kingdom", Level = "k") %>%
                dplyr::select(Tax_level, Level, ID, Taxon = Kingdom, Bit, Length, Ident, Cons)
        }

        long <- list(k, p, c, o, f, g, s)
        long <- dplyr::bind_rows(long[!is.na(long)]) %>% dplyr::ungroup()

    }

    # Assign final taxonomy with each method

    # just the Hit with best Bit score

    if(any(Assignation_method == c("TopHit"))){

        TopHit <-
            long %>%
            dplyr::group_by(Tax_level, Level, ID) %>%
            dplyr::slice(which.max(Bit)) %>%
            dplyr::mutate(Method = "TopHit") %>%
            dplyr::relocate(Method, .before = Tax_level)


    }



    # Hit with best Bit score and in case of ties choose the one with
    # the best Consensus score

    if(any(Assignation_method == c("TopHitPlus"))){

        TopHitPlus <-
            long %>%
            dplyr::arrange(dplyr::desc(Level), ID, -Bit, -Cons) %>%
            dplyr::group_by(Tax_level, Level,ID) %>%
            dplyr::slice(1) %>%
            dplyr::mutate(Method = "TopHitPlus") %>%
            dplyr::relocate(Method, .before = Tax_level)

    }


    # Hit with the best consensus score, in case of ties choose the one with the best
    # Bit score
    if(any(Assignation_method == c("TopN"))){

        TopN <-
            long %>%
            dplyr::arrange(dplyr::desc(Level), ID, -Cons, -Bit) %>%
            dplyr::group_by(Tax_level, Level,ID) %>%
            dplyr::slice(1) %>%
            dplyr::mutate(Method = "TopN") %>%
            dplyr::relocate(Method, .before = Tax_level)
    }






    # Prepare the data and compute consensus for each taxonomic level for TopNPlus
    # method
    # We need to recompute everything because this method implies prefiltering
    # of some of the hits.

    if(any(Assignation_method == c("TopNPlus"))){


        k <- p <- c <- o <- f <- g <- s <- NA

        if(!is.na(Species)){
        s <-
            tmp %>%
            dplyr::group_by(ID) %>%
            dplyr::filter((max(Ident)-Ident) <= Ident_max_diff_s) %>%
            dplyr::filter(Ident >= Ident_min_s) %>%
            dplyr::filter(Length > Length_min & E < E_max) %>%
            dplyr::mutate(Nb_tot_hits = dplyr::n()) %>%
            dplyr::group_by(ID, Species) %>%
            dplyr::summarise(Bit = max(Bit),
                      Length = Length[which.max(Bit)],
                      Ident = Ident[which.max(Bit)],
                      Nb_tot_hits = mean(Nb_tot_hits),
                      Cons = round(dplyr::n() *100/ Nb_tot_hits, 2),
                      .groups = "drop_last") %>%
            dplyr::mutate(Tax_level = "Species", Level = "s") %>%
            dplyr::select(Tax_level, Level, ID, Taxon = Species, Bit, Length, Ident, Cons)
        }


        if(!is.na(Genus)){
        g <-
            tmp %>%
            dplyr::group_by(ID) %>%
            dplyr::filter((max(Ident)-Ident) <= Ident_max_diff_g) %>%
            dplyr::filter(Ident >= Ident_min_g) %>%
            dplyr::filter(Length > Length_min & E < E_max) %>%
            dplyr::mutate(Nb_tot_hits = dplyr::n()) %>%
            dplyr:: group_by(ID, Genus) %>%
            dplyr::summarise(Bit = max(Bit),
                      Length = Length[which.max(Bit)],
                      Ident = Ident[which.max(Bit)],
                      Nb_tot_hits = mean(Nb_tot_hits),
                      Cons = round(dplyr::n() *100/ Nb_tot_hits, 2),
                      .groups = "drop_last") %>%
            dplyr::mutate(Tax_level = "Genus", Level = "g") %>%
            dplyr::select(Tax_level, Level, ID, Taxon = Genus, Bit, Length, Ident, Cons)
        }

        if(!is.na(Family)){
        f <-
            tmp %>%
            dplyr::group_by(ID) %>%
            dplyr::filter((max(Ident)-Ident) <= Ident_max_diff_f) %>%
            dplyr::filter(Ident >= Ident_min_f) %>%
            dplyr::filter(Length > Length_min & E < E_max) %>%
            dplyr::mutate(Nb_tot_hits = dplyr::n()) %>%
            dplyr::group_by(ID, Family) %>%
            dplyr::summarise(Bit = max(Bit),
                      Length = Length[which.max(Bit)],
                      Ident = Ident[which.max(Bit)],
                      Nb_tot_hits = mean(Nb_tot_hits),
                      Cons = round(dplyr::n() *100/ Nb_tot_hits, 2),
                      .groups = "drop_last") %>%
            dplyr::mutate(Tax_level = "Family", Level = "f") %>%
            dplyr::select(Tax_level, Level, ID, Taxon = Family, Bit, Length, Ident, Cons)
        }

        if(!is.na(Order)){
        o <-
            tmp %>%
            dplyr::group_by(ID) %>%
            dplyr::filter((max(Ident)-Ident) <= Ident_max_diff_o) %>%
            dplyr::filter(Ident >= Ident_min_o) %>%
            dplyr::filter(Length > Length_min & E < E_max) %>%
            dplyr::mutate(Nb_tot_hits = dplyr::n()) %>%
            dplyr::group_by(ID, Order) %>%
            dplyr::summarise(Bit = max(Bit),
                      Length = Length[which.max(Bit)],
                      Ident = Ident[which.max(Bit)],
                      Nb_tot_hits = mean(Nb_tot_hits),
                      Cons = round(dplyr::n() *100/ Nb_tot_hits, 2),
                      .groups = "drop_last") %>%
            dplyr::mutate(Tax_level = "Order", Level = "o") %>%
            dplyr::select(Tax_level, Level, ID, Taxon = Order, Bit, Length, Ident, Cons)
        }

        if(!is.na(Class)){
            c <-
                tmp %>%
                dplyr::group_by(ID) %>%
                dplyr::filter((max(Ident)-Ident) <= Ident_max_diff_o) %>%
                dplyr::filter(Ident >= Ident_min_o) %>%
                dplyr::filter(Length > Length_min & E < E_max) %>%
                dplyr::mutate(Nb_tot_hits = dplyr::n()) %>%
                dplyr::group_by(ID, Class) %>%
                dplyr::summarise(Bit = max(Bit),
                                 Length = Length[which.max(Bit)],
                                 Ident = Ident[which.max(Bit)],
                                 Nb_tot_hits = mean(Nb_tot_hits),
                                 Cons = round(dplyr::n() *100/ Nb_tot_hits, 2),
                                 .groups = "drop_last") %>%
                dplyr::mutate(Tax_level = "Class", Level = "c") %>%
                dplyr::select(Tax_level, Level, ID, Taxon = Class, Bit, Length, Ident, Cons)
        }

        if(!is.na(Phylum)){
            p <-
                tmp %>%
                dplyr::group_by(ID) %>%
                dplyr::filter((max(Ident)-Ident) <= Ident_max_diff_o) %>%
                dplyr::filter(Ident >= Ident_min_o) %>%
                dplyr::filter(Length > Length_min & E < E_max) %>%
                dplyr::mutate(Nb_tot_hits = dplyr::n()) %>%
                dplyr::group_by(ID, Phylum) %>%
                dplyr::summarise(Bit = max(Bit),
                                 Length = Length[which.max(Bit)],
                                 Ident = Ident[which.max(Bit)],
                                 Nb_tot_hits = mean(Nb_tot_hits),
                                 Cons = round(dplyr::n() *100/ Nb_tot_hits, 2),
                                 .groups = "drop_last") %>%
                dplyr::mutate(Tax_level = "Phylum", Level = "p") %>%
                dplyr::select(Tax_level, Level, ID, Taxon = Phylum, Bit, Length, Ident, Cons)
        }

        if(!is.na(Kingdom)){
            k <-
                tmp %>%
                dplyr::group_by(ID) %>%
                dplyr::filter((max(Ident)-Ident) <= Ident_max_diff_o) %>%
                dplyr::filter(Ident >= Ident_min_o) %>%
                dplyr::filter(Length > Length_min & E < E_max) %>%
                dplyr::mutate(Nb_tot_hits = dplyr::n()) %>%
                dplyr::group_by(ID, Kingdom) %>%
                dplyr::summarise(Bit = max(Bit),
                                 Length = Length[which.max(Bit)],
                                 Ident = Ident[which.max(Bit)],
                                 Nb_tot_hits = mean(Nb_tot_hits),
                                 Cons = round(dplyr::n() *100/ Nb_tot_hits, 2),
                                 .groups = "drop_last") %>%
                dplyr::mutate(Tax_level = "Kingdom", Level = "k") %>%
                dplyr::select(Tax_level, Level, ID, Taxon = Kingdom, Bit, Length, Ident, Cons)
        }

        long <- list(k, p, c, o, f, g, s)
        long <- dplyr::bind_rows(long[!is.na(long)]) %>% dplyr::ungroup()


        TopNPlus <-
            long %>%
            dplyr::arrange(dplyr::desc(Level), ID, -Cons, -Bit) %>%
            dplyr::group_by(Tax_level, Level,ID) %>%
            dplyr::slice(1) %>%
            dplyr::mutate(Method = "TopNPlus") %>%
            dplyr::relocate(Method, .before = Tax_level)



        # NAs for taxa that did not match the threshold

        # retrieve IDs that have been droped either because they didn't match certain
        # thresholds in TopNPlus method or because they had missing values
        # in Bit, Ident,...(no blast hit)

        TopNPlus<-
            TopNPlus %>%
            dplyr::full_join(IDs, by = c("Tax_level", "Level", "ID")) %>%
            dplyr::mutate(Method = replace(Method, which(is.na(Method)), "TopNPlus"))

    }

    long <- list(TopHit, TopHitPlus, TopN, TopNPlus)
    long <- dplyr::bind_rows(long[!is.na(long)]) %>% dplyr::ungroup()


    # Homogenise the missing values
    long[which(long$Bit == 0), c("Bit", "Length", "Ident", "Cons", "Taxon")] <- NA
    long[is.na(long$Bit), c("Bit", "Length", "Ident", "Cons", "Taxon")] <- NA
    long[which(long$Taxon == ""), "Taxon"] <- NA



    # Merge with TRUE taxonomies

    if(! is.null(taxo)){

        # Read the file, add column names and clean the species
        if(is.character(taxo)){
            reftaxo <- unique(data.table::fread(taxo, sep = "\t", header = FALSE))
        } else {
            reftaxo <- taxo
        }
        colnames(reftaxo) <- c("Feature.ID", "Taxon")
        reftaxo <- cbind(reftaxo, split_taxonomy(reftaxo$Taxon, clean = TRUE))
        reftaxo$Taxon <- NULL

        # Transform into long format
        reftaxo <-
            reftaxo %>%
            dplyr::select(Feature.ID:Species) %>%
            tidyr::pivot_longer(cols = -Feature.ID, names_to = "Tax_level",
                         values_to = "Taxon_true") %>%
            dplyr::distinct()

        # merge with the taxonomic assignments results
        long <-
        long %>%
            dplyr::left_join(reftaxo,  by = c("Tax_level" = "Tax_level",
                                       "ID" = "Feature.ID"))


    }

    long %>%
        dplyr::mutate(Tax_level = factor(Tax_level, levels = taxonomic_levels)) %>%
        dplyr::mutate(Level = factor(Level, levels = tolower(substr(levels(Tax_level), 1, 1)))) %>%
        dplyr::arrange(Method, Tax_level) %>%
        return()

}






# /*
# ----------------------------- pivot_assign_taxonomy ---------------------------
# */






#' Pivot outputs of `assign_taxonomy` from a long to a wide format
#'
#' This function will simply reorganize the output of the function [`assign_taxonomy()`]
#' in order to obtain one row for each `Method` and each query `ID` with the
#' results of each taxonomic level side by side instead of being on top of each other.
#' For example the species level output will be placed in columns
#' `s_Taxon_true`, `s_Taxon`, `s_Ident`,  `s_Cons`, `s_Length`, `s_Bit`.
#' The other taxonomic levels results will be placed in columns with similar names
#' with a different prefix : `g_` for genus, `f_` for family, `o_` for order,
#' `c_` for class, `p_` for phylum and `k_` for kingdom.
#'
#' @param assigned A data.frame produced by the function [`assign_taxonomy()`]
#'
#' @return A data.frame with the same data as the input but reorganized in a wide format
#' @export
#'
#' @examples
#'
#' # Fake top 15 blast hits for a sequence corresponding to Scirpoides holoschoenus
#' # and inspired by a real life example where the 4 assignment methods give
#' # different results
#'
#' example_top15_blast_hits <-
#'     data.frame(
#'         TaxID_query = "ID001",
#'         Species_true = "Scirpoides holoschoenus",
#'
#'         Bit_score = c(1075,704,682,651,630,630,549,549,529,525,525,520,520,518,512),
#'         Length = c(589,607,608,595,597,597,568,568,570,541,571,616,609,569,569),
#'         Identity = c(99.660,88.138,87.500,86.891,86.265,86.265,84.683,84.507,83.860,
#'                      84.843,83.713,82.630,82.594,83.480,83.304),
#'         E_value = 10^-14,
#'
#'         Phylum = "Streptophyta",
#'         Class = "Magnoliopsida",
#'         Order = "Poales",
#'         Family = "Cyperaceae",
#'         Genus = c("Scirpoides","Erioscirpus","Erioscirpus","Dracoscirpoides",
#'                   "Dracoscirpoides","Dracoscirpoides", "Cyperus","Cyperus","Cyperus",
#'                   "Cyperus","Cyperus","Cyperus","Cyperus","Cyperus","Cyperus"),
#'         Species = c("Scirpoides holoschoenus","Erioscirpus microstachyus",
#'                     "Erioscirpus comosus", "Dracoscirpoides ficinioides",
#'                     "Dracoscirpoides falsa", "Dracoscirpoides falsa",
#'                     "Cyperus pulchellus", "Cyperus leucocephalus"," Cyperus haspan",
#'                     "Cyperus spiralis", "Cyperus flaccidus", "Cyperus sp.",
#'                     "Cyperus diffusus", "Cyperus tenuispica", "Cyperus haspan")
#'     )
#'
#' # corresponding true taxonomy
#' example_taxonomy <-
#' data.frame(ID = "ID001",
#' Taxo = "k__Viridiplantae; p__Streptophyta; c__Magnoliopsida; o__Poales; f__Cyperaceae; g__Scirpoides; s__holoschoenus")
#'
#'
#' # Assign the taxonomy with mostly default options and in "long" format
#'
#' tax_pred <-
#'     assign_taxonomy(df = example_top15_blast_hits,
#'                     taxo = example_taxonomy,
#'                     Order = NA, Family = NA) # ignore the family and genus levels
#' tax_pred
#'
#'
#' # Present the same data in a different, wider, format
#'
#' pivot_assign_taxonomy(tax_pred)
#'

#'
pivot_assign_taxonomy <- function(assigned){

    # add empty columns if ever any of these columns are missing from the output
    # these columns will be removed at the end
    if(all(colnames(assigned) != "Taxon_true")) {assigned$Taxon_true <- NA}
    if(all(colnames(assigned) != "Taxon")) {assigned$Taxon_true <- NA}
    if(all(colnames(assigned) != "Bit")) {assigned$Bit <- NA}
    if(all(colnames(assigned) != "Length")) {assigned$Length <- NA}
    if(all(colnames(assigned) != "Ident")) {assigned$Ident <- NA}
    if(all(colnames(assigned) != "Cons")) {assigned$Cons <- NA}

    # removing levels will cause unusable results
    if(all(colnames(assigned) != "Level")) {assigned$Level <- NA}

    # Removing ID or Method might result into column lists
    # (more than one possible value per row)
    if(all(colnames(assigned) != "Method")) {assigned$Method <- ""}
    if(all(colnames(assigned) != "ID")) {assigned$ID <- ""}

    assigned %>%
        dplyr::select(-Tax_level) %>%
        tidyr::pivot_wider(
            id_cols = c(Method, ID),
            names_from = Level,
            values_from = c("Taxon_true", "Taxon", "Bit",
                            "Length", "Ident", "Cons"),
            names_glue = "{Level}_{.value}",
            names_sort = TRUE) %>%
        dplyr::select(dplyr::any_of(
            c("Method", "ID",
              # "k_Taxon_true", "p_Taxon_true", "c_Taxon_true",
              # "o_Taxon_true", "f_Taxon_true", "g_Taxon_true", "s_Taxon_true",
              "k_Taxon_true", "k_Taxon", "k_Bit", "k_Length", "k_Ident", "k_Cons",
              "p_Taxon_true", "p_Taxon", "p_Bit", "p_Length", "p_Ident", "p_Cons",
              "c_Taxon_true", "c_Taxon", "c_Bit", "c_Length", "c_Ident", "c_Cons",
              "o_Taxon_true", "o_Taxon", "o_Bit", "o_Length", "o_Ident", "o_Cons",
              "f_Taxon_true", "f_Taxon", "f_Bit", "f_Length", "f_Ident", "f_Cons",
              "g_Taxon_true", "g_Taxon", "g_Bit", "g_Length", "g_Ident", "g_Cons",
              "s_Taxon_true", "s_Taxon", "s_Bit", "s_Length", "s_Ident", "s_Cons"
            ))) %>%
        # Remove empty columns
        dplyr::select(tidyselect::where(function(x) any(!is.na(x)))) %>%
        dplyr::arrange(ID, Method) %>%
        return()
}









