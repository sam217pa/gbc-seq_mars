#' La fonction suivante permet donc de ne garder que les séquences dont la
#' qualité médiane est supérieure à 40.
keep_clean_only <- function(data, qual_ = 40) {

    get_dirty_seq <- function(data)
        group_by(data, name) %>%
            summarise(median = median(qual)) %>%
            filter(median < qual_) %>%
            select(name) %>% unlist()

    data %>%
        filter(!(name %in% get_dirty_seq(data)))
}

#' Détermine la polarité de la transition
#' @param reference: la base de référence
#' @param experimental: la base séquencée
find_polarity <- function(reference, experimental) {
    is_w <- function(base) { ifelse(base == "A" | base == "T", TRUE, FALSE)}
    is_s <- function(base) { ifelse(base == "C" | base == "G", TRUE, FALSE)}

    if      (is_w(reference) & is_w(experimental)) "ww"
    else if (is_w(reference) & is_s(experimental)) "ws"
    else if (is_s(reference) & is_w(experimental)) "sw"
    else if (is_s(reference) & is_s(experimental)) "ss"
    else ""
}

#' La fonction suivante renvoit un `tbl_df` avec le nom, le type de mutant, le
#' genotype.
#'
#' @param data: le jeu de donnée snp,
#' @param mutant: crée un filtre par mutant. Selon les mutants, on ne recherche
#'   pas les mêmes séquences.
make_genotype <- function(data, mutant_ = c("ws", "sw", "s", "w"), qual_ = 30) {
    weak_or_strong <- function(base) {
        is_w <- function(base) { ifelse(base == "A" | base == "T", TRUE, FALSE)}
        is_s <- function(base) { ifelse(base == "C" | base == "G", TRUE, FALSE)}

        if (is_w(base)) "W"
        else if (is_s(base)) "S"
        else if (base == "N") "N"
        else stop(base, " is not a DNA base")
    }


    data %>%
        keep_clean_only() %>%
        filter(cons == "x" | cons == "X") %>% # ne garde que les positions de SNP
        filter(mutant %in% mutant_) %>%       # filtre selon le mutant
        rowwise() %>%
        mutate(seqb = ifelse(qual < qual_, "N", seqb )) %>% # convertit les bases d'une qualité inférieure à qual_ en N
        mutate(geno = weak_or_strong(seqb)) %>% # détermine le génotype de la base séquencée
        ungroup() %>% group_by(name, mutant) %>%
        summarise(geno = toString(geno) %>% gsub(", ", "", .)) # pretty print
}

#' La fonction suivante permet de compter le nombre d'anomalies, autrement dit
#' de SNP inattendus dans les traces de conversion. Renvoit une `tbl_df`
#'
#' @param data : le jeu de donné snp
#' @param kmer : le kmer à recenser dans le génotype.
count_kmer <- function(data, kmer) {

    ## closure around find_xxx, which count the number of occurence of the kmer in
    ## the string.
    count_anomaly <- function(string, kmer_) {
        find_xxx <- function(string) {
            function(kmer_) {
                XxX <- gregexpr(kmer_, string) %>% unlist() # compte le nombre d'occurences du kmer
                if   (XxX[1] > 0) length(XxX)
                ## if gregexpr find something, give its length (as integer)
                else 0L # 0 otherwise.
            }
        }
        find_xxx(string)(kmer_) # count anomaly uses find_xxx to count kmer in string.
    }

    ## applique la fonction count_anomaly par ligne, et ne conserve que les
    ## séquences qui montrent une ou plus d'anomalies.
    data %>%
        rowwise() %>%
        mutate(anom = count_anomaly(geno, kmer)) %>%
        filter(anom > 0)
}

#' Compte les doublets des séquences propres.
#' Un wrapper autour de make_genotype et count_kmer pour l'analyse interactive.
count_doublet <- function(data, mut, kmer) {
    data %>%
        keep_clean_only() %>%
        make_genotype(mut) %>%
        count_kmer(kmer) %>% {
            if (nrow( . ) > 0) select(., anom) %>% sum()
            else 0
        }
}

get_marker_only <- function(data, snp_only = TRUE) {
    ## permet de ne conserver que les marqueurs dans le jeu de donnée
    snp_without_N <-  data %>% filter(cons != "N", cons != "-", name != "psw21")

    if (snp_only) filter(snp_without_N, cons == "X" | cons == "x")
    else snp_without_N
}

#' Attribue les bons niveaux de facteurs aux donneurs.
set_mutant_factor <- function(data) {
    mutate(data, mutant = factor(mutant, levels = c("s", "w", "sw", "ws")))
}

#' attribue les bons niveaux de facteurs aux différentes séquences. Choisies de
#' façon à ce que la base expérimentale apparaisse en rouge.
set_seq_factor <- function(data) {
    mutate(data,
           seq = factor(seq, levels = c("seqb", "refb", "snpb"),
                        labels = c("Expérimental", "Référence", "Donneur")))
}

#' attribue les bons niveaux de facteurs aux différentes bases, de façon à ce
#' que S et W clusterent.
set_base_factor <- function(data) {
    mutate(data, base = factor(base, levels = c("G", "C", "T", "A")))
}

get_gc_content <- function(data, mut=NULL, tidy = TRUE) {
    clean_seq <- function(x) toString(x) %>% gsub(", ", "", . )

    get_gc_rate <- function(dna_seq) {

    }
    if (is.null(mut)) mut <- c("ws", "sw", "w", "s")

    gc_content <- data %>% filter(mutant %in% mut) %>% group_by(name, mutant) %>%
        summarise(exp = seqinr::GC(seqb),
                  snp = seqinr::GC(snpb),
                  ref = seqinr::GC(refb))
    if (tidy) gc_content %>% tidyr::gather("type", "GC", 3:5)
    else gc_content
}
