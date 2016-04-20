#' ---
#' title: Phruscling SNP et joyeusetés associées…
#' author: Samuel BARRETO
#' date:
#' output:
#'   tufte::tufte_html:
#'     highlight: tango
#'     toc: true
#'     toc_depth: 2
#' ---

#' # Introduction
#'
#' Cette page correspond à l'analyse des résultats du SNP calling effectué avec
#' mon petit programme `phruscle` ("phred + muscle").
#'
#' Pour rappel, ce programme utilisait phred sur les fichiers de séquençage pour
#' déterminer les bases. Les séquences obtenues étaient ensuite alignées à un
#' alignement de référence composé de la séquence sauvage et de la séquence du
#' gène synthétique.
#'
#' J'ai ensuite analysé les sorties de phred qui recense la qualité des bases
#' appelées. J'ai pu ainsi attribuer un score de qualité à toutes les bases non
#' chimériques, _ie_ avec une correspondance dans l'alignement de référence
#' **et** dans le fichier de séquençage.

#+ setup, include=FALSE
library(knitr)
library(tufte)
opts_chunk$set(cache = FALSE, dev = 'png', include = TRUE,
               echo = TRUE, warning = FALSE, error = FALSE,
               message = FALSE)


#+ setup2, include = TRUE
library(dplyr)
library(ggplot2)
library(readr)
library(viridis)
library(cowplot)
library(extrafont)
#
fte_theme <- function() {
    ## inspiré de http://minimaxir.com/2015/02/ggplot-tutorial/
    ## Generate the colors for the chart procedurally with RColorBrewer
    palette <- RColorBrewer::brewer.pal("Greys", n=9)
    color.background = palette[2]
    color.grid.major = palette[3]
    color.axis.text = palette[6]
    color.axis.title = palette[7]
    color.title = palette[9]
    ## Begin construction of chart
    theme_bw(base_size=9) +
        ## Set the entire chart region to a light gray color
        theme(
            text = element_text(family = "Ubuntu Light"),
            panel.background=element_rect(
                fill=alpha(color.background, 1), color=color.background),
            plot.background=element_rect(
                fill=alpha(color.background, 0.5), color=color.background),
            panel.border=element_rect(color=color.background),
            panel.grid.major=element_line(
                color=color.grid.major,size=.25, linetype = "dotted"),
            panel.grid.minor=element_blank(),
            axis.ticks=element_blank(),
            strip.text.y = element_text(angle = 0, size = 8),
            legend.position="none", # Format the legend, but hide by default
            legend.background = element_blank(),
            legend.text = element_text(size=9,color=color.axis.title),
            legend.key = element_rect(
                fill = scales::alpha("gray", 0),
                colour = scales::alpha("gray", 0)),
            ## Set title and axis labels, and format these and tick marks
            plot.title = element_text(
                color = color.title, size = 14, face = "bold", vjust = 1.25,
                hjust = 0, family = "Ubuntu"),
            axis.text.x = element_text(size = 7,color = color.axis.text),
            axis.text.y = element_text(size = 7,color = color.axis.text),
            axis.title.x = element_text(
                size = 8, color = color.axis.title, vjust = 0, hjust = 0.9),
            axis.title.y = element_text(
                size = 8, color = color.axis.title, vjust = 0.9, angle = 0 ),
            ## Plot margins
            plot.margin = unit(c(0.35, 0.2, 0.3, 0.35), "cm")
        )
}
#
theme_set(theme_bw() + fte_theme())
#
legend_position <- function(x=NULL, y=NULL) {
    custom_legend <- theme(
        legend.margin = unit(-0.3,"lines"),
        legend.key = element_rect(fill = scales::alpha("gray", 0),
                                  colour = scales::alpha("gray", 0)))
    if (is.null(x) & is.null(y)) {
        custom_legend + theme(legend.position = "bottom")
    } else {
        custom_legend + theme(legend.position = c(x, y))
    }
}
#
data_location <- "../../data/phruscle_snpcall.csv"
#
snp <- read_csv(data_location, col_types = "ccciciccid", trim_ws = TRUE) %>%
    mutate(
        name = gsub("-1073.+$", "", name),
        base = toupper(base)#,
    )
#
find_mutant <- function(name) {
    if      (grepl("ws", name)) "ws"
    else if (grepl("sw", name)) "sw"
    else if (grepl("W", name )) "w"
    else if (grepl("S", name )) "s"
    else ""
}
#
# neat little trick to reduce time of rowwise application of find_mutant.
snp <- snp %>%
    group_by(name) %>%
    summarise(count = n()) %>%
    rowwise() %>%
    mutate(mutant=find_mutant(name)) %>%
    inner_join(snp, .) %>%
    filter(mutant != "")

#' J'ai donc obtenu un tableau de donnée de la forme suivante, où :
#'
#' - `cons`: un code maison pour déterminer le _consensus_. `.` si les trois
#'   bases sont identiques, `x` si la base séquencée est la base introduite, `X`
#'   si la base séquencée est la base sauvage, `N` si c'est encore autre chose.
#' - `name`: le nom expérimental de la séquence
#' - `refb`: la base de référence
#' - `refp`: la position sur la référence.
#' - `seqb`: le reverse complement de la base séquencée.
#' - `seqp`: la position sur la séquence expérimentale.
#' - `snpb`: la base sur le gène synthétique introduit.
#' - `base`: la base séquenceé.
#' - `qual`: sa qualité.
#' - `phase`: la position sur le spectrogramme.
#' - `mutant`: ordinal, la manip.
#'
#' Vincent appelle ça joliment un alignement pairwise en colonne.

head(snp)

#' # Observations générales
#'
#' ## Nombres de séquences ?
n_distinct(snp$name)

#' On a donc bien toutes nos séquences.

snp %>%
    group_by(mutant) %>%
    summarise(count = n()) %>%
    knitr::kable(align = "c")

snp %>%
    group_by(mutant, name) %>%
    summarise(count = n()) %>%
    ggplot(data = ., aes(x = count )) +
    geom_histogram(binwidth = 1) +
    labs(x = "Nombre de positions", y = "",
         title = "Distribution du nombre\nde positions par séquence")

#' J'ai voulu déterminer s'il existe des bases dans les sorties de phruscle qui
#' ne correspondent pas aux sorties de phd. Autrement dit si la base `base` est
#' bien toujours la même que la base `seqb`.

FALSE %in% (snp$base == snp$seqb)

## filter(snp, cons == "X")

#' Plutôt bon signe ! Phruscle a fonctionné comme il faut, _ie_ toutes les bases
#' qu'on a sorties dans nos alignements ont une correspondance dans le fichier
#' phred !.
#'

#' # Analyse de la qualité des séquences

#' J'ai regardé globalement la distribution de la qualité des bases :

#+ qual0, fig.margin=TRUE
snp %>%
    ggplot(data = ., aes(x = qual )) +
    geom_histogram(binwidth = 1)

#' J'ai regardé séquence par séquence les variation de qualité le long de la
#' séquence :

#+ qual1, fig.fullwidth=TRUE, fig.width=15
snp %>%
    ggplot(data = ., aes(x = refp, y = qual, color = qual )) +
    geom_point(alpha = 0.1, size = 0.1) +
    geom_line(aes(group = name), alpha = 1/20, size = 0.1 ) +
    scale_color_viridis()

#+ qual2, fig.margin=TRUE, fig.cap="Position des bases de qualité < 30"
snp %>%
    filter(qual < 30) %>%
    ggplot(data = ., aes(x = refp, y = qual )) +
    geom_point(alpha = 0.1, size = 0.1)

#' Certaines séquences ne sont pas des plus propres. On peut certainement les
#' virer pour ne garder que les bonnes séquences. Comme pour le cas de la
#' séquence psw21.
#'
#' J'ai choisi pour statistique la qualité *médiane*, la qualité *moyenne* étant
#' trop influencée par les baisses de qualité dans les positions terminales. De
#' plus, la médiane a une distribution plus centrée que la moyenne.

#+ median_qual, fig.margin = TRUE
snp %>%
    group_by(name) %>%
    summarise(median = median(qual)) %>%
    qplot(data = . , median, binwidth = 0.1) +
    geom_vline(xintercept = 40, color = "red")

#+ mean_qual, fig.margin = TRUE
snp %>%
    group_by(name) %>%
    summarise(mean = mean(qual)) %>%
    qplot(data = . , mean, binwidth = 0.1) +
    geom_vline(xintercept = 40, color = "red")

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

#' Comme attendu, la qualité est moindre en fin de run, et plutôt bonne au
#' début^[On séquence de gauche à droite, la cassette de résistance est en 3'].
#' Ce sont des séquences trimmées par phd. On n'est pas obligé de trimmer, mais
#' ça permet de ne garder que les séquences qui vallent le coup. Peut-être
#' refaire les analyses avec les séquences non trimmées.

#' # Distribution de la longueur des séquences.

#' Comme attendu, la longueur des séquences est plus élevée avec les
#' constructions SW et WS, puisqu'on a placé l'amorce légérèment en amont de la
#' construction.
snp %>%
    group_by(name, mutant) %>%
    summarise(len = max(refp) - min(refp)) %>%
    ggplot(data = ., aes(x = len )) +
    geom_histogram(binwidth = 1) +
    facet_grid(mutant ~ .)

#+ snp_distrib0, fig.margin=TRUE
snp %>%
    filter(cons == "x" | cons == "X") %>%
    ## filter(refp < 750 & refp > 40) %>%
    ggplot(aes(x = refp)) +
    geom_histogram(binwidth = 1)


find_polarity <- function(reference, experimental) {
    is_w <- function(base) { ifelse(base == "A" | base == "T", TRUE, FALSE)}
    is_s <- function(base) { ifelse(base == "C" | base == "G", TRUE, FALSE)}

    if      (is_w(reference) & is_w(experimental)) "ww"
    else if (is_w(reference) & is_s(experimental)) "ws"
    else if (is_s(reference) & is_w(experimental)) "sw"
    else if (is_s(reference) & is_s(experimental)) "ss"
    else ""
}

#' La fonction suivante permet de créer les alignement via ggplot2. Les clones
#' sont triés par longueur de trace de conversion par la fonction embarquée
#' sort_by_tract_length.
plot_align <- function(data, mutant_) {

    sort_by_tract_length <- function(mut) {
        data %>%
            filter(mutant %in% mut) %>%
            group_by(name, mutant) %>%
            filter(cons == "x") %>%           # dernière position correspondant au donneur
            summarise(len = min(refp)) %>%
            ungroup() %>%
            arrange(len) %>%
            {.$name}
    }


    is_inside_conv <- function(data, mut) {
        data %>%
            filter(mutant %in% mut, cons == "x") %>%
            group_by(name) %>%
            summarise(bascul = min(refp)) %>%
            inner_join(data, ., by = "name") %>%
            keep_clean_only() %>%
            ## garde les positions après le point de switch qui correspondent au génotype parental.
            filter(refp > bascul, cons == "X" | cons == "x") %>%
            mutate(tract = TRUE) %>%
            select(name, bascul, refp, tract) %>%
            left_join(data, .)
    }

    is_restoration <- function(data, mut) {
        ## détermine si le SNP est une restauration du génotype parental. Autrement
        ## dit les SNP qui sont dans la trace de conversion mais qui sont X, *ie*
        ## des SNPs parentaux.
        data %>%
            filter(mutant %in% mut, cons == "X", tract == TRUE)
    }

    ## default color
    red    <- RColorBrewer::brewer.pal(n = 4, "Set1")[1]
    blue   <- RColorBrewer::brewer.pal(n = 4, "Set1")[2]
    green  <- RColorBrewer::brewer.pal(n = 4, "Set1")[3]
    violet <- RColorBrewer::brewer.pal(n = 4, "Set1")[4]

    data <-
        is_inside_conv(data, mutant_) %>%
        group_by(name, mutant) %>%
        filter(mutant %in% mutant_, cons == "x" | cons == "X" ) %>%
        rowwise() %>%
        mutate(sens = find_polarity(refb, seqb)) %>%
        ungroup() %>%
        keep_clean_only()

    data$tract[is.na(data$tract)] <- FALSE

    ggplot(data = data,
           aes(x = factor(refp),
               y = factor(name, levels = sort_by_tract_length(mut=mutant_)))) +
        geom_point(aes(color = tract, size = qual, alpha=qual)) +
        ## représente les cas complexes
        geom_point(data = is_restoration(data, mutant_),
                   aes(alpha = qual), color = green, size = 5) +
        ## représente la séquence du donneur
        geom_text(aes(label = snpb, x = factor(refp), y = -3), color = blue,
                  vjust = -0.5, size = 2, family = "Ubuntu Light") +
        annotate("text", x = 1.2, y = -0.8, label = "Donneur : ", color = blue, size = 2 ) +
        ## représente la séquence du receveur.
        geom_text(aes(label = refb, x = factor(refp),
                      y = length(sort_by_tract_length(mutant_)) + 10),
                  color = red, vjust = 9, size = 2, family = "Ubuntu Light") +
        annotate("text", x = 1.2, y = length(sort_by_tract_length(mutant_)) + 5,
                 label = "Receveur : ", color = red, size = 2 ) +
        scale_color_brewer(palette = "Set1", guide = "none") +
        scale_alpha( range=c(1/5, 0.8), guide=FALSE ) +
        scale_size(range = c(0.5, 2), breaks = c(10, 50),
                   labels = c("Faible", "Forte")) +
        labs(x = "", y = "", size = "Qualité",
             title = paste("Alignement pour la manip", toupper(mutant_))) +
        theme(legend.direction = "horizontal",
              ## legend.position  = c(0.8, -0.2),
              legend.margin = unit(0,"lines"),
              ## legend.key = element_rect(fill = scales::alpha("gray", 0),
              ##                           colour = scales::alpha("gray", 0)),
              panel.grid.major.y = element_line(size = 0.1, linetype = "dotted")) +
        legend_position(0.8, 0.98)
}


plot_qual <- function(data, mutant_) {
    red    <- RColorBrewer::brewer.pal(n = 4, "Set1")[1]
    blue   <- RColorBrewer::brewer.pal(n = 4, "Set1")[2]
    data %>%
        filter(mutant == mutant_) %>%
        ggplot(aes(x = refp, y = qual)) +
        geom_point(aes(color = qual), size = 1/10, alpha = 1/20) +
        geom_smooth(se = FALSE, color = red, linetype = "dotted") +
        scale_color_continuous(high = blue, low = red) +
        labs(x = "Position sur la référence", y = "")
}

plot_align_qual <- function(data, mutant_) {
    plot_grid(
        data %>% plot_align(mutant_),
        data %>% plot_qual(mutant_),
        ncol = 1, rel_heights = c(8.5, 1.5),
        align = 'v', labels = c("  1", "  2")
    )
}


#+ ws_tract, fig.fullwidth=TRUE, fig.width=8, fig.height=10, cache=FALSE
snp %>% plot_align_qual("ws")

#+ sw_tract, fig.fullwidth=TRUE, fig.width=8, fig.height=10, cache=FALSE
snp %>% plot_align_qual("sw")

#+ s_tract, fig.fullwidth=TRUE, fig.width=8, fig.height=10, cache=FALSE
snp %>% plot_align_qual("s")

#+ w_tract, fig.fullwidth=TRUE, fig.width=8, fig.height=10, cache=FALSE
snp %>% plot_align_qual("w")

#+ align1, fig.fullwidth=TRUE, fig.width=21, fig.height=21, cache=FALSE
plot_grid(
    snp %>% plot_align_qual("s"),
    snp %>% plot_align_qual("w"),
    snp %>% plot_align_qual("ws"),
    snp %>% plot_align_qual("sw"),
    labels = c("A", "B", "C", "D"), ncol = 2
)
## interactive only
## save_as_a3("test.pdf")

#' # Comparaison des longueurs de trace de conversion
#'
#' On peut s'attendre sous l'hypothèse gBGC à ce que les traces de conversions
#' soient plus longues lorsqu'on transforme par un plasmide porteur de S que
#' lorsqu'on transforme par un plasmide porteur de W.
#'

#+ len_distrib, fig.margin=TRUE
snp %>%
    filter(cons == "x") %>%
    group_by(name, mutant) %>%
    summarise(len = max(refp) - min(refp) ) %>%
    ggplot(aes(x = factor(mutant, levels = c("s", "w", "sw", "ws")),
               y = len)) +
    geom_point(alpha = 0.2) +
    stat_summary(fun.y = "mean", geom = "point", color = "red", shape = "|",
                 size = 7) +
    labs(x = "mutant", y = "Longueur\ndu tract") +
    coord_flip()

snp %>%
    filter(cons == "x") %>%
    group_by(name, mutant) %>%
    summarise(len = max(refp) - min(refp) ) %>%
    group_by(mutant) %>%
    summarise(mean_len = mean(len))

#+ len_distrib2, fig.margin=TRUE
snp %>%
    filter(cons == "x") %>%
    group_by(name, mutant) %>%
    summarise(len = max(refp) - min(refp) ) %>%
    ggplot(data = ., aes(x = len, fill = mutant )) +
    geom_histogram(binwidth = 5, position = "dodge") +
    facet_grid(mutant ~ .) +
    scale_fill_brewer(palette = "Set1") +
    labs(x = "Longueur de trace de conversion",
         y = "",
         title = "Distribution des\nlongueurs de trace de conversion")

#' Pour l'instant on ne voit quand même pas grand chose de flagrant. Il
#' semblerait même que les plasmides `S` introduisent des traces de conversion
#' plus courtes, même si rien de quantifié pour l'instant.

#' # Comptage des alternances SNP SS et WW
#'
#' ## Comptage des doublets
#'
#' En fin de trace de conversion, on s'attend à une alternance entre le génotype
#' du donneur et le génotype du receveur. En fin de trace de conversion, on aura
#' donc des doublets (GC)(GC) ou (AT)(AT), Sous l'hypothèse gBGC, les traces de
#' conversions se terminent préférentiellement sur des SNPs GC. On aurait donc
#' plus de doublets (GC)(GC).
#'
#' On pourrait simplement compter les doublets en question par construction SW
#' et WS, mais il existe des cas intermédiaires où le génotype parental est
#' restauré. Il faut donc discrimer les doublets en fin de trace de conversion
#' des doublets intermédiaires correspondant aux cas complexes. En fait, ces cas
#' complexes se manifestent par des *triplets* (GC)(GC)(GC) ou (AT)(AT)(AT) dans
#' la trace de conversion, où le SNP médian correspond au génotype parental.

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

ws_et_sw = c("ws", "sw")
w_et_s = c("w", "s")

#' Dans les manips ws et sw, on a
#' `r snp %>% count_doublet(mut = ws_et_sw, "WW")` doublets WW contre
#' `r snp %>% count_doublet(mut = ws_et_sw, "SS")` doublets SS.
#'
#' Si on sépare manip par manip, on a
#' `r snp %>% count_doublet(mut = "ws", "WW")` doublets WW dans la manip WS
#' contre `r snp %>% count_doublet(mut = "ws", "SS")` doublets SS, et
#' `r snp %>% count_doublet(mut = "sw", "WW")` doublets WW dans la manip sw
#' contre `r snp %>% count_doublet(mut = "sw", "SS")` SS.
#'
#' Si on compte les triplets, autrement dits les cas complexes, on a
#' `r snp %>% count_doublet(mut = ws_et_sw, "SSS")` triplets SSS dans les manips
#' alternants, et `r snp %>% count_doublet(mut = ws_et_sw, "WWW")` WWW.
#'
#' Dans les manips sans alternance, on a
#' `r snp %>% count_doublet(mut = w_et_s, "WSW")` WSW et
#' `r snp %>% count_doublet(mut = w_et_s, "SWS")` SWS.
#'

#' On peut représenter ça de façon graphique par la méthode suivante.
#'

get_marker_only <- function(data, snp_only = TRUE) {
    ## permet de ne conserver que les marqueurs dans le jeu de donnée
    snp_without_N <-  data %>% filter(cons != "N", cons != "-", name != "psw21")

    if (snp_only) filter(snp_without_N, cons == "X" | cons == "x")
    else snp_without_N
}

set_mutant_factor <- function(data) {
    ## attribue les bons niveaux de facteurs aux donneurs s.
    mutate(data, mutant = factor(mutant, levels = c("s", "w", "sw", "ws")))
}

set_seq_factor <- function(data) {
    ## attribue les bons niveaux de facteurs aux différentes séquences. Choisies
    ## de façon à ce que la base expérimentale apparaisse en rouge.
    mutate(data, seq = factor(seq, levels = c("seqb", "refb", "snpb"),
                              labels = c("Expérimental", "Référence", "Donneur")))
}

set_base_factor <- function(data) {
    ## attribue les bons niveaux de facteurs aux différentes bases, de façon à ce
    ## que S et W clusterent.
    mutate(data, base = factor(base, levels = c("G", "C", "T", "A")))
}

#+ bascul1, fig.fullwidth=TRUE, fig.width=15
snp %>%
    keep_clean_only() %>%
    filter(cons == "x") %>%
    group_by(name) %>%
    summarise(bascul = min(refp)) %>%
    inner_join(snp, .) %>%
    filter(refp == bascul) %>%
    rowwise() %>% mutate(sens = find_polarity(refb, seqb)) %>%
    set_mutant_factor() %>%
    ggplot(aes(x = refp, fill = sens)) +
    geom_histogram(binwidth  = 10, position = "dodge") +
    facet_grid( ~mutant) +
    scale_fill_brewer(palette = "Set1", labels = c("vers AT", "vers GC")) +
    labs(x = "Position sur la référence", y = "", color = "", fill = "") +
    legend_position(0.9, 0.6)

## /*
## #' Si on compte le nombre de positions WWW, on en a 3 dans la manip sw :
## snp %>% keep_clean_only() %>% make_genotype(c("ws", "sw"), qual_ = 0) %>% count_kmer("WWW") #%>% knitr::kable()

## #' Si on compte le nombre de marqueurs SSS, on en a 9 sur les deux manips :
## snp %>% make_genotype(c("ws", "sw"), qual_ = 0) %>% count_kmer("SSS") %>% select(anom) %>% sum()
## #' Si on compte seulement les doublets manips par manips :
## snp %>% make_genotype(c("ws")) %>% count_kmer("WWW") #%>% select(anom) %>% sum()
## snp %>% make_genotype(c("ws")) %>% count_kmer("SS") %>% select(anom) %>% sum()
## snp %>% make_genotype(c("sw")) %>% count_kmer("SS") %>% select(anom) %>% sum()
## snp %>% make_genotype(c("sw")) %>% count_kmer("WW") %>% select(anom) %>% sum()

## #' Les mêmes comptages pour les manips précédentes confirment ce qu'on savait déjà :
## snp %>% make_genotype("w") %>% count_kmer("WSW") %>% knitr::kable()
## snp %>% make_genotype("s") %>% count_kmer("WSW")
##                                         # peanuts.

## #' Les mêmes comptages sans filtrer sur la qualité, on retrouve les 6 positions
## #' que Laurent voyait déjà sur les alignements:
## snp %>% make_genotype("w", qual_ = 0) %>% count_kmer("WSW") %>% knitr::kable()
## snp %>% make_genotype("s", qual_ = 0) %>% count_kmer("SWS") %>% knitr::kable()

## snp %>% count_doublet("w", "WSW")
## snp %>% count_doublet("s", "SWS")

## */
#'
#'
#' # Quantification du rapport en base
#'
#' Franck se disait qu'il serait peut-être judicieux de comparer globalement le
#' taux de GC global de nos séquences, sachant qu'on introduit des SNPs AT et GC
#' dans toutes nos séquences, en tout cas pour la manip SW / WS. Au vu des
#' graphiques suivants, on peut déjà se rendre compte d'une chose : le taux de
#' GC ne se comporte pas comme attendu sous l'hypothèse d'une conversion biaisée
#' dans les manips SW et WS.
#'
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


snp %>%
    get_gc_content() %>%
    filter(name != "psw21") %>%
    mutate(mutant = factor(mutant, levels = c("s", "w", "sw", "ws"))) %>%
    ## group_by(mutant, type) %>%
    ## summarise(gc_mean = mean(GC)) %>%
    ggplot(aes(y = GC,
               x = factor(type, levels = c("snp", "exp", "ref"),
                          labels = c("Donneur", "Expérimental", "Receveur")))) +
    geom_jitter(aes(color = mutant), width = 1/4, alpha = 1/3, size = 1) +
    geom_line(aes(group = name, color = mutant), alpha = 1/20) +
    stat_summary(fun.y = mean, geom = "point", color = "black", shape = "¦",
                 size = 6) +
    facet_grid(mutant ~ .) +
    labs(y = "% GC", x = "", title = "Comparaison des taux de GC") +
    scale_color_brewer(palette = "Set1") +
    coord_flip()

#' # Enquête de néomutation
#'

is_inside_conv <- function(data, mut) {
  data %>%
      filter(mutant %in% mut, cons == "x") %>%
      group_by(name) %>%
      summarise(bascul = min(refp)) %>%
      inner_join(data, ., by = "name") %>%
      keep_clean_only() %>%
      ## garde les positions après le point de switch qui correspondent au génotype parental.
      filter(refp > bascul, cons == "X" | cons == "x") %>%
      mutate(tract = TRUE) %>%
      select(name, bascul, refp, tract) %>%
      left_join(data, .)
}

## wip
snp %>%
    is_inside_conv(c("ws", "sw", "s", "w")) %>%
    filter(name == "pS10")

## /wip



sessionInfo()
