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
library(ggthemes)
library(cowplot)
library(extrafont)
library(gcbiasr)

set_gcbiasr_theme()

source("load_data.R")

snp <- read_phruscle("../../data/phruscle_snpcall.csv")

source("analysis_function.R")

#' J'ai donc obtenu un tableau de donnée de la forme suivante, où :
#'
#' - `cons`: un code maison pour déterminer le _consensus_. `.` si les trois
#'   bases sont identiques, `x` si la base séquencée est la base introduite, `X`
#'   si la base séquencée est la base sauvage, `N` si c'est encore autre chose.
#' - `name`: le nom expérimental de la séquence
#' - `refb`: la base de référence
#' - `expb`: la base séquencée.
#' - `snpb`: la base sur le gène synthétique introduit.
#' - `refp`: la position sur la référence.
#' - `expp`: la position sur la séquence expérimentale.
#' - `qual`: sa qualité.
#' - `mutant`: ordinal, la manip.
#' - `lastmp`: last marker position. 
#' - `switchp`: la position de bascule. 
#' - `switchb`: la base à la position de bascule. 
#' - `inconv`: logical, TRUE si la base est dans la trace de conversion.
#' - `isrestor`: logical, TRUE si la base est dans la trace de conversion, est
#' un marqueur et correspond à l'haplotype du donneur.
#' Vincent appelle ça joliment un alignement pairwise en colonne.

head(snp)

#' # Observations générales
#'
#' ## Nombres de séquences ?
n_distinct(snp$name)

snp %>%
    select(mutant, name) %>%
    group_by(mutant) %>%
    summarise(count = n_distinct(name)) %>%
    knitr::kable()

#' On a donc bien toutes nos séquences.

snp %>%
    group_by(mutant) %>%
    summarise(count = n()) %>%
    knitr::kable(align = "c")

#+ pos_number1, fig.margin=TRUE
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
    facet_grid(mutant ~ .) +
    labs(x = "longueur", y = "")

#+ snp_distrib0, fig.margin=TRUE
snp %>%
    filter(cons == "x" | cons == "X") %>%
    ## filter(refp < 750 & refp > 40) %>%
    ggplot(aes(x = refp)) +
    geom_histogram(binwidth = 1)
snp

## "cuicui"
## filter(snp, cons == "x") %>% summarise(min = min(refp))

## snp %>% plot_align("ws")
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
    snp %>% plot_align("s"),
    snp %>% plot_align("w"),
    snp %>% plot_align("ws"),
    snp %>% plot_align("sw"),
    labels = c("A", "B", "C", "D"), ncol = 2
)


## interactive only
## save_as_a3("plot_align_moire.pdf")

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
    summarise(len = min(refp) ) %>%
    ggplot(aes(x = factor(mutant, levels = c("s", "w", "sw", "ws")), y = len)) +
    geom_point(alpha = 0.2) +
    stat_summary(fun.y = "mean", geom = "point", color = "red", shape = "|",
                 size = 7) +
    labs(x = "mutant", y = "Longueur\ndu tract") +
    coord_flip()

snp %>%
    filter(cons == "x") %>%
    group_by(name, mutant) %>%
    summarise(len = min(refp) ) %>%
    group_by(mutant) %>%
    summarise(mean_len = mean(len)) %>%
    knitr::kable(align = "c")

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

## /* wip
snp %>%
    ## get_marker_only() %>%
    filter(cons == "x") %>%
    mutate(trans = paste0(refb, seqb)) %>%
    group_by(mutant, trans) %>%
    summarise(count = n()) %>%
    ggplot(aes(x = trans, y = count, color = trans )) +
    geom_point() +
    facet_grid( mutant ~ .) +
    coord_flip()


## */

sessionInfo()
