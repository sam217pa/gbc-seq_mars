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
snp <- read_phruscle("../../data/phruscle_snpcall.csv")

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

#+ pos_number1, fig.margin=TRUE, cache = TRUE
plot_grid(
    snp %>%
    group_by(mutant, name) %>%
    summarise(count = n()) %>%
    ggplot(aes(x = count, fill = mutant)) +
    geom_histogram(binwidth = 1) +
    scale_fill_brewer(palette = "Set1") +
    labs(x = "Nombre de positions", y = "",
         title = "Distribution de la longueur des séquences"),
    snp %>%
    keep_clean_only() %>%
    group_by(mutant, name) %>%
    summarise(count = n())  %>%
    ggplot(aes(x = count, fill = mutant)) +
    geom_histogram(binwidth = 1) +
    labs(x = "Nombre de positions", y = "",
         title = "Distribution de la longueur des séquences filtrées") +
    ## coord_cartesian(xlim = c(660, 800)) +
    scale_fill_brewer(palette = "Set1"),
    ncol = 1
)

#' J'ai voulu déterminer s'il existe des bases dans les sorties de phruscle qui
#' ne correspondent pas aux sorties de phd. Autrement dit si la base `base` est
#' bien toujours la même que la base `expb`.

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

#+ qual1, fig.fullwidth=TRUE, fig.width=15, cache = TRUE
quality_control(snp)

#' Comme attendu, la qualité est moindre en fin de run, et plutôt bonne au
#' début^[On séquence de gauche à droite, la cassette de résistance est en 3'].
#' Ce sont des séquences trimmées par phd. On n'est pas obligé de trimmer, mais
#' ça permet de ne garder que les séquences qui vallent le coup. Peut-être
#' refaire les analyses avec les séquences non trimmées.

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

#' Après avoir appliqué tout ça, on a perdu `r n_distinct(keep_clean_only(snp, qual = 0)$name) -  n_distinct(keep_clean_only(snp)$name)` séquences
#'
#'
#' # Distribution de la longueur des séquences.
#'
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

#' # Analyses visuelles des traces de conversion

#+ ws_tract, fig.fullwidth=TRUE, fig.width=8, fig.height=10
plot_align(snp, "ws")

#+ sw_tract, fig.fullwidth=TRUE, fig.width=8, fig.height=10
plot_align(snp, "sw")

#+ s_tract, fig.fullwidth=TRUE, fig.width=8, fig.height=10
plot_align(snp, "s")

#+ w_tract, fig.fullwidth=TRUE, fig.width=8, fig.height=10
plot_align(snp, "w")

#+ align1, fig.fullwidth=TRUE, fig.width=21, fig.height=21
plot_grid(
    snp %>% plot_align("s") + theme(legend.position = "bottom")
   ,snp %>% plot_align("w") + theme(legend.position = "bottom")
   ,snp %>% plot_align("ws") + theme(legend.position = "bottom")
   ,snp %>% plot_align("sw") + theme(legend.position = "bottom")
   ,labels = c("A", "B", "C", "D")
   ,ncol = 4
)

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

#+ bascul1, fig.fullwidth=TRUE, fig.width=15
plot_conv_len(snp)

#' Pour l'instant on ne voit quand même pas grand chose de flagrant. Il
#' semblerait même que les plasmides `S` introduisent des traces de conversion
#' plus courtes, même si rien de quantifié pour l'instant.
#'
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

#+ doublet, fig.margin = TRUE
knitr::kable(count_last_snp(snp))
plot(count_last_snp(snp))

#+ restor, fig.margin = TRUE
knitr::kable(count_restor(snp))
plot(count_restor(snp))

#'


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

get_gc_content(data = snp) %>% plot()

#' # Enquête de néomutation
#'

devtools::session_info()
