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
      panel.background=element_rect(fill=alpha(color.background, 1), color=color.background),
      plot.background=element_rect(fill=alpha(color.background, 0.5), color=color.background),
      panel.border=element_rect(color=color.background),
      ## panel.border=element_rect(size = 1, color=color.grid.major),
      panel.grid.major=element_line(color=color.grid.major,size=.25, linetype = "dotted"),
      panel.grid.minor=element_blank(),
      axis.ticks=element_blank(),
      strip.text.y = element_text(angle = 0, size = 8),
      ## Format the legend, but hide by default
      legend.position="none",
      legend.background = element_blank(),
      legend.text = element_text(size=9,color=color.axis.title),
      legend.key = element_rect(fill = scales::alpha("gray", 0),
                                colour = scales::alpha("gray", 0)),
      ## Set title and axis labels, and format these and tick marks
      plot.title = element_text(color = color.title, size = 14, face = "bold", vjust = 1.25, hjust = 0, family = "Ubuntu"),
      axis.text.x = element_text(size = 7,color = color.axis.text),
      axis.text.y = element_text(size = 7,color = color.axis.text),
      axis.title.x = element_text(size = 8, color = color.axis.title, vjust = 0, hjust = 0.9),
      axis.title.y = element_text(size = 8, color = color.axis.title, vjust = 0.9, angle = 0 ),
      ## Plot margins
      plot.margin = unit(c(0.35, 0.2, 0.3, 0.35), "cm")
    )
}
#
theme_set(theme_bw() + fte_theme())
#
legend_position <- function(x=NULL, y=NULL) {
  custom_legend <- theme(legend.margin = unit(-0.3,"lines"),
                         legend.key = element_rect(fill = scales::alpha("gray", 0),
                                                   colour = scales::alpha("gray", 0)))
  if (is.null(x) & is.null(y)) custom_legend + theme(legend.position = "bottom")
  else custom_legend + theme(legend.position = c(x, y))
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

summary(snp$qual)

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

filter(snp, cons == "X")

#' Plutôt bon signe ! Phruscle a fonctionné comme il faut, _ie_ toutes les bases
#' qu'on a sorties dans nos alignements ont une correspondance dans le fichier
#' phred !.
#'

#' # Analyse de la qualité des séquences
#'

#' On regarde globalement la distribution de la qualité des bases :

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
  geom_line(aes(group = name), alpha = 0.1, size = 0.1 ) +
  scale_color_viridis()

#+ qual2, fig.margin=TRUE, fig.cap="Position des bases de qualité < 30"
snp %>%
  filter(qual < 30) %>%
  ggplot(data = ., aes(x = refp, y = qual )) +
  geom_point(alpha = 0.1, size = 0.1)

#' Comme attendu, la qualité est moindre en fin de run, et plutôt bonne au
#' début^[On séquence de droite à gauche, la cassette de résistance est à la
#' position avant -1]. Ce sont des séquences trimmées par phd. On n'est pas
#' obligé de trimmer. Peut-être refaire les analyses avec les séquences non
#' trimmées.
#'

#' # Distribution de la longueur des séquences.

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

#' Une fonction pour faire un alignement des positions d'intérêt, sans avoir à
#' répéter le code deux fois pour sw et ws. J'ai _hardcodé_ les limites de la
#' référence basées sur le graphique dans la marge.

find_polarity <- function(reference, experimental) {
  is_w <- function(base) { ifelse(base == "A" | base == "T", TRUE, FALSE)}
  is_s <- function(base) { ifelse(base == "C" | base == "G", TRUE, FALSE)}

  if      (is_w(reference) & is_w(experimental)) "ww"
  else if (is_w(reference) & is_s(experimental)) "ws"
  else if (is_s(reference) & is_w(experimental)) "sw"
  else if (is_s(reference) & is_s(experimental)) "ss"
  else ""
}

plot_align <- function(data, mutant_) {

  sort_by_tract_length <- function(mut) {
    data %>%
      filter(mutant == mut) %>%
      group_by(name, mutant) %>%
      filter(cons == "x") %>%           # dernière position correspondant au donneur
      summarise(len = min(refp)) %>%
      ungroup() %>%
      arrange(len) %>%
      {.$name}
  }
  ## print(sort_by_tract_length(mutant_) %>% length())
  red <- RColorBrewer::brewer.pal(n = 4, "Set1")[1]
  blue <- RColorBrewer::brewer.pal(n = 4, "Set1")[2]
  green <- RColorBrewer::brewer.pal(n = 4, "Set1")[3]
  violet <- RColorBrewer::brewer.pal(n = 4, "Set1")[4]

  data %>%
    group_by(name, mutant) %>%
    filter(#refp < 750 & refp > 57,
           mutant == mutant_,
           cons == "x" | cons == "X" ) %>%
    rowwise() %>%
    mutate(sens =  find_polarity(refb, seqb)) %>%
    ggplot(aes(x = factor(refp),
               y = factor(name, levels = sort_by_tract_length(mut=mutant_)))) +
    geom_point(aes(color = sens, size = qual, alpha=qual)) +
    ## add the donneur seq
    geom_text(aes(label = snpb, x = factor(refp), y = -3), color = red,
              vjust = -0.5, size = 2, family = "Ubuntu Light") +
    ## add the ref seq
    geom_text(aes(label = refb, x = factor(refp), y = length(sort_by_tract_length(mutant_)) + 10),
              color = blue, vjust = 6.5, size = 2, family = "Ubuntu Light") +
    scale_color_brewer(palette = "Set1", guide = "none") +
    scale_alpha( range=c(1/5, 0.8), guide=FALSE ) +
    scale_size(range = c(1, 3), breaks = c(10, 50),
               labels = c("Faible", "Forte")) +
    labs(x = "Position sur la référence", y = "", size = "Qualité",
         title = paste("Alignement pour la manip", toupper(mutant_))) +
    theme(legend.direction = "horizontal",
          ## legend.position  = c(0.8, -0.2),
          legend.margin = unit(0,"lines"),
          ## legend.key = element_rect(fill = scales::alpha("gray", 0),
          ##                           colour = scales::alpha("gray", 0)),
      panel.grid.major.y = element_line(size = 0.1, linetype = "dotted")) +
    legend_position(0.8, 0.98)
}

snp %>% plot_align("ws")

## snp %>% plot_align("w")
snp
snp %>%
  filter(cons == "x" | cons == "X")

plot_qual <- function(data, mutant_) {
  data %>%
    filter(mutant == mutant_) %>%
    ggplot(aes(x = refp, y = qual)) +
    geom_point(size = 0.1, alpha = 1/20) +
    scale_color_viridis(end = 0.95) +
    labs(x = "Position sur la référence",
         y = "") +
    ## coord_cartesian(xlim = c(50, 730)) +
    geom_smooth(se = FALSE, color = "red", linetype = "dotted")
}

plot_align_qual <- function(data, mutant_) {
  plot_grid(
    data %>% plot_align(mutant_),
    data %>% plot_qual(mutant_),
    ncol = 1, rel_heights = c(8.5, 1.5),
    align = 'v', labels = c("  1", "  2")
  )
}


n_distinct(snp %>% filter(mutant == "ws") %>% group_by(name) %>% .$name )

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
  labels = c("A", "B", "C", "D"), ncol = 4
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
  scale_fill_viridis(discrete = TRUE) +
  labs(x = "Longueur de trace de conversion",
       y = "",
       title = "Distribution des\nlongueurs de trace de conversion")

#' Pour l'instant on ne voit quand même pas grand chose de flagrant. Il
#' semblerait même que les plasmides `S` introduisent des traces de conversion
#' plus courtes, même si rien de quantifié pour l'instant.
#'

#' Comptage des alternances SNP SS et WW
#'
#' Pour la manip WS et SW, le but est de compter les positions où on voit des alternances de WWW ou SSS.

#' Certaines séquences ne sont pas des plus propres. On peut certainement les
#' virer pour ne garder que les bonnes séquences. Comme pour le cas de la
#' séquence psw21.

#+ mean_qual, fig.margin = TRUE
snp %>%
  group_by(name) %>%
  summarise(median = median(qual)) %>%
  qplot(data = . , median, binwidth = 1) +
  geom_vline(xintercept = 40, color = "red")

#' Dans la fonction suivante, je me suis donc dit qu'il était possible de ne
#' garder que les séquences dont la qualité médianne est supérieure à 40.
keep_clean_only <- function(data, qual_ = 40) {

  get_dirty_seq <- function(data)
    group_by(data, name) %>%
    summarise(median = median(qual)) %>%
    filter(median < qual_) %>%
    select(name) %>% unlist()

  data %>%
    filter(!(name %in% get_dirty_seq(data)))
}

#' La fonction suivante renvoit un `tbl_df` avec le nom, le type de mutant, le
#' genotype.
#' @param data: le jeu de donnée snp,
#' @param mutant: crée un filtre par mutant. Selon les mutants, on ne recherche
#'   pas les mêmes séquences.
make_genotype <- function(data, mutant_, qual_ = 30) {
  weak_or_strong <- function(base) {
    is_w <- function(base) { ifelse(base == "A" | base == "T", TRUE, FALSE)}
    is_s <- function(base) { ifelse(base == "C" | base == "G", TRUE, FALSE)}

    if (is_w(base)) "W"
    else if (is_s(base)) "S"
    else if (base == "N") "N"
    else stop(base, " is not a DNA base")
  }


  data %>%
    filter(cons == "x" | cons == "X", name != "psw21") %>%
    filter(mutant %in% mutant_) %>%
    rowwise() %>%
    # convertit les bases d'une qualité inférieure à qual_ en N
    mutate(seqb = ifelse(qual < qual_, "N", seqb )) %>%
    mutate(geno = weak_or_strong(seqb)) %>%
    ungroup() %>% group_by(name, mutant) %>%
    summarise(geno = toString(geno) %>% gsub(", ", "", .))
}

#' La fonction suivante permet de compter le nombre d'anomalies, autrement dit
#' de SNP inattendus dans les traces de conversion. Renvoit une `tbl_df`
#' @param data : le jeu de donné snp
#' @param kmer : le kmer à recenser dans le génotype.
count_kmer <- function(data, kmer) {
  count_anomaly <- function(string, kmer_) {
    find_xxx <- function(string) {
      function(kmer_) {
        XxX <- gregexpr(kmer_, string) %>% unlist()
        xXx <- gregexpr(kmer_, string) %>% unlist()
        counter <- 0
        if      (XxX[1] > 0) counter <- counter + length(XxX)
        else if (xXx[1] > 0) counter <- counter + length(xXx)
        counter
      }
    }
    find_xxx(string)(kmer_)
  }

  data %>%
    rowwise() %>%
    mutate(anom = count_anomaly(geno, kmer)) %>%
    filter(anom > 0)
}

#' Compte les doublets des séquences propres.
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
snp %>% count_doublet(mut = ws_et_sw, "WW")
snp %>% count_doublet(mut = ws_et_sw, "SS")

#' Si on compte le nombre de positions WWW, on en a 3 dans la manip sw :
snp %>% keep_clean_only() %>% make_genotype(c("ws", "sw"), qual_ = 0) %>% count_kmer("WWW") #%>% knitr::kable()
#' Si on compte le nombre de marqueurs SSS, on en a 9 sur les deux manips :
snp %>% make_genotype(c("ws", "sw"), qual_ = 0) %>% count_kmer("SSS") %>% select(anom) %>% sum()
#' Si on compte seulement les doublets manips par manips :
snp %>% make_genotype(c("ws")) %>% count_kmer("WWW") #%>% select(anom) %>% sum()
snp %>% make_genotype(c("ws")) %>% count_kmer("SS") %>% select(anom) %>% sum()
snp %>% make_genotype(c("sw")) %>% count_kmer("SS") %>% select(anom) %>% sum()
snp %>% make_genotype(c("sw")) %>% count_kmer("WW") %>% select(anom) %>% sum()

#' Les mêmes comptages pour les manips précédentes confirment ce qu'on savait déjà :
snp %>% make_genotype("w") %>% count_kmer("WSW") %>% knitr::kable()
snp %>% make_genotype("s") %>% count_kmer("WSW")
                                        # peanuts.

#' Les mêmes comptages sans filtrer sur la qualité, on retrouve les 6 positions
#' que Laurent voyait déjà sur les alignements:
snp %>% make_genotype("w", qual_ = 0) %>% count_kmer("WSW") %>% knitr::kable()
snp %>% make_genotype("s", qual_ = 0) %>% count_kmer("SWS") %>% knitr::kable()

snp %>% count_doublet("w", "WSW")
snp %>% count_doublet("s", "SWS")

#' # Quantification du rapport en base
#'
#' Franck se disait qu'il serait peut-être judicieux de comparer globalement le
#' taux de GC global de nos séquences, sachant qu'on introduit des SNPs AT et GC
#' dans toutes nos séquences, en tout cas pour la manip SW / WS.
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

snp %>%
  get_gc_content() %>%
  filter(name != "psw21") %>%
  mutate(mutant = factor(mutant, levels = c("s", "w", "sw", "ws"))) %>%
  group_by(mutant, type) %>%
  ## summarise(gc_mean = mean(GC)) %>%
  ggplot(
    aes(y = GC,
        x = factor(type, levels = c("snp", "exp", "ref"),
                   labels = c("Donneur", "Expérimental", "Receveur")))) +
  geom_jitter(aes(color = mutant), width = 1/3, alpha = 0.5) +
  stat_summary(fun.y = mean, geom = "point", color = "black", shape = "¦",
               size = 6) +
  facet_grid(mutant ~ .) +
  labs(y = "% GC", x = "", title = "Comparaison des taux de GC") +
  scale_color_brewer(palette = "Set1") +
  coord_flip() +
  theme(text = element_text(family = "Palatino"))

snp %>%
  keep_clean_only() %>%
  get_gc_content()  %>%
  set_mutant_factor() %>%
  filter(name != "psw21") %>%
  mutate(type = factor(type, levels = c("ref", "exp", "snp"),
                       labels = c("Receveur", "Expérimental", "Donneur"))) %>%
  ggplot(aes(x = type, y = GC, color = mutant)) +
  ## geom_point() +
  geom_line(aes(group = name), alpha = 0.3) +
  scale_color_brewer(palette = "Set1") +
  labs(x = "", y = "%GC", color = "Donneur :") +
  legend_position() +
  theme(text = element_text(family = "Ubuntu Light"))

## fonts()

snp %>%
  set_mutant_factor() %>%
  get_marker_only(TRUE) %>%
  select(name, mutant, refp, seqb, refb, snpb) %>%
  tidyr::gather(key = seq, "base", 4:6) %>%
  group_by(mutant, seq, base) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  set_seq_factor() %>%
  set_base_factor() %>%
  ggplot(aes(x = base, y = count, color = seq, shape = mutant )) +
  geom_point() +
  facet_grid(mutant ~ .) +
  scale_shape_discrete(guide = FALSE) +
  scale_color_brewer(palette = "Set1") +
  labs(y = "", x = "", shape = "Séquence :", color = "Donneur :") +
  coord_flip() +
  legend_position()

sessionInfo()
