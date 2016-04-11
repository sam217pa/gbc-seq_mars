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
opts_chunk$set(cache = TRUE, dev = 'png', include = TRUE,
               echo = TRUE, warning = FALSE, error = FALSE,
               message = FALSE)

fte_theme <- function() {
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
    theme(panel.background=element_rect(fill=color.background, color=color.background),
          plot.background=element_rect(fill=color.background, color=color.background),
          panel.border=element_rect(color=color.background),
          panel.grid.major=element_line(color=color.grid.major,size=.25),
          panel.grid.minor=element_blank(),
          axis.ticks=element_blank(),
          ## Format the legend, but hide by default
          legend.position="none",
          legend.background = element_rect(fill = scales::alpha(color.background, 0.3)),
          legend.text = element_text(size=7,color=color.axis.title),
          ## Set title and axis labels, and format these and tick marks
          plot.title=element_text(color=color.title, size=14,
                                  face = "bold", vjust=1.25, hjust = 0),
          axis.text.x=element_text(size=7,color=color.axis.text),
          axis.text.y=element_text(size=7,color=color.axis.text),
          axis.title.x=element_text(size=8, color=color.axis.title, vjust=0, hjust = 0.9),
          axis.title.y=element_text(size=8, color=color.axis.title, vjust = 0.9, angle = 0 ),
          ## Plot margins
          plot.margin = unit(c(0.35, 0.2, 0.3, 0.35), "cm"))
}

legend_position <- function(x=NULL, y=NULL) {
  if (is.null(x) & is.null(y)) theme(legend.position = "bottom")
  else theme(legend.position = c(x, y))
}

#+ setup2, include = TRUE
library(dplyr)
library(ggplot2)
library(readr)
library(viridis)
## library(purrr)
## library(ggthemes)
library(cowplot)
## library(gridExtra)

theme_set(theme_bw() + fte_theme())


data_location <- "../../data/phruscle_snpcall.csv"

snp <- read_csv(data_location) %>%
  mutate(
    name = gsub("-1073.+$", "", name),
    base = toupper(base)#,
    ## mutant = ifelse(grepl("ws", name), "ws",
    ##          ifelse(grepl("sw", name), "sw",
    ##          ifelse(grepl("W", name), "w",
    ##          ifelse(grepl("S", name), "s", ""))))
  )

find_mutant <- function(name) {
  if      (grepl("ws", name)) "ws"
  else if (grepl("sw", name)) "sw"
  else if (grepl("W", name )) "w"
  else if (grepl("S", name )) "s"
  else ""
}

snp <- snp %>%
  group_by(name) %>%
  summarise(count = n()) %>%
  rowwise() %>%
  mutate(mutant=find_mutant(name)) %>%
  inner_join(snp, .)

#'
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

head(snp)

#' # Observations générales
#'
#' ## Nombres de séquences ?
                                        #
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
#' répéter le code deux fois pour sw et ws. J'ai ̉hardcodé' les limites de la
#' référence basé sur le graphique dans la marge.

plot_align <- function(data, mutant_) {

  sort_by_tract_length <- function(mut) {
    data %>%
      filter(mutant == mut) %>%
      group_by(name, mutant) %>%
      filter(cons == "x") %>%
      summarise(len = max(refp) - min(refp)) %>%
      ungroup() %>%
      arrange(len) %>%
      {.$name}
  }

  find_polarity <- function(reference, experimental) {
    is_w <- function(base) { ifelse(base == "A" | base == "T", TRUE, FALSE)}
    is_s <- function(base) { ifelse(base == "C" | base == "G", TRUE, FALSE)}

    if      (is_w(reference) & is_w(experimental)) "ww"
    else if (is_w(reference) & is_s(experimental)) "ws"
    else if (is_s(reference) & is_w(experimental)) "sw"
    else if (is_s(reference) & is_s(experimental)) "ss"
    else ""
  }

  data %>%
    group_by(name, mutant) %>%
    filter(
      refp < 750 & refp > 57,
      mutant == mutant_,
      cons == "x" | cons == "X"
    ) %>%
    rowwise() %>%
    mutate(sens =  find_polarity(refb, seqb)) %>%
    ggplot(aes(x = factor(refp),
               y = factor(name, levels = sort_by_tract_length(mut=mutant_)))) +
    geom_point(aes(color = sens, size = qual, alpha=qual)) +
    scale_color_viridis(
      discrete = TRUE,
      end=0.95
      ## labels = c("Indel", "S vers S", "S vers W", "W vers S", "W vers W")
    ) +
    scale_alpha(
      range=c(1/5, 1),
      guide=FALSE
    ) +
    scale_size(
      trans = "reverse",
      range = c(1, 3),
      breaks = c(10, 50),
      labels = c("10", "50")) +
    labs(
      x = "Position sur la référence",
      y = "",
      color = "Remplacement",
      size = "Qualité",
      title = paste("Alignement pour la manip", toupper(mutant_))
    ) +
    theme(
      legend.position = "top",
      legend.margin=unit(-0.3,"lines"),
      legend.key = element_rect(
        fill = scales::alpha("gray", 0),
        colour = scales::alpha("gray", 0)),
      panel.grid.major.y = element_line(size = 0.1, linetype = "dotted")
    )
}

plot_qual <- function(data, mutant_) {
  data %>%
    filter(mutant == mutant_) %>%
    ggplot(aes(x=refp, y=qual, color=qual)) +
    geom_point(size=0.1, alpha=0.1) +
    scale_color_viridis(end=0.95) +
    labs(x="Position sur la référence",
         y="") +
    coord_cartesian(xlim=c(50, 730)) +
    geom_smooth(se = FALSE, color=viridis(1), linetype="dotted")
}

plot_align_qual <- function(data, mutant_) {
  plot_grid(
    data %>% plot_align(mutant_),
    data %>% plot_qual(mutant_),
    ncol=1, rel_heights = c(8.5, 1.5),
    align='v', labels=c("  1", "  2")
  )
}

#+ ws_tract, fig.fullwidth=TRUE, fig.width=8, fig.height=10
snp %>% plot_align_qual("ws")

#+ sw_tract, fig.fullwidth=TRUE, fig.width=8, fig.height=10
snp %>% plot_align_qual("sw")

#+ s_tract, fig.fullwidth=TRUE, fig.width=8, fig.height=10
snp %>% plot_align_qual("s")

#+ w_tract, fig.fullwidth=TRUE, fig.width=8, fig.height=10
snp %>% plot_align_qual("w")

## snp %>%
##   filter(cons != ".") %>%
##   rowwise() %>%
##   mutate(sens=find_polarity(refb, seqb)) %>%
##   filter(sens != "") %>%
##   print(n=100)

#+ align1, fig.fullwidth=TRUE, fig.width=21, fig.height=21
plot_grid(
  snp %>% plot_align_qual("s"),
  snp %>% plot_align_qual("w"),
  snp %>% plot_align_qual("ws"),
  snp %>% plot_align_qual("sw"),
  labels=c("A", "B", "C", "D"), ncol=4
)

## save_as_a3("test.png")

## snp %>%
##   group_by(mutant) %>%
##   summarise(count = n())
