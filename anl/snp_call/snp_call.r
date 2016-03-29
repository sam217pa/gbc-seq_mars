#' ---
#' title: Analyse globale des SNP
#' author: Samuel BARRETO
#' date:
#' output:
#'   tufte::tufte_html:
#'     highlight: tango
#'     toc: true
#'     toc_depth: 2
#' ---


#' # Lecture des données et nettoyage
#'
#' Il faut dans un premier temps lire les données, les nettoyer et associer à
#' chaque séquence la catégorie de plasmide transformeur.
#'
#' La table des identifiants est faite par le script python.
#'

#+ setup, include=FALSE
library(knitr)
library(tufte)
opts_chunk$set(cache = FALSE, dev = 'png', include = TRUE,
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
          legend.background = element_rect(fill=color.background),
          legend.text = element_text(size=7,color=color.axis.title),
          ## Set title and axis labels, and format these and tick marks
          plot.title=element_text(color=color.title, size=14,
                                  face = "bold", vjust=1.25, hjust = 0),
          axis.text.x=element_text(size=7,color=color.axis.text),
          axis.text.y=element_text(size=7,color=color.axis.text),
          axis.title.x=element_text(size=8, color=color.axis.title, vjust=0, hjust = 0.8),
          axis.title.y=element_text(size=8, color=color.axis.title, vjust = 0.9, angle = 0 ),
          ## Plot margins
          plot.margin = unit(c(0.35, 0.2, 0.3, 0.35), "cm"))
}

#+ setup2, include = TRUE
library(dplyr)
library(ggplot2)
library(readr)
library(viridis)
library(purrr)
library(ggthemes)

theme_set(theme_bw() + fte_theme())

data_location <- "../../data/snp_call/snp_calling.dat"
id_table_loc <- "../../data/id_table.dat"

#' J'ai relié les données par un `inner_join` sur la base du nom de
#' l'identifiant de la séquence GATC.

snp <- inner_join(
  x = read_delim(file = data_location, delim = " "),
  y = read_delim(file = id_table_loc, delim = " "),
  by=c("read_name" = "id")) %>%
  select(
    -match,
    -subject_name,
    -index_of_subject,
    read = read_name,
    refb = s_base,
    readb = q_base,
    -s_qual,
    base_q = q_qual,
    refp = offset_on_subject,
    readp = offset_on_read,
    -length_of_subject,
    -length_of_snp,
    readfp = start_match_of_read,
    readlp = end_match_of_read,
    dir = match_direction) %>%
  mutate(name = gsub("-1073bis", "", name))

snppos_sw <- filter(snp, mutant == "sw") %>%
  select(refp) %>%
  distinct() %>%
  unlist() %>%
  as.vector()
snppos_ws <- filter(snp, mutant == "ws") %>%
  select(refp) %>%
  distinct() %>%
  unlist() %>%
  as.vector()


#' # Observations générales
#'
#' ## Nombres de reads par transformation

snp %>%
  group_by(read, mutant) %>%
  summarise(count = n()) %>%
  group_by(mutant) %>%
  summarise(count = n()) %>%
  knitr::kable(col.names = c("Transformant", "Nombre de read"), align = "c")

#' ## Nombre de SNP par transformation

snp %>%
  group_by(mutant) %>%
  summarise(count = n()) %>%
  knitr::kable(col.names = c("Transformant", "Nombre de SNP"), align = "c")

#' ## Distribution du nombre de SNP

snp %>%
  group_by(read, mutant) %>%
  summarise(count = n()) %>%
  ggplot(aes(x = count)) +
  geom_density(adjust = 0.4, fill = "gray", alpha = 0.5) +
  labs(x = "Nombre de SNP",
       y = "Densite",
       title = "Distribution du nombre de SNP" )


#' ## Distribution de la qualité des SNP

#+ spark, cache=TRUE, fig.width=15, fig.height=7, fig.fullwidth=TRUE
snp %>%
  ggplot(aes(x = refp, y = base_q)) +
  geom_point(alpha = 1/5, size = 0.1) +
  geom_line(alpha = 1) +
  facet_wrap( ~ name) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size= 6),
        axis.text = element_text(size = 5),
        panel.margin = unit(0.1, "lines"))

#' # Répartition des SNP sur les séquences
#'
#' J'ai voulu regarder globalement la distribution des SNP sur les
#' séquences.

#+ repar, cache=FALSE, fig.margin=TRUE
snp %>%
  ggplot(aes(x = refp, fill = mutant )) +
  geom_histogram(binwidth = 10) +
  facet_grid(mutant~.) +
  labs(x = "Position", y = "")

#+ repar2, cache=FALSE, fig.margin=FALSE
snp %>%
  ggplot(aes(x = refp, fill = mutant )) +
  geom_density(aes(group = mutant), adjust = 1/3, alpha = 1/2) +
  geom_vline(xintercept = snppos_sw, alpha = 0.1) +
  labs(x = "Position",
       y = "Densite",
       title = "Distribution des SNPs sur les reads",
       fill = "Transformant") +
  theme(panel.grid.major.x = element_blank(),
        legend.position = c(0.8, 0.8))

#' ## Définition des transitions

#' J'ai définit une fonction qui permet de déterminer le type de
#' mutation, de strong à weak (`sw`), de weak à strong (`ws`) ou autre.

## #' une fonction qui définit le type de transition, de la référence
## #' à la séquence observée.
find_transition <- function(ref, base) {
  stopifnot(typeof(ref) == "character", typeof(base) == "character")

  is_W <- function(base) ifelse(base == "A" || base == "T", TRUE, FALSE)
  is_S <- function(base) ifelse(base == "C" || base == "G", TRUE, FALSE)

  if (is_W(ref) & is_S(base)) "ws"
  else if (is_S(ref) & is_W(base)) "sw"
  else if (is_S(ref) & is_S(base)) "ss"
  else if (is_W(ref) & is_W(base)) "ww"
  else stop("Base ", ref, " or base ", base, " is not a DNA base.",
            call. = FALSE)
}

#' J'ai donc regardé si on avait des SNPs improbables, *ie* des
#' mutations de W à W (`ww`) ou de S à S (`ss`).
snp %>%
  rowwise() %>%
  mutate(trans = find_transition(refb, readb)) %>%
  filter(trans %in% c("ss", "ww")) %>%
  print()

#' Il n'y en a pas. Plutôt une bonne nouvelle !
#'
#' ## Distribution par type de substitutions
#'
#' J'ai voulu voir si les différents types de substitutions ont des
#' répartitions différentes.

snp %>%
  rowwise() %>%
  mutate(trans = find_transition(refb, readb)) %>%
  ggplot(aes(x = refp, fill = trans)) +
  geom_density(adjust = 1/3, alpha = 1/2) +
  geom_vline(xintercept = snppos_sw, alpha = 0.1) +
  theme(panel.grid.major.x = element_blank()) +
  labs(x = "Position",
       y = "Densité",
       fill = "Transition",
       title = "Distribution des substitutions") +
  theme(legend.position = c(0.8, 0.8))

snp %>%
  rowwise() %>%
  mutate(trans = find_transition(refb, readb)) %>%
  ggplot(aes(x = refp )) +
  geom_histogram(aes(fill = trans, color = trans), binwidth = 10,
                 position = "dodge" ) +
  labs(x = "Position", y = "",
       title = "Distribution des transitions",
       fill = "Transition",
       color = "Transition") +
  theme(legend.position = c(0.75, 0.85))

#' # Trace de conversion
#'
#' ## Distribution de la longueur
#'
#' J'ai ensuite regardé comment étaient réparties les longeurs de
#' trace de conversion sur les séquences.
#'
snp %>%
  group_by(read) %>%
  summarise(max = max(refp)) %>%
  ggplot(aes(x = max)) +
  geom_histogram(binwidth = 10) +
  labs(x = "Longueur de trace de conversion",
       y = "",
       title = "Distribution de la longueur\nde trace de conversion" )

#' À première vue, le dernier SNP est souvent dans la partie 3' du
#' gène. Cette distribution là n'est pas celle qu'on avait avec les
#' séquences W et S.
#'
#' J'ai donc voulu voir s'il y avait une différence de longueur de
#' trace de conversion entre les séquences de SW et les séquences WS,
#' dans un premier temps en fréquence, puis en densité.

snp %>%
  group_by(read, mutant) %>%
  summarise(max = max(refp)) %>%
  ggplot(aes(x = max)) +
  geom_histogram(aes(fill = mutant), position = "dodge") +
  labs(fill = "Transformant", x = "Longueur de conversion tract", y = "") +
  ggtitle( "Distribution de la\nlongueur de trace de conversion")

snp %>%
  group_by(read, mutant) %>%
  summarise(max = max(refp)) %>%
  ggplot(aes(x = max, fill = mutant, color = mutant)) +
  geom_density(adjust = 1/3, alpha = 0.5) +
  labs(y = "Densite", x = "Position du point de basculement",
       title = "Position du point de basculement") +
  theme(legend.position = c(0.2, 0.81))

#' J'ai voulu voir si le point de switch était plus souvent à un SNP
#' strong que weak.

snp %>%
  group_by(read, mutant) %>%
  summarise(max = max(refp)) %>%
  inner_join(x = snp, y = .) %>%
  filter(refp == max) %>%
  rowwise() %>%
  mutate(trans = find_transition(refb, readb)) %>%
  ggplot(aes(x = max)) +
  geom_histogram(aes(fill = trans, color = trans), position = "dodge", binwidth = 10) +
  theme(legend.position = c(0.2, 0.8)) +
  labs(x = "Position du point de basculement",
       y = "",
       fill = "Mutation",
       color = "Mutation",
       title = "Type de mutation au point de switch")

#' # En(-)quête de néo-mutations…
#'
#' J'ai voulu voir si on retrouvait les néomutations du séquençage
#' précédent dans la trace de conversion. Le problème c'est que cette
#' fois la distribution du nombre de SNP à une position donnée ne
#' permet pas de discriminer facilement les positions de SNP calibré
#' avec les néo-mutations. J'ai donc décidé de ne regarder que les
#' positions où le nombre de SNP est de 1, en supposant qu'elles
#' incluent les positions de néo-mutation.
#'

snp %>%
  group_by(refp) %>%
  summarise(count = n()) %>%
  filter(count == 1) %>%
  print()

#' En clair, il n'y a pas de position où on ne retrouve qu'un
#' SNP. Donc si néo-mutation il y a, elle occasionne dans au moins
#' deux séquences différentes un même SNP. Ce qui est peu probable.

#+ nombresnp, fig.margin = TRUE
snp %>%
  group_by(refp) %>%
  summarise(count = n()) %>%
  ## filter(count < 5) %>%
  ggplot(aes(x = count)) +
  geom_histogram(binwidth = 1) +
  geom_vline(xintercept = 5, color = "red") +
  labs(x = "Nombre de SNP", y = "", title = "Distribution du nombre de SNP par position")

#' En fait aucune position de SNP calibrée ne génère moins de 5
#' mutations. On peut donc en conclure qu'il n'y a pas de néomutations
#' dans cette manip.
