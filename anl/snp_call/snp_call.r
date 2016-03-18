#' ---
#' title: "Rapide analyse des SNP"
#' author: "Samuel BARRETO"
#' date: "31 janvier 2016"
#' output:
#'   html_document:
#'     highlight: tango
#'     theme: journal
#'     toc: yes
#' ---

#' # Lecture des données
#'
#' Les données sont les tableurs de snp calling que GATC nous a fait
#' gratuitement. Ce script n'est que préliminaire, je referai le snp
#' calling comme il faut par la suite.

library(dplyr)
library(ggplot2)
library(readr)
library(viridis)

#' Lecture des données et combinaison dans un seul tableur. 
ws <- read_delim("../../data/csv/1582203.SNP.csv", delim = ";")
sw <- read_delim("../../data/csv/1582443.SNP.csv", delim = ";")
ws$type <- "ws"
sw$type <- "sw"

## #' une fonction qui définit le type de transition, de la référence à la séquence observée. 
find_transition <- function(ref, base) {
  stopifnot(typeof(ref) == "character",
            typeof(base) == "character")

  is_W <- function(base) ifelse(base == "A" || base == "T", TRUE, FALSE)
  is_S <- function(base) ifelse(base == "C" || base == "G", TRUE, FALSE)

  if (is_W(ref) & is_S(base)) "ws"
  else if (is_S(ref) & is_W(base)) "sw"
  else if (is_S(ref) & is_S(base)) "ss"
  else if (is_W(ref) & is_W(base)) "ww"
  else stop("Base ", ref, " or base ", base, " is not a DNA base.", call. = FALSE)
}

data <- rbind(ws, sw) # lie les deux jeux de données
colnames(data) <-  c("ref", "pos", "ref.base", "alt.base", "qual",
                     "query", "qpos", "qlen", "type") # change les noms de colonne
data <- mutate(data, query = gsub("-1073bis", "", query)) # supprime le nom de l'amorce

## ggplot default theme
theme_set(theme_minimal(base_size = 10, base_family = "Courier"))

#' Je regarde la distribution des SNP sur les séquences
#+ distrisnp
data %>%
  ggplot(aes(x = pos, fill = type)) +
  geom_histogram(binwidth = 10) +
  facet_grid(type~.) +
  xlab("Position sur la reference") +
  ylab("Nombre de SNP")

#' Je fais rapidement un petit alignement des séquences
#'
data %>%
  ggplot(data = ., aes(x = pos, y  = query )) +
  ## geom_text(aes(label = alt.base)) +
  geom_point(aes(color = alt.base), alpha = 0.4)

data %>%
  group_by(query) %>%
  summarise(end.tc = max(pos)) %>%
  ggplot(aes(x = end.tc)) +
  geom_histogram()

data %>% select(query) %>% unique()

#' Il n'y a que 127 séquences analysées dans ce tableau. Ça ne va pas.
#' On devrait s'attendre à 192 séquences normalement. Il manque du
#' signal. Je vais faire l'analyse moi-même.
