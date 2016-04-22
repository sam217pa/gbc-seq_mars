#' La fonction suivante permet de créer les alignement via ggplot2. Les clones
#' sont triés par longueur de trace de conversion par la fonction embarquée
#' sort_by_tract_length.
plot_align <- function(data, mutant_) {
                                        # FIXME : why NA in seq ?
    sort_by_tract_length <- function(mut) {
        ## combine les données pour les clones avec conversion et les clones
        ## sans, pour lesquels on ne peut pas déterminer de longueur de trace de
        ## conversion.
        rbind(
            ## ceux là convertissent, donc …
            data %>%  filter(mutant == mut) %>%
              group_by(name, mutant) %>%
            ## on peut filtrer sur le dernier point de conversion
              filter(cons == "x") %>%
              summarise(len = min(refp)) %>%
              ungroup(),
            ## mais ceux là non.
            data %>% filter(mutant == mut) %>%
              group_by(name, mutant) %>%
            ## donc on filtre sur ceux qui ne basculent pas.
              filter(!("x" %in% cons)) %>%
              summarise(len = max(refp)) %>%
              ungroup()
        ) %>%
            ## trie par longueur de trace
            arrange(len) %>%
            ## convertit en vecteur
            select(name) %>% unlist() %>%  as.vector()
    }

    is_inside_conv <- function(data, mut) {
        ## détermine si une position est après ou avant le point de bascule.
        data %>%
            filter(mutant %in% mut, cons == "x") %>%
            group_by(name) %>%
            summarise(bascul = min(refp)) %>%
            inner_join(data, ., by = "name") %>%
            ## garde les positions après le point de switch qui correspondent au
            ## génotype parental.
            filter(refp > bascul, cons == "X" | cons == "x") %>%
            mutate(tract = TRUE) %>%
            select(name, bascul, refp, tract) %>%
            left_join(data, .)
    }

    is_restoration <- function(data, mut) {
        ## détermine si le SNP est une restauration du génotype parental. Autrement
        ## dit les SNP qui sont dans la trace de conversion mais qui sont X, *ie*
        ## des SNPs parentaux.
        data %>% filter(mutant %in% mut, cons == "X", tract == TRUE)
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
        mutate(name = factor(name, levels = sort_by_tract_length(mutant_))) %>%
        keep_clean_only()

    data$tract[is.na(data$tract)] <- FALSE

    first_snp <- filter(snp, cons == "X") %>% summarise(min = min(refp))
    last_snp <- filter(snp, cons == "x") %>% summarise(max = max(refp))

    ggplot(data = data, aes(x = refp, y = name)) +
        geom_point(aes(color = tract, size = qual, alpha=qual)) +
        ## représente les cas complexes
        geom_point(data = is_restoration(data, mutant_),
                   aes(alpha = qual), color = green, size = 5) +
        ## représente la séquence du donneur
        geom_text(aes(label = snpb, x = refp, y = -3), color = blue,
                  vjust = -0.5, size = 2, family = "Ubuntu Light") +
        annotate("text", x = 1.2, y = -0.8, label = "Donneur : ", color = blue, size = 2 ) +
        ## représente la séquence du receveur.
        geom_text(aes(label = refb, x = refp,
                      y = length(sort_by_tract_length(mutant_)) + 10),
                  color = red, vjust = 9, size = 2, family = "Ubuntu Light") +
        ## j'ai tenté de représenter les séquences par des segments de couleur,
        ## mais je trouve que ça perturbe la lecture. On a un effet moiré pas du
        ## tout attendu, qui n'apporte rien.
        ##
        ## geom_segment(data =
        ##                  filter(data, cons == "x") %>%
        ##                  group_by(name) %>%
        ##                  summarise(max = max(refp), min = min(refp)),
        ##              aes(x = max, xend = min, y = name, yend = name),
        ##              color = blue, alpha = 0.2) +
        ## geom_segment(data =
        ##                  filter(data, cons == "X") %>%
        ##                  group_by(name) %>%
        ##                  summarise(max = max(refp), min = min(refp)),
        ##              aes(x = max, xend = min, y = name, yend = name),
        ##              color = red, alpha = 0.2) +
        annotate("text", x = 1.2, y = length(sort_by_tract_length(mutant_)) + 5,
                 label = "Receveur : ", color = red, size = 2 ) +
        scale_color_brewer(palette = "Set1", guide = "none") +
        scale_alpha( range=c(1/5, 0.8), guide=FALSE ) +
        scale_size(range = c(0.5, 2), breaks = c(10, 50),
                   labels = c("Faible", "Forte")) +
        scale_x_continuous(breaks = extended_range_breaks()(data$refp)) +
        coord_cartesian(xlim = c(first_snp, last_snp)) +
        labs(x = "", y = "", size = "Qualité",
             title = paste("Alignement pour la manip", toupper(mutant_))) +
        theme(legend.direction = "horizontal",
              legend.margin = unit(0,"lines"),
              panel.grid.major.y = element_line(size = 0.1, linetype = "dotted")) +
        legend_position(0.8, 0.98)
}

#' Représente la distribution de la qualité par position.
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

#' Représente l'alignement superposé avec la qualité des bases.
plot_align_qual <- function(data, mutant_) {
    plot_grid(
        data %>% plot_align(mutant_),
        data %>% plot_qual(mutant_),
        ncol = 1, rel_heights = c(8.5, 1.5),
        align = 'v', labels = c("  1", "  2")
    )
}
