data_location <- "../../data/phruscle_snpcall.csv"

snp <- read_csv(data_location, col_types = "ccciciccid", trim_ws = TRUE) %>%
    mutate(
        name = gsub("-1073.+$", "", name),
        base = toupper(base)#,
    )

find_mutant <- function(name) {
    if      (grepl("ws", name)) "ws"
    else if (grepl("sw", name)) "sw"
    else if (grepl("W", name )) "w"
    else if (grepl("S", name )) "s"
    else ""
}

# neat little trick to reduce time of rowwise application of find_mutant.
snp <- snp %>%
    group_by(name) %>%
    summarise(count = n()) %>%
    rowwise() %>%
    mutate(mutant=find_mutant(name)) %>%
    inner_join(snp, .) %>%
    filter(mutant != "")
