# phruscle_call

This dir contains the R script and its outcomes as html to analyse the results
of the basecalling by phred and muscle, _ie_ `phruscle basecall` and `phrusce
table`.

The R script is turned into a beautiful tufte style formatted html page by :

```sh
R --vanilla -e "rmarkdown::render('phruscle_call.R')"
```
