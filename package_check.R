#!/usr/bin/env Rscript
# update.packages(ask = FALSE, checkBuilt = TRUE)

c("dplyr", "plyr", "ggplot2", "data.table", "devtools") -> chk_pkgs

suppressPackageStartupMessages(
  sapply(chk_pkgs, require, character.only=TRUE, quietly=FALSE, warn.conflicts=FALSE)
) -> ret

missing_pkgs <- sort(names(which(ret == FALSE)))

if (length(missing_pkgs) > 0) {
  warning("The following packages are not installed: %s", 
          paste0(sprintf("  - %s", missing_pkgs), collapse="\n"))
}

quit(save=FALSE, status=length(names) == 0, runLast = FALSE)