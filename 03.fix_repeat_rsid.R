library(tidyverse)
library(vroom)

setwd("/project/richards/guillaume.butler-laporte/bin/hisat-folder/hisat_index/")

vars<-vroom("snp151Common.txt.gz", col_names=FALSE, col_types = cols(.default = "c"))

vars_fixed<-vars %>%
  filter(str_detect(X2, "^chr[0-9XYM]*$") | str_detect(X2, "^chr6_GL000")) %>%
  group_by(X5) %>%
  mutate(instance=1:n()) %>%
  mutate(X5=paste0(X5, ".", instance)) %>%
  ungroup() %>%
  dplyr::select(-instance)

vroom_write(vars_fixed, "snp151Common_fixed.txt.gz", col_names=FALSE)

