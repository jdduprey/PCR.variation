
debug.rep.tally <- PCR.proportions %>%
  unite(bio, Hash, col=ID, sep='.', remove=TRUE) %>%
  group_by(ID) %>%
  add_tally()

table(debug.rep.tally$n)


debug.long <-by.tech.table %>%
  unite(bio, Hash, col=ID, sep='.', remove=TRUE) %>%
  group_by(ID) %>% 
  add_count()

table(debug.long$n)
