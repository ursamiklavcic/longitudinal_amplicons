# Import csv into R matrix 

otutab = read_csv('data/csv_files/otutab_filt.csv') %>%
  column_to_rownames('Group') %>%
  as.matrix()


