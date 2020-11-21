chemical <-
  readr::read_table(
    "data-raw/poisox.txt",
    col_names = c("age", "sex", "prior", "poisox", "after", "mortal")
  ) %>%
  dplyr::mutate(
    sex = factor(sex),
    poisox = factor(poisox),
    blood = after - prior,
    mortal = factor(mortal)
  ) %>%
  dplyr::select(-after, -prior)

usethis::use_data(chemical, overwrite = TRUE)
