#source("code/genomics_functions.R")

count.df <- function(df, gcol, scol) {
  df %>% group_by(!!as.name(gcol)) %>% summarize(!!scol := n()) %>% as.data.frame()
}

first <- dplyr::first

index.r_join <- function(a, b) {
  a <- as.data.frame(a)
  if (!"index" %in% colnames(a)) {
    a$index <- 1:nrow(a)
  }
  
  b <- as.data.frame(b)
  
  right_join(a, b, by="index")
}

strsplit.n <- function(x, split, n=1) {
	unlist(purrr::map(strsplit(x, split), n))
}