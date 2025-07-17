library(dplyr)
library(stringr)

# take gff and expand to dataframe with attribute values as separate columns
gff2df <- function(fname, hasAlias=T) {
  # open fimo gff output and label columns
  res <- read.table(fname)
  colnames(res) <- c('seqname', 'source', 'feature', 'start', 'end',
                     'score', 'strand', 'frame', 'attributes')
  
  # parse attribute column to separate columns
  attrs <- c('Name', 'Alias', 'pvalue', 'qvalue', 'sequence')
  for (a in attrs) {
    res[[a]] <- str_extract(res$attributes, paste("(?<=",a,"=)[^;]+",sep=''))
  }
  
  # if alias flag set (for Hv13 results) set Alias to gene name from ID
  if (!hasAlias) {
    res$Alias <- str_extract(res$Name, "^[^.]+")
  } 
  
  # drop attribute columns and set numeric 
  res <- res %>% select(-attributes) %>% mutate_at(c('start', 'end', 'score', 'pvalue', 'qvalue'), as.numeric)
  
  # return df
  return(res)
}

# take dataframe and reformat attributes to one column for gff format
df2gff <- function(d, fname) {
  d$attributes = apply(d[, c('Name', 'Alias', 'pvalue', 'qvalue', 'sequence')], 1, function(row) {
    paste0(names(d)[9:13], "=", row, ";", collapse="")
  })
  d <- d %>% select(-c(Name, Alias, pvalue, qvalue, sequence))
  write.table(d, file=fname, sep='\t', col.name=F,
              row.names=F, quote=F)
} 
