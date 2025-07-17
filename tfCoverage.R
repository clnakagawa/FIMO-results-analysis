library(data.table)
library(ggplot2)
source("C:/Users/cnakagawa/OneDrive - Cystic Fibrosis Foundation/Desktop/tfbs data/parseGff.R")

# read in files 
# files <- list.files(pattern = 'calu3_reg[12345678]_fimo.gff')
# files <- lapply(files, function(x) gff2df(x, hasAlias = F))
# allReg <- rbindlist(files)
# head(allReg)

# helper function for row merging
mergeRanges <- function(dt) {
  setorder(dt, start, end)  # Sort by start and end
  
  merged <- dt[1]  # Initialize with the first range
  for (i in 2:nrow(dt)) {
    last <- merged[.N]
    current <- dt[i]
    # If overlapping or contiguous
    if (current$start <= last$end) {
      merged[.N, end := max(last$end, current$end)]  # Extend the range
    } else {
      merged <- rbind(merged, current)  # Add new non-overlapping range
    }
  }
  return(merged)
}

# helper function for overlap in the table with given range
getOverlap <- function(s, e, rdt) {
  testOverlap <- rdt[end >= s & start <= e]
  testOverlap[, rCov := (abs(e - start) 
                         + abs(end - s) 
                         - abs(e - end)
                         - abs(start - s)) / 2]
  sum(testOverlap$rCov) / (e-s)
}

# function to repeat process for any data.table of fimo results gff 
gffCov <- function(gdt, qCut, bSize, startLim = 0, endLim = Inf, title = '') {
  # filter
  gdt <- gdt[qvalue < qCut]
  gdt <- gdt[start <= endLim & end >= startLim]
  
  # merge
  gdt <- mergeRanges(gdt)
  
  # full range
  if (startLim == 0) {
    frange <- c(min(gdt$start), max(gdt$end))
  }
  else {
    frange <- c(startLim, endLim)
  }
  
  
  # get average coverage
  gCov <- getOverlap(frange[1], frange[2], gdt)
  
  # get bins
  rBins <- data.table(
    pos = head(seq(frange[1], frange[2], by = bSize), -2) 
  )
  rBins[, binCov := as.numeric(lapply(pos, function(x) getOverlap(x, x+bSize, gdt)))]
  
  # plot bins
  plt <- ggplot(rBins, aes(pos, binCov)) + 
    geom_bar(stat = 'identity', just = 0, width = bSize, fill = "lightblue", color = "black") +
    scale_x_continuous(name = 'chr position', 
                       breaks = seq(frange[1], frange[2], by = bSize * 2),
                       expand = c(0, 2)) +
    scale_y_continuous(name = '% coverage',
                       labels = seq(0, 100, by = 20),
                       breaks = seq(0, 1, by = 0.2), 
                       expand = c(0, 0), 
                       limits = c(0, 1.05)) +
    geom_hline(yintercept = gCov, 
               color = "red", 
               linetype = "dotted", 
               size = 1) +
    theme(axis.ticks.x=element_blank(), axis.text.x=element_blank(),
          plot.title=element_text(size=7)) +
    ggtitle(title)
  
  # return plot object
  return(plt)
}
