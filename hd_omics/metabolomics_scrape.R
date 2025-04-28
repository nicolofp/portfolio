url = "https://www.metabolomicsworkbench.org/data/DRCCStudySummary.php?Mode=StudySummary&SortBy=release_date&AscDesc=desc&ResultsPerPage=100"

library(rvest)
library(xml2)
library(dplyr)
library(data.table)
library(GEOquery)

met_scrape = lapply(1:34, function(i){
  print(paste0("Page ",i))
  url = paste0("https://www.metabolomicsworkbench.org/data/DRCCStudySummary.php?",
               "Mode=StudySummary&SortBy=release_date&AscDesc=desc&ResultsPerPage=100&",
               "CurrentResultsPageNum=",i,"&CurrentStartResultsNum=",(i-1)*100+1,
               "&CurrentEndResultsNum=",i*100)
  page = read_html(url)
  tab = page %>% 
    html_element(xpath = "/html/body/div[2]/div/div[1]/div/div/table[2]") %>% 
    html_table(fill = T) %>%
    as.data.table()
  return(tab)
})
met_scrape = rbindlist(met_scrape)


View(met_scrape[Species == "Homo sapiens"][order(Samples,na.last = T)])
