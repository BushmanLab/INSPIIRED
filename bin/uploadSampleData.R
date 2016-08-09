# uploadSampleData.R [path to SQLite db]  [path to csv data file]

library('RSQLite')
options(useFancyQuotes = F)

args <- commandArgs(trailingOnly = TRUE)
dbConn <- dbConnect(RSQLite::SQLite(), args[1])

dataTable <- read.table(args[2], sep=',', header=T, stringsAsFactors=F, strip.white=T)
data <- apply(dataTable, 1, function(x) split(unname(unlist(x)), names(x)))

for (row in data)
{
   rowValues <- unlist(row)
 
   for (i in 1:length(rowValues)){
      if ((grepl('[A-Za-z]', rowValues[i]) || grepl('\\d\\d\\d\\d\\-\\d\\d\\-\\d\\d', rowValues[i])) && rowValues[i] != 'NULL') rowValues[i] <- dQuote(rowValues[i]) 
   }

   fields <- gsub('([A-Za-z]+)\\.([A-Za-z]+)', '\\1 \\2', names(row)) 
   fields <- paste(dQuote(fields), collapse=", ")  
   values <- paste(rowValues, collapse=", ")

   command <- paste0('INSERT INTO "gtsp" (', fields, ') VALUES (', values, ')')
   message(command)
   r <- dbSendQuery(dbConn, command)
}

dbClearResult(dbListResults(dbConn)[[1]])
q()
