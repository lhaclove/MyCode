this R script used for PAST analysize. Tang has published this method and R packages, but it need a lot of memory while dealing with large
dateset ( 770k snp for me). so i changed 1). doParallel to doSNOWsome (with has a progress bar), and 2). %dopar% to foreach, and 3). 
some rm() to remove the unused variable. and this 770k snp will cause about 60gb ram and will take about 20 hour to finish.
