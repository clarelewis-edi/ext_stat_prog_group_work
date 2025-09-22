# Clare Lewis (s2879721), Grace  (s), Luke  (s)
# brief description of each person's contribution and the % they did

setwd("/Users/clarelewis/Documents/github/ext_stat_prog_group_work") ## comment out of submitted
a <- scan("shakespeare.txt",what="character",skip=83,nlines=196043-83,
          fileEncoding="UTF-8")

# 4B

# look for all uppercase words (excluding 'I', 'O', and 'A')
# there are some instances of these single letter words with punctuation, we want to include these

# WE WANT TO DO PUNCTUATION REMOVAL BEFORE THIS STEP TO INCLUDE INSTANCES OF 'I,' AS AN EXAMPLE
uppercase_indices <- which(a == toupper(a) & a != 'A' & a != 'I' & a != 'O')# include all of the options after we remove punctuation
non_uppers <- a[-uppercase_indices]

# return to this after we remove the punctuation
# uppers <- unique(uppers[grep('I',uppers,fixed=TRUE)])
