# Clare Lewis (s2879721), Grace  (s), Luke  (s)
# brief description of each person's contribution and the % they did

setwd("/Users/clarelewis/Documents/github/ext_stat_prog_group_work")
a <- scan("shakespeare.txt",what="character",skip=83,nlines=196043-83,
          fileEncoding="UTF-8")

#### Look at the shakespeare text ####

length(a)


#### 4(c) ####

# View all words that contain a hyphen or an underscore in "a" 
a_punc_hyph = grep("-", a, fixed = T )
a[a_punc_hyph]
length(a_punc_hyph)

a_punc_under = grep("_", a, fixed = T )
a[a_punc_under]
length(a_punc_under)

# Remove all cases of hyphens and underscores in each word in the text
a_c <- gsub(pattern = "[_-]", replacement = "", a)

#### NOTE: Consider splitting hypthenated words into two seperate words. May yield better predictive power and more sensensical prompts




#### 4(d) ####
=======
# 4B

# look for all uppercase words (excluding 'I', 'O', and 'A')
# there are some instances of these single letter words with punctuation, we want to include these

# WE WANT TO DO PUNCTUATION REMOVAL BEFORE THIS STEP TO INCLUDE INSTANCES OF 'I,' AS AN EXAMPLE
uppercase_indices <- which(a == toupper(a) & a != 'A' & a != 'I' & a != 'O')# include all of the options after we remove punctuation
non_uppers <- a[-uppercase_indices]

# return to this after we remove the punctuation
# uppers <- unique(uppers[grep('I',uppers,fixed=TRUE)])
>>>>>>> ba76274733b876dcbf8990ee7a2f85585df46eeb
