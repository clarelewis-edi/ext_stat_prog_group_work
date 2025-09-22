# Clare Lewis (s2879721), Grace  (s), Luke  (s)
# brief description of each person's contribution and the % they did

setwd("C:\\Users\\Luke Egan\\Desktop\\Extended Statistical Programming\\ext_stat_prog_group_work") ## comment out of submitted
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
