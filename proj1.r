# Clare Lewis (s2879721), Grace  (s), Luke  (s)
# brief description of each person's contribution and the % they did

setwd("/Users/clarelewis/Documents/github/ext_stat_prog_group_work")
a <- scan("shakespeare.txt",what="character",skip=83,nlines=196043-83,
          fileEncoding="UTF-8")

#### Look at the shakespeare text ####

#### 4-a ####
# Remove stage directions (indicated by'[_' and '_]')

open_direction <- grep('\\[', a)  
direction <- length(open_direction)

direction_length <- rep(0, d) 
for(i in 1:d) {direction_length[i] <- grep('\\]', a[open_direction[i]:(open_direction[i]+100)])[1]}  

unclosed_bracket <- which(is.na(direction_length))
direction_length[1338] <- 15
direction_length[2026] <- 0

close_direction <- rep(0, d)  
for(i in 1:d){close_direction[i] <- (open_direction[i] + direction_length[i] -1)}  

##which(duplicate[close_direction]) 

direction_words <- rep(0, sum(direction_length)) 
for(i in 1:d){direction_words[(sum(direction_length[0:(i-1)])+1):(sum(direction_length[0:i]))] <- (open_direction[i]:close_direction[i])}

a <- a[-a[direction_words]]

a <- gsub('\\[_', '', a)

# 4B

# look for all uppercase words (excluding 'I', 'O', and 'A')
# there are some instances of these single letter words with punctuation, we want to include these

# WE WANT TO DO PUNCTUATION REMOVAL BEFORE THIS STEP TO INCLUDE INSTANCES OF 'I,' AS AN EXAMPLE
uppercase_indices <- which(a == toupper(a) & a != 'A' & a != 'I' & a != 'O')# include all of the options after we remove punctuation
non_uppers <- a[-uppercase_indices]

# return to this after we remove the punctuation
# uppers <- unique(uppers[grep('I',uppers,fixed=TRUE)])


#### 4(c) ####

# View all words that contain a hyphen or an underscore in "a" 
a_punc_hyph = grep("-", a, fixed = T )
a[a_punc_hyph]
length(a_punc_hyph)

a_punc_under = grep("_", a, fixed = T )
a[a_punc_under]
length(a_punc_under)

# Remove all cases of hyphens and underscores in each word in the text
a <- gsub(pattern = "[_-]", replacement = "", a)

#### NOTE: Consider splitting hypthenated words into two seperate words. May yield better predictive power and more sensensical prompts




#### 4(d) ####

split_punc = function(text_vec, punc_vec){
  
  i_punc = grep(punc_vec, text_vec, fixed = T)
  
}

#### 4-f ####
#Make all the text lower case

a <- tolower(a)


