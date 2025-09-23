# Clare Lewis (s2879721), Grace  (s), Luke  (s)
# brief description of each person's contribution and the % they did

##setwd("/Users/clarelewis/Documents/github/ext_stat_prog_group_work")
a <- scan("shakespeare.txt",what="character",skip=83,nlines=196043-83,
          fileEncoding="UTF-8")

#### Look at the shakespeare text ####

#### 4-a ####
# Remove stage directions (indicated by'[_' and '_]')

open_direction <- grep('\\[', a)  
direction <- length(open_direction)

direction_length <- rep(0, direction) 
for(i in 1:direction) {direction_length[i] <- grep('\\]', a[open_direction[i]:(open_direction[i]+100)])[1]}  

unclosed_bracket <- which(is.na(direction_length))
direction_length[1338] <- 15
direction_length[2026] <- 0

close_direction <- rep(0, direction)  
for(i in 1:direction){close_direction[i] <- (open_direction[i] + direction_length[i] -1)}  

##which(duplicate[close_direction]) 

direction_words <- rep(0, sum(direction_length)) 
for(i in 1:direction){direction_words[(sum(direction_length[0:(i-1)])+1):(sum(direction_length[0:i]))] <- (open_direction[i]:close_direction[i])}

a <- a[-direction_words]

a <- gsub('\\[_', '', a)

#### 4(d) ####

split_punc = function(text_vec, punc_vec){
  
  # Correct the format of the punctuation vector to be used in the "grep" function
  punc_vec = paste(punc_vec, collapse = "|")
  
  # Find the indices of all words in the text vector that contain punctuation marks given in the punctuation vector
  i_punc = grep(punc_vec, text_vec)
  
  # Create a holding vector for the elements of the text vector seperated from their corresponding punctuation, such that the length of this vector is the length of text + number of punc cases
  # Note that this assumes that all words dont have two cases of punctuation in them. Sound assumption
  text_vec_extra = rep(0, length(i_punc)+length(text_vec))
  
  # Indices in the holding vector which can be used as a reference to put the seperated punctuation as new elements
  i_extra_elements = i_punc+1:length(i_punc)-1
  
  # Take out the last character of words with puncuation using the "substr" function and put them in their corresponding element in the new text vector
  # Note: This assumes that all punctuation marks appear at the end of words. Sound assumption for the text file we are working with
  text_vec_extra[i_extra_elements+1] = substr(text_vec[i_punc], nchar(text_vec[i_punc]), nchar(text_vec[i_punc]))
  
  # Create a copy of the text vector with all punctuation removed
  text_vec_no_punc = gsub(pattern = punc_vec, replacement = "", text_vec)
  
  # Puts the elements/words without punctuation in their corresponding position in the new word vector
  text_vec_extra[-i_extra_elements-1] = text_vec_no_punc
  
  # Return the new text vector with punctuation seperated from all elements and are noew considered there own unique element in the text vector
  return(text_vec_extra)
}

# 4B

# look for all uppercase words (excluding 'I', 'O', and 'A')
# there are some instances of these single letter words with punctuation, we want to include these

# WE WANT TO DO PUNCTUATION REMOVAL BEFORE THIS STEP TO INCLUDE INSTANCES OF 'I,' AS AN EXAMPLE
uppercase_indices <- which(a == toupper(a) & a != 'A' & a != 'I' & a != 'O')# include all of the options after we remove punctuation
a <- a[-uppercase_indices]

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

#### NOTE: Consider splitting hyphenated words into two separate words. May yield better predictive power and more sensensical prompts



#### 4(e) ####
punc_vec = c(",","\\.",";","!",":","\\?")
a = split_punc(a, punc_vec)

#### 4-f ####
#Make all the text lower case

a <- tolower(a)


#### 5 ####
#a  
#Created a vector of unique words in the text
words <- unique(a)

#b
# Indicates which location in the words vector each word in the full text corresponds to
index <- match(a, words)

#c
# Counts the occurences of each word
word_occurences <- tabulate(index)

#d
# Creates a vector of the 1000 most common words by using a for loop and rank to locate and order these words
b<- rep(0, 1000)
for(i in 1:1000){b[i] <- words[which((rank(word_occurences, ties.method = "first")) == (length(word_occurences) + 1 - i))]}

#### 6 ####
# a
mlag <- 4
token <- match(a,b)
token

# b
#### test
d <- c(1:20)
R <- matrix(NA,length(d)-mlag,mlag+1)
for (i in 1:(length(d)-mlag)) {
  for (j in 1:(mlag+1)) {
    R[i,j] <- d[i + j - 1]
  }
}

## real
M <- matrix(NA,length(a)-mlag,mlag+1)
for (i in 1:(length(a)-mlag)) {
  for (j in 1:(mlag+1)) {
    M[i,j] <- token[i + j - 1]
  }
}

#### 7

key <- c("never", "die", ",", "but")

key_tokens <- rep(0,length(key))

for (i in 1:length(key)) {
  key_tokens[i] <- grep(paste0("^", key[i], "$"),b)
}
full_sentence <- key_tokens

repeat {
  for (i in 1:mlag) {
    matching_row_index <- which(apply(M[, i:mlag, drop = FALSE], 1, function(row) all(row == key_token)))
    if (length(matching_row_index) > 0) break
  }
  next_word <- M[matching_row_index,mlag+1]
  key_tokens <- c(key_tokens[2:mlag],sample(na.omit(next_word),1))
  full_sentence <- append(full_sentence,key_tokens[mlag])
  
  if (key_tokens[mlag] == match(".",b))
    break
}
output <- cat(b[full_sentence], sep=" ")
output



#### this is the instructions and code provided (we didnt use this code yet)

next.word <- function(key,M,M1,w=rep(1,ncol(M)-1)) {
  
}
where key is the word sequence for which the next word is to be generated, M is as defined above, M1 is
the vector of word tokens for the whole text and w is the vector of mixture weights (which actually don’t
                                                                                     need to be normalized).
The function should return a token for the next word, generated according to the
model described above. It should be able to deal appropriately with any length of key: using reduced order
versions of the model for short keys, and using only data from the end of key if it is too long

The crucial part of the function is (repeatedly) finding the rows of M that match key (or its reduced versions).
Suppose the current key is to be matched to columns mc:mlag of M. Now compute

ii <- colSums(!(t(M[,mc:mlag,drop=FALSE])==key))
If ii[j]=0 and is finite (see ?is.finite) then row j of M contains a match


