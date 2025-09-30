# Clare Lewis (s2879721), Grace  (s), Luke  (s)
# brief description of each person's contribution and the % they did

##setwd("/Users/clarelewis/Documents/github/ext_stat_prog_group_work")
#setwd("C:/Users/Grace Sheahan/Documents/Extended Statistical Programming/Assessment 1") ## Setting working directory to folder for Assessment 1

a <- scan("shakespeare.txt",what="character",skip=83,nlines=196043-83,fileEncoding="UTF-8") # Scanning the complete works of shakespear in under the name 'a'

#### Look at the Shakespeare text ####

#### 4-a ####

# Removing stage directions (indicated by'[_' and '_]')

open_direction <- grep('\\[_', a)  # Identifying opening index of stage directions by the occurrence of [_ 
direction <- length(open_direction) # The total number of stage directions identified

direction_length <- rep(0, direction) 
for(i in 1:direction) {direction_length[i] <- grep("\\_]|\\.]|’]", a[open_direction[i]:(open_direction[i]+100)])[1]} # Vector defines the length of directions using the position of the next closed bracket relative to the corresponding opening.

unclosed_direction <- which(is.na(direction_length)) # Locating any directions for which a close was not found

# unclosed_direction

# a[open_direction[unclosed_direction]:(open_direction[unclosed_direction]+20)] # Examined this to identify the causes of the issue - Is not a direction

open_direction <- open_direction[-unclosed_direction[1]] # Removing this location as there is not actually a direction there (per manual check)
direction_length <- direction_length[-unclosed_direction[1]] # Removing the corresponding (NA) length to the removed opening location
direction <- direction-1 # Accounting for the removal of the index of the inaccurately located direction

close_direction <- rep(0, direction)  
for(i in 1:direction){close_direction[i] <- (open_direction[i] + direction_length[i] -1)}  # The close location is the length of the direction after the direction opening (minus 1 as the length counts both the starting and ending values)

direction_words <- rep(0, sum(direction_length)) 
for(i in 1:direction){direction_words[(sum(direction_length[0:(i-1)])+1):(sum(direction_length[0:i]))] <- (open_direction[i]:close_direction[i])} # Contains the indices of all words which are stage directions (all words between open and close of a direction inclusive)

a <- a[-direction_words] # Removes all direction words from a

a <- gsub('\\[_', '', a) # Removing the [_ which was identified as an error above using manual checking


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

a_punc_under = grep("_", a, fixed = T )

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
# 
#### 7 - Old Code ####
# 
# key <- c("never", "die", ",", "but")
# 
# key_tokens <- rep(0,length(key))
# 
# for (i in 1:length(key)) {
#   key_tokens[i] <- grep(paste0("^", key[i], "$"),b)
# }
# full_sentence <- key_tokens
# 
# repeat {
#   for (i in 1:mlag) {
#     matching_row_index <- which(apply(M[, i:mlag, drop = FALSE], 1, function(row) all(row == key_token)))
#     if (length(matching_row_index) > 0) break
#   }
#   next_word <- M[matching_row_index,mlag+1]
#   key_tokens <- c(key_tokens[2:mlag],sample(na.omit(next_word),1))
#   full_sentence <- append(full_sentence,key_tokens[mlag])
#   
#   if (key_tokens[mlag] == match(".",b))
#     break
# }
# output <- cat(b[full_sentence], sep=" ")
# output
# 
# 
# 
# #### this is the instructions and code provided (we didnt use this code yet)
# 
# next.word <- function(key,M,M1,w=rep(1,ncol(M)-1)) {
#   
# }
# where key is the word sequence for which the next word is to be generated, M is as defined above, M1 is
# the vector of word tokens for the whole text and w is the vector of mixture weights (which actually don’t
#                                                                                      need to be normalized).
# The function should return a token for the next word, generated according to the
# model described above. It should be able to deal appropriately with any length of key: using reduced order
# versions of the model for short keys, and using only data from the end of key if it is too long
# 
# The crucial part of the function is (repeatedly) finding the rows of M that match key (or its reduced versions).
# Suppose the current key is to be matched to columns mc:mlag of M. Now compute
# 
# ii <- colSums(!(t(M[,mc:mlag,drop=FALSE])==key))
# If ii[j]=0 and is finite (see ?is.finite) then row j of M contains a match
# 
# 



#### Section 7#####

mlag <- 4   # Setting mlag (temporarily)
# 
M <- matrix(NA,length(a)-mlag,mlag+1)
for (i in 1:(length(a)-mlag)) {
  for (j in 1:(mlag+1)) {
    M[i,j] <- token[i + j - 1]
  }
}

M1 <- token

key <-c()

next.word <- function(key,M,M1,w=rep(1,ncol(M)-1)) { ## Returns NA for key = c()??
  if(length(key) < mlag){key <- append(key, rep("", (mlag-length(key))), 0)}  # Ensures the key vector is at least as long as mlag before beginning - inserts blanks at the start of the vector if this is not already the case
  weights = c()
  next_token = c()
  
  for(i in 1:mlag){  #Starting with 1, as 1:mlag will check for matches of mlag entried. the for loop will reduce matching requirements by 1 (2:mlag) if no mlag length matches found, and will continue to reduce by 1 until a match is found
    token_key <- match(key[(length(key)-mlag+i):length(key)],b)  # Defines a vector of length mlag of the tokens of the words in the key vector. If the key is longer than mlag it only takes the last mlag entries
    ii <- colSums(!(t(M[,i:mlag,drop=FALSE]) == token_key))  ## Code provided by SW; this will be equal to zero only if rows match completely (!() means that matches will be marked FALSE, and vice versa. As sych when summing matches will be adding 0, so any columns without full matches will be non-zero (and non-finite if NAs))
    
    full_matches <- which(is.finite(ii) & ii == 0) # Returns the indices of the rows for which there is no NAs (is.finite()) and there is an exact match (== 0) 
    
    next_word_tokens <- na.omit(M[full_matches,mlag+1]) # Creates a vector of the tokens for all the words which are preceded by an exact match 
    weights <- append(weights, rep(w[i],length(next_word_tokens)))
    next_token <- append(next_token, next_word_tokens)
  }
  if(length(next_token)==0){
    next_word <- b[sample(na.omit(M1),1)]
    # Checks if there are tokens in nwt; if this is the case we will sample from one of these as these are the best match
    # sampled_word <- next_word_tokens[na.omit(sample(length(next_word_tokens)),1)] # Taking a sample token from the vector. We do this by sampling the length of the vector and taking the corresponding index of the vector - this works even when the vector only has a single entry (whereas a direct sample would take a sample from 1:(token number) instead of the entry itself)
    # next_word <- b[sampled_word] # Finding the word in b which corresponds to the sample word
    
    return(next_word) # Once the next word has been defined, the function prints this and stops
  }
  else{
    next_word <- b[sample(next_token, size = 1, prob = weights)]
    return(next_word)
  }
}

# next.word(key_token,M,M1,w = c(1000,0,0,1))

#### 8 ####

starter.word.token <- function(start_word){ # Defining a starter word token function with a parameter start_word
  punc_vec <- append(punc_vec, c('?', '.')) # Adding ? and . to punc_vec, as the current versions contained in this vector (\\. and \\?) will not be recognized as matches to themselves in b
  punc_tokens <- na.omit(b[punc_vec]) # Finding the tokens of any punctuation in b
  
  if (missing(start_word)) {  ## Checking if a starter word is specified, and begins defining action in the case it is not
    start_token <- sample(na.omit(M1[! M1 %in% punc_tokens]), 1) #  A starter token is randomly sampled from the text token vector M1 (excluding punctuation)
#    start_word <- b[start_token]  # Defines the word corresponding to the sampled token
    return(start_token)
  }
  start_word <- tolower(start_word) # Converts to lower case to ensure match is found in b if possible
  
  if (start_word %in% punc_vec){ # Checks if the provided starter word is actually punctuation
    print("Error: This is punctuation, not a valid start word; A random token has been generated instead.") # Error message and notification to user of randomly sampled replacement
    start_token <- sample(na.omit(M1[! M1 %in% punc_tokens]), 1) # Randomly samples word_token as provided starter not suitable
  }
  
  else if (!start_word %in% b){
    print("Error: This is not a recognized common word; A random token has been generated instead.") # Error message and notification to user of randomly sampled replacement
    start_token <- sample(na.omit(M1[! M1 %in% punc_tokens]), 1) # Randomly samples word_token as provided starter word is not tokenized
  }
  
  else {start_token <- match(start_word, b)} # Finds corresponding token if the starter word is appropriate (word in b)
  return(start_token)
}

#### 9 ####

sentence <- c()

#punc_tokens <- na.omit(b[punc_vec])

generate.sentence <- function(sentence_prompt){ # Defining a generate sentence function with a parameter sentence prompt
  
  if(missing(sentence_prompt)){sentence_prompt <- b[starter.word.token()] # Accounts for no prompt being provided by generating a starter word (using starter.word.token() to generate a token and taking the corresponding b value)
  }
  
  sentence_prompt <- unlist(strsplit(sentence_prompt, split = " ")) # If the sentence prompt is several words long, it is redefined as a vector split into its individual words
  
  repeat{
    next_word <- next.word(sentence_prompt, M, M1, w=c(1000, 100, 10, 1)) # Runs next_word, to find the appropriate next word
    if(next_word %in% punc_vec & sentence_prompt[length(sentence_prompt)] %in% punc_vec){ # Checks if the last and next word are both punctuation, if so the next word is discarded
    } else
    {sentence_prompt <- append(sentence_prompt, next_word)} # The next word is added to the end of the sentence prompt vector
    
    if(sentence_prompt[length(sentence_prompt)] == '.'){ # Checks if the final value in the sentence prompt is a full stop
      break # Ends the repeat as the sentence has ended
    }
  }
  full_sentence <- paste(sentence_prompt, collapse = " ") # Concatenates the vector of sentence words into a single string, separating the words with a space
  return(full_sentence) # Prints the full sentence
}