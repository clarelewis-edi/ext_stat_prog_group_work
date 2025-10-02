# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- 
# Clare Lewis (s2879721), Grace Sheahan (s2898645), Luke Egan (s2837709)
# brief description of each person's contribution and the % they did
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- 

# setwd("/Users/clarelewis/Documents/github/ext_stat_prog_group_work")
# setwd("C:/Users/Grace Sheahan/Documents/Extended Statistical Programming/Assessment 1") ## Setting working directory to folder for Assessment 1

# --- Introduction ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ 

# This project aims to produce a rudimentary 'small language model', with which to simulate sentences written by Shakespeare.
# The model has been fed a copy of the complete works of Shakespeare; which it uses to find and return a likely follow up to the latest word based on a prompt or random starting point.
# In mathematical terms, the code uses a Markov model and randomly samples next words based on weighted probabilities of all potential follow ups.
 
# The provided code includes steps taken to 'clean' the text until it is fit for purpose, and the prior definitions of several variables and functions which are later used in the final function.
# The culmination of the code is a function 'generate.sentence()', which meets the objective of simulating a Shakespearean sentence, and does so with optional user prompts.
 
# The complete works of Shakespeare used to build this model were taken from: https://www.gutenberg.org/cache/epub/100/pg100.txt

# --- Section 1 -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# We created a GitHub repository for our group project at the following location; https://github.com/clarelewis-edi/ext_stat_prog_group_work through which we could all collaborate on the project file

# --- Section 2 -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# The text file was saved from the provided location. A copy of this file was added to our repository.

# --- Section 3 -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Change the path to point our local repositories, and read in the text file under the name 'a'

# setwd("/Users/clarelewis/Documents/github/ext_stat_prog_group_work")
# setwd("C:/Users/Grace Sheahan/Documents/Extended Statistical Programming/Assessment 1") 
# setwd("C:/Users/Luke Egan/Desktop/Extended Statistical Programming")

a <- scan("shakespeare.txt",what="character",skip=83,nlines=196043-83,fileEncoding="UTF-8") # The conditions limit to only Shakespeare's works in the file, cutting out unnecessary information

# --- Section 4 -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Complete the necessary pre-processing steps to clean and correctly format the text before building of the model could begin

# a) Remove stage directions (indicated by'[_' and '_]')
open_direction <- grep('\\[_', a)  # The locations of the start of each direction (indicated by '[_')
direction <- length(open_direction) # Number of directions

direction_length <- rep(0, direction) 
# This loop finds the next close of direction corresponding to each opening and assigns the distance to this as the length of the direction i
for (i in 1:direction) {
  direction_length[i] <- grep('\\_]|\\.]', a[open_direction[i]:(open_direction[i] + 100)])[1] # '.]' was also included in the search as there were a number of cases in which this is used instead of '_]'
}

unclosed_direction <- which(is.na(direction_length)) # Check for any errors or missing data in direction_length

# The unclosed direction vector, and subsequently the corresponding directions were examined manually to identify the issues
# (The following lines are commented out to avoid printing of unnecessary information when running the full code)
# unclosed_direction

# a[open_direction[unclosed_direction[1]]:(open_direction[unclosed_direction[1]+20)]
# a[open_direction[unclosed_direction[2]]:(open_direction[unclosed_direction[2]+20)]

direction_length[unclosed_direction[2]] <- 10 # Found by manual inspection

# The first input of unclosed_direction was found to be an error usage of '[_' and not actually a direction.
# The input/value corresponding to this error is removed from open_direction, direction_length, and direction
open_direction <- open_direction[-unclosed_direction[1]]
direction_length <- direction_length[-unclosed_direction[1]]
direction <- direction - 1


close_direction <- rep(0, direction)  
for (i in 1:direction) {
  close_direction[i] <- (open_direction[i] + direction_length[i] - 1)
}  

direction_words <- rep(0, sum(direction_length)) 
# This loop assigns the location of the close of the corresponding directions.
for (i in 1:direction) {
  direction_words[(sum(direction_length[0:(i - 1)]) + 1):(sum(direction_length[0:i]))] <- (open_direction[i]:close_direction[i]) # The '-1' accounts for the fact that the length includes the open_location
}

a <- a[-direction_words] # Redefine the text to remove all directions based on the locations of their words

a <- gsub('\\[_', '', a) # Remove the '[_' from the erroneously indicated direction (as found from the unclosed_bracket inspection)

# b) Look for and remove cases where the entire word is uppercase (excluding 'I', 'O', and 'A' because these are words themselves)
# Note: This will remove cases of single letter words immediately preceded or followed by punctuation (ex. 'I,'), however the occurrence of these cases were not significant enough to impact the subsequent 'b' vector

uppercase_indices <- which(a == toupper(a) & !(a %in% c('A', 'I', 'O'))) 
a <- a[-uppercase_indices]


# c) Find and remove all cases of hyphens and underscores in each word in the text
# Note: We explored the option of keeping hyphenated words as is, however, none of these words appeared in the vector 'b'
# Note: We chose to remove the hyphen rather than split into two words because we felt this kept the contextual integrity of the words in the text

a <- gsub(pattern = '_|-', replacement = '', a) 

# d) Create a function that will split all punctuation into their own entries
# The inputs for this function: text_vec, punc_vec
# - text_vec: a vector containing the text strings from which punctuation will be split
# - punc_vec: a vector containing all of the potential punctuation marks that may be attached to words in the text vector
# The output for this function: text_vec_split
# - text_vec_split: text_vec with the occurrences of punctuation from punc_vec split into separate elements

split_punc <- function(text_vec, punc_vec) {
  punc_vec <- paste(punc_vec, collapse = '|') # Correct the format of the punctuation vector to be used in the "grep" function
  
  i_punc <- grep(punc_vec, text_vec) # Locate all the indices of punctuation occurrences in the text vector
  
  text_vec_split <- rep(0, length(i_punc) + length(text_vec)) # The length of the new vector will be that of the current vector plus the number of punctuation occurrences
  
  i_extra_elements <- i_punc + 1:length(i_punc) - 1 # Locations of the extra element from splitting the punctuation occurrences; these will be offset by 1 each time to account for the other punctuation having already been inserted ahead of it
  
  text_vec_split[i_extra_elements + 1] <- substr(text_vec[i_punc], nchar(text_vec[i_punc]), nchar(text_vec[i_punc])) # Take out the last character of words containing punctuation using the "substr" function and put them in their corresponding element in the new text vector
  
  text_vec_no_punc <- gsub(pattern = punc_vec, replacement = "", text_vec)
  
  text_vec_split[-i_extra_elements - 1] <- text_vec_no_punc
  
  return(text_vec_split)
}


# e) Implement the split_punc vector defined above to separate punctuation in 'a'
punc_vec <- c(",", "\\.", ";", "!", ":", "\\?")
a <- split_punc(a, punc_vec) # 'a' is being passed as the text_vec and punc_vec is defined above

# f) 
a <- tolower(a) # Convert all the text lower case

# --- Section 5 -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#### questionably all the comments in 5 can be removed

# a) Create a vector of unique words in the text

words <- unique(a)

# b) Indicate which location in the words vector each word in the full text corresponds to

index <- match(a, words)

# c) Count the occurrences of each word

word_occurences <- tabulate(index)

# d) Create a vector of the 1000 most common words by using a for loop and rank to locate and order these words

b <- rep(0, 1000)
for (i in 1:1000) {
  b[i] <- words[which((rank(word_occurences, ties.method = "first")) == (length(word_occurences) + 1 - i))]
}

# --- Section 6 -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#### is mlag just to define the longest text string that our functions will look at?
# mlag is a value that will determine the width of the matrix M, mlag can be any non-zero number, for this exercise we assigned it the value of 4, but that value can be changed
############### changing mlag produces an error
mlag <- 4 # This value can be changed by the user to change the output results of the model
M1 <- match(a,b) # M1 is a vector of length 'a' with token values representing the index of 'b' where that given word shows up, NA values exist everywhere that a word in 'a' does not exist in 'b'

# b) Create a matrix M that has dimensions a-mlag*mlag+1

M <- matrix(NA, length(a) - mlag, mlag + 1)

# Insert values into M such that row 1 contains values M1[1:mlag+1], row 2 contains M1[2:mlag+2], row 3 contains M1[3:mlag+3] etc.
for (i in 1:(length(a) - mlag)) {
  for (j in 1:(mlag + 1)) {
    M[i,j] <- M1[i + j - 1]
  }
}


# --- Section 7 -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# The inputs for this function: key, M, M1, w
# - key:
# - M:
# - M1:
# - w:
# The output for this function: next_word
# - next_word

#key <- c()
############### changing mlag produces an error in next.word, also running just next.word gives an error if key is not defined (cant remember how it's supposed to handle this)
next.word <- function(key, M, M1, w = rep(1, ncol(M) - 1)) { 
  if (length(key) < mlag) {
    key <- append(key, rep("", (mlag - length(key))), 0)
  }  # Ensures the key vector is at least as long as mlag before beginning - inserts blanks at the start of the vector if this is not already the case
  
  weights = c()
  next_token = c()
  
  for (i in 1:mlag) {
    #Starting with 1, as 1:mlag will check for matches of mlag entried. the for loop will reduce matching requirements by 1 (2:mlag) if no mlag length matches found, and will continue to reduce by 1 until a match is found
    token_key <- match(key[(length(key) - mlag + i):length(key)], b)  # Defines a vector of length mlag of the tokens of the words in the key vector. If the key is longer than mlag it only takes the last mlag entries
    ii <- colSums(!(t(M[ ,i:mlag, drop=FALSE]) == token_key))  ## Code provided by SW; this will be equal to zero only if rows match completely (!() means that matches will be marked FALSE, and vice versa. As sych when summing matches will be adding 0, so any columns without full matches will be non-zero (and non-finite if NAs))
    
    full_matches <- which(is.finite(ii) & ii == 0) # Returns the indices of the rows for which there is no NAs (is.finite()) and there is an exact match (== 0) 
    
    next_word_tokens <- na.omit(M[full_matches, mlag + 1]) # Creates a vector of the tokens for all the words which are preceded by an exact match 
    weights <- append(weights, rep(w[i]/length(next_word_tokens),length(next_word_tokens)))
    next_token <- append(next_token, next_word_tokens)
  }

  if (length(next_token) == 0) {
    next_word <- b[sample(na.omit(M1), 1)]
    # Checks if there are tokens in nwt; if this is the case we will sample from one of these as these are the best match
    # sampled_word <- next_word_tokens[na.omit(sample(length(next_word_tokens)),1)] # Taking a sample token from the vector. We do this by sampling the length of the vector and taking the corresponding index of the vector - this works even when the vector only has a single entry (whereas a direct sample would take a sample from 1:(token number) instead of the entry itself)
    # next_word <- b[sampled_word] # Finding the word in b which corresponds to the sample word
    return(next_word) # Once the next word has been defined, the function prints this and stops
  }
  else {
    next_word <- b[sample(next_token, size = 1, prob = weights)]
    return(next_word)
  }
}

# --- Section 8 -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
########## since we skip this function as long as there's not a missing value do we want to keep the print statements?

# The inputs for this function: start_word
# - start_word:
# The output for this function: start_token
# - next_word:

starter.word.token <- function(start_word) { ###Defining a starter word token function
  punc_vec <- append(punc_vec, c('?', '.'))
  punc_tokens <- na.omit(b[punc_vec])
  
  if (missing(start_word)) {  ## If a starter word is not specified a starter token is randomly sampled from the text (excluding punctuation)
    start_token <- sample(na.omit(M1[! M1 %in% punc_tokens]), 1)
    start_word <- b[start_token]
    return(start_token)
  }
  start_word <- tolower(start_word)
  
  if (start_word %in% punc_vec){
    print("Error: This is punctuation, not a valid start word; A random token has been generated instead.")
    start_token <- sample(na.omit(M1[! M1 %in% punc_tokens]), 1)
  }
  
  else if (!start_word %in% b){
    print("Error: This is not a recognized common word; A random token has been generated instead.")
    start_token <- sample(na.omit(M1[! M1 %in% punc_tokens]), 1)
  }
  
  else {start_token <- match(start_word, b)}
  return(start_token)
}

# --- Section 9 -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# The inputs for this function: sentence_prompt
# - sentence_prompt:
# The output for this function: full_sentence
# - full_sentence:

#punc_tokens <- na.omit(b[punc_vec]) do we need this line?

generate.sentence <- function(sentence_prompt) { # Defining a generate sentence function with a parameter sentence prompt
  
  if(missing(sentence_prompt)){
    sentence_prompt <- b[starter.word.token()] # Accounts for no prompt being provided by generating a starter word (using starter.word.token() to generate a token and taking the corresponding b value)
  }
  
  sentence_prompt <- unlist(strsplit(sentence_prompt, split = " ")) # If the sentence prompt is several words long, it is redefined as a vector split into its individual words
  
  repeat{
    next_word <- next.word(sentence_prompt, M, M1, w = rep(1, ncol(M) - 1)) # Runs next_word, to find the appropriate next word
    if (next_word %in% punc_vec & sentence_prompt[length(sentence_prompt)] %in% punc_vec) { # Checks if the last and next word are both punctuation, if so the next word is discarded
    } 
    
    else {
      sentence_prompt <- append(sentence_prompt, next_word)
    } # The next word is added to the end of the sentence prompt vector
    
    if (sentence_prompt[length(sentence_prompt)] == '.') { # Checks if the final value in the sentence prompt is a full stop
      break # Ends the repeat as the sentence has ended
    }
  }
  
  full_sentence <- paste(sentence_prompt, collapse = " ") # Concatenates the vector of sentence words into a single string, separating the words with a space
  return(full_sentence) # Prints the full sentence
}

generate.sentence()
generate.sentence('Blood')
generate.sentence('Simon')
generate.sentence('Today I went')
