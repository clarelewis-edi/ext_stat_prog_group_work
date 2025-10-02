# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- 
# Clare Lewis (s2879721), Grace Sheahan (s2898645), Luke Egan (s2837709)
# Clare: Completed pre-processing steps 4b and 6, we all collaborated on the later sections, working in person together but I lead the code writing for section 7
# Grace: Completed pre-processing steps 4a and 5, we all collaborated on the later sections, working in person together but I lead the code writing for section 9
# Luke: Completed pre-processing steps 4c-4f, we all collaborated on the later sections, working in person together but I lead the code writing for section 8
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- 

# Setting working directory to folder for Assessment 1
# setwd("/Users/clarelewis/Documents/github/ext_stat_prog_group_work")
# setwd("C:/Users/Grace Sheahan/Documents/Extended Statistical Programming/Assessment 1")
# setwd("C:/Users/Luke Egan/Desktop/Extended Statistical Programming")

# --- Introduction ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ 

# This project aims to produce a rudimentary 'small language model', with which to simulate sentences written by Shakespeare.
# The model has been fed a copy of the complete works of Shakespeare; which it uses to find and return a likely follow up to the latest word based on a prompt or random starting point.
# In mathematical terms, the code uses a Markov model and randomly samples next words based on weighted probabilities of all potential follow ups.
 
# The provided code includes steps taken to 'clean' the text until it is fit for purpose, and the defining of several variables and functions which are later used in the final function.
# The culmination of the code is a function 'generate.sentence()', which meets the objective of simulating a Shakespearean sentence, and does so with optional user prompts.
 
# The complete works of Shakespeare used to build this model were taken from: https://www.gutenberg.org/cache/epub/100/pg100.txt

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# We created a GitHub repository for our group project at the following location; https://github.com/clarelewis-edi/ext_stat_prog_group_work through which we could all collaborate on the project file

# The text file was saved from the provided location. A copy of this file was added to our repository.

# Change the path to point our local repositories, and read in the text file under the name 'a'

# setwd("/Users/clarelewis/Documents/github/ext_stat_prog_group_work")
# setwd("C:/Users/Grace Sheahan/Documents/Extended Statistical Programming/Assessment 1") 
# setwd("C:/Users/Luke Egan/Desktop/Extended Statistical Programming")

a <- scan("shakespeare.txt",what="character",skip=83,nlines=196043-83,fileEncoding="UTF-8") # The conditions limit to only Shakespeare's works in the file, cutting out unnecessary information

# --- Pre-processing -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# The necessary pre-processing steps to clean and correctly format the text before building of the model could begin

# Removing stage directions (indicated by'[' and ']')
open_direction <- grep('\\[', a)  # The locations of the start of each direction (indicated by '[')
direction <- length(open_direction) # Number of directions

direction_length <- rep(0, direction) 
# Finds the next close of direction corresponding to each opening and assigns the distance to this as the length of the direction i
for (i in 1:direction) {
  direction_length[i] <- grep('\\]', a[open_direction[i]:(open_direction[i] + 100)])[1] # The end of the direction is identified as the next occurrence of ']'
}

unclosed_direction <- which(is.na(direction_length)) # Checking for any errors or missing data in direction_length

# The unclosed direction vector, and subsequently the corresponding directions were examined manually to identify the issues
# (The following lines are commented out to avoid printing of unnecessary information when running the full code)
# unclosed_direction

# a[open_direction[unclosed_direction[1]]:(open_direction[unclosed_direction[1]+20)]
# a[open_direction[unclosed_direction[2]]:(open_direction[unclosed_direction[2]+20)]

direction_length[unclosed_direction[1]] <- 15 # Found by manual inspection

# The second input of unclosed_direction was found to be an error usage of '[_' and not actually a direction.
# The input/value corresponding to this error is removed from open_direction, direction_length, and direction
open_direction <- open_direction[-unclosed_direction[2]]
direction_length <- direction_length[-unclosed_direction[2]]
direction <- direction - 1


close_direction <- rep(0, direction)  
for (i in 1:direction) {
  close_direction[i] <- (open_direction[i] + direction_length[i] - 1)
}  

direction_words <- rep(0, sum(direction_length)) 
# Assigns the location of the close of the corresponding directions.
for (i in 1:direction) {
  direction_words[(sum(direction_length[0:(i - 1)]) + 1):(sum(direction_length[0:i]))] <- (open_direction[i]:close_direction[i]) # The '-1' accounts for the fact that the length includes the open_location
}

a <- a[-direction_words] # Redefining the text to remove all directions based on the locations of their words

a <- gsub('\\[_', '', a) # Removed the '[_' from the erroneously indicated direction (as found from the unclosed_bracket inspection)

# Look for and remove cases where the entire word is uppercase (excluding 'I', 'O', and 'A' because these are words themselves)
# Note: This will remove cases of single letter words immediately preceded or followed by punctuation (ex. 'I,'), however the occurrence of these cases were not significant enough to impact the subsequent 'b' vector
uppercase_indices <- which(a == toupper(a) & !(a %in% c('A', 'I', 'O'))) 
a <- a[-uppercase_indices]


# Finds and removes all cases of hyphens and underscores in each word in the text
# Note: We explored the option of keeping hyphenated words as is, however, none of these words appeared in the vector 'b'
# Note: We chose to remove the hyphen rather than split into two words because we felt this kept the contextual integrity of the words in the text
a <- gsub(pattern = '_|-', replacement = '', a) 

# This creates a function that will split all punctuation into their own entries
# The inputs for this function: text_vec, punc_vec
# - text_vec: a vector containing the text strings from which punctuation will be split
# - punc_vec: a vector containing all of the potential punctuation marks that may be attached to words in the text vector
# The output for this function: text_vec_split
# - text_vec_split: text_vec with the occurrences of punctuation from punc_vec split into separate elements

split_punc <- function(text_vec, punc_vec) {
  punc_vec <- paste(punc_vec, collapse = '|') # Correct the format of the punctuation vector to be used in the "grep" function
  i_punc <- grep(punc_vec, text_vec) # Locate all the indices of punctuation occurrences in the text vector
  
  text_vec_split <- rep(0, length(i_punc) + length(text_vec)) # The length of the new vector will be that of the current vector plus the number of punctuation occurrences
  
  i_extra_elements <- i_punc + 1:length(i_punc) - 1 # Locations of the extra element from splitting the punctuation occurrences; these will be offset by 1 each time to account for the other punctuation having already been split ahead of it
  text_vec_split[i_extra_elements + 1] <- substr(text_vec[i_punc], nchar(text_vec[i_punc]), nchar(text_vec[i_punc])) # Take out the last character of words containing punctuation and put them in their corresponding element in the new text vector
  
  text_vec_no_punc <- gsub(pattern = punc_vec, replacement = "", text_vec)
  text_vec_split[-i_extra_elements - 1] <- text_vec_no_punc
  
  return(text_vec_split)
}


# Implement the split_punc vector defined above to separate punctuation in 'a'
punc_vec <- c(",", "\\.", ";", "!", ":", "\\?", "—") # The last item in this list is an emdash, not a hyphen which is why it is included with the punctuation
a <- split_punc(a, punc_vec) # 'a' is being passed as the text_vec and punc_vec is defined above
# Due to instances of punctuation ending a quote (ex. 'this is a quote.') the end quotation mark is being split and considered a word.
# This occurs enough times for this "word" to show up in 'b', we have decided to manually remove this "word" from 'a' to help our model
a <- gsub("’", "", a)

a <- tolower(a) # Convert all the text lower case
words <- unique(a) # Unique words in the text
index <- match(a, words) # Location in the words vector each word in the full text corresponds to
word_occurences <- tabulate(index) # Number of occurrences of each word

# Create a vector of the 1000 most common words by using a for loop and rank to locate and order these words
b <- rep(0, 1000)
for (i in 1:1000) {
  b[i] <- words[which((rank(word_occurences, ties.method = "first")) == (length(word_occurences) + 1 - i))]
}

# mlag determines the maximum number of words our model will use to predict the next word, this can be any non-zero positive integer, reasonably this is restricted by the length of your text vector
# mlag will determine the width of the matrix M, we set it to 4, but that value can be changed
mlag <- 4
M1 <- match(a,b) # M1 is a vector of length 'a' with token values representing the index of 'b' where that given word shows up, NA values represent words in 'a' that do not exist in 'b'

# Create a matrix M that has dimensions (length(a)-mlag)x(mlag+1)

M <- matrix(NA, length(a) - mlag, mlag + 1)

# Insert values into M such that row i contains the next mlag+1 tokens starting from the token for the i-th word of the text
# Ex. For mlag = 4, row 1 would contain the first five tokens of the text, row 2 would contain the second token to the sixth token and so on for the length of M
# The rows represent all of the mlag+1 tokens corresponding to the sequences of words in the text
for (i in 1:(length(a) - mlag)) {
  for (j in 1:(mlag + 1)) {
    M[i,j] <- M1[i + j - 1]
  }
}

# --- Building the Model -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# NEXT WORD FUNCTION
# A function that takes a (vector of) string(s) and predicts a word to follow
# The inputs for this function: key, M, M1, w
# - key: The string of words that is used to generate the next word; the length of this is not restricted or defined
# - M: The matrix defined above, this will be used to find matches of the key in the text, in order to identify possible next words
# - M1: The token vector of the text
# - w: Vector of mixture weights of length mlag which are used to assign probability values to the potential next words based on the length of the matched string.
#      The mixture weight is assigned based on the number of directly preceding words that are matched to produce a next word result
# The output for this function: next_word
# - next_word: The function's prediction of a word that should succeed the key

key <- c()

next.word <- function(key, M, M1, w = rep(1, ncol(M) - 1)) { 
  # Ensures the key vector is at least as long as mlag before beginning - inserts blanks at the start of the vector if this is not already the case
  if (length(key) < mlag) {
    key <- append(key, rep("", (mlag - length(key))), 0)
  }  
  
  # Initialise holding vectors
  weights <- c()
  next_token <- c()
  
  # Checks M for exact matches to the last mlag-i values of the key, adding successive words of matches to the vector of potential next words, with corresponding mixture weights being added to the weight vector
  for (mc in 1:mlag) {
    token_key <- match(key[(length(key) - mlag + mc):length(key)], b)  # Vector of the tokens corresponding to the last mlag values of the key
    ii <- colSums(!(t(M[, mc:mlag, drop=FALSE]) == token_key))  # Checks for matches of the token key in M (returning 0 for matches)
    full_matches <- which(is.finite(ii) & ii == 0) # Row indices of M that contain a match
    
    next_word_tokens <- na.omit(M[full_matches, mlag + 1]) # Vector of the tokens for all the words which are preceded by an exact match 
    weights <- append(weights, rep(w[mc]/length(next_word_tokens), length(next_word_tokens))) # For each entry in next_word_tokens, this vector contains the corresponding mixture weight
    next_token <- append(next_token, next_word_tokens)
  }
  
# If no matches were found, a random word is selected from a token in M1, otherwise a word will be sampled from the vector of matches based on their mixture weights
  if (length(next_token) == 0) {
    next_word <- b[sample(na.omit(M1), 1)]
    
    return(next_word)
  }
  else {
    next_word <- b[sample(next_token, size = 1, prob = weights)]
    return(next_word)
  }
}

# STARTER WORD TOKEN FUNCTION
# A function that provides a token representing a common word in the text
# The inputs for this function: start_word
# - start_word: A word for which the token will be found - this is optional
# The output for this function: start_token
# - start_token: A token of a word in vector b

starter.word.token <- function(start_word) {
  punc_vec <- append(punc_vec, c('?', '.')) # Adding these marks as their current format in the vector won't be matched properly
  punc_tokens <- as.numeric(na.omit(match(punc_vec, b)))

  # If no input is given a random, non-punctuation token is selected from M1
  if (missing(start_word)) { 
    start_token <- sample(na.omit(M1[! M1 %in% punc_tokens]), 1)
    start_word <- b[start_token]
    return(start_token)
  }
  start_word <- tolower(start_word) # Match is case-sensitive, so standardising to lowercase
  
  # If punctuation is provided, a random word token is selected from M1 instead, as punctuation is not an appropriate sentence starter
  if (start_word %in% punc_vec){
    print("Error: This is punctuation, not a valid start word; A random token has been generated instead.")
    start_token <- sample(na.omit(M1[! M1 %in% punc_tokens]), 1)
  }
  
  # If the word provided is not one of the common words, a random word token is selected from M1
  else if (!start_word %in% b){
    print("Error: This is not a recognized common word; A random token has been generated instead.")
    start_token <- sample(na.omit(M1[! M1 %in% punc_tokens]), 1)
  }
  
  else {start_token <- match(start_word, b)}
  return(start_token)
}

# GENERATE SENTENCE FUNCTION
# A function that generates a full sentence based on the (optional) prompt, stopping when a full stop is reached
# The inputs for this function: sentence_prompt
# - sentence_prompt: An optional input that may contain a user specified string.
# The output for this function: full_sentence
# - full_sentence: The string of words predicted based on the prompt

generate.sentence <- function(sentence_prompt) {
  punc_vec <- append(punc_vec, c('?', '.')) # Adding these marks as their current format in the vector won't be matched properly

  # The sentence will be started by a random common word if no prompt is provided (using the starter word function)
  if (missing(sentence_prompt)) {
    sentence_prompt <- b[starter.word.token()]
  }
  
  # If punctuation is provided, a random word is assigned as the prompt instead, as punctuation is not an appropriate sentence starter
  if (sentence_prompt %in% punc_vec) {
    sentence_prompt <- b[starter.word.token()]
  }
  
  # If the sentence prompt is several words long, it is redefined as a vector split into its individual words
  sentence_prompt <- unlist(strsplit(sentence_prompt, split = " "))
  
  # Continuously use the next word function to add words to the sentence until a full stop is reached
  repeat {
    next_word <- next.word(sentence_prompt, M, M1, w = rep(1, ncol(M) - 1))
    
    # Adds the next word to the end of the sentence prompt vector unless the last and next word are both punctuation
    if (!(next_word %in% punc_vec & sentence_prompt[length(sentence_prompt)] %in% punc_vec)) { 
      sentence_prompt <- append(sentence_prompt, next_word)
    }
    
    # Ends repeat loop when a full stop is reached as the sentence has ended
    if (sentence_prompt[length(sentence_prompt)] == '.') {
      break 
    }
  }
  
  full_sentence <- paste(sentence_prompt, collapse = " ") # Concatenates the vector of sentence words into a single string, separating the words with a space
  return(full_sentence)
}

# Example uses of the generate sentence function
generate.sentence()
# generate.sentence('If')
# generate.sentence('Simon')
# generate.sentence('Today I went')
