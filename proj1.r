# Clare Lewis (s2879721), Grace  (s), Luke  (s)
# brief description of each person's contribution and the % they did

##setwd("/Users/clarelewis/Documents/github/ext_stat_prog_group_work")
#setwd("C:/Users/Grace Sheahan/Documents/Extended Statistical Programming/Assessment 1") ## Setting working directory to folder for Assessment 1

a <- scan("shakespeare.txt",what="character",skip=83,nlines=196043-83,fileEncoding="UTF-8") # Scanning the complete works of shakespear in under the name 'a'

#### Look at the Shakespeare text ####

# Remove stage directions (indicated by'[_' and '_]')
open_direction <- grep('\\[_', a)  
direction <- length(open_direction)

direction_length <- rep(0, direction) 
for(i in 1:direction) {direction_length[i] <- grep('\\_]|\\.]', a[open_direction[i]:(open_direction[i]+100)])[1]}

unclosed_direction <- which(is.na(direction_length))

# unclosed_direction

# a[open_direction[unclosed_direction[1]]:(open_direction[unclosed_direction[1]+20)]
# a[open_direction[unclosed_direction[2]]:(open_direction[unclosed_direction[2]+20)]
direction_length[unclosed_direction[2]] <- 10

open_direction <- open_direction[-unclosed_direction[1]]
direction_length <- direction_length[-unclosed_direction[1]]
direction <- direction-1

close_direction <- rep(0, direction)  
for(i in 1:direction){close_direction[i] <- (open_direction[i] + direction_length[i] -1)}  

direction_words <- rep(0, sum(direction_length)) 
for(i in 1:direction){direction_words[(sum(direction_length[0:(i-1)])+1):(sum(direction_length[0:i]))] <- (open_direction[i]:close_direction[i])}

a <- a[-direction_words]

a <- gsub('\\[_', '', a)

#### 4B ####
# Look for and remove cases where the entire word is uppercase (excluding 'I', 'O', and 'A')

uppercase_indices <- which(a == toupper(a) & !(a %in% c('A', 'I', 'O')))
a <- a[-uppercase_indices]

#### 4C ####
# Find and remove all words that contain a hyphen or an underscore 

a_punc_hyph = grep("-", a, fixed = T )
a_punc_under = grep("_", a, fixed = T )

# Remove all cases of hyphens and underscores in each word in the text
a <- gsub(pattern = "[_-]", replacement = "", a)

#### NOTE: Consider splitting hyphenated words into two separate words. May yield better predictive power and prompts

#### 4D ####
# Create function that will split all punctuation into their own entries

split_punc <- function(text_vec, punc_vec){
  # Correct the format of the punctuation vector to be used in the "grep" function
  punc_vec <- paste(punc_vec, collapse = "|")
  
  # Find the indices of all words in the text vector that contain punctuation marks given in the punctuation vector
  i_punc <- grep(punc_vec, text_vec)
  
  # Create a holding vector for the elements of the text vector separated from their corresponding punctuation, such that the length of this vector is the length of text + number of punc cases
  # Note that this assumes that all words don't have two cases of punctuation in them. Sound assumption
  text_vec_extra <- rep(0, length(i_punc) + length(text_vec))
  
  # Indices in the holding vector which can be used as a reference to put the separated punctuation as new elements
  i_extra_elements <- i_punc + 1:length(i_punc) - 1
  
  # Take out the last character of words with punctuation using the "substr" function and put them in their corresponding element in the new text vector
  # Note: This assumes that all punctuation marks appear at the end of words. Sound assumption for the text file we are working with
  text_vec_extra[i_extra_elements + 1] <- substr(text_vec[i_punc], nchar(text_vec[i_punc]), nchar(text_vec[i_punc]))
  
  # Create a copy of the text vector with all punctuation removed
  text_vec_no_punc <- gsub(pattern = punc_vec, replacement = "", text_vec)
  
  # Puts the elements/words without punctuation in their corresponding position in the new word vector
  text_vec_extra[-i_extra_elements - 1] <- text_vec_no_punc
  
  # Return the new text vector with punctuation separated from all elements and are now considered there own unique element in the text vector
  return(text_vec_extra)
}

#### 4E ####
# Implement the function to separate punctuation

punc_vec <- c(",", "\\.", ";", "!", ":", "\\?")
a <- split_punc(a, punc_vec)

#### 4F ####
# Convert all the text lower case

a <- tolower(a)

#### 5A ####
#Created a vector of unique words in the text

words <- unique(a)

#### 5B ####
# Indicates which location in the words vector each word in the full text corresponds to

index <- match(a, words)

#### 5C ####
# Counts the occurrences of each word

word_occurences <- tabulate(index)

#### 5D ####
# Creates a vector of the 1000 most common words by using a for loop and rank to locate and order these words

b <- rep(0, 1000)
for (i in 1:1000) {
  b[i] <- words[which((rank(word_occurences, ties.method = "first")) == (length(word_occurences) + 1 - i))]
}

#### 6A ####

mlag <- 4 # this value can be changed by the user to change the output results of the model
token <- match(a,b)

#### 6B ####

M <- matrix(NA, length(a) - mlag, mlag + 1)
for (i in 1:(length(a) - mlag)) {
  for (j in 1:(mlag + 1)) {
    M[i,j] <- token[i + j - 1]
  }
}


#### 7 #####

M1 <- token
#key <-c() I think we dont need this command but need to test
w = rep(1, ncol(M) - 1) # can be changed by the user

next.word <- function(key ,M, M1, w) { 
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
    weights <- append(weights, rep(w[i],length(next_word_tokens)))
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

#### 8 ####

starter.word.token <- function(start_word){ ###Defining a starter word token function
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

#### 9 ####

#punc_tokens <- na.omit(b[punc_vec]) do we need this line?

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