# Vector of names
names <- c("Ben", "Alice", "John", "Stephan", "Anna", "Mary",
           "Clark", "Mackenzie", "Eleanor")

# How many names did we write?
number_of_names <- length(names)

# How long is each of the names?
length_of_names <- nchar(names)

# What is the mean name length in our sample?
mean <- mean(length_of_names)

# What is the standard deviation of lengths?
std <- sd(length_of_names)

# What is the standard error?
# The standard error by definition is equal to
# the standard deviation, divided by square root of N)
se <- std/sqrt(number_of_names)
