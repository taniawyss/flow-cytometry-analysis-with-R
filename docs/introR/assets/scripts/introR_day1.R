# Introduction to R for flow cytometry analysis - September 2023

# --- Type directly within Console:
1 + 1
x   <- 128.5
x
abs(-11)
?p.adjust

# --- Create a Project
# in existing directory

# ---- Working directory
getwd()
setwd("/export/scratch/twyss/SIB_training/flowCyt_2023/")

# To list the objects or variables that are in your workspace, type
ls()
# To remove (delete) an object from the workspace, use function rm():
rm(x)
# To remove (delete) all objects from the workspace, type
rm(list=ls())

# ---- Create a script
# First Steps and commands, ex. 1  
w <- 3
h <- 0.5
area <- w * h 
area

# Spaces and new lines
1 + 1 +
  + 2

# ---- Workspace

ls()
# To remove (delete) an object from the workspace, use function rm:
rm(x)
# To remove (delete) all objects from the workspace, type
rm(list=ls())


# Install packages hosted on CRAN: use a function from the utils package:
install.packages("stringi") # stringi is a package for character string manipulations

install.packages("rmarkdown") # rmarkdown package that allows to create pdf reports (see day 2)

# Install packages hosted on bioconductor: first install the BiocManager package that is available on CRAN:
install.packages("BiocManager")

## !! takes time, run during break
# Then use the install() function from the BiocManager package
# Install flowCore:
BiocManager::install("flowCore")
BiocManager::install("flowCore", force = TRUE)

# Install ggcyto:
BiocManager::install("ggcyto")
install("ggcyto")

# Load the packages:
library(limma)  
library(DESeq2)  
library(MASS)  
library(ggplot2)

# ---- R and package version
# Prints the currently used R version
R.version.string    

# Print version information about R and all attached or loaded packages
sessionInfo() 

# Print the version of a specific package:
packageVersion("stringi")

# Navigate:
# tab key
# up and down arrows
# ctrl-l to erase console window
# Ctrl-1 and Ctrl-2 to jump between the script and the console windows

# incomplete statement (+) (in the console, will wait), or Esc
1 + 1+
  + 2


# ---- Let's practice 3

# 1) Assign the values 6.7 and 56.3 to variables "a" and "b", respectively.
a <- 6.7
b <- 56.3

# 2) Calculate (2*a)/b + (a*b) and assign the result to variable "x". Display the content of "x".
x <- (2*a)/b + a*b
x

# 3) Find out how to compute the square root of variables. Compute the square roots of "a" and "b" and of the ratio "a/b".
sqrt(a) #using function sqrt()
b^0.5 # power 0.5 is the square root
(a/b)**0.5 # another way of specifying power

# 4) a) Calculate the logarithm to the base 2 of "x".
#   b) Calculate the natural logarithm of "x".
log2(x) # Function specifically for Log 2. Alternatively: log(x, base=2)
log(x)  # If we don't specify the base, default is the natural logarithm. 


# ---- Let's practice 4

# 1) Create two vectors, "vector_a" and "vector_b", containing the values from −5 to 5 and from 10 down to 0, respectively.
vector_a <- -5:5
vector_b <- seq(10,0) # alternatively: vector_b <- c(10,9,8,7,6,5,4,3,2,1,0)

# 2) Calculate the (elementwise) sum, difference and product between the elements of "vector_a" and "vector_b". 
vector_a + vector_b #sum
vector_a - vector_b #difference
vector_a * vector_b #product

# 3) a) Calculate the sum of elements in "vector_a"
#    b) Calculate the overall sum of elements in both  "vector_a" and "vector_b".
sum(vector_a)
sum(vector_a, vector_b) # alternatively : sum(vector_a + vector_b)

# 4) Identify the smallest and the largest value among both "vector_a" and "vector_b".
min(vector_a, vector_b)
max(vector_a, vector_b)

# 5) Compute the overall mean of the values among both "vector_a" and "vector_b"
mean( c( vector_a, vector_b) ) # mean() works only on a single vector, unlike sum, min and max! 
# Concatenate both vectors (using c() ) before computing the mean

# ---- Let's practice 5

# 1)Install and load the package MASS (or other CRAN packages).
# install.packages("MASS")
library(MASS)

# 2) The following command line loads the bacteria data.frame present in the MASS package. Execute it:
data(bacteria)
?bacteria

# 3) What are the names of the columns of the bacteria data.frame ?
names(bacteria)

# 4) Use the [ ] , to select in bacteria rows 100 to 119 in the column "ap".
bacteria[ 100:119 , "ap" ]

# 5) Use $ to get the column "week" and check how many missing values it has.
sum(is.na(bacteria$week))

# Optional : 6) use comparison operators to count how many rows correspond to a “placebo” treatment (“trt” column).
sum(bacteria$trt == "placebo")




# ---- Let's practice 6
# A clinical dataset from patients with lung cancer is available in the file clinical_data2.csv. 
# Let's explore the dataset to see  what it contains.

# 1) Open a new script file in R studio, comment it and save it.
# 2) Have look at the csv file in R studio's file explorer. What do you need to check in order to be able to  read in the file correctly?

# 3) Read the file into R, assign its content to object "clinical_data2". Examine the object.
# Adapt the path to the path in your own system!
clinical_data2 <- read.csv("course_datasets/clinical_data2.csv")

# 4) How many observations and variables does the dataset have?
dim(clinical_data2)

# 5) What is the structure of the dataset? What are the names and classes  of the variables?
str(clinical_data2)

# 6) Which variables appear to be categorical? Convert them to factors.

clinical_data2$gender <- factor(clinical_data2$gender)
clinical_data2$stage <- factor(clinical_data2$stage, levels = c("I","II","III","IV"))
clinical_data2$treatment_status <- factor(clinical_data2$treatment_status)
clinical_data2$response_to_treatment <- factor(clinical_data2$response_to_treatment,levels = c("PD","SD","PR","CR"))

# 7) Get the summary statistics of "clinical_data2"
summary(clinical_data2)


# ---- Let's practice 6bis

# 8) Use the function table() to compute the number of samples in  different patient groups.
# Hint : try some of the example in the help(table) page.

# a) How many samples are included of each gender (male, female)? 
table(clinical_data2$gender)

# b) How many samples are included per level of response to treatment (PD, SD, PR, CR)? 
table(clinical_data2$response_to_treatment)

# c) Make a 2x2 table gender and level of response to treatment.
table(clinical_data2[,c("gender","response_to_treatment")])

# 9) Isolate the samples from male patients using  subset(). 
# Compute a summary statistics just for the weights of the  subset. 
# Then do the same for the samples from female patients. 
# Export the data of each subgroup to a csv file.

# Isolate the samples from male patients
clinical_data2_male <- subset(clinical_data2, gender=="male")

# Compute a summary statistics just for the weights of the  subset
summary(clinical_data2_male$weight)

# Export the data to a csv file.
write.table(clinical_data2_male,file = "course_datasets/clinical_data2_male.csv",
            quote=FALSE, 
            sep=",",
            row.names=FALSE)

# Isolate the samples from female patients
clinical_data2_female <- subset(clinical_data2, gender=="female")

# Compute a summary statistics just for the weights of the subset
summary(clinical_data2_female$weight)

# Export the data to a csv file.
write.table(clinical_data2_female,file = "course_datasets/clinical_data2_female.csv",
            quote=FALSE, 
            sep=",",
            row.names=FALSE)

# 10) Compute the means and standard deviations for male and female patient weights using tapply(). 
# Then do the same by level of response to treatment.

# by gender
tapply(clinical_data2$weight, clinical_data2$gender, mean) # mean: by gender
tapply(clinical_data2$weight, clinical_data2$gender, sd) # standard deviation by gender

# by response to treatment
tapply(clinical_data2$weight, clinical_data2$response_to_treatment, mean) # mean: by level of response to treatment
tapply(clinical_data2$weight, clinical_data2$response_to_treatment, sd) # standard deviation by level of response to treatment



