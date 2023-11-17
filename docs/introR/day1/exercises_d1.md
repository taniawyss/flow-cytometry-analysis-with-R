
In this section, you will find the R code that we will use during the course. We will explain the code and output during correction of the exercises.

<!-- This is commented text -->

<!-- ## Source of data -->


## First exploration of R using RStudio

Four windows are displayed by default within RStudio. One of the windows corresponds to the R console.

<figure>
  <img src="../../assets/images/Rstudio.png" width="800"/>
</figure>


Type the following commands within the console (bottom left window in RStudio) at the prompt (">"), followed by the "Enter" key after each one to view the output printed on the console.

```r
1 + 1
```
The first command (1 + 1) prints "2" in the console. 
```r
x	<- 128.5
```
The second command does not print anything in the console, but a new variable called "x" and that contains the value 128.5 is created and listed in the Workspace (top right window in RStudio).
```r
x
```
The third command prints the value stored within the "x" variable in the console.
```r
abs(-11)
```
The fourth command, with the use of the abs() function, prints the absolute value of -11 in the console.
```r
?p.adjust
```
Finally, the fifth command opens the help page for the p.adjust function (bottom right "Help" window in RStudio).


## Working directory

To manipulate data within R, we first need to import it. R needs a way to locate the files within the hard drive or system. Therefore, we can specify the working directory, i.e. the location where R will look for files.

!!! warning
    To run the code below with setwd() make sure you put within the quotes a path that exists within your system.

```r

# To see what is the current working directory, use the function:
getwd()
# [1] "C:/Users/twyss/Documents/Rcourse"

# To change the working directory to any existing folder on your hard drive or system, use setwd() and the file path within quotes, e.g.
setwd("D:/R_exercises/")

``` 

## Workspace - Environment and history

Once a value has been assigned to a named variable, as we did assigning 128.5 to x above, the variable is saved and listed within the Workspace, which is displayed in one of the RStudio windows. 

Explore your workspace using the command line:

```r
# To list the objects or variables that are in your workspace, type
ls()
# To remove (delete) an object from the workspace, use function rm():
rm(x)
# To remove (delete) all objects from the workspace, type
rm(list=ls())

```

## Let's practice 2 - Create a script

R scripts allow you to save all code for further use or reference. For big projects, it is essential to create an R script.
To create a script, go to File > New File > R Script. Save it with file name "ex1.R" or any other that is suitable for you.
Add a comment symbol (#, the pound or hash sign) at the beginning of the first line.

Type or paste the following code. Look at the script (before running it).

Can you understand each line? What do you expect it to print to the console? Next, run the script and explore RStudio features such as the Workspace (Environment). Run the script line by line. Try both the "Run" button and the keyboard shortcut.  Watch variables appear in the Environment window (top right).
Watch what is printed to the console (bottom left window). Does it match your expectation?

```r

# First Steps and commands, ex. 1  
w <- 3
h <- 0.5
area <- w * h 
area
```

## Packages

When R is installed for the first time, a set of "base" packages is installed along the R software. The list of available packages can be viewed in the package Explorer window within RStudio (bottom right "Packages" window). Each package is a bundle of functions designed and created by an author to perform specific, usually related tasks. When working with "non-standard" data types, eg in bioinformatics or flow cytometry analysis, packages with bioinformatics-related functions need to be installed by the user. 

Common repositories for packages are CRAN and Bioconductor.

Install packages from [CRAN](https://cran.r-project.org/web/packages/available_packages_by_name.html) with the install.packages() function.

To install packages hosted on [Bioconductor](https://bioconductor.org/), we need 2 steps. First, we install a package called BiocManager, that will allow us to have access to the install() function to download Bioconductor packages.

```r
# Install packages hosted on CRAN: use a function from the utils package:
install.packages("stringi") # stringi is a package for character string manipulations

install.packages("rmarkdown") # rmarkdown package that allows to create pdf reports (see day 2)

# Install packages hosted on bioconductor: first install the BiocManager package that is available on CRAN:
install.packages("BiocManager")

# Then use the install() function from the BiocManager package
# !! This takes time to complete, run it during coffee or lunch break!
# Install flowCore:
BiocManager::install("flowCore")

# Install ggcyto:
BiocManager::install("ggcyto")

```

Once a package is installed, its content and functions need to be made accessible to R. library() loads the package for the current session.
It is good practice to load all needed packages at the top of a script.

```r
# My Script

library(limma)  
library(DESeq2)  
library(MASS)  
library(ggplot2)

# Here my data analysis begins
```

If you run the above code, what is the output on the Console? What does it mean?

??? done "Answer"
  	Packages such as limma or DESeq2 are not installed as base packages. They are hosted on Bioconductor and provide functions for RNAseq or microarray data analysis. The error message indicates that these packages were not installed and need to be installed before being able to load them.


## R version and session information

R is constantly upgraded by developers, which release a new version of R about every 6 months. Along with R upgrades, packages also get upgrades. From one version to the other of a package, it may happen that the default parameters of functions change. Therefore, it is important to always have in mind which current version of R and packages have been used for any analysis. Print the current R version and versions of attached or loaded packages using:

```r
# Prints the currently used R version
R.version.string	

# Print version information about R and all attached or loaded packages
sessionInfo() 

# Print the version of a specific package:
packageVersion("stringi")
```

## Let's practice - 3

1) Assign the values 6.7 and 56.3 to variables a and b, respectively.

2) Calculate (2*a)/b + (a*b) and assign the result to variable x. Display the content of x.

3) Find out how to compute the square root of variables. Compute the square roots of a and b and of the ratio a/b.

4) a) Calculate the logarithm to the base 2 of x (i.e., log2 x).
   b) Calculate the natural logarithm of x (i.e., loge x).

??? done "Answer"
	```r
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
    ```

## Let's practice - 4

1) Create two vectors, vector_a and vector_b, containing values from −5 to 5 and from 10 down to 0, respectively.

2) Calculate the (element-wise) sum, difference and product between the elements of vector_a and vector_b.

3) a) Calculate the sum of elements in vector_a.
   b) Calculate the overall sum of elements in both vector_a and vector_b.

4) a) Identify the smallest and the largest value in vector_a 
   b) among both vector_a and vector_b.

5) Compute the overall mean of the values among both vector_a and vector_b.

Hint: Each task in exercises 1-5 can be performed in a single statement per  vector (the minimum and maximum count as 2 tasks)

??? done "Answer"
	```r
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
    ```

## Let's practice - 5

1) In your script, write the command to load the package "MASS".

2) Write the following command to load the bacteria data set from the package MASS:
data(bacteria) # loads the bacteria data set (from MASS)

Execute the command. Check: You should have a variable named "bacteria" in your Environment/Workspace.

3) What are the names of the columns of the bacteria data.frame ?

4) Use [ ] to select rows 100 to 119 of the column “ap” .

5) Use $ to get the column "week" and check how many missing values it has.

Optional : 6) Count how many rows correspond to a “placebo” treatment (“trt”  column) using the comparison operator "==".

??? done "Answer"
	```r
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
    ```

## Let's practice - 6

A clinical dataset from patients with lung cancer is available in the file clinical_data2.csv.

The clinical_data2.csv file was generated from the clinical_data.csv source file using the following code:

```r
## !! Adapt the path to the path in your own system if you wish to import the data available in file "clinical_data.csv"
clinical_data <- read.csv("course_datasets/clinical_data.csv")

# View the format of the data:
head(clinical_data)

# Number of rows and columns:
dim(clinical_data)

# Column names:
colnames(clinical_data)

# Structure of the data:
str(clinical_data)

# Convert the gender to a factor and re-order the disease stage:
clinical_data$gender <- factor(clinical_data$gender)
clinical_data$stage <- factor(clinical_data$stage, levels = c("I","II","III","IV"))

str(clinical_data)

# Obtain a summary for each variable:
summary(clinical_data)

# View the data in rown number 2:
clinical_data[2,]

# View the data in column named "age":    
clinical_data[,"age"]

# Check the element number 30 of the disease stage column
clinical_data$stage[30]

# View the data corresponding to patients with stage II disease
subset(clinical_data, stage=="II")

# View the data corresponding to patients with stage II disease and that are female:
subset(clinical_data, stage=="II" & gender=="female") 

# View the data corresponding to patients with stage I or II disease and that are female:
subset(clinical_data, (stage=="I" | stage=="II") & gender=="female")

tapply(X=clinical_data$age, INDEX=clinical_data$stage, FUN=min)

# Add a new patient by concatenating rows and assign the result to a new variable called "clinical_updated":
clinical_updated <- rbind(clinical_data,
                          data.frame(sample_id = "LC02", 
                                     collection_date = "18.02.2021",
                                     age=71,
                                     gender= "female",
                                     stage="I"))

# Create a new treated variable:
treated <- rep( c("yes","no"), nrow(clinical_data)/2)
clinical_mod <- cbind(treated, clinical_data)

# Remove the first column and view the format:
clinical_orig <- clinical_mod[,-1] 
head(clinical_orig)

clinical_orig <- clinical_mod[,2:dim(clinical_mod)[2]] 

# Export to a new csv file:
write.table(clinical_updated, file="course_datasets/clinical_updated.csv",
            quote=FALSE, sep=",",row.names=FALSE)
```

Let's explore the dataset to see  what it contains.
1) Optional: Open a new script file in R studio, comment it and save it.

2) Have a look at the csv file in R studio's file explorer. What do you need to check in order to be able to read in the file correctly?

3) Read the file into R, assign its content to object "clinical_data2". Examine the object.

4) How many observations and variables does the dataset have?

5) What is the structure of the dataset? What are the names and classes of the variables?

6) Which variables appear to be categorical? Convert them to factors.

7) Get the summary statistics of "clinical_data2"

??? done "Answer"
	```r
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
    ```

## Let's practice 6bis

8) Use the function table() to compute the number of samples in  different patient groups. a) How many samples are included of each gender (male, female)? b) How many samples are included per level of response to treatment (PD, SD, PR, CR)? c) Make a 2x2 table gender and level of response to treatment.
Hint : try some of the example in the help(table) page.

9) Isolate the samples from male patients using subset(). Compute a summary statistics just for the weights of the  subset. Then do the same for the samples from female patients. Export the data of each subgroup to a csv file.

10) Compute the means and standard deviations for male and  female patient weights using tapply(). Then do the same by level of response to treatment.

??? done "Answer"
	```r
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
    ```


**End of Day 1, good job!** :sparkle:

<!--
## Feedback :sparkle:
-->









