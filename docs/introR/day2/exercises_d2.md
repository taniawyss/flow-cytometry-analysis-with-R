
In this section, you will find the R code that we will use during the course. We will explain the code and output during correction of the exercises.

R offers many options for creating graphics and figures. Two common ways to create graphics within R use "base" functions such as plot(), or functions of the package [ggplot2](https://ggplot2.tidyverse.org/). 

So many possible chart types exist, and some are better than others at representing different data types. The online tool ["From Data to Viz"](https://www.data-to-viz.com/#explore) provides a nice decision tree to suggest a set of potentially appropriate graphics to visualize the dataset according to the input data format and type. Feel free to explore the decision tree! Most of the graphics suggested by "Data to Viz" are generated with ggplot2 functions and other packages providing extensions for ggplot2 functions. 

<!-- This is commented text -->

In these exercises, we will create graphics using the "base" R functions, i.e. functions of the graphics package that is installed at the same time as R is installed.

## Let's practice 7

Import the clinical data from the file clinical_data_mod.csv. This files contains the same data as clinical_data.csv and in addition, one more column.

1) Run str() to check your data frame: did it load correctly?

2) Convert gender and response_to_treatment to factor variables.

3) Plot a histogram of patient weight and customize it with colours, labels, title.

4) Make a scatter plot of height against patient weights using the function plot(). 

Function arguments:

* use solid circles as plotting symbol

* add a title

* customize the axis labels (“Weight [kg]”, “Height [m]”)

* color the points by gender.

Add a legend for the gender. Fit a trend line using the function abline().

5) Create a new column called “BMI” and compute the BMI of patients from their weight and height 

6) Make boxplots of BMI from patients with different responses to treatment. Customize with title, labels, colors. Add points to the boxplots to show the individual values.

7) Optional: Repeat 6 with stage, instead of response to treatment. Hint: what kind of variable is stage?


??? done "Answer"
	```r
    clinical_data <- read.table("course_datasets/clinical_data_mod.csv", header=TRUE, sep=",") # define classes for columns

    str(clinical_data)


    # 2) Convert gender and response_to_treatment to factor variables
    # define the order of factor levels
    clinical_data$gender <- factor(clinical_data$gender)
    clinical_data$response_to_treatment <- factor(clinical_data$response_to_treatment, levels = c("PD","SD","PR","CR"))

    # 3) Plot an histogram of patient weight and customize it with colours, labels, title and represent the density line on top.
    hist(clinical_data$weight,
     freq=FALSE, breaks=8,
     main="Patient Weight",
     col="orange" ,
     xlab="Weight [kg]")
    # lines(density(clinical_data$weight), col='blue')
    # Note: freq=FALSE makes the histogram density based, which makes it scale well with the density line

    # 4) Make a scatter plot of height against patient weights using the function plot().
    #    Function arguments:
    #     - use solid circles as plotting symbol
    #     - add a title
    #     - customise the axis labels  (“Weight [kg]”, “Height [m]”)
    #     - colour the points by gender. 
    # Add a legend for the gender.

    plot(clinical_data$weight,clinical_data$height, 
     pch=19,
     main="Weight vs Height in Patients",
     xlab="Weight [kg]", ylab="Height [m]",
     col=c("orange", "blue")[clinical_data$gender]
     )

    legend("bottomright",
       legend=levels(clinical_data$gender),
       col=c("orange","blue"),
       pch=19)

    abline(lm(clinical_data$height ~ clinical_data$weight),
       col="black", lwd=1.5)

    # 5) Compute the BMI = Weight / Height^2
    clinical_data$BMI <- clinical_data$weight / (clinical_data$height^2)

    # 6) Make boxplots of BMIs from patients with different responses to treatment. Customize with title, labels, colors.
    boxplot(BMI ~ response_to_treatment, data= clinical_data,
        col=c("red", "orange","green","blue"),
        main="BMI by level of Response to treatment",
        xlab="Response to treatment", ylab="BMI"
    )
    points(BMI ~ response_to_treatment, data= clinical_data)

    # 7) Optional: Repeat 6 with BMI and stage, instead of weight and gender.
    clinical_data$stage <- factor(clinical_data$stage, levels = c("I","II","III","IV"))
    boxplot(BMI ~ stage, data= clinical_data,
        col=c("blue", "green","orange","red"),
        main="Patient BMI by Stage",
        xlab="Stage", ylab="BMI"
    )
    points(BMI ~ stage, data= clinical_data)
    
    getwd() # check where you are. If you didn't change anything, you will be in the folder with the .Rproj file (rproject root)
    ```

## Let's practice 8

1) Create a multi-panel figure with the four graphics (3, 4, 6 and 7 from the previous exercise) on one page, exporting the figure to a pdf file with paper size A4. Set width and height arguments in the call to pdf() to make it look nice.

2) Optional: Export the histogram (3 from previous exercise) to a png file. Set the width and height arguments in the call to png() to make it look nice.

??? done "Answer"
	```r
	  # 1) Make a multi-panel figure with the four graphics on one page, exporting the figure to a png file.
    # Set width and height arguments in the call to png() to make it look nice.

    pdf("clinical_data_plots.pdf", width=7, height=7, paper="a4") 

    # we want to put 4 plots on the same panel -> 2 rows and 2 columns
    par(mfrow=c(2,2))

    # Plot 1
    hist(clinical_data$weight,
     freq=FALSE, breaks=8,
     main="Patient Weight",
     col="orange" ,
     xlab="Weight [kg]")
    lines(density(clinical_data$weight), col='blue')

    # Plot 2
    plot(clinical_data$weight,clinical_data$height, 
     pch=19,
     main="Weight vs Height in Patients",
     xlab="Weight [kg]", ylab="Height [m]",
     col=c("orange", "blue")[clinical_data$gender]
    )

    legend("bottomright",
       legend=levels(clinical_data$gender),
       col=c("orange","blue"),
       pch=19)

    abline(lm(clinical_data$height ~ clinical_data$weight),
       col="black", lwd=1.5)

    # Plot 3
    boxplot(BMI ~ response_to_treatment, data= clinical_data,
        col=c("red", "orange","green","blue"),
        main="BMI by level of Response to treatment",
        xlab="Response to treatment", ylab="BMI"
    )
    points(BMI ~ response_to_treatment, data= clinical_data)


    # Plot 4
    boxplot(BMI ~ stage, data= clinical_data,
        col=c("blue", "green","orange","red"),
        main="Patient BMI by Stage",
        xlab="Stage", ylab="BMI"
    )
    points(BMI ~ stage, data= clinical_data)

    dev.off()

    # 2) Optional: Export the histogram (3 from previous exercise) to a png file. 
    # Set width and height arguments in the call to png() to make it look nice.

    # png: width and height are in pixels by default

    png("hist_weight.png", width=800, height=600)
    hist(clinical_data$weight,
     freq=FALSE, breaks=8,
     main="Patient Weight",
     col="orange" ,
     xlab="Weight [kg]")
    lines(density(clinical_data$weight), col='blue')
    dev.off()

	```

**Bonus** Create plots with [ggplot2](https://ggplot2.tidyverse.org/).

ggplot2 is a package that allows you to create graphics by adding customization "layers" one after the other. You first provide the data, tell ggplot2 how to map variables to aesthetics (i.e. the variables you want to have on the x and y coordinates), what plot type to use, and it provides some default colors for the categorical variables, and automatic legend positioning. Everything can then be further customized by the user by adding additional layers using ggplot2 functions.

```r
library(ggplot2)

# Import data:
clinical_data <- read.table("course_datasets/clinical_data_mod.csv", header=TRUE, sep=",") # define classes for columns
str(clinical_data)

# Convert gender  to factor variables
clinical_data$gender <- factor(clinical_data$gender)

# Compute the BMI = Weight / Height^2
clinical_data$BMI <- clinical_data$weight / (clinical_data$height^2)

# Simple boxplot of BMI vs stage
ggplot(data=clinical_data, aes(x=stage, y=BMI)) + 
  geom_boxplot()

# Boxplot of BMI vs stage, coloring according to stage
# With legend:
ggplot(data=clinical_data, aes(x=stage, y=BMI, color=stage)) + 
  geom_boxplot() 

# Boxplot of BMI vs stage, coloring according to stage, flipping orientation
ggplot(data=clinical_data, aes(x=stage, y=BMI, color=stage)) + 
  geom_boxplot() + 
  coord_flip()  

# Without legend:
ggplot(data=clinical_data, aes(x=stage, y=BMI, color=stage)) + 
  geom_boxplot() + 
  theme(legend.position = "none")

# Scatter plot of weight vs height, coloring by gender
ggplot(data=clinical_data, aes(x=weight, y=height, color=gender)) +
  geom_point()

# Compare plot generated earlier with base R to the one with ggplot2:
# Reproduce Plot 2 of exercise 8
plot(clinical_data$weight,clinical_data$height, 
     pch=19,
     main="Weight vs Height in Patients",
     xlab="Weight [kg]", ylab="Height [m]",
     col=c("orange", "blue")[clinical_data$gender])
legend("bottomright",
       legend=levels(clinical_data$gender),
       col=c("orange","blue"),
       pch=19)
abline(lm(clinical_data$height ~ clinical_data$weight),
       col="black", lwd=1.5)

# Same plot with ggplot2
# Scatter plot of weight vs height, coloring by gender
ggplot(data=clinical_data, aes(x=weight, y=height, color=gender)) +
  geom_point() +
  scale_color_manual(values=c("female"="orange", "male"="blue")) + 
  ggtitle("Weight vs Height in Patients") +
  xlab("Weight [kg]") + ylab("Height [m]") +
  theme_bw() +
  geom_smooth(method = "lm",
              formula = y ~ x,
              se=TRUE) # display confidence interval around smoothed curve

# Create a separate plot for males and females:
# Use Facet wrap: separate the plots according to a categorical (factor) variable
ggplot(data=clinical_data, aes(x=weight, y=height, color=gender)) +
  geom_point() +
  scale_color_manual(values=c("female"="orange", "male"="blue")) + 
  ggtitle("Weight vs Height in Patients") +
  xlab("Weight [kg]") + ylab("Height [m]") +
  theme_bw() +
  geom_smooth(method = "lm",
              formula = y ~ x,
              se=TRUE) +
  facet_wrap(~gender)  
  
  
### Multi panel figures with ggplot2
# save each plot in an object, create 2 plots that are the same except for the color of the dots
p1 <- ggplot(data=clinical_data, aes(x=weight, y=height, color=gender)) +
  geom_point() +
  scale_color_manual(values=c("female"="orange", "male"="blue")) + 
  ggtitle("Weight vs Height in Patients") +
  xlab("Weight [kg]") + ylab("Height [m]") +
  theme_bw() 

p2 <- ggplot(data=clinical_data, aes(x=weight, y=height, color=gender)) +
  geom_point() +
  scale_color_manual(values=c("female"="turquoise", "male"="plum")) + 
  ggtitle("Weight vs Height in Patients") +
  xlab("Weight [kg]") + ylab("Height [m]") +
  theme_bw() 

# install the cowplot package to arrange several plots on a single page
# install.packages("cowplot")
library(cowplot)
plot_grid(p1, p2, nrow=1)  

```

## Let's practice 9


In this exercise we will use a 36-color spectral flow cytometry dataset from a study performed in the context of Covid-19 research. Only a subset from 4 healthy donors will be used. For each healthy donor, there are three time points, as indicated in the FCS file names (day 0, day 7, day 14). Data was downloaded through the Flow Repository database, [FR-FCM-Z3WR](https://flowrepository.org/id/FR-FCM-Z3WR). FCS files were pre-gated on live CD3+CD19-T cells in FlowJo.

Create a new script in which you will:

1) Import the FCS files (located in course_datasets/FR_FCM_Z3WR/). Do not transform or truncate the values.

2) Create a data frame with the list of channels and corresponding antigens, and plot it. Hint: get the antigens from the parameters of one of the flowFrame in the set

3) Convert the channel names in the expression matrices to the corresponding antigen names (where applicable)

4) Add a new column to the phenotypic data with the time point of the sample. Plot the phenotypic data

5) Create a bivariate density plot showing «FSC-H» againts «HLA-DR» for all samples from day 0. Apply a flowJo inverse hyperbolic sine scale to the y axis («HLA-DR»)

??? done "Answer"
	```r
    # load packages
    library(flowCore)
    library(ggcyto)

    # 1) Import the FCS files (course_datasets/FR_FCM_Z3WR/). 
    # Do not transform or truncate the values. 

    # path to the directory with the fcs files
    fcs.dir<- file.path( "course_datasets/FR_FCM_Z3WR/")

    # read fcs files into a flowSet
    fcs_data <- read.flowSet(path=fcs.dir, pattern="*.fcs", transformation = FALSE, truncate_max_range = FALSE) 

    #2) Create a data frame with the list of channels and corresponding antigens, and plot it.
    #Hint: get the antigens from the parameters of one of the flowFrame in the set
    channels <- colnames(fcs_data)
    antigen <- pData(parameters(fcs_data[[1]]))$desc
    panel <- data.frame(channels = channels, antigen= antigen)

    # View the panel in the console
    panel

    #3) Convert the channel names in the expression matrices to the 
    # corresponding antigen names (where applicable)

    colnames(fcs_data)[!is.na(antigen)] <- antigen[!is.na(antigen)] 

    # check
    head(exprs(fcs_data[[1]])[,c(5:10)])

    #4) Add a new column to the phenotypic data with the time point of the sample
    # check sample names
    sampleNames(fcs_data)
    # [1] "0E1F8E_0.fcs"  "0E1F8E_14.fcs" "0E1F8E_7.fcs"  "180E1A_0.fcs"  "180E1A_14.fcs" "180E1A_7.fcs" 
    # [7] "1A9B20_0.fcs"  "1A9B20_14.fcs" "1A9B20_7.fcs"  "61BBAD_0.fcs"  "61BBAD_14.fcs" "61BBAD_7.fcs" 

    # add column with time point
    pData(fcs_data)$time_point <- rep(c("D0","D14","D7"),4)

    # View the phenotypic data
    pData(fcs_data)

    # 5) Create a bivariate density plot showing "FSC-H" against "HLA-DR" for all samples from day 0. 
    # Apply a flowJO inverse hyperbolic sine scale to the y axis ("HLA-DR")

    # split by time point 
    fcs_data.split <- split(fcs_data, pData(fcs_data)$time_point)

    # create the bivariate density plot
    autoplot(fcs_data.split$D0, x="FSC-H",y="HLA-DR", bins = 64) + 
    scale_x_flowjo_biexp() + 
    scale_y_flowjo_fasinh()

    
	```

## Let's practice - R markdown

If a data analysis project involves many steps and generation of various plots, one of the easy and very practical ways to bundle and organize all steps of analysis together is to use R markdown files to generate PDF or html reports. These reports both display the R code used as well as the output generated, such as graphics, tables, statistical test results, ...

The difference between an R script and an R markdown file, is that the code is organized within chunks in the R markdown file. In between the chunks, the user can write text that contains information about the analysis.

To create an R markdown file, go to File > New File > R markdown. Add a name.
This will create a new file that already has some example content. As you can see, the R code is organized in chunks highlighed in grey, with details written as free text in between the chunks. 
We can see that the pound sign (#) is used outside of the R code chunks. In this case, the # symbol does not correspond to a comment, but will indicate header levels for the titles and subtitles within your final document obtained after report generation.

Once the Rmd is ready, the report can be generated by hitting the "Knit" button at the top of the window. 
<figure>
  <img src="../../assets/images/example_Rmd.png" width="800"/>
</figure>


The example Rmd generates the following html report (saved in the same folder as the Rmd file by default), that shows both the code and the resulting output:
<figure>
  <img src="../../assets/images/example_html.png" width="800"/>
</figure>


You can find a short video that introduces some of the principles of R markdown on [Youtube](https://www.youtube.com/watch?v=2YZSDGGoQzQ), from the beginning up to minute 23:30. Starting at minute 23:30, this video also introduces ggplot2.


If you would like to practice creating your own R markdown, modify the one that is generated with the example content when you select File > New File > R Markdown.

1) Create a new Rmd file with the following options at the top (in the top YAML instructions within the 2 dash sequences "- - -")

* Title: «Let’s practice»
* Author: your name
* Select the «use current date when rendering object» option
* Default output format: HTML

2) We will repeat exercise 7, but this time by creating a report:
Within an R code chunk, import the data from the file clinical_data_mod.csv, convert gender and response to treatment to factors and compute the BMI of patients in a separate code chunk called “prepare_data”. Change the chunk options so that code will not appear in the output.

Then, create a new code chunk for each plot. Make sure the plot is centered.
Add a header (preceded by the # symbols outside of the code chunks) before each plot with some suggestive plot title.

3) Save the Rmd file and produce the html document by «knitting» it.

[Download solution Rmd file](../assets/pdf/R_practice_Rmd_solution.Rmd){: .md-button }

For tweaking your reports, such as chosing different output formats, or hiding or showing the code within the report, we recommend that you consult the R markdown documentation provided in this [Definite guide eBook](https://bookdown.org/yihui/rmarkdown/).

Another useful resource is [RStudio's R Markdown tutorial](https://rmarkdown.rstudio.com/lesson-1.html). 


## Let's practice - R markdown bis

Create an R markdown file and knit to an html report the flow cytometry data analysis performed in Exercise number 9. 
Proceed similarly, creating a code chunk where you load the libraries, import and pre-process the data, and create a separate code chunk for every plot.

[Download solution Rmd file](../assets/pdf/R_practice_Rmdbis_solution.Rmd){: .md-button }


## Let's practice - Introduction to statistics with R.

The SIB course "First Steps with R in life sciences" provides additional material to perform statistics with R.

The slides introducing statistics with R can be found [here](https://github.com/sib-swiss/first-steps-with-R-training/blob/master/slides/First-steps-with-R_day2_afternoon.pdf) on the [course material github page](https://github.com/sib-swiss/first-steps-with-R-training). 

The R code to run Wilcoxon or T tests can be found [here (ex.9)](https://github.com/sib-swiss/first-steps-with-R-training/blob/master/solutions/R_practice9_solution.R) using source data that can be found [here](https://github.com/sib-swiss/first-steps-with-R-training/tree/master/course_datasets). The R code to run a linear model can be found [here (ex.10)](https://github.com/sib-swiss/first-steps-with-R-training/blob/master/solutions/R_practice10_solution.R). Feel free to try it out!


**End of Day 2, good job!** :sparkle:





