# Day 2

setwd("/export/scratch/twyss/SIB_training/flowCyt_2023/")

# Plotting - the basics
x <- 1:100
y <- log(x) + (x/100)^5
plot(x,y) # equivalent to plot(x, y, type="p")

# Adding elements to plots:
x <- seq(0,100, by=10)
y <- log(x) + (x/100)^5

plot(x, y, ylim=c(0,7))
lines(x, y+1)
points(x ,y-1, type="b")

# Customizing plots
plot(x, y, type="l", col="red",  lwd=7)
lines(x, y+1, col="blue",  lty="dashed")
points(x, y-1, type="b",
        pch=19)

# pch: possible shapes
help(points)

colors() # list possible color names

# Customizing plots further:
x <- seq(0, 100, length.out=10)
y <- log(x) + (x/100)^5
plot(x,y, type="l", col="red", 
     ylim=c(1,7),  
     xlab="The variable x", main ="x vs. y" )
lines(x, y+1, lwd=3, lty="dashed", col="blue")
points(x, y-1, type="b", pch=15)

# Add a legend:
legend(x="bottomright",  
       legend=c("red line","blue line", "black line"),  
       lty=c(1,2,1), 
       pch=c(NA,NA,19),
       col=c("red", "blue", "black"),  bg="gray90")

# Use data available within packages to play:
data()
data(package = .packages(all.available = TRUE))
data(iris) #data sets in standard packages
?iris
head(iris)

# Generate random numbers that are normally distributed
rnorm(10) # mean=0, sd=1  by default
rnorm(10, mean=10, sd=2) #customized mean and sd

# Histograms
x <- rnorm(10000) # generate 10k numbers
?hist
hist(x, breaks=20,  freq=FALSE,  
     main="Hist",  col="pink")


# Boxplots
library(MASS)
data(Melanoma) #Data from MASS package. 205  Denmark with malignant melanoma
# Check what is in the data frame:
head(Melanoma)

# A simple boxplot
?boxplot
boxplot(Melanoma$thickness)
boxplot(Melanoma$thickness,
         ylab="Tumour thickness (mm)",  col="white")

# Check the class of the $status column:
str(Melanoma)

# Convert status to factor:
?boxplot
Melanoma$status <- factor(Melanoma$status)
levels(Melanoma$status)

boxplot(Melanoma$thickness ~ Melanoma$status)

# More boxplots:

# Provide several vectors or data subsets:
boxplot(Melanoma$thickness[Melanoma$status=="1"],  
        Melanoma$thickness[Melanoma$status=="2"],  
        Melanoma$thickness[Melanoma$status=="3"],
        main="Thickness of melanoma per patient status", xlab="status", 
        ylab="Tumour thickness",  names=c("1","2","3"))
points(Melanoma$status, Melanoma$thickness,  col="blue",pch=19) # x, y then plotting parameters


# Only plot the thickness for 2 status groups:
boxplot(Melanoma$thickness[Melanoma$status=="1"],  
        Melanoma$thickness[Melanoma$status=="3"],
        main="Thickness of melanoma per patient status", xlab="status", 
        ylab="Tumour thickness",  names=c("1","3"))

# Provide a formula:
#  ?~.

boxplot(Melanoma$thickness ~ Melanoma$status,
        main="Thickness of melanoma per patient status", xlab="status", 
        ylab="Tumour thickness",  names=c("1","2","3"))
points(thickness ~ status, data=Melanoma, # y divided into the levels of x
       col="blue", pch=19)


# Abline
?airquality
data(airquality) # Daily measurements, New York,
plot(airquality$Wind, airquality$Ozone, pch=20,  
     xlab= "Wind (mph)", ylab="Ozone (ppb)")
abline(h=60, col="red", lty="dashed")
abline(v=seq(3,21,3), col="grey", lty="dotdash")
legend("topright", "Maximum allowable ozone concentration",  col="red", lty="dashed")

# Fitting a trend line:
plot(airquality$Wind, airquality$Ozone, pch=20,  xlab= "Wind (mph)", ylab="Ozone (ppb)")
?lm
abline(lm(airquality$Ozone ~ airquality$Wind),  col=2, lwd=2)
legend("topright", legend= c("measures","fitted line"),  
       pch= c(20, NA), 
       lty = c(0, 1), 
       lwd=c(NA, 2),
       col = c(1, 2), bg = "gray90")

# Scatter plot pairs
head(iris) 
pairs(iris[,1:4], main="Edgar Anderson's Iris Data",  
      col=c("red", "green3", "blue")[iris$Species])

# Use bg because pch=21 is a solid circle
pairs(iris[,1:4], main="Edgar Anderson's Iris Data",  
      pch=21, bg=c("red", "green3", "blue")[iris$Species]) 



# colors:
# setosa in red
# versicolor in green
# virginica in blue
# Why? Each of the 3 listed colors gets attributed to one of the factors
# following the levels order.
levels(iris$Species)
View(cbind(iris, 
           c("red", "green3", "blue")[iris$Species]))


# ---- Let's practice 7

clinical_data <- read.table("course_datasets/clinical_data_mod.csv", header=TRUE, sep=",") # define classes for columns

str(clinical_data)


# 2) Convert gender and response_to_treatment to factor variables
# define the order of factor levels
clinical_data$gender <- factor(clinical_data$gender)
clinical_data$response_to_treatment <- factor(clinical_data$response_to_treatment, 
                                              levels = c("PD","SD","PR","CR"))

# 3) Plot a histogram of patient weight and customize it with colors, labels,
# title
hist(clinical_data$weight,
     freq=FALSE, breaks=8,
     main="Patient Weight",
     col="orange" ,
     xlab="Weight [kg]")
# Note: freq=FALSE makes the histogram density based, which makes it scale 
# well with the density line
# represent the density line on top of the histogram
lines(density(clinical_data$weight), col='blue') # smoothed 


# 4) Make a scatter plot of height against patient weights using the function plot().
#    Function arguments:
#     - use solid circles as plotting symbol
#     - add a title
#     - customize the axis labels  (“Weight [kg]”, “Height [m]”)
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
        main="BMI by Revel of Response to treatment",
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

# --- Permanent graphical parameter changes
?par
# Set ploting color and pch (symbol shape)
par(col="red", pch=15)
plot(1:10, 2:11)

# Set margins:
par(mar=c(5.1,4.1,4.1,2.1))	#set margins in lines

# Normal margins
par(mar=c(5.1,4.1,4.1,2.1))
plot(1:10)

# Wide margins
par(mar=c(8.1,8.1,8.1,8.1))
plot(1:10)

par(mar=c(5.1,4.1,4.1,2.1))

# Multi-panel figures using par(mfrow=c())
par(mfrow=c(1,2),col="firebrick", pch=19, 
    mar=c(5.1,4.1,4.1,2.1))	#1x2 plot array
x <- seq(-100, 100, 0.1)
plot(x, y=x^2, ylim=c(-10000, 10000), main="quadratic")
plot(x, y=x^3, ylim=c(-10000, 10000), main="cubic")

plot(1,1)
plot(1:10)

par(mfrow=c(1,1))

# Current settings of par()
par()
par()$pch
par()$col

# Resetting par
# par() is automatically reset to defaults when you:
# Restart R or close/switch Rstudio projects
# Run dev.off(), which closes the most recent plot/plotting device
# Run graphics.off(), which closes plots/plotting devices
# In RStudio, clear all plots using the broom icon

# Saving figures to files
# The new file with the plot is saved inside the working directory
pdf(file="quadratic_cubic.pdf", width=7, height=4,  paper="a4")
par(mfrow=c(1,2),col="firebrick", pch=19)
x <- seq(-100, 100, 0.1)
plot(x, y=x^2, ylim=c(-10000, 10000), main="quadratic")
plot(x, y=x^3, ylim=c(-10000, 10000), main="cubic")
dev.off()

# --- Let's practice 8

# 1) Create a multi-panel figure with the four graphics on one page, exporting the figure to a png file.
# Set width and height arguments in the call to png() to make it look nice.

pdf("clinical_data_plots.pdf", # width=7, height=7, 
    paper="a4") 

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
        main="BMI by Revel of Response to treatment",
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

# Bonus ggplot2
library(ggplot2)

# Boxplot of BMI vs stage, coloring according to stage
# With legend:
ggplot(data=clinical_data, aes(x=stage, y=BMI, color=stage)) + 
  geom_boxplot() 

# Boxplot of BMI vs stage, coloring according to stage, flipping orientation
ggplot(data=clinical_data, aes(x=stage, y=BMI, color=stage)) + 
  geom_boxplot() + 
  coord_flip()  

# Without legend:
ggplot(data=clinical_data, aes(x=BMI, y=stage, color=stage)) + 
  geom_boxplot() + 
  coord_flip() + 
  theme(legend.position = "none")

# Scatter plot of weight vs height, coloring by gender
ggplot(data=clinical_data, aes(x=weight, y=height, color=gender)) +
  geom_point()


# Compare plot generated earlier with base R with ggplot2:
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

# Customizing with ggplot2
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
# Facet wrap: separate the plots according to a categorical (factor) variable
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
p1 <- ggplot(data=clinical_data, aes(x=weight, y=height, color=gender)) +
  geom_point() +
  scale_color_manual(values=c("female"="orange", "male"="blue")) + 
  ggtitle("Weight vs Height in Patients") +
  xlab("Weight [kg]") + ylab("Height [m]") +
  theme_bw() 

# Create a separate plot for males and females:
# Facet wrap: separate the plots according to a categorical (factor) variable
p2 <- ggplot(data=clinical_data, aes(x=weight, y=height, color=gender)) +
  geom_point() +
  scale_color_manual(values=c("female"="turquoise", "male"="plum")) + 
  ggtitle("Weight vs Height in Patients") +
  xlab("Weight [kg]") + ylab("Height [m]") +
  theme_bw() 

# install.packages("cowplot")
library(cowplot)
plot_grid(p1, p2, nrow=1)





