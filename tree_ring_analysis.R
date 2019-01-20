## Tree-Ring Width Analysis
## Data from Mary Ranch, Santa Clara, California. 
## Blue oak trees widths (1932-1976). Available at:
## ftp://ftp.ncdc.noaa.gov/pub/data/paleo/treering/measurements/northamerica/usa/ca645.rwl
##
## Presented as Final Project: Functional Data Analysis, Master in Mathematical 
## Ingeneering. Universidad Carlos III de Madrid, January, 2019.
##
## Author: Harold Antonio Hernández Roig (hahernan@est-econ.uc3m.es)
##
## This work mimics Section 4. Data Illustrations / 4.1 Stringing of
## Tree-Ring Widhts in [1].
##
## References:
## [1] K Chen, H G Müller, and J Wang. Stringing High-Dimensional Data for 
## Functional Analysis. Journal of the American Statistical Association, 
## 106(493):275–284, 2011.

library(fdapace) # Functional Principal Component Analysis (FPCA) via 
# Principal Analysis by Conditional Estimation (PACE)
# Includes Stringing functionality!
# Usage (vignette): 
# https://cran.r-project.org/web/packages/fdapace/vignettes/fdapaceVig.html


library(dplR)   # Dendrochronology (or tree-ring dating) Program Library in R
# Usage (vignette): 
# https://cran.r-project.org/web/packages/dplR/vignettes/intro-dplR.pdf

## Loading data: 
# data = read.rwl('ca645.rwl') # Error in read.tucson(fname, ...) : precision unknown in series B2707A
data = read.rwl('ca645_2.rwl') # (series B2707A deleted due to "unknown precision")

dim(data) # gives years and number of series
colnames(data) # the series IDs
rwl.report(data) 

plot(data, plot.type="spag") 

## Observations from 1932-1976
years = 1932:1976
data2 = data[as.character(years),]

plot(data2, plot.type="spag") 

width_gain = colSums(data2,na.rm = TRUE)

## Divide each tree measure by total width gained
for (i in 1:42) {
  data2[,i] = data2[,i]/width_gain[i]
}

## Plot tree-rings after scalling...
require(ggplot2)
require(reshape2)

final_data = data2
final_data$year = row.names(data2)  
df = melt(final_data, id= 'year')
names(df) = c('year', 'Tree', 'Width')
ggplot() + 
  geom_line(data = df, aes(x = year, y = Width, color = Tree, group = Tree), size = 0.5) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

## Stringing
stringingfit = Stringing(X = t(as.matrix(data2)))

data4 = matrix(0,45,42)
for (i in 1:42) {
  data4[,i] = unlist(stringingfit$Ly[i])
}
data4 = data.frame(data4)
colnames(data4) = colnames(data2)

## Plot tree-rings after stringing...
final_data = data4
final_data$id = 1:45
df = melt(final_data, id= 'id')
names(df) = c('id', 'Tree', 'Width')
ggplot(df, aes(id, Width)) + 
  geom_line(data = df, aes(x = id, y = Width, color = Tree, group = Tree), size = 0.5)

## Plot of Stringin-function:
df = data.frame(Stringed_Order = stringingfit$StringingOrder,
                Year = years)
ggplot(data=df, aes(x=Year, y= Stringed_Order, group=1)) +
  geom_line()+
  geom_point()

## Load Precipitations 
# (file generated at: https://cefa.dri.edu/Westmap/Westmap_home.php)
# 5-months period total precipitations from Dec. previous to Apr. current Year.

prec_sc = read.delim("precip_Santa_ClaraCA_1932_1976.txt", header = TRUE, sep = "", dec = ".")

df_prec = data.frame(prec_sc)
ggplot(data=df_prec, aes(x=Year, y= Total_Precipitation, group=1)) +
  geom_line()+
  geom_point()

## Now both graphs in the same plot and different y-axis

df_prec$Stringed_Order = df$Stringed_Order

p <- ggplot(df_prec, aes(x = Year))
p <- p + geom_line(aes(y = Total_Precipitation, colour = "Precipitation (dashed)"), linetype = "dashed")  

# adding stringed order
p <- p + geom_line(aes(y = Stringed_Order, colour = "Stringed Order"))

# now adding the secondary axis
p <- p + scale_y_continuous(sec.axis = sec_axis(~.*1, name = "Rainy season total precipitation"))

# modifying colours and theme options
p <- p + scale_colour_manual(values = c("blue", "red"))
p <- p + labs(y = "Stringed Order",
              x = "Year",
              colour = "Comparison of Curves")
p <- p + theme(legend.position = c(0.1, 0.9)) 
p


## Running PPCA on dense-stringed-dataset
L3 <- MakeFPCAInputs(tVec=1:45, yVec = t(as.matrix(data4)), na.rm = TRUE)
FPCAdense <- FPCA(L3$Ly, L3$Lt)

# Plot the FPCA object
plot(FPCAdense)

# Variance explained by PC
FPCAdense$cumFVE # 8 PC needed to explain 95% of variability!

fitted_dense = fitted(FPCAdense, K = 8)
IDs = colnames(data2)

# Plot of Fitted-Curves:
par(mfrow=c(2,2),oma = c(0, 0, 2, 0))
CreatePathPlot(FPCAdense, subset = c(5), K = 8, main = IDs[5], showObs = TRUE, ylim=c(0, 0.06)) ; grid()
CreatePathPlot(FPCAdense, subset = c(10), K = 8, main = IDs[10], showObs = TRUE, ylim=c(0, 0.06)) ; grid()
CreatePathPlot(FPCAdense, subset = c(30), K = 8, main = IDs[30], showObs = TRUE, ylim=c(0, 0.06)) ; grid()
CreatePathPlot(FPCAdense, subset = c(41), K = 8, main = IDs[41], showObs = TRUE, ylim=c(0, 0.06)) ; grid()
mtext("Fitted Tree-Rings (K = 8)", outer = TRUE, cex = 1)


## Running PPCA on dense-stringed-dataset with smoothing!
FPCAsmoothed <- FPCA(L3$Ly, L3$Lt, optns = list(methodMuCovEst = 'smooth', userBwCov = 2))

# Plot the FPCA object
plot(FPCAsmoothed)
# Variance explained by PC
FPCAsmoothed$cumFVE # 2 PC needed to explain 95% of variability!

# Plot of Fitted-Curves:
par(mfrow=c(2,2),oma = c(0, 0, 2, 0))
CreatePathPlot(FPCAsmoothed, subset = c(5), K = 4, main = IDs[5], showObs = TRUE, ylim=c(0, 0.06)) ; grid()
CreatePathPlot(FPCAsmoothed, subset = c(10), K = 4, main = IDs[10], showObs = TRUE, ylim=c(0, 0.06)) ; grid()
CreatePathPlot(FPCAsmoothed, subset = c(30), K = 4, main = IDs[30], showObs = TRUE, ylim=c(0, 0.06)) ; grid()
CreatePathPlot(FPCAsmoothed, subset = c(41), K = 4, main = IDs[41], showObs = TRUE, ylim=c(0, 0.06)) ; grid()
mtext("Fitted Tree-Rings (Smoothed & K = 4)", outer = TRUE, cex = 1)
