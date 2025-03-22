
#2021. 1. 20 
# Staitistical Analysis of Microbiome data R code

ng <- layout(matrix(c(1,3,2,3),2,2, byrow=TRUE), widths=c(5,2),height=c(3,4))
layout.show(ng) 
ng1 <- layout(matrix(c(1,2,3,4),2,2, byrow=TRUE), widths=c(5,2),height=c(3,4))
layout.show(ng1) 

group_s <- ifelse(iris$Species %in% "setosa",1, + ifelse(iris$Species %in% "versicolor",2,3))
re <- gsub(pattern = "\\$", replacement = ".", "metacharacters$uses$in$regular$expressions")

tax[1:3]
#[1] Root;p__Proteobacteria;c__Betaproteobacteria;o__Neisseriales;
#f__Neisseriaceae;g__Neisseria
#[2] Root;p__Bacteroidetes;c__Flavobacteria;o__Flavobacteriales;
#f__Flavobacteriaceae;g__Capnocytophaga
#[3] Root;p__Actinobacteria;c__Actinobacteria;o__Actinomycetales;
#f__Actinomycetaceae;g__Actinomyces

tax_1 <- gsub(".+__", "", tax)
tax_1[1:3]

grep("[wd]", c("Sepal.Length", "Sepal.Width", "Petal.Length",  Petal.Width","Species"))
grep("[wd]", c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width","Species"), value = TRUE)
grep("Width", c("Sepal.Length", "Sepal.Width", "Petal.Length","Petal.Width","Species"), value = TRUE,fixed =TRUE) 

fecal <- grep("drySt", colnames(tab))

4.2 Introduction to the dplyr Package
The dplyr package provides a set of functions for efficiently manipulating datasets
in R.
# Select all columns except female
head(select(tab, -female))


# Select all columns except those from female to prog (inclusive)
head(select(tab, -(female:prog )))

#select(), such as starts_with(), ends_with(), matches(), contains(), and one_of ().

# select all columns that start with the character string "s"
head(select(tab, starts_with("s")))

#Selecting rows using filter()
# Filter the rows for students with reading score greater than or equal 70.
filter(tab, read >= 70)

#Filter the rows for students with both reading and math scores greater than or equal 70
 filter(tab, read >= 70, math >= 70)

#Create new columns using mutate()
#Calculate average read and write scores
 head(mutate(tab, avg_read = sum(read)/n()))

#To keep only the new variables, use transmute()
> head(transmute(tab,avg_read = sum(read)/n()))

#Create new columns using mutate() and pipe operator
> tab %>% mutate(avg_read = sum(read/n())) %>%

#To collapses a data frame to a single row.
summarise(tab, avg_read = mean(read, na.rm = TRUE))

#Create summaries of the data frame using summarise() and pipe operator
tab %>% summarise(avg_read = mean(read), + min_read = min(read),+ max_read = max(read),+ n = n())

Grouping Observations Using group_by()
#First group by gender, and then get the summary statistics of reading by g
ender

by_gender <- group_by(tab, female)
read_by_gender <- summarise(by_gender,+ n = n(),  avg_read = mean(read, na.rm = TRUE), min_read = min(read,na.rm = TRUE),
max_read = max(read,na.rm = TRUE)) 


#The algorithm of ggplot2
A layer is composed of four parts: data and aesthetic mapping, a statistical
transformation (stat), a geometric object (geom), and a position adjustment
(Wickham 2010).
In ggplot2, there are two main high-level functions to create a plot: qplot(), short for
quick plot, and ggplot().
With qplot(), a plot is created in a way of all at once; with
ggplot(), a plot is created piece-by-piece and layer functions.

library(ggplot2)
 ggplot() +
 layer(	data = iris, mapping = aes(x = Sepal.Width, y = Sepal.Length),
	geom = "point", stat = "identity", position = "identity") +
 scale_y_continuous() +
scale_x_continuous() +
coord_cartesian()
 #vs
plot(iris$Sepal.Width, iris$Sepal.Length,cex=0.7,pch=16)

qplot(Sepal.Width, Sepal.Length, data = iris, geom = c("point", "smooth"), method = "lm", log = "xy")

library(ggplot2)
 p <- ggplot(iris, aes(x=Sepal.Width, y=Sepal.Length))
 # Sepal.Width and Sepal.Length are columns in iris dataframe
p1 <- p + geom_point()
ggplot(iris, aes(x=Sepal.Width, y=Sepal.Length)) + geom_point()

#Add smoothing geom (layer2)
 p2 <- p1 + geom_smooth(method="lm")
summary(p2)

p3 <- ggplot(iris, aes(x=Sepal.Width, y=Sepal.Length)) +
 #Add scatterplot geom (layer1)
 geom_point(col="blue", size=3) +
 #Add smoothing geom (layer2)
 geom_smooth(method="lm",col="red",size=2)

p4 <- ggplot(iris, aes(x=Sepal.Width, y=Sepal.Length)) +
 #Add scatterplot geom (layer1)
 geom_point(aes(col=Species), size=3) +
 #Add smoothing geom (layer2)
geom_smooth(method="lm",col="red",size=2)

p5 <- p4 + coord_cartesian(xlim=c(2.2,4.2), ylim=c(4, 7)) # zooms in

p6 <- p5 + labs(title="Sepal width vs sepal length", subtitle="Using iris dataset", y="Length of Sepal", x="Width of Sepal")
print(p6)#Or plot(p6)

#Add Title and Labels using ggtitle(), xlab() and ylab()
 p7 <-p5 + ggtitle("Sepal width vs sepal length", subtitle="Using iris dataset") + ylab("Length of Sepal") + xlab("Width of Sepal")
 print(p7)

#Spliting plots by rows
 ggplot(iris, aes(x=Sepal.Width, y=Sepal.Length)) +
 geom_point(aes(col=Species), size=3) +
 geom_smooth(method="lm",col="red",size=2) +
 coord_cartesian(xlim=c(2.2,4.2), ylim=c(4, 7)) +
 # Add Facet Grid
 facet_grid(Species ~.)

#Spliting plots by columns
 ggplot(iris, aes(x=Sepal.Width, y=Sepal.Length)) +
 geom_point(aes(col=Species), size=3) +
 geom_smooth(method="lm",col="red",size=2) +
 coord_cartesian(xlim=c(2.2,4.2), ylim=c(4, 7)) +
 # Add Facet Grid
 facet_grid(.~ Species)

#Spliting plots by columns containing all the data
> ggplot(iris, aes(x=Sepal.Width, y=Sepal.Length)) +
 geom_point(aes(col=Species), size=3) +
 geom_smooth(method="lm",col="red",size=2) +
 coord_cartesian(xlim=c(2.2,4.2), ylim=c(4, 7)) +
 # Add Facet Grid
facet_grid(.~Species, margin=TRUE)

The main difference between wrap faceting and grid faceting is that with wrap
faceting, it is possible to choose the number of rows and columns in the grid.

#Facet Wrap
 #Splitting plots by columns
 ggplot(iris, aes(x=Sepal.Width, y=Sepal.Length)) +
 geom_point(aes(col=Species), size=3) +
 geom_smooth(method="lm",col="red",size=2) +
 coord_cartesian(xlim=c(2.2,4.2), ylim=c(4, 7)) +
 #Add Facet Wrap
 facet_wrap(~ Species, nrow=2)


#5.2.2 Diversity Data for ALS Study
##load abundance table
 abund_table=read.csv("ALSG93AGenus.csv",row.names=1,check.names=FALSE)
 abund_table_t<-t(abund_table)

library(vegan)
 #use the diversity function (vegan package) to calculate Shannon index
 #make a data frame of Shannon index
 H<-diversity(abund_table_t, "shannon")
 df_H<-data.frame(sample=names(H),value=H,measure=rep("Shannon",length(H)))

#Obtain grouping information from sample data
 df_H$Group <- with(df_H,
 ifelse(as.factor(sample)%in% c("A11-28F","A12-28F","A13-28F","A14-28F","A15-28F","A16-28F"),c("G93m1"),
 ifelse(as.factor(sample)%in% c("A21-28F","A22-28F","A23-28F","A24-28F","A25-28F","A26-28F"),c("WTm1"),
 ifelse(as.factor(sample)%in% c("C11-28F","C12-28F","C13-28F"),c("G93m4"),
 ifelse(as.factor(sample)%in% c("C21-28F","C22-28F","C23-28F"),c("WTm4"),
 ifelse(as.factor(sample)%in% c("B11-28F","B12-28F","B13-28F","B14-28F","B15-28F","D11-28F","D12-28F","D13-28F","D14-28F"),c("BUm3to3.5"),
c("NOBUm3to3.5")))))))

The whole data set include sample data from months 1, 3, 3.5 and 4. 
We interest in comparisons of treatment versus not treatment during 3 to 3.5 months. 
So we subset the data from 3 to 3.5 months as below:

library(dplyr)
 df_H_G6 <- select(df_H, Group,value)
 df_H_G93BUm3 <- filter(df_H_G6,Group=="BUm3to3.5"|Group=="NOBUm3to3.5")

5.2.3 Calculating Power or Sample Size Using R Function power.t.test()
where, n is the number of sample size per group, delta is true difference in
means, sd is the standard deviation, sig.level is the significance level (Type I error
probability), power is the power of test (1 minus Type II error probability), type is
the type of t test, and alternative is one-or-two sided test.

aggregate(formula = value ~Group, data = df_H_G93BUm3, FUN = var)
n1 <- 9;  n2 <-7;  s1<-sqrt(0.02892) ; s2<-sqrt(0.04349)
 s=sqrt((n1-1)*s1^2+(n2-1)*s2^2)/(n1+n2-2)

power.t.test(n=2:10,delta=2.504-2.205,sd=0.05012)

plot_richness(physeq, measures = c("Chao1", "Shannon"), x = "Group", color= "Group")
plot_bar(physeq, x=“Sample”, y=“Abundance”, fill=NULL, title=NULL,facet_grid=NULL)

plot_bar(physeq, fill="Group")  # fill different colors for Vdr−/− and WT groups to which each sample belongs.

TopNGenus <- names(sort(taxa_sums(physeq), TRUE)[1:5])
Top5Genus <- prune_taxa(TopNGenus,physeq)
plot_bar(Top5Genus, fill="Group", facet_grid=~Group)

plot_heatmap(Top5Genus)
(p <- plot_heatmap(Top5Genus, "NMDS", "bray"))  # NMDS ordination method and Bray-Curtis distance method to plot the top five genera.
plot_heatmap(Top5Genus, "PCoA", "bray")  # heatmap uses PCoA ordination on the default Bray-Curtis distance (fig 7.4)

plot_network(g, physeq=NULL, type=“samples”, color=“Group”, shape=“Group”)
plot_net(physeq, distance = “bray”, type = “samples”, maxdist = 0.7, color =NULL, shape = NULL)
ig <- make_network(physeq, max.dist=0.8)
plot_network(ig, physeq, color="Group", shape="Group")  

plot_tree(physeq, method = “sampledodge”, color = NULL, shape = NULL, ladderize = FALSE)

library(GUniFrac)
data(throat.otu.tab)  # The smokers data set
data(throat.tree)
data(throat.meta)

install.packages("BiocManager")
BiocManager::install("phyloseq")
#McMurdie PJ, Holmes S (2013).
 “phyloseq: An R package for reproducible interactive analysis and graphics of microbiome census data.” 
PLoS ONE, 8(4), e61217.

library(phyloseq)
 #Convert the data to phyloseq format
 OTU = otu_table(as.matrix(throat.otu.tab), taxa_are_rows = FALSE)
 SAM = sample_data(throat.meta)
 TRE <-throat.tree
 physeq<-merge_phyloseq(phyloseq(OTU),SAM,TRE) )
ntaxa(physeq) 

physeq = prune_taxa(taxa_names(physeq)[1:50], physeq) # prune_taxa() function prunes just the first 50 OTUs in the data set.
plot_tree(physeq, ladderize = "left", color = "SmokingStatus") # map color to smoking status variable
plot_tree(physeq, ladderize = "left", color = "SmokingStatus",shape ="SmokingStatus")
plot_tree(physeq,color = "SmokingStatus", shape = "SmokingStatus", ladderize = "left") + coord_polar(theta = "y")

abund_table_norm <- decostand(abund_table, "normalize")
 bc_dist<- vegdist(abund_table_norm , method = "bray")
cluster_single <- hclust (bc_dist, method = 'single') ; plot(cluster_single)
cluster_complete <- hclust (bc_dist, method = 'complete') ;plot(cluster_average)
cluster_ward <- hclust (bc_dist, method = 'ward.D2') ; plot(cluster_ward)

par (mfrow = c(2,2))
> plot(cluster_single)
> plot(cluster_complete)
> plot(cluster_average)
> plot(cluster_ward)


stand_abund_table <- decostand(abund_table, method = "total")
 PCA <-rda(stand_abund_table)
The display option “species” is the vegan package label for OTUs/taxa.
biplot(PCA, display = 'species')

sum (apply (stand_abund_table, 2, var)) #Total variation is a sum of variations of each genus in analyzed matrix.
biplot(PCA, display = 'species') # draw the diagrams using the function biplot(). The display option “species” is the vegan package label for OTUs/taxa.
ordiplot(PCA, display = "sites", type = "text") #to draw both genus and sample scores as centroids as below

##########################################################################3
"cleanplot.pca" <- function(res.pca, ax1=1, ax2=2, point=FALSE,
ahead=0.07, cex=0.7)
{
# A function to draw two biplots (scaling 1 and scaling 2) from an object
# of class "rda" (PCA or RDA result from vegan's rda() function)
#
# License: GPL-2
# Authors: Francois Gillet & Daniel Borcard, 24 August 2012

require("vegan")
par(mfrow=c(1,2))
p <- length(res.pca$CA$eig)
# Scaling 1: "species" scores scaled to relative eigenvalues
sit.sc1 <- scores(res.pca, display="wa", scaling=1, choices=c(1:p))
spe.sc1 <- scores(res.pca, display="sp", scaling=1, choices=c(1:p))
plot(res.pca, choices=c(ax1, ax2), display=c("wa", "sp"), type="n",
main="PCA - scaling 1", scaling=1)
if (point)
{
points(sit.sc1[,ax1], sit.sc1[,ax2], pch=20)
text(res.pca, display="wa", choices=c(ax1, ax2), cex=cex,
pos=3, scaling=1)
}
else
{
text(res.pca, display="wa", choices=c(ax1, ax2), cex=cex,
scaling=1)
}
text(res.pca, display="sp", choices=c(ax1, ax2), cex=cex, pos=4,
col="red", scaling=1)
arrows(0, 0, spe.sc1[,ax1], spe.sc1[,ax2], length=ahead, angle=20,
col="red")
pcacircle(res.pca)
# Scaling 2: site scores scaled to relative eigenvalues
sit.sc2 <- scores(res.pca, display="wa", choices=c(1:p))
spe.sc2 <- scores(res.pca, display="sp", choices=c(1:p))
plot(res.pca, choices=c(ax1,ax2), display=c("wa","sp"), type="n",
main="PCA - scaling 2")
if (point) {
points(sit.sc2[,ax1], sit.sc2[,ax2], pch=20)
text(res.pca, display="wa", choices=c(ax1 ,ax2), cex=cex, pos=3)
}
else
{
text(res.pca, display="wa", choices=c(ax1, ax2), cex=cex)
}
text(res.pca, display="sp", choices=c(ax1, ax2), cex=cex, pos=4,
col="red")
arrows(0, 0, spe.sc2[,ax1], spe.sc2[,ax2], length=ahead,
angle=20, col="red")
}

"pcacircle" <- function (pca)
{
# Draws a circle of equilibrium contribution on a PCA plot
# generated from a vegan analysis.
# vegan uses special constants for its outputs, hence
# the 'const' valuebelow.
eigenv <- pca$CA$eig
p <- length(eigenv)
n <- nrow(pca$CA$u)
tot <- sum(eigenv)
const <- ((n - 1) * tot)^0.25
radius <- (2/p)^0.5
radius <- radius * const
symbols(0, 0, circles=radius, inches=FALSE, add=TRUE, fg=2)
}
#########################################################################3

library(GUniFrac)
data(throat.otu.tab)
data(throat.meta)
library(dplyr)
throat_meta <- select(throat.meta, SmokingStatus, Age, Sex, PackYears)


#7.4.7 Constrained Analysis of Principal Coordinates (CAP)
throat_cap <- capscale(throat.otu.tab ~ SmokingStatus + Sex + PackYears + Condition(Age), throat_meta, dist="bray")
throat_cap
anova(throat_cap)
#The model is statistically significant with p-value of 0.004.

plot(throat_cap)  # it is not informative
groups <- throat.meta$SmokingStatus
groups

plot(throat_cap, type="n")
points(throat_cap, col=as.numeric(as.factor(groups)),pch=as.numeric(as.factor(groups)))
ordispider(throat_cap, groups, lty=2, col="grey", label=T)
ordiellipse(throat_cap, groups, lty=2, col="grey", label=F)

set.seed (123)
anova(throat_cap, step=1000)
anova(throat_cap, by="axis", step=1000)
 #We can see that the first axis is significant with p-value of 0.001 and second axis is marginally significant with p-value of 0.073
anova(throat_cap, by="terms", step=1000)
#both SmokingStatus and Sex are significant with p-values = 0.003, and 0.035, respectively.

#let’s do forward selection using the function ordistep() in vegan package
step_forward <- ordistep(capscale(throat.otu.tab ~ 1, data=throat_meta),scope=formula(throat_cap), direction="forward", pstep=1000)
#The selected variable is Smoking Status. The final model is fitted below.
cap_final<- capscale(throat.otu.tab ~ SmokingStatus, throat_meta, dist="bray")


