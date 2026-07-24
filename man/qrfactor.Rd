\name{qrfactor}
\alias{qrfactor}
\alias{qrfactor.default}
\alias{print.qrfactor}
\alias{summary.qrfactor}
\alias{plot.qrfactor}
\title{Perfoming Factor Analysis with GIS data}
\description{The anovagis function performs  Factor Analysis on GIS data. }
\usage{
qrfactor(source,layer='',var=NULL,type='',p="Yes",scale="sd",t='',nf=2,m=NULL,f=NULL,...)
\method{qrfactor}{default}(source,layer='',var=NULL,type='',p="Yes",scale="sd",t='',nf=2,m=NULL,f=NULL,...)
\method{print}{qrfactor}(x,...)
\method{summary}{qrfactor}(object,...)
\method{plot}{qrfactor}(x,factors=c(1,2),type="loading",plot="",
cex="",pch=15,pos=3,main="",xlim="optimise",
ylim="optimise",abline=TRUE,legend="topright",legendvalues=c(100),
values=FALSE,nfactors=3,rowname=TRUE,par=c(1,2),...)
}
\arguments{
\item{source}{Folder path of the layer. Please quote the full folder path with forward slash "/".  You can use R object as a source but you must set the layer parameter to "nofile"; see below}
\item{layer}{ The layer qrfactor in the folder that you want to work with. It is the file name of qrfactor. This is case sensitive, please. In case you want to use non spatial data such as ".csv", ".txt", "dat" or ".tab" insert the full file name as layer. In case of using R object as a source set "layer" parameter to "nofile"  }
\item{var}{The attributes or variables of the layer. In case of using non spatial data such as ".csv", ".txt", "dat" or ".tab"  var are variables or column names }
\item{type}{Types of plots 'mds'for multidimensional scale, 'coordinate' for principal coordinate analyse. Or The type of results one wants to plot. It takes "scores", "loadings", pca or eigenvectors. The deault is loadings. }
\item{p}{Determine whether prediction must be done:"Yes". The scores are appended to the GIS data  }
\item{t}{The list of variables that one wants to transform eg. transform=c("gold","diamond")}
\item{scale}{scale the data:"sd","pca","data". The default is "sd" that is the scaled data divided by the standard deviation. It can lso take "log" or "sqrt" and use the default "sd" for normal distribution transformation  }
\item{m}{the the match field: the common variable on both the table and spatial data. This name must be identical to both sets of data}
\item{f}{The full path of csv file and the name of csv eg. C:/Users/owusu/Documents/Rpackages/qrfactor14/inst/external/farms.csv}
\item{x}{an object of class \code{"qrfactor"}, i.e., a fitted model. }
\item{object}{an object of class \code{"qrfactor"}, i.e., a fitted model. }
\item{plot}{The type of plots one desires. It takes "all" for all the 3 plots or "q" for q plot  or "r" for r plot or 'qr' for both q and r plots}
\item{factors}{list of factors one wants to plot. The default is factors=c(1,2). Please do not forget "c" in the list. }
\item{cex}{A numerical value giving the amount by which plotting text and symbols should be magnified relative to the default. It also accepts a vector of values which are recycled eg cex=c("gold")}
\item{nfactors}{The number of factors to extract}
\item{pch}{Either an integer specifying a symbol or a single character to be used as the default in plotting points.}
\item{pos}{The position of text labels}
\item{main}{Main title of the graph}
\item{xlim}{x-coordinates of the axis eg xlim=c(-1.5,1.5)}
\item{ylim}{y-coordinates of the axis eg ylim=c(-1.5,1.5)}
\item{abline}{the intercept and slope, single values of straight lines through the current plot. eg. abline(-0.5,0.5)}
\item{legend}{position of legend: it takes topright,topleft, bottomright,bottomleft, top, left, bottom, right }
\item{legendvalues}{The values of the legend}
\item{values}{Incase one wants to label the graph with another variables. eg. values=c("gold")}
\item{nf}{The number of factors to extract}
\item{rowname}{rownames of the data}
\item{par}{the layout setteing in a form of list}
\item{...}{any other parameter can be added}
}
\value{
Objects of the class that basically list its elements
\item{data}{ Original data for the model. All records must be numeric. It also accepts continous data}
\item{gisdata}{ GIS data for the model incase you use shape files}
\item{x.standard}{ it is the scale matrix of the original data}
\item{correlation}{ The correlation matrix for the data }
\item{eigen.value}{ eigen value of correlation matrix of the data}
\item{eigen.vector}{ eigen vector of correlation matrix of the data}
\item{diagonal.matrix}{ diagonal matrix of eigen vector }
\item{pca}{ pca loadings}
\item{pcascores}{ PCA scores }
\item{r.loading}{ R-mode loadings}
\item{q.loading}{ Q-mode loadings }
\item{loadings}{ combined loadings of R and Q on the same axis}
\item{q.scores}{ computed Q-mode scores }
\item{scores}{ combined R-mode and Q-mode scores on the same axis }
\item{rownames}{ row names of the loadings}
\item{variables}{ variables names of the loadings, of the original data}
}
\author{George Owusu}
\examples{
\dontrun{
#apply qrfactor to csv data

csv= system.file("external", "Africanfreshwater.csv", package = "qrfactor") #list the csv file
var=c( "Domestic", "Industry",   "Agricultur", "Resources",  "withdrawal","perCapitaW")
mod0=qrfactor(csv,var=var)
plot(mod0,rowname="COUNTRY")

#apply qrfactor on shapefile
source<- system.file("external", package = "qrfactor")
layer="Africanfreshwater"
mod1=qrfactor(source,layer,var=var)
plot(mod1,rowname="COUNTRY")

#apply qrfactor on imported spatial data into R
gisdata <- na.omit(readOGR(source, layer))
mod2=qrfactor(gisdata,var=var)

#join CSV data and shapefile
mod3=qrfactor(source,layer,var=var,m="COUNTRY",f=csv)
mod5=qrfactor(mod3$gisdata,var=var,m="COUNTRY",f=csv) #multiple join

par(mfrow=c(1,2))
plot(mod2,rowname="COUNTRY",cex=c("means"),legend="topleft",values=c("cluster"),pch=23)
#plot(mod2,cex=c("means"),type="cluster")# cluster analyses
plot(mod2,type="map")#plots several maps
#plot(mod2,type="diagnose")#plots histograms and qqplots
}

}
\references{Bivand, R. S., Pebesma, E. J., Gomez-Rubio, V. (2008) Applied Spatial Data Analysis with R. Springer
Kabacoff, I. R. (2011) R in Action. Data Analysis and Graphics with R. Manning Publications Co
}
\keyword{Factor Analysis}