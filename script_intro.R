require(pcadapt)
require(RcppRoll)
require(MASS)
require(shiny)
require(data.table)
require(Rcpp)
require(shapes)
sourceCpp("~/Documents/thesis/git/adaptive-introgression/aiUtils.cpp")
setwd("~/Documents/thesis/git/Introgression/populus/Original_data_set_ch6_12_15/")

#### pcadapt file: individuals in columns ###
filename <- "comt.chr06.snp.full.final.pcadapt"

#### population file: individuals in rows ###
popfile <- "populus.txt"

################## rs file ##################
#rsfile <- "tmp.map"

################ Parameters #################
K <- 1
ploidy <- 2
min.maf <- 0.001
adm.pop <- "4"
pop.1 <- "1"
pop.2 <- "3"
window.size <- 15000

if ((K %% 2)==1){
  k <- K + 1
} else {
  k <- K
}
#############################################
geno.matrix <- as.matrix(fread(filename))
pop <- as.character(read.table(popfile, header = TRUE)$POP)
na.snp.vec <- apply(geno.matrix, MARGIN = 2, FUN = function(x){sum(x == 9)})
na.ind.vec <- apply(geno.matrix, MARGIN = 1, FUN = function(x){sum(x == 9)})
#rsinfo <- read.table(rsfile,header = FALSE)
nIND <- length(pop)

imptd.geno <- impute_geno(geno.matrix)
#geno.matrix <- t(geno.matrix)
geno.matrix <- t(imptd.geno)

########## Compute global scores ############
x <- pcadapt(imptd.geno, K = k, ploidy = ploidy, min.maf = min.maf)

ss <- NULL
ss$u <- x$scores[,1:k]
ss$v <- x$zscores[x$maf>=min.maf,1:k]/sqrt(length(x$zscores[,1]))
ss$d <- x$singular.values[1:k]
global.scores <- x$scores
D <- array(0,dim=c(k,k))
diag(D) <- ss$d

############ Process geno matrix ############
geno.matrix <- geno.matrix[,x$maf>=min.maf]
geno.matrix[which(geno.matrix==9)] <- NA 
reduced.maf <- x$maf[x$maf>=min.maf]
geno.scale <- scale(geno.matrix,center=T,scale=sqrt(2*(reduced.maf*(1-reduced.maf))))

############# Compute prediction ############
pred <- (cbind(ss$u[,K]*ss$d[K]))%*%(t(ss$v[,K]))
geno.scale[which(is.na(geno.matrix))] <- pred[which(is.na(geno.matrix))]

############# Compute residuals #############
residuals <- geno.scale - pred

############# Compute statistics ############
mean.stat <- apply(residuals[pop==adm.pop,],MARGIN=2,FUN=function(x){mean(x,na.rm=TRUE)})
null <- cov.rob(mean.stat)
norm.stat <- (mean.stat-null$center)/sqrt(null$cov)
enlarged.stat <- c(rep(norm.stat[1],window.size-1),norm.stat)

########## Averaging over windows ###########
final.stat <- roll_mean(enlarged.stat^2,n = window.size,by = 1)

############### Launch Shiny ################
ui <- fluidPage(
  sliderInput(inputId="i",label="PC to be displayed on the x-axis",value=1,min=1,max=k),
  sliderInput(inputId="j",label="PC to be displayed on the y-axis",value=2,min=1,max=k),
  sliderInput("range",label="Range:",min = 1, max = ncol(residuals), value = c(1,ncol(residuals))),
  fluidRow(
    splitLayout(cellWidths = c("50%","50%"), plotOutput(outputId="scores"), plotOutput("window.scores")),
    plotOutput(outputId = "stat")
  )
)

server <- function(input,output){
  output$scores <- renderPlot({
    plot(global.scores[,input$i],global.scores[,input$j],xlab=paste0("PC",input$i),ylab=paste0("PC",input$j),main="Global PCA");
    points(global.scores[pop==pop.1,input$i],global.scores[pop==pop.1,input$j],col="coral2",pch=19);
    points(global.scores[pop==adm.pop,input$i],global.scores[pop==adm.pop,input$j],col="skyblue2",pch=19);
    points(global.scores[pop==pop.2,input$i],global.scores[pop==pop.2,input$j],col="yellowgreen",pch=19)
  })
  output$window.scores <- renderPlot({
    small.x <- reactive({
      window <- (min(input$range)):(max(input$range));
      aux <- geno.scale[,window]%*%ss$v[window,]%*%solve(D);
      proc.aux <- procOPA(global.scores[pop!=adm.pop,],aux[pop!=adm.pop,],scale = FALSE)
      proj.aux <- aux[pop==adm.pop,]%*%t(proc.aux$R)
      aux[pop==adm.pop,] <- proj.aux
      n.adm.IND <- length(which(pop==adm.pop))
      
      for (p in 1:k){
        scale.1 <- sqrt(sum(global.scores[,p]^2))
        scale.2 <- sqrt(sum(aux[,p]^2))
        aux[,p] <- aux[,p]*scale.1/scale.2
      }
      list(scores = aux)
    });
    plot(small.x()$scores[,input$i],small.x()$scores[,input$j],xlab=paste0("PC",input$i),ylab=paste0("PC",input$j),main="Local PCA");
    points(small.x()$scores[pop==pop.1,input$i],small.x()$scores[pop==pop.1,input$j],col="coral2",pch=19);
    points(small.x()$scores[pop==adm.pop,input$i],small.x()$scores[pop==adm.pop,input$j],col="skyblue2",pch=19);
    points(small.x()$scores[pop==pop.2,input$i],small.x()$scores[pop==pop.2,input$j],col="yellowgreen",pch=19);
#     plot(global.scores[,input$i],global.scores[,input$j],xlab=paste0("PC",input$i),ylab=paste0("PC",input$j),main="Global PCA",xlim=c(-0.1,0.1),ylim=c(-0.2,0.2));
#     points(global.scores[pop=="AEN",input$i],global.scores[pop=="AEN",input$j],col="coral2",pch=19);
#     points(global.scores[pop=="CEM",input$i],global.scores[pop=="CEM",input$j],col="skyblue2",pch=19);
#     points(global.scores[pop=="HG",input$i],global.scores[pop=="HG",input$j],col="yellowgreen",pch=19);
#     arrows(aux[pop=="CEM",input$i],aux[pop=="CEM",input$j],global.scores[pop=="CEM",input$i],global.scores[pop=="CEM",input$j],code=1,length=0.15)
  })
  output$stat <- renderPlot({
    plot(seq(1,ncol(residuals),by=20),final.stat[seq(1,ncol(residuals),by=20)],main="Window of interest",xlab="SNP index",ylab="Stat",cex=0.2,pch=19);
    abline(v=min(input$range),col="red",cex=5);
    abline(v=max(input$range),col="red",cex=5);
  })
}

shinyApp(ui=ui,server=server)

