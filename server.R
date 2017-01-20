server <- function(input,output){
  output$scores <- renderPlot({
    plot(global.scores[,input$i],global.scores[,input$j],xlab=paste0("PC",input$i),ylab=paste0("PC",input$j),main="Global PCA");
    points(global.scores[pop=="AEN",input$i],global.scores[pop=="AEN",input$j],col="coral2",pch=19);
    points(global.scores[pop=="CEM",input$i],global.scores[pop=="CEM",input$j],col="skyblue2",pch=19);
    points(global.scores[pop=="HG",input$i],global.scores[pop=="HG",input$j],col="yellowgreen",pch=19)
  })
  output$window.scores <- renderPlot({
    small.x <- reactive({
      window <- (input$min.window):(input$max.window);
      aux <- geno.scale[,window]%*%ss$v[window,]%*%solve(D);
      small.mean.ind <- apply(aux[pop!=adm.pop,],MARGIN=2,FUN = function(X){mean(X,na.rm=TRUE)})
      for (p in (1:ncol(aux))){
        aux[,p] <- aux[,p]*mean.ind[p]/small.mean.ind[p]
      }
      list(scores = aux)
    });
    plot(small.x()$scores[,input$i],small.x()$scores[,input$j],xlab=paste0("PC",input$i),ylab=paste0("PC",input$j),main="Local PCA");
    points(small.x()$scores[pop=="AEN",input$i],small.x()$scores[pop=="AEN",input$j],col="coral2",pch=19);
    points(small.x()$scores[pop=="CEM",input$i],small.x()$scores[pop=="CEM",input$j],col="skyblue2",pch=19);
    points(small.x()$scores[pop=="HG",input$i],small.x()$scores[pop=="HG",input$j],col="yellowgreen",pch=19)
  })
  output$stat <- renderPlot({
    plot(seq(1,ncol(residuals),by=20),final.stat[seq(1,ncol(residuals),by=20)],main="Window of interest",xlab="SNP index",ylab="Stat",cex=0.2,pch=19);
    abline(v=input$min.window,col="red",cex=5);
    abline(v=input$min.window,col="red",cex=5);
    abline(v=input$max.window,col="red",cex=5);
  })
}

shinyApp(ui=ui,server=server)


