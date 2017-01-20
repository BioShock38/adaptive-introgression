ui <- fluidPage(
  sliderInput(inputId="i",label="PC to be displayed on the x-axis",value=1,min=1,max=k),
  sliderInput(inputId="j",label="PC to be displayed on the y-axis",value=2,min=1,max=k),
  sliderInput(inputId="min.window",label="Window min. value",value=1,min=1,max=(ncol(residuals)-1)),
  sliderInput(inputId="max.window",label="Window max. value",value=ncol(residuals),min=2,max=ncol(residuals)),
  
  fluidRow(
    splitLayout(cellWidths = c("50%","50%"), plotOutput(outputId="scores"), plotOutput("window.scores")),
    plotOutput(outputId = "stat")
  )
)