library(shiny)
library(ggplot2)
library(ape)
library(cowplot)
library(gridExtra)
library(ggtree)
library(dplyr)
library(tidyr)

# Define UI for app that draws a histogram ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Phylogenetic Random Effects"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Input: Slider for the number of bins ----
      sliderInput(inputId = "treeseed",
                  label = "Tree seed",
                  min = 1,
                  max = 100,
                  value = 1),
      
      sliderInput(inputId = "simseed",
                  label = "Simulation seed",
                  min = 1,
                  max = 100,
                  value = 1),
      
      sliderInput(inputId = "nspp",
                  label = "Number of Species",
                  min = 2,
                  max = 100,
                  value = 3),
      
      sliderInput(inputId = "nrep",
                  label = "Number of repeated measures",
                  min = 1,
                  max = 100,
                  value = 1),
      
      sliderInput(inputId = "resid.sd",
                  label = "Residual SD",
                  min = 0,
                  max = 50,
                  value = 0),
      
      sliderInput(inputId = "phylo.int.sd",
                  label = "Phylogenetic Intercept SD",
                  min = 0,
                  max = 50,
                  value = 0),
      
      sliderInput(inputId = "int.sd",
                  label = "Intercept SD",
                  min = 0,
                  max = 50,
                  value = 0),
      
    sliderInput(inputId = "phylo.slope.sd",
                label = "Phylogenetic Slope SD",
                min = 0,
                max = 50,
                value = 0),
    
    sliderInput(inputId = "slope.sd",
                label = "Slope SD",
                min = 0,
                max = 50,
                value = 0)
    
  ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Histogram ----
      plotOutput(outputId = "distPlot")
      
    )
  )
)

# Define server logic required to draw a histogram ----
server <- function(input, output) {
  
  # Histogram of the Old Faithful Geyser Data ----
  # with requested number of bins
  # This expression that generates a histogram is wrapped in a call
  # to renderPlot to indicate that:
  #
  # 1. It is "reactive" and therefore should be automatically
  #    re-executed when inputs (input$bins) change
  # 2. Its output type is a plot
  dat <- reactive({
    set.seed(input$treeseed)
    nspp <- input$nspp
    nrep <- input$nrep
    
    phy <- rtree(n = nspp)
    Vphy <- ape::vcv(phy)
    physd.y <- input$phylo.int.sd
    physd.x <- input$phylo.slope.sd
    
    sd.y <- input$int.sd
    sd.x <- input$slope.sd
      
    sd.resid <- input$resid.sd
    phycormat <- matrix(c(1,1,1,1),nrow=2)
    physdvec <- c(physd.y, physd.x)
    phyvarmat <- physdvec %*% t(physdvec)
    phycovmat <- phyvarmat * phycormat
    
    physigma <- kronecker(phycovmat, Vphy) 
    
    betas <- c(0,0)
    
    set.seed(input$simseed)
    
    b_phy <- MASS::mvrnorm(n=1, mu=rep(betas,each=nspp), Sigma=physigma)
    
    Y.phy <- rep(head(b_phy,nspp), each = nrep)
    X.phy <- rep(tail(b_phy,nspp), each = nrep) 
    
    cormat <- matrix(c(1,1,1,1),nrow=2)
    sdvec <- c(sd.y, sd.x)
    varmat <- sdvec %*% t(sdvec)
    covmat <- varmat * cormat
    
    sigma_b <- kronecker(covmat, diag(nspp)) 
    
    b <- MASS::mvrnorm(n=1, mu=rep(betas,each=nspp), Sigma=sigma_b)
    
    Y.int <- rep(head(b,nspp), each = nrep)
    X.int <- rep(tail(b,nspp), each = nrep)
    
    print(data.frame(Y.int,X.int))
    
    Y.e <- rnorm(nspp*nrep, sd=sd.resid)
    X.e <- rnorm(nspp*nrep, sd=1)
    
    Y.phyint <- Y.phy + Y.int
    X.phyint <- X.phy + X.int
    X.phyint.e <- X.phyint + X.e
    Y.phyint.e <- Y.phyint + Y.e
    Y <- Y.phyint.e + X.phyint.e
    
    dat <- data.frame(sp = rep(rownames(Vphy), each = nrep)
                      , Y.phy
                      , Y.phyint
                      , Y.phyint.e
                      , Y
                      , X.phy
                      , X.phyint.e
    )
    # x    <- rnorm(1000, mean=input$int.sd, sd=input$phylo.int.sd)
    # dd <- data.frame(x)
    dd2 <- list(phy,dat)
  })
  
  output$distPlot <- renderPlot({
    dd <- dat()
    phy <- dd[[1]]
    Vphy <- ape::vcv(phy)
    dat <- dd[[2]]
    # bins <- seq(min(x), max(x), length.out = input$phylo.int.sd + 1)
    
#     hist(x, breaks = bins, col = "#75AADB", border = "white",
#          xlab = "Waiting time to next eruption (in mins)",
#          main = "Histogram of waiting times")
    
    nspp <- input$nspp
    nrep <- input$nrep
#     
    dat <- dat %>% mutate(sp = factor(sp,levels=rownames(Vphy)), obs=sp)
    
    dat2 <- (dat
             %>% gather(key = sd, value, -c(sp,obs,X.phyint.e, X.phy, Y))
             %>% mutate(sd = factor(sd,levels=c("Y.phy","Y.phyint","Y.phyint.e","Y")))
    )
    
    # print(dat2)
    
    
    dat2 <- dat2 %>% mutate(ll = paste(sp,1:(nspp*nrep),sep="_"))
    
    gg <- (ggplot(dat2, aes(x=sd, y=value, group=ll, colour=obs))
           + geom_point()
           + geom_line(alpha=0.5)
           + theme_bw()
           + theme(legend.position = "none")
    )
    
    gg2 <- (ggplot(dat, aes(y=Y.phy, x="0", color=obs))
            + geom_point()
            + theme_bw()
    )
    
    gg_tree <- (ggtree(phy)+ geom_tiplab())
    
    gg_all <- plot_grid(gg_tree, gg, gg2, nrow=1)
    gg_all
#     gg2
  })
  
}

# Create Shiny app ----
shinyApp(ui = ui, server = server)