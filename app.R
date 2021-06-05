#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)
library(mvtnorm)
library(lme4)

#setwd("/Users/syavuz/Box/BDB/Project2/2-level/")


read <- function(file, cycle, N) {
    all <- cbind(read.csv(file), cycle = cycle) %>%
        slice_sample(n = N, replace = F) %>%
        mutate(cid = row_number()) %>%
        pivot_longer(time0:time11, 'time', names_prefix = 'time', values_to = 'y',
                     names_transform = list(time = as.integer))
}

N=200
set.seed(12345)
data <- rbind(read('ecls_11_imputed.csv', 1, N), read('ecls_98_imputed.csv', 2, N)) %>%
    mutate(id = as.integer(factor(paste(cycle, 'X', cid))))
fit <- lmer(y ~ time + I(time ^ 2) + (Male + SES) * (time + I(time ^ 2)) + (1 + time | id), data, subset = cycle == 2,
            control = lmerControl('bobyqa'))
beta.0 <- fit@beta
Sigma <- as.data.frame(VarCorr(fit))$vcov
sigma.1 <- sqrt(Sigma[4])
sigma.2 <- matrix(Sigma[c(1, 3, 3, 2)], nrow = 2)
beta <- beta.0

# beta+abs(t(t(c(.2,.5,.8))) %*% beta)

#for(ses in 1:2){
#t2, male, t2:male
# t=1

#beta[0] <- 50
# beta[3] <- beta[3] + c(-0.8, 0, 0.5)[t]
# 
# beta[5] <- beta[5] + beta[5]
# beta[8] <- beta[8] + beta[8]*2
# beta[9] <- beta[9] + beta[9]*1

#beta[6] <- beta[6] + c(-0.8, 0, 0.5)[t]

#beta[7] <- beta[7] + c(-0.8, 0, 0.5)[t]
#beta[9] <- beta[9] + c(-0.8, 0, 0.5)[t]
#beta[5] <- beta[5]  + beta[5]*c(-5,.5)[ses]

#beta[-1] <- beta[-1] + abs(beta[-1]*.5) 

re <- rmvnorm(N, sigma = sigma.2)
x <- model.matrix(~ time + I(time ^ 2) + (Male + SES) * (time + I(time ^ 2)), data = filter(data, cycle == 1))
summary(fit)
fit2 <- lmer(y ~ time + I(time ^ 2) + (Male + SES) * (time + I(time ^ 2)) + (1 + time | id), data, subset = cycle == 1,
             control = lmerControl('bobyqa'))
data$cycle <-   factor(data$cycle, labels =  c("Current", "Historical"))



# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Historical and Current Cycle BDB project"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            sliderInput("beta1",
                        "Intercept",
                        min = round(beta[1]*.2,2),
                        max = round(beta[1]*1.8,2),
                        value = beta[1]),

        sliderInput("beta2",
                    "Linear Time",
                    min = round(beta[2]*.2,2),
                    max = round(beta[2]*1.8,2),
                    value = beta[2]),
        sliderInput("beta3",
                    "Quadratic Time",
                    min = beta[3]-.8,
                    max = beta[3]+.5,
                    value = beta[3]),
        sliderInput("beta4",
                    "Male",
                    min = round(beta[4]-2,2),
                    max = round(beta[4]+2,2),
                    value = beta[4]),
        sliderInput("beta5",
                    "SES",
                    min = round(beta[5]*.2,2),
                    max = round(beta[5]*1.8,2),
                    value = beta[5]),
   
    ),
    
        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("distPlot"),
           textOutput("selected_var")
        )
    )
)

#beta1 <- seq(beta[1]*.2, beta[1]*1.8, length.out = 12)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$distPlot <- renderPlot({
        # generate bins based on input$bins from ui.R
        beta1 <- seq(input$beta1, input$beta1, length.out = 1)
        beta.0[1] <- beta1
        
        beta2 <- seq(input$beta2, input$beta2, length.out = 1)
        beta.0[2] <- beta2
        
        beta3 <- seq(input$beta3, input$beta3, length.out = 1)
        beta.0[3] <- beta3
        
        beta4 <- seq(input$beta4, input$beta4, length.out = 1)
        beta.0[4] <- beta4
        
        beta5 <- seq(input$beta5, input$beta5, length.out = 1)
        beta.0[5] <- beta5
        
        #the following line might have a problem
        y <- pmax(pmin(x %*% beta.0 + rowSums(x[, 1:2] * matrix(rep(re, each = 6), ncol = 2)) + rnorm(N * 6, sd = sigma.1), 205), 0)
        data$y[data$cycle == "Current"] <- y

        fit2 <- lmer(y ~ time + I(time ^ 2) + (Male + SES) * (time + I(time ^ 2)) + (1 + time | id), data, subset = cycle == "Current",
                     control = lmerControl('bobyqa'))
        betas <- data.frame(rbind(beta,beta.0,fit2@beta))
        rownames(betas) <- c("Historical cycles","Generating coefficients", "Current cycle estimates")
        colnames(betas) <- rownames(summary(fit2)$coefficients)
        print(betas)
        # draw the histogram with the specified number of bins
        ggplot(subset(data), aes(x = time, y = y))+ 
            theme_bw() + geom_jitter(size = 0.01) +
            facet_grid(vars(cycle)) +
            theme(plot.title = element_text(hjust = 0.5), text = element_text(size=18))  + scale_x_continuous(breaks = c(0,1,2,3,7,11)) +
            stat_smooth(method = 'loess', formula = y ~ x, size = 2, se = FALSE, color = "red", linetype = "dashed") +
            stat_smooth( method = 'lm', formula = y ~ x + I(x^2), size = .05, se = FALSE, mapping = aes(x = time, y = y, group = ChildID))
        },height = 800,width = 1200)
  
}

# Run the application 
shinyApp(ui = ui, server = server)

