library(shiny)
library(plotly)

# Create a sample dataframe with 100 points
set.seed(123)  # For reproducibility
n_points <- 100

df <- data.frame(
  x = runif(n_points, 1, 10),
  y = runif(n_points, 1, 20)
  # x=as.numeric(mydata_ow[[2]]$POSIX.time-first(mydata_ow[[2]]$POSIX.time)),
  # y=mydata_ow[[2]]$CH4dry_ppb
)



# Create an interactive scatter plot with "box select" mode
p <- plot_ly(data = df, x = ~x, y = ~y, type = 'scatter', mode = 'markers+select',
             selectedpoints = list(marker = list(size = 10)))

# Define a Shiny app
ui <- fluidPage(
  plotlyOutput("scatterplot"),
  actionButton("saveSelectedData", "Save Selected Data")
)

server <- function(input, output, session) {
  # Create a list to store selected points
  selected_points <- list()

  # Handle selected points using box select
  observeEvent(event_data("plotly_selected"), {
    event_info <- event_data("plotly_selected")
    if (!is.null(event_info)) {
      selected_indices <- event_info$pointNumber + 1  # Adjust for 1-based indexing
      selected_points <<- c(selected_points, selected_indices)
    }
  })

  # Save the selected points as a new dataframe
  observeEvent(input$saveSelectedData, {
    if (length(selected_points) > 0) {
      selected_data <- df[unlist(selected_points), ]
      assign("selected_data", selected_data, envir = .GlobalEnv)
    }
  })

  output$scatterplot <- renderPlotly({
    p
  })
}

shinyApp(ui, server)


ggplot()+
  geom_point(data = df, aes(x,y))+
  geom_point(data = selected_data, aes(x,y, colour = "selected data"))+
  theme_article()
