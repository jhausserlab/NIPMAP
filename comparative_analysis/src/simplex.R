library(plotly)

## Create a color map assigned to long vs short survivors for the simplex
create_color_map <- function(long_survivors4000, pca3D) {
  site_groups <- rep(1, length(pca3D[1,]))
  value_to_set <- 2
  
  for (i in seq_along(long_survivors4000)) {
    start_index <- 100 * (long_survivors4000[i] - 1) + 1
    end_index <- 100 * long_survivors4000[i]
    site_groups[start_index:end_index] <- value_to_set
  }
  
  map_site_group_to_color <- function(site_group) {
    if (site_group == 1) {
      return("rgba(200, 0, 0, 0.2)")
    } else if (site_group == 2) {
      return("green")
    }
  }
  
  color <- sapply(site_groups, map_site_group_to_color)
  return(color)
}


## Plot simplex with long vs short survivors
plot_simplex <- function(pca3D, Archs_3D, color, custom_nichesLabels, colNiches.hex) {
  # Define simplex vertices coordinates
  simplex_x <- c(Archs_3D[1, 1], Archs_3D[2, 1], Archs_3D[3, 1], Archs_3D[4, 1], Archs_3D[1, 1])
  simplex_y <- c(Archs_3D[1, 2], Archs_3D[2, 2], Archs_3D[3, 2], Archs_3D[4, 2], Archs_3D[1, 2])
  simplex_z <- c(Archs_3D[1, 3], Archs_3D[2, 3], Archs_3D[3, 3], Archs_3D[4, 3], Archs_3D[1, 3])
  simplex_x <- c(simplex_x, Archs_3D[c(1, 3), 1], Archs_3D[c(2, 4), 1])
  simplex_y <- c(simplex_y, Archs_3D[c(1, 3), 2], Archs_3D[c(2, 4), 2])
  simplex_z <- c(simplex_z, Archs_3D[c(1, 3), 3], Archs_3D[c(2, 4), 3])
  
  plotly::plot_ly(x = pca3D[1, ],
                  y = pca3D[2, ],
                  z = pca3D[3, ],
                  type = "scatter3d",
                  mode = "markers",
                  marker = list(symbol = "triangle", size = 4, color = color),
                  name = "sites",
                  mode = "text") %>%
    add_trace(x = Archs_3D[, 1],
              y = Archs_3D[, 2],
              z = Archs_3D[, 3],
              type = "scatter3d",
              mode = "markers+text",
              text = custom_nichesLabels,
              textposition = c('top right', 'bottom right', 'top left', 'top right'),
              textfont = list(color = '#000000', size = 16),
              showlegend = TRUE,
              name = "niches",
              marker = list(color = ~colNiches.hex, symbol = "star-diamond", size = 12),
              inherit = FALSE) %>%
    add_trace(x = simplex_x,
              y = simplex_y,
              z = simplex_z,
              type = "scatter3d",
              mode = "markers+lines",
              name = "Simplex",
              line = list(color = "blue", width = 2)) %>%
    layout(scene = list(xaxis = list(title = "PC1"),
                        yaxis = list(title = "PC2"),
                        zaxis = list(title = "PC3"))
    )
}

