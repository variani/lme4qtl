#--------------------
# Function ggimage
#--------------------

#' @export
ggimage <- function(x, limits, 
  low = "#132B43", high = "#56B1F7", 
  #low = "#d9d9d9", high = "#000000",
  #low = "#023858", high = "#3182bd",  
  ...) 
{
  x <- as(x, "sparseMatrix")
  x <- as(x, "dgTMatrix")
  data <- data.frame(x = x@x, j = x@j + 1, i = x@i + 1)
  
  if (missing("limits")) {
    limits <- c(min(data$x, 0, na.rm=TRUE), max(data$x, 2, na.rm=TRUE))
  }
  
  nrow <- nrow(x)
  ncol <- ncol(x)
  
  # levelplot(x ~ j + i, data)
  ggplot(data, aes(j, i, fill = x)) + geom_tile() + coord_equal(ratio = 1) +
    scale_fill_continuous(limits = limits, low = low, high = high, na.value = NA) +
    scale_x_continuous(expand = c(0, 0), limits = c(0.5, ncol + 0.5)) +
    scale_y_reverse(expand = c(0, 0), limits = c(nrow + 0.5, 0.5)) +
    theme_bw() + theme(legend.position = "bottom") +
    labs(y = "Row", x = "Column", fill = "")
}


