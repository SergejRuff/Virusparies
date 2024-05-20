#' @title select_theme
#'
#' @description This function selects different themes for plots.
#'
#' @param theme_choice A character indicating which theme to use. Options include "minimal", "classic", "light", "dark", and "void".
#'
#' @return A ggplot2 theme object to be used in plots.
#'
#' @details This function is an internal utility function used within the package.
#' It selects the theme to use for the different plots in this package.
#'
#' @importFrom ggplot2 theme_minimal theme_classic theme_light theme_dark theme_void theme_gray
#' theme_bw theme_linedraw theme_test
#'
#' @keywords internal
select_theme <- function(theme_choice) {
  theme_selected <- theme_minimal()
  if (theme_choice == "classic") {
    theme_selected <- theme_classic()
  } else if (theme_choice == "light") {
    theme_selected <- theme_light()
  } else if (theme_choice == "dark") {
    theme_selected <- theme_dark()
  } else if (theme_choice == "void") {
    theme_selected <- theme_void()
  } else if (theme_choice == "grey" || theme_choice == "gray") {
    theme_selected <- theme_gray()
  } else if (theme_choice == "bw") {
    theme_selected <- theme_bw()
  } else if (theme_choice == "linedraw") {
    theme_selected <- theme_linedraw()
  } else if (theme_choice == "test") {
    theme_selected <- theme_test()
  }
  return(theme_selected)
}
