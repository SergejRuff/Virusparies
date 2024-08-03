#' @title select_theme
#'
#' @description This function selects different themes for plots.
#'
#' @param theme_choice A character indicating which theme to use. Options include "minimal", "classic", "light", "dark", and "void".
#'
#' @return A ggplot2 theme object to be used in plots.
#' @author Sergej Ruff
#' @details This function is an internal utility function used within the package.
#' It selects the theme to use for the different plots in this package.
#'
#' @importFrom ggplot2 theme_minimal theme_classic theme_light theme_dark theme_void theme_gray
#' theme_bw theme_linedraw theme_test theme
#'
#' @noRd
select_theme <- function(theme_choice) {

  add_dotted_grid_lines <- grepl("_dotted", theme_choice)
  theme_choice <- sub("_dotted$", "", theme_choice)

  base_theme <- theme_minimal()

  if (theme_choice == "classic") {
    base_theme <- theme_classic()
  } else if (theme_choice == "light") {
    base_theme <- theme_light()
  } else if (theme_choice == "dark") {
    return(theme_dark())
  } else if (theme_choice == "void") {
    base_theme <- theme_void()
  } else if (theme_choice == "grey" || theme_choice == "gray") {
    base_theme <- theme_gray()
  } else if (theme_choice == "bw") {
    base_theme <- theme_bw()
  } else if (theme_choice == "linedraw") {
    base_theme <- theme_linedraw()
  } else if (theme_choice == "test") {
    base_theme <- theme_test()
  }

  if (add_dotted_grid_lines) {
    # Add dotted grid lines
    base_theme <- base_theme +
      theme(
        panel.grid.major = element_line(linewidth = 0.5, linetype = "dotted", colour = "grey50"),  # Dotted major grid lines
        panel.grid.minor = element_line(linewidth = 0.25, linetype = "dotted", colour = "grey75")   # Dotted minor grid lines
      )
  }

  return(base_theme)
}
