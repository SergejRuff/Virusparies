
#' keep colours constistent for each virusgroup
#'
#' @param vh_file vhfile
#' @param groupby x_column
#'
#' @importFrom utils data
#' @keywords internal
consistentColourPalette <- function(vh_file = vh_file, groupby = "best_query") {

  unique_phyla <- c(
    Taleaviricota = "#800000",
    Peploviricota = "#9A6324",
    Uroviricota = "#ffd8b1",
    Hofneiviricota = "#469990",
    Phixviricota = "#000075",
    Cossaviricota = "#e6194B",
    Cressdnaviricota = "#f58231",
    Saleviricota = "#ffe119",
    Ambiviricota = "#bfef45",
    Duplornaviricota = "#3cb44b",
    Kitrinoviricota = "#42d4f4",
    Lenarviricota = "#4363d8",
    Negarnaviricota = "#911eb4",
    Pisuviricota = "#f032e6",
    Artverviricota = "#808000",
    Nucleocytoviricota = "#fffac8",
    Preplasmiviricota = "#aaffc3",
    Dividoviricota = "#dcbeff",
    unclassified = "#a9a9a9"
  )

  # Assuming ICTV_data is properly loaded and formatted
  ICTV_data <- ICTV_data  # Replace with your actual ICTV_data

  # Split the data frame by the 'Phylum' column
  ICTV_data_split <- split(ICTV_data, ICTV_data$Phylum)

  # Initialize family colors vector
  family_colors_vector <- character(length = length(unique(ICTV_data$Family)))
  names(family_colors_vector) <- unique(ICTV_data$Family)

  # Iterate over each sublist of ICTV_data_split
  for (phylum_name in names(ICTV_data_split)) {
    # Get the color corresponding to the current phylum
    phylum_color <- unique_phyla[phylum_name]

    # Get unique families in the sublist
    unique_families <- unique(ICTV_data_split[[phylum_name]]$Family)

    # Map each family to its corresponding color
    family_colors <- rep(phylum_color, length(unique_families))

    # Assign colors to the virus family names
    family_colors_vector[unique_families] <- family_colors
  }

  # After assigning the colors, remove "viridae" from the names of family_colors_vector
  names(family_colors_vector) <- sub("viridae$", "", names(family_colors_vector), ignore.case = TRUE)
  names(family_colors_vector) <- gsub("^(Ortho|Pseudo)", "", names(family_colors_vector), ignore.case = TRUE)

  # Get unique families in vh_file[groupby]
  unique_families_in_data <- unique(vh_file[[groupby]])

  # Initialize a vector to store matched colors
  matched_colors <- rep(NA_character_, length(unique_families_in_data))

  # Iterate over the colors in family_colors_vector
  for (i in seq_along(family_colors_vector)) {
    color <- family_colors_vector[i]

    # Check for partial matches between color names and family names
    partial_matches <- grepl(names(color), unique_families_in_data, ignore.case = TRUE)

    # Assign the color to matched family names
    matched_colors[partial_matches] <- color
  }

  # Handle unmatched virus groups (NA values)
  unmatched <- is.na(matched_colors)
  matched_colors[unmatched] <- unique_phyla["unclassified"]

  # Filter out NA values (no matches) and corresponding family names
  matched_families <- unique_families_in_data[!is.na(matched_colors)]
  matched_colors <- matched_colors[!is.na(matched_colors)]

  # Handle empty strings represented as ""
  empty_strings <- matched_colors == ""
  matched_colors[empty_strings] <- unique_phyla["unclassified"]

  # Create matched_vector with matched colors and their names
  matched_vector <- c(matched_colors)
  names(matched_vector) <- matched_families

  return(matched_vector)
}

