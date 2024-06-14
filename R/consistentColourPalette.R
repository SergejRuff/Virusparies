
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

  phylum_mapping <- list()

  for (phylum_name in names(ICTV_data_split)) {
    phylum_color <- unique_phyla[phylum_name]
    unique_families <- unique(ICTV_data_split[[phylum_name]]$Family)
    family_colors <- rep(phylum_color, length(unique_families))
    family_colors_vector[unique_families] <- family_colors

    phylum_mapping[unique_families] <- phylum_name
  }

  names(family_colors_vector) <- sub("viridae$", "", names(family_colors_vector), ignore.case = TRUE)
  names(family_colors_vector)[names(family_colors_vector) == "Pseudo"] <- "PseudoPseudo"
  names(family_colors_vector) <- gsub("^(Allo|Ortho|Pseudo)", "", names(family_colors_vector), ignore.case = TRUE)

  names(phylum_mapping ) <- sub("viridae$", "", names(phylum_mapping ), ignore.case = TRUE)
  names(phylum_mapping )[names(phylum_mapping ) == "Pseudo"] <- "PseudoPseudo"
  names(phylum_mapping ) <- gsub("^(Allo|Ortho|Pseudo)", "", names(phylum_mapping ), ignore.case = TRUE)

  # Get unique families in vh_file[groupby]
  unique_families_in_data <- unique(vh_file[[groupby]])

  # Initialize a vector to store matched colors
  matched_colors <- rep(NA_character_, length(unique_families_in_data))

  # Iterate over the colors in family_colors_vector
  for (i in seq_along(family_colors_vector)) {
    color <- family_colors_vector[i]


    # Check for partial matches between color names and family names
    partial_matches <- grepl(names(color), unique_families_in_data, ignore.case = TRUE)
    #print(partial_matches)

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


  # Iterate over the colors in family_colors_vector
  legend_labels <- matched_vector
  for (i in seq_along(phylum_mapping)) {
    color <- names(phylum_mapping)[i]
    partial_matches <- grepl(color, names(legend_labels), ignore.case = TRUE)
    if (!any(is.na(partial_matches)) && any(partial_matches) && !is.na(phylum_mapping[[color]])) {
      legend_labels[partial_matches] <- phylum_mapping[[color]]
    }
  }

  legend_labels[legend_labels == "#a9a9a9"]<- "unclassified"

  # Extract unique names and corresponding values
  unique_labels <- unique(legend_labels)
  labels <- matched_vector[match(unique_labels, legend_labels)]

  # Create a named vector
  names(labels) <- unique_labels





  return(list(legend_labels = legend_labels, labels = labels))
}

