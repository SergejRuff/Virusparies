#' @title VhgRunsTable: Generate a gt summary table of unique runs for each virus group
#'
#' @description VhgRunsTable generates a summary table of unique runs for each virus group
#' based on the input data set.
#'
#' @param vh_file A data frame containing the Virushunter hittables results.
#' @param groupby (optional): A character specifying the column containing the groups (default: "best_query").
#' Note: Gatherer hittables do not have a "best_query" column. Please provide an appropriate column for grouping.
#' @param taxa_rank (optional): When `groupby` is set to "ViralRefSeq_taxonomy", specify the taxonomic rank to group your data by.
#' Supported ranks are:
#' - "Subphylum"
#' - "Class"
#' - "Subclass"
#' - "Order"
#' - "Suborder"
#' - "Family" (default)
#' - "Subfamily"
#' - "Genus" (including Subgenus)
#' @param cut (optional): A numeric value representing the cutoff for the refseq E-value (default: 1e-5).
#' Removes rows in file with values larger than cutoff value in "ViralRefSeq_E" column.
#' @param title (optional): The title of the plot (default: "Summary of unique runs by virus group").
#' @param title_align (optional): A character vector specifying the alignment of title (and subtitle) text.
#' Possible values are "left" (default), "center", or "right".
#' @param names_ (optional): A vector of length 3 containing column names (default: c("Virus Group","Number of Unique SRA Runs","SRAs Found")).
#' @param align (optional): A character vector specifying the alignment of text in the table columns.
#' Possible values are "left" (default), "center", or "right".
#' @param subtitle (optional): A character specifying the subtitle of the plot (default: NULL).
#' @param data_row.pad (optional): Numeric value specifying the row padding (default: 6).
#' @param column_colour (optional): character specifying the background colour for the column header (default: "dodgerblue4").

#' @param title_size (optional): The size of the title text (default: 26).
#' @param subtitle_size (optional): Numeric specifying the size of the subtitle text (default: 14).
#' @param title_weight (optional): Character or numeric value specifying title font weight.
#' The weight of the font can be modified thorough a text-based option such as "normal", "bold" (default),
#' "lighter", "bolder", or, a numeric value between 1 and 1000, inclusive.
#'
#' @param title_colour (optional): A character specifying the color for the title text (default: "dodgerblue4").
#' @param table_font_size  (optional): Numeric value specifying table font size. This will change font
#' size for the column header and for all values in each cell (default: 14).
#'
#' @param cell_colour (optional): Character specifying cell colour (default: "grey90").
#'
#' @param col_everyrow (optional): Bool value specifying if every row or every second row of the table
#' should be filled with the colour from the cell_colour argument. col_everyrow = TRUE colors every
#' row and col_everyrow = FALSE (default) colors every second row.
#'
#' @return A formatted gt table summarizing unique runs for each virus group
#'
#' @details
#' VhgRunsTable calculates the number of unique runs for each virus group from the input data set.
#' It takes VirusHunter hittables as input and determines how many runs (SRA runs or local FASTQ files) found a specific virus family.
#'
#' A graphical table is returned with two columns. The first column contains the name of the virus group,
#' while the second column contains all IDs of SRA runs or local files that found that virus group.
#'
#' @examples
#' path <- system.file("extdata", "virushunter.tsv", package = "Virusparies")
#' vh_file <- ImportVirusTable(path)
#'
#' # example 1: generate table with defaul arguments
#' table <- VhgRunsTable(vh_file,cut = 1e-5)
#'
#' table
#'
#' # example 2: generate table with custom arguments
#'
#' table_2 <- VhgRunsTable(vh_file,title = "test",title_align="right",
#' names_ = c("column_1","column_2","column_3"),align = "right",subtitle="subtitlele")
#'
#' table_2
#'
#' # example 3: virusgatherer example
#' # import gatherer files
#' path2 <- system.file("extdata", "virusgatherer.tsv", package = "Virusparies")
#' vg_file <- ImportVirusTable(path2)
#'
#'
#' table_3 <- VhgRunsTable(vg_file,groupby = "ViralRefSeq_taxonomy")
#'
#' table_3
#'
#' @seealso \code{\link{VhgRunsBarplot}}
#'
#' @importFrom dplyr group_by summarise n_distinct across any_of
#' @importFrom gt gt
#' @importFrom rlang .data as_string ensym
#' @export
VhgRunsTable <- function(vh_file,
                         groupby = "best_query",
                         taxa_rank = "Family",
                         cut = 1e-5,
                         title="Summary of unique runs by virus group",
                         title_align = "left",
                         names_=NULL,
                         align = "left",
                         subtitle =NULL,
                         data_row.pad=6,
                         column_colour="dodgerblue4",
                         title_size = 26,
                         subtitle_size=14,
                         title_weight="bold",
                         title_colour = "dodgerblue4",
                         table_font_size = 14,
                         cell_colour="grey90",
                         col_everyrow=FALSE){


  if (is_file_empty(vh_file)) {
    #message("Skipping VhgBoxplot generation due to empty data.")
    return(invisible(NULL))  # Return invisible(NULL) to stop further execution
  }

  # Convert non-quoted column names to strings
  groupby <- rlang::as_string(rlang::ensym(groupby))
  taxa_rank <- rlang::as_string(rlang::ensym(taxa_rank))

  if(groupby == "ViralRefSeq_taxonomy"){

    vh_file <- VhgPreprocessTaxa(vh_file,taxa_rank)

  }

  if(is.null(names_)){
    names_ <- c("Virus Group","Number of Unique Runs","Runs Found")
  }

  vh_file <- vh_file[vh_file$ViralRefSeq_E < cut,]

  if (is_file_empty(vh_file)) {
    #message("Skipping VhgBoxplot generation due to empty data.")
    return(invisible(NULL))  # Return invisible(NULL) to stop further execution
  }

  which_runs <- vh_file %>%
    group_by(.data[[groupby]]) %>%
    summarise(unique_SRA_run = n_distinct(across(any_of(c("SRA_run", "run_id")))),
              SRAs_found = paste(unique(across(any_of(c("SRA_run", "run_id")))), collapse = ", "))

  # Remove "c()", "\", and quotation marks
  which_runs$SRAs_found <- gsub("[c()\\\\\"]", "", which_runs$SRAs_found)

  which_runs <- setNames(which_runs, names_)

  row_number <- nrow(which_runs)

  # creat a table
  which_runs_table <- which_runs %>%
    gt()%>%
    tab_header(
      title = title
    )%>%
    cols_align(
      align = align
    )%>%
    opt_align_table_header(align = title_align)

  if(!is.null(subtitle)){
    which_runs_table <- which_runs_table %>%
      tab_header(
      title = title,
      subtitle = subtitle # Add the subtitlele
    )%>%
      tab_options(
        heading.subtitle.font.size = px(subtitle_size)

      )
  }

  which_runs_table <- which_runs_table %>%
    tab_options(
      data_row.padding = px(data_row.pad),
      column_labels.background.color = column_colour,
      heading.title.font.size = px(title_size)

    )%>%
    tab_style(
      style = cell_text(
        color = title_colour,
        weight = title_weight
      ),
      locations = cells_title(groups = "title")
    ) %>%
    tab_options(table.font.size = px(table_font_size))

  if(col_everyrow==FALSE){
    which_runs_table <- which_runs_table%>%
      tab_style(
        style = cell_fill(color = cell_colour),
        locations = cells_body(rows = seq(1,row_number,2))
      )
  }

  if(col_everyrow==TRUE){
    which_runs_table <- which_runs_table%>%
      tab_style(
        style = cell_fill(color = cell_colour),
        locations = cells_body(rows = seq(1,row_number))
      )
  }


  return(which_runs_table)

}
