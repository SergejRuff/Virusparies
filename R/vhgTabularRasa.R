#' @title VhgTabularRasa: Function for Generating Custom Formatted Graphical Tables
#'
#' @description This function creates a formatted table using the `gt` package.
#'
#' @param file A data frame.
#' @param title (optional): a custom title.
#' Default is "Graphical Table".
#' @param title_align (optional): a character vector specifying the alignment of title (and subtile) text.
#' Possible values are "left", "center", or "right". Default is "left".
#' @param names_ (optional): a vector containing column names, with a length matching the number
#' of columns in the file argument. Default is names(file)
#' @param align (optional): a character vector specifying the alignment of text in the table columns.
#' Possible values are "left", "center", or "right". Default is "left".
#' @param subtit (optional): a character vector specifying the subtitle. Default is NULL (no subtitle).
#' @param data_row.pad (optional): numeric value specifying the row padding. Default is 6.
#' @param column_colour (optional): character specifying the background colour for the column header.
#' Default is dodgerblue4.
#' @param title_size (optional): numeric value specifying title size. Default is 26.
#' @param subtitle_size (optional): numeric value specifying subtitle size. Default is 14.
#' @param title_weight (optional): character or numeric value specifying title font weight.
#' The weight of the font can be modified thorough a text-based option such as "normal", "bold",
#' "lighter", "bolder", or, a numeric value between 1 and 1000, inclusive.Default value is "bold".
#' @param title_colour (optional): character specifying title colour. Default is dodgerblue4.
#' @param table_font_size  (optional): numeric value specifying table font size. This will change font
#' size for the column header and for all values in each cell. Default is 14.
#' @param cell_colour (optional): character specifying cell colour.Default is grey90.
#' @param col_everyrow (optional): bool value specifying if every row or every second row of the table
#' should be filled with the colour from the cell_colour argument. col_everyrow = TRUE colours every
#' row and col_everyrow = FALSE colours every second row. Default is FALSE.
#'
#' @return Returns a `gt` table object formatted according to the specified parameters.
#'
#' @details This function creates a formatted table using the `gt` package, based on input data with specified column names.
#' It is particularly useful for generating tables that cannot be produced with `vhRunsTable`,
#' when the input data does not originate from the `vhRunsBarplot` functions.
#'
#' The `VhgTabularRasa` function allows users to generate tables styled like Virusparies tables using their own input data.
#' Additionally, users can create custom tables by adjusting the parameters within this function.
#'
#' `VhgTabularRasa` just like `vhRunsTable` empowers users to tailor their tables to suit their needs.
#' With its customizable features, users can effortlessly modify titles, subtitles, and column names.
#' By default, column names are derived from the dataframe's structure using `names(file)`.
#' If the dataframe lacks column names, an error message is triggered. However, users can
#' supply their own column names via the `names_` argument, requiring a vector of names matching the dataframe's column count.
#'
#' `VhgTabularRasa` takes dataframes as Input. Passing any other object type results in an error message.
#' Users can fine-tune their tables with options to adjust text attributes, alignment, and size, as well as row padding
#' and color schemes for titles, subtitles, columns, and backgrounds.
#' If that is not enough, `VhgTabularRasa` returns an gt tables object, which
#' can be further manipulated with the `gt` package available on CRAN.
#' @author Sergej Ruff
#' @examples
#' path <- system.file("extdata", "virushunter.tsv", package = "Virusparies")
#' vh_file <- ImportVirusTable(path)
#'
#' # plot 1: plot boxplot for "identity"
#' identity <- VhgBoxplot(vh_file,y_column = "ViralRefSeq_ident")
#'
#' # generate table
#' VhgTabularRasa(identity$summary_stats)
#'
#'
#' # example 2: plot part of Vh_file (could be any other table)
#' # using first 10 rows of SRA_run,num_hits,bestquery,ViralRefSeq_E and Identity col.
#' vh_file_part <- vh_file[c(1:10),c(1,7,9,10,11)]
#'
#' VhgTabularRasa(vh_file_part,title = "first 10 rows of vh_file",subtit =
#' "example for any table",names_ = c("Runs","Number of Contigs","Best Query Result",
#' "Refrence E-Value","Refrence Identity"))
#'
#'
#'
#' @references
#'
#' Iannone R, Cheng J, Schloerke B, Hughes E, Lauer A, Seo J(2024). _gt: Easily Create Presentation-Ready DisplayTables_.
#' R package version 0.10.1,\url{https://CRAN.R-project.org/package=gt}.
#'
#' @seealso
#' VirusHunterGatherer is available here: \url{https://github.com/lauberlab/VirusHunterGatherer}.
#'
#' @import gt
#' @importFrom stats setNames
#' @export
VhgTabularRasa <- function(file,title="Graphical Table",
                               title_align = "left",names_=NULL,align = "left",subtit =NULL,
                               data_row.pad=6,column_colour="dodgerblue4",title_size = 26,subtitle_size=14,
                               title_weight="bold",title_colour = "dodgerblue4",table_font_size = 14,
                               cell_colour="grey90",col_everyrow=FALSE){

  # check file and arguments
  check_is_dataframe(file)
  #is_file_empty(file)
  if (is_file_empty(file)) {
    #message("Skipping VhgBoxplot generation due to empty data.")
    return(invisible(NULL))  # Return invisible(NULL) to stop further execution
  }
  has_columnnames(file)
  arg_character(title_align)
  arg_character(align)
  arg_character(column_colour)
  arg_character(title_colour)
  arg_character(cell_colour)
  arg_logical(col_everyrow)

  if(is.null(names_)){
    names_ <- names(file)
  }

  file <- setNames(file, names_)

  # Create gt table
  gt_table <- file %>%
    gt()%>%
    tab_header(
      title = title
    )%>%
    cols_align(
      align = align
    )%>%
    opt_align_table_header(align = title_align)

  if(!is.null(subtit)){
    gt_table <- gt_table %>%
      tab_header(
        title = title,
        subtitle = subtit # Add the subtitle
      )%>%
      tab_options(
        heading.subtitle.font.size = px(subtitle_size)

      )
  }

  gt_table <- gt_table %>%
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

  row_number <- nrow(gt_table$`_data`)

  if(col_everyrow==FALSE){
    gt_table <- gt_table%>%
      tab_style(
        style = cell_fill(color = cell_colour),
        locations = cells_body(rows = seq(1,row_number,2))
      )
  }

  if(col_everyrow==TRUE){
    gt_table <- gt_table%>%
      tab_style(
        style = cell_fill(color = cell_colour),
        locations = cells_body(rows = seq(1,row_number))
      )
  }


  # Print the gt table
  #print(gt_table)

  # Return the original summary statistics data frame
  return(gt_table)
}


