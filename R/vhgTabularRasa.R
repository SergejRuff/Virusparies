#' @title vhgTabularRasa: Function for Generating Custom Formatted Graphical Tables
#'
#' @description This function creates a formatted table using the `gt` package.
#'
#' @param file A data frame containing summary statistics.
#' @param title (optional): a custom title.
#' Default is "Graphical Table".
#' @param title_align (optional): a character vector specifying the alignment of title (and subtile) text.
#' Possible values are "left", "center", or "right". Default is "left".
#' @param names_ (optional): a vector of length 3 containing column names.
#' Default is names(file)
#' @param align (optional): a character vector specifying the alignment of text in the table columns.
#' Possible values are "left", "center", or "right". Default is "left".
#' @param subtit (optional): a character vector specifying the subtitle. Default is NULL.
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
#' It is particularly useful for generating tables that cannot be produced with `vhEvalBoxTable` or `vhRunsTable`,
#' when the input data does not originate from the `vhRunsBarplot` or `VhgBoxplot` functions.
#'
#' The `vhgTabularRasa` function allows users to generate tables styled like Virusparies tables using their own input data.
#' Additionally, users can create custom tables by adjusting the parameters within this function.
#' @examples
#' path <- system.file("extdata", "virushunter.tsv", package = "Virusparies")
#' vh_file <- importVirusTable(path)
#'
#'# plot 1: plot boxplot for "identity"
#' #identity <- VhgBoxplot(vh_file,eval_vs_iden="identity",cut = 1e-5)
#'
#' # generate table
#' #iden_table <- vhEvalBoxTable(identity$summary_stats)
#'
#' #iden_table
#'
#' @seealso \code{\link{VhgBoxplot}}
#' @import gt
#' @export
vhgTabularRasa <- function(file,title="Graphical Table",
                               title_align = "left",names_=NULL,align = "left",subtit =NULL,
                               data_row.pad=6,column_colour="dodgerblue4",title_size = 26,subtitle_size=14,
                               title_weight="bold",title_colour = "dodgerblue4",table_font_size = 14,
                               cell_colour="grey90",col_everyrow=FALSE){


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
  print(gt_table)

  # Return the original summary statistics data frame
  return(gt_table)
}


