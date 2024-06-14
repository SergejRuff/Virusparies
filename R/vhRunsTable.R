#' @title vhRunsTable: Generate a gt summary table of unique runs for each virus group
#'
#' @description This function generates a summary table of unique runs for each virus group
#' based on the input dataset.
#'
#' @param vh_file A data frame containing the Virushunter Hittables Results
#' @param groupby (optional) A string indicating the column used for grouping the data.
#' "best_query" and "ViralRefSeq_taxonomy" can be used. Default is "best_query".
#' Note: Gatherer hittables do not have a "best_query" column. Please provide an appropriate column for grouping.
#' @param cut The significance cutoff value for E-values (default: 1e-5). Removes rows in vh_file
#' with values larger than cutoff value in ViralRefSeq_E column.
#' @param title (optional): a title for the summary title.
#' Default is "Summary of unique runs by virus group"
#' @param title_align (optional): a character vector specifying the alignment of title (and subtile) text.
#' Possible values are "left", "center", or "right". Default is "left".
#' @param names_ (optional): a vector of length 3 containing column names.
#' Default is c("Virus Group","Number of Unique SRA Runs","SRAs Found")
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
#' @return A formatted gt table summarizing unique runs for each virus group
#'
#' @details
#' vhRunsTable calculates the number of unique runs for each virus group from the input dataset.
#' It takes VirusHunter hittables as input and determines how many runs (SRA runs or local FastQ files) found a specific virus family.
#'
#' A graphical table is returned with two columns. The first column contains the name of the virus group,
#' while the second column contains all IDs of SRA runs or local files that found that virus group.
#'
#' @examples
#' path <- system.file("extdata", "virushunter.tsv", package = "Virusparies")
#' vh_file <- importVirusTable(path)
#'
#' # example 1: generate table with defaul arguments
#' table <- vhRunsTable(vh_file,cut = 1e-5)
#'
#' table
#'
#' # example 2: generate table with custom arguments
#'
#' table_2 <- vhRunsTable(vh_file,title = "test",title_align="right",
#' names_ = c("column_1","column_2","column_3"),align = "right",subtit="subtitle")
#'
#' table_2
#'
#' @seealso \code{\link{vhRunsBarplot}}
#'
#' @importFrom dplyr group_by summarise n_distinct across any_of
#' @importFrom gt gt
#' @importFrom rlang .data
#' @export
vhRunsTable <- function(vh_file,groupby = "best_query",cut = 1e-5,title="Summary of unique runs by virus group",
                        title_align = "left",names_=NULL,align = "left",subtit =NULL,
                        data_row.pad=6,column_colour="dodgerblue4",title_size = 26,subtitle_size=14,
                        title_weight="bold",title_colour = "dodgerblue4",table_font_size = 14,
                        cell_colour="grey90",col_everyrow=FALSE){


  if(groupby == "ViralRefSeq_taxonomy"){

    vh_file <- taxonomy_group_preprocess(vh_file)

  }

  if(is.null(names_)){
    names_ <- c("Group","Number of Unique SRA Runs","SRAs Found")
  }

  vh_file <- vh_file[vh_file$ViralRefSeq_E < cut,]

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

  if(!is.null(subtit)){
    which_runs_table <- which_runs_table %>%
      tab_header(
      title = title,
      subtitle = subtit # Add the subtitle
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
