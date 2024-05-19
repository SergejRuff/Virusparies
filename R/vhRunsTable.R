#' @title Generate a gt summary table of unique runs for each virus group
#'
#' @description This function generates a summary table of unique runs for each virus group
#' based on the input dataset.
#'
#' @param vh_file A data frame containing the Virushunter Hittables Results
#' @param cut The significance cutoff value for E-values (default: 1e-5). Removes rows in vh_file
#' with values larger than cutoff value in ViralRefSeq_E column.
#' @param title (optional): a title for the summary title.
#' Default is "Summary Table of Unique Runs for Each Virus Group"
#' @param title_align (optional): a character vector specifying the alignment of title (and subttile) text.
#' Possible values are "left", "center", or "right". Default is "left".
#' @param names_ (optional): a vector of length 3 containing column names.
#' Default is c("Virus Group","Number of Unique SRA Runs","SRAs Found")
#' @param align (optional): a character vector specifying the alignment of text in the table columns.
#' Possible values are "left", "center", or "right". Default is "left".
#' @param subtit (optional): a character vector specifying the subtitle. Default is NULL.
#' @param data_row.pad (optional): numeric value specifying the row padding. Default is 6.
#' @param column_colour (optional): character specifying the background colour for the column header.
#' @param title_size (optional): numeric value specifying title size. Defaulit is 26.
#' @param subtitle_size (optional): numeric value specifying subtitle size. Defaulit is 14.
#' @param title_weight (optional): character specifying title weight. Default value is "bold".
#' @param title_colour (optional): character specifying title colour. Default is NULL.
#'
#' @return A formatted gt table summarizing unique runs for each virus group
#'
#' @details This function is an internal utility function used within the package.
#' It calculates the number of unique runs for each virus group from the input dataset.
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
vhRunsTable <- function(vh_file,cut = 1e-5,title="Summary Table of Unique Runs for Each Virus Group",
                        title_align = "left",names_=NULL,align = "left",subtit =NULL,
                        data_row.pad=6,column_colour=NULL,title_size = 26,subtitle_size=14,
                        title_weight="bold",title_colour=NULL){

  if(is.null(names_)){
    names_ <- c("Virus Group","Number of Unique SRA Runs","SRAs Found")
  }

  vh_file <- vh_file[vh_file$ViralRefSeq_E < cut,]

  which_runs <- vh_file %>%
    group_by(.data$best_query) %>%
    summarise(unique_SRA_run = n_distinct(across(any_of(c("SRA_run", "run_id")))),
              SRAs_found = paste(unique(across(any_of(c("SRA_run", "run_id")))), collapse = ", "))

  # Remove "c()", "\", and quotation marks
  which_runs$SRAs_found <- gsub("[c()\\\\\"]", "", which_runs$SRAs_found)

  # creat a table
  which_runs_table <- which_runs %>%
    gt()%>%
    tab_header(
      title = title
    )%>%
    cols_label(
      best_query = names_[[1]],
      unique_SRA_run = names_[[2]],
      SRAs_found = names_[[3]]
    ) %>%
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
    )


  return(which_runs_table)

}
