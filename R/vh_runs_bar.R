#' Generate a bar plot showing the number of datasets with unique runs found for each Virus group.
#'
#' @param vh_file Dataframe containing Virushunter hittables results
#' @param generate_table Logical indicating whether to generate a table along with the plot
#' (default: FALSE)
#'
#' @return A list containing the bar plot and optionally the generated table and processed data
#' @import ggplot2
#' @importFrom rlang .data
#' @export
vh_runs_bar <- function(vh_file,generate_table=FALSE){


  # preprocess data for plot
  sample_run <- preprocess_runs_bar(vh_file)

  if(generate_table == TRUE){

    # creat a table
    which_runs_table <- generate_run_bar_table(vh_file)
  }

  #####################
  ### Generate Plot ###
  ### #################

  run_bar <- ggplot(data = sample_run,aes(x=reorder(.data$best_query,.data$unique_SRA_run,FUN=max),y=.data$unique_SRA_run,fill= .data$best_query))+
    geom_bar(stat = "identity")+
    labs(title = "Number of datasets with hits found for each Virus group",
         x="Virus family found in query",
         y="Number of datesets with hits",
         subtitle = paste0("total number of datasets with hits: ",sum(sample_run$unique_SRA_run)))+
    theme_linedraw()+
    theme(legend.position = "bottom",
          axis.text.y = element_text(size = 10),
          axis.text.x = element_text(size = 10),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 10),
          legend.key.size = unit(1.5, "lines"),
          legend.title = element_text(size = 12, face = "bold"))+
    guides(fill=guide_legend(title="virus family"))+
    geom_text(aes(label = .data$res),
              hjust = -0.1,
              color = "black",
              position = position_dodge(0.9),
              size = 3.5)+
    coord_flip()+
    scale_y_continuous(expand = c(0, 0), limits = c(0, max(sample_run$unique_SRA_run)+max(sample_run$unique_SRA_run)/10))+
    theme(
      # This is the new default font in the plot
      text = element_text( size = 8, color = "black"),
      plot.title = element_text(

        size = 16,
        face = "bold",
        color = "#2a475e"),
      # Statistical annotations below the main title
      plot.subtitle = element_text(

        size = 12,
        face = "bold",
        color="#1b2838"
      ))

  plot(run_bar)

  if(generate_table == TRUE){

    return(list(run_bar=run_bar,which_runs_table=which_runs_table,
                sample_run=sample_run))
  }

  if(generate_table == FALSE){

    return(list(run_bar=run_bar,
                sample_run=sample_run))
  }


}





