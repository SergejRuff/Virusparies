#' @title vhEvalIdenBoxplot: Generate boxplots comparing E-values or identity for each virus group
#'
#' @description
#'  This function generates boxplots comparing either E-values or identity
#' for each virus group from VirusHunter Hittables results.
#'
#' @param vh_file A data frame containing VirusHunter Hittable results.
#' @param cut The significance cutoff value for E-values (default: 1e-5)
#' @param eval_vs_iden Specifies whether to generate boxplots for E-values ("evalue")
#' or identity ("identity")
#' @param theme_choice A character indicating the ggplot2 theme to apply. Options include "minimal",
#'  "classic", "light", "dark", "void", "grey" (or "gray"), "bw", "linedraw", and "test".
#'  Default is "minimal".
#'
#' @return A list containing the generated boxplot, summary statistics, and outliers
#'
#' @details This function generates boxplots comparing either E-values or identity for each virus
#' group from the VirusHunters Hittable.
#' It also calculates summary statistics and  identifies outliers for further analysis.
#' The user can specify whether to generate boxplots for E-values or identity by setting
#' the parameter `eval_vs_iden`.
#'
#' @examples
#' path <- system.file("extdata", "virushunter.tsv", package = "Virusparies")
#' vh_file <- importVirusTable(path)
#'
#' # plot 1: plot boxplot for "evalue"
#' eval <- vhEvalIdenBoxplot(vh_file,eval_vs_iden="evalue",cut = 1e-5)
#'
#'
#' # plot 2: plot boxplot for "identity"
#' identity <- vhEvalIdenBoxplot(vh_file,eval_vs_iden="identity",cut = 1e-5)
#'
#' @import ggplot2
#' @importFrom rlang .data
#' @export
vhEvalIdenBoxplot <- function(vh_file,eval_vs_iden="evalue",cut = 1e-5,theme_choice = "minimal"){

  if(eval_vs_iden=="evalue"){
    # define a cut off fot evalue significance
    cutoff <- -log10(cut)

  }

  # Apply the selected theme
  theme_selected <- select_theme(theme_choice)




  ########################
  ### generate boxplot ###
  ########################

  if(eval_vs_iden=="evalue"){

    message("boxplot plotting RefSeq evalues for each virus group")
    print(paste0("using the following cut off: ", cutoff))



    boxp <- ggplot(vh_file,aes(x=reorder(.data$best_query,-log10(.data$ViralRefSeq_E),FUN=median),
                               y=-log10(.data$ViralRefSeq_E),fill=.data$best_query))+
      geom_boxplot(staplewidth = 0.4)+
      labs(x="virus family query",
           y="-log10 of viral Reference E-values",
           title="Boxplot plotting viral Refrence E-Values for each virus family",
           subtitle = paste0("red line shows viral Refrence E-values under user-defined threshold: ",10^(-cutoff)," (-log10 scale: ",cutoff,")"))+
      geom_hline(aes(yintercept=cutoff), colour="#990000")+
      theme_selected+
      theme(legend.position = "bottom")+
      guides(fill=guide_legend(title="virus family"))+
      coord_flip()+
      theme(
        # This is the new default font in the plot
        text = element_text( size = 8, color = "black"),
        plot.title = element_text(
          size = 16,
          face = "bold",
          color = "#2a475e"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.key.size = unit(1.5, "lines"),
        legend.title = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 10)
      )

  }

  if(eval_vs_iden=="identity"){

    message("boxplot plotting identity for each virus group")

    boxp <- ggplot(vh_file,aes(x=reorder(.data$best_query,.data$ViralRefSeq_ident,FUN=median),
                                          y=.data$ViralRefSeq_ident,fill=.data$best_query))+
      geom_boxplot(staplewidth = 0.4)+
      labs(x="virus family query",
           y="Viral Reference Identity in %",
           title="Boxplot plotting viral Refrence Identity for each virus family")+
      theme_selected+
      theme(legend.position = "bottom")+
      guides(fill=guide_legend(title="virus family"))+
      coord_flip()+
      theme(
        # This is the new default font in the plot
        text = element_text( size = 8, color = "black"),
        plot.title = element_text(
          size = 16,
          face = "bold",
          color = "#2a475e"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.key.size = unit(1.5, "lines"),
        legend.title = element_text(size = 12, face = "bold")
      )


  }


  if(eval_vs_iden=="evalue"){

    print("generating summary stats in dataframe and gt-table outlier extraction for evalue-boxplots")

    summary_stats <- vh_sum_stat_evavlue_boxplot(vh_file,cutoff)
    outlier <- find_outlier_eval_box(vh_file)


  }



  if(eval_vs_iden=="identity"){

    print("generating summary stats in dataframe and gt-table for identity-boxplots")

    summary_stats <- summary_stats_identity(vh_file)


  }




  plot(boxp)

  if(eval_vs_iden=="evalue"){
    return(list(boxp=boxp,summary_stats=summary_stats,outlier=outlier))
    #return(list(boxp=boxp,outlier=outlier))
  }

  if(eval_vs_iden=="identity"){
    return(list(boxp=boxp,summary_stats=summary_stats))
  }


}







