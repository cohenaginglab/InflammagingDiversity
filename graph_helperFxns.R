#########################
## Helper functions for generating the graphics
##
############################

## Function for mixing two colors (biplot)
mix_colors <- function(red_intensity, blue_intensity) {
  red_intensity <- max(min(red_intensity, 1), 0)
  blue_intensity <- max(min(blue_intensity, 1), 0)
  rgb(red = red_intensity, green = 0, blue = blue_intensity, maxColorValue = 1)
}

## function for looking up color where x is the biomarker
## and cols is the data.frame with color and biomarker mapping
## will return gray if the biomarker doesn't exist in the mapping
colorLookup <- function(x, cols){
  bm <- x  
  l <- which(cols$bm == bm)
  if(length(l)>0){ 
    return(cols$color[l])
  }else{ 
    return(cGray)
  }
}

# Determine significance level
sign_level <- function(p) {
  if (p < 0.001) return(3)
  else if (p < 0.01) return(2)
  else if (p < 0.05) return(1)
  else return(0)
}


# Function to obtain color according to comorbidity
outcome.col <- function(outcome) {
  return(outcome_colors[outcome])
}

# Function to determine line type according to cohort
seg.type <- function(measure) {
  return(line_types[measure])
}

# Function for determining the significance symbol
sign.symbol <- function(x) {
  symbols <- c("", "*", "**", "***")
  return(symbols[x + 1])
}


# Function to draw the graph frame with an arrow for extended UCIs
plot_frame <- function(data_subset, title, add_legend = FALSE, x_axis_title = "Odds Ratio", 
                       lw=2, figLetter, c_names, intraSpace=1, interSpace=1, start_y=1, isTsim=F) {
  num_rows <- nrow(data_subset) + length(unique(data_subset$Outcome)) - 1
  
  plot(1, type = "n", xlim = common_xlim, ylim = c(1, num_rows),
       xlab = "", ylab = "", xaxt = "n", yaxt = "n", main = title, 
       cex.main = 2, cex.lab = 1.7) 
  
  abline(v = 1, col = "black", lwd = 2, lty = 2)
  
  current_y <- start_y
  
  for (outcome in ordered_outcomes) {
    if (!outcome %in% unique(data_subset$Outcome)) next
    outcome_data <- data_subset[data_subset$Outcome == outcome,]
    for (measure in selected_measures) {
      if (!measure %in% unique(outcome_data$Measure)){
        #current_y <- current_y + intraSpace
        next
      }
      row <- outcome_data[outcome_data$Measure == measure,]
      
      # Draw the segment to the UCI value or to 6 if an arrow is required
      segments(x0 = row$LCI, y0 = current_y, x1 = min(row$UCI, 6.6), y1 = current_y,
               col = outcome.col(row$Outcome), lty = seg.type(row$Measure), lwd = lw)
      
      points(row$Effect, current_y, pch = point_types[row$Measure], col = outcome.col(row$Outcome), cex = 1, lwd=1.8)
      
      # Display effect and meaning at the end of the segment or next to the arrow
      effect_text <- paste(format(row$Effect, digits = 2, nsmall = 2), sign.symbol(row$Sign_level))
      text_x_pos <- min(row$UCI, 6.6) + 0.1
      text(x = text_x_pos, y = current_y, labels = effect_text, pos = 4, cex = 1.4, offset = 0.5)
      
      current_y <- current_y + intraSpace
    }
    
  }
  
  text(x = common_xlim[1]-0.1, y = num_rows, labels = figLetter, pos = 4, col = "black", font = 2,cex = 1.75)
  title(xlab=x_axis_title, mgp=c(2.3,1,0), cex.lab=1.6)
  axis(side = 1, at = seq(0, 6, by = 1), las = 1, cex.axis=1.2)
  
  if (add_legend) {
    reversed_outcomes <- rev(ordered_outcomes)
    cohort_names <- c_names
    line_types <- c(1,2) # solide, pointillé, pointillé, pointillé
    legend_names <- c(reversed_outcomes, cohort_names)
    legend_colors <- c(outcome_colors[reversed_outcomes], rep("black", length(cohort_names)))
    legend_line_types <- c(rep(1, length(reversed_outcomes)), line_types)
    legend_point_types <- c(rep(NA, length(reversed_outcomes)), point_types[2], point_types[1])
    
    legend(x = par("usr")[2] + 0, y = par("usr")[4], 
           legend = legend_names,
           col = legend_colors, 
           lty = legend_line_types, seg.len=3,
           pch = legend_point_types,
           cex = 1.4, 
           lwd = c(rep(lw, length(reversed_outcomes)),1.5,1.5),
           xjust = 1, yjust = 1,
           x.intersp = 0.2, 
           y.intersp = 1) 
  }
}
