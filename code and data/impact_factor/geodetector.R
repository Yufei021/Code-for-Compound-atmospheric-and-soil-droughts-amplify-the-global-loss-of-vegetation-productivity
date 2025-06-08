library(GD)
library(openxlsx)
library(readxl)
setwd("G:/paper01/VPD&SM_CDE(M03)/GIMMS_NIRv/impact_factors/cor_re/SM")

# Convert tibble format to data.frame
data <- read_excel("SM_loss.xlsx", range = cell_cols(1:11))
# View(data)
# str(data)

# Convert tibble to data.frame
data <- as.data.frame(data) 

# Before conducting Geographical Detector analysis, continuous variables need to be discretized.
# The optimal discretization method and number of categories will be selected for better modeling.
# Available discretization methods include: "equal", "natural", "quantile", "geometric", and "sd".
# The optidisc() function automatically selects the best discretization method and number of categories.
# For multiple variables, including continuous ones:
discmethod <- c("equal", "natural", "quantile", "geometric", "sd")
discitv <- c(3:7)
dataFin <- data
data.continuous <- dataFin[, c(1:11)]  # Discretize continuous variables in columns 1 to 11

# Discretize the data
odc1 <- optidisc(GPP_loss ~ ., data = data.continuous, discmethod, discitv)
dim(data.continuous)

# Assign the discretized values back to the original data
dataFin[, c(1:10)] <- data.continuous

# Perform single-factor geographical detector analysis
gd <- gd(GPP_loss ~ ., data = dataFin[, c(11, 1:10)])
gd

str(gd)

gd_result <- gd$Factor

colnames(gd_result)


# plot the results--------------------------------------------------------------
# Set font to Times New Roman
windowsFonts(TimesNewRoman = windowsFont("Times New Roman"))

library(ggplot2)
ggplot(gd_result, aes(x = reorder(variable, qv), y = qv)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  geom_text(aes(label = round(qv, 3)), 
            hjust = -0.1, 
            size = 4, 
            family = "TimesNewRoman") +
  labs(
    # title = "Explanatory power of each factor for GPP loss (q value)",
    # x = "Factor",
    y = "Q value") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.23))) +
  theme_minimal(base_size = 14) +
  theme(text = element_text(family = "TimesNewRoman"),
        plot.title = element_text(hjust = 0.5),
        aspect.ratio = 9 / 6,  # 👈 Set aspect ratio (height/width)
        panel.grid.major = element_blank(),       # Remove major grid lines
        panel.grid.minor = element_blank(),       # Remove minor grid lines
        panel.border = element_blank(),           # Remove all borders
        axis.line = element_line(color = "black", linewidth = 0.6),  # Add bottom and left axis lines
        axis.line.y.right = element_blank(),      # Avoid duplicate right axis line
        axis.line.x.top = element_blank(),        # Avoid duplicate top axis line
        
        # set axis tick text color to black
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black")
  )


















