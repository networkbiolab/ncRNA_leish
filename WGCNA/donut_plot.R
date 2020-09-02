library(ggplot2)

# Create test data.
data <- data.frame(
  category=c("Upregulated", "Downregulated"),
  count=c(297,190)
)

# Compute percentages
data$fraction <- data$count / sum(data$count)

# Compute the cumulative percentages (top of each rectangle)
data$ymax <- cumsum(data$fraction)

# Compute the bottom of each rectangle
data$ymin <- c(0, head(data$ymax, n=-1))

# Compute label position
data$labelPosition <- (data$ymax + data$ymin) / 2

# Compute a good label
data$label <- paste0(data$category, ":\n", data$count)

pal<-c("#377eb8", "#e41a1c")
# Make the plot
ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
  geom_rect(colour="black") +
  geom_label( x=3.5, aes(y=labelPosition, label=label), size=5.5) +
  scale_fill_manual(values = pal)+
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void() +
  theme(legend.position = "none")