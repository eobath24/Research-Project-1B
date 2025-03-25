library(ggplot2)
library(dplyr)

install.packages("cowplot")
library(cowplot)

read_length_df <- read.csv("~/shared-team/Ella_RP1/greenbeard/lengths.csv")

summary_df <- ddply(read_length_df, "Strain", summarise, grp.mean=mean(Length))

total.length.plot <- ggplot(read_length_df, aes(x=Length, fill=Strain, color=Strain)) +
geom_histogram(binwidth=100, alpha=0.5, position="dodge") +
geom_vline(data=summary_df, aes(xintercept=grp.mean, color=Strain), linetype="dashed", linewidth =0.6) +
scale_x_continuous(breaks = c(0, 15000, 30000, 45000, 60000, 75000, 90000, 105000, 120000)) +
scale_y_continuous(limit= c(0, 30000) , breaks = c(0, 5000, 10000, 15000, 20000, 25000, 30000)) +
labs(x = "Read length (bp)", y = "Count") +
theme_bw()

total.length.plot

ten_kb.length.plot <- ggplot(read_length_df, aes(x=Length, fill=Strain, color=Strain)) +
geom_histogram(binwidth=100, alpha=0.5, position = "dodge") +
geom_vline(data=summary_df, aes(xintercept=grp.mean, color=Strain), linetype="dashed", linewidth=0.5, show.legend = TRUE) +
scale_x_continuous(limit = c(0,10000), breaks = c(0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000)) +
scale_y_continuous(limit= c(0, 30000) , breaks = c(0, 5000, 10000, 15000, 20000, 25000, 30000)) +
labs(x = "Read length (bp)", y = "Count") +
theme_bw()

ten_kb.length.plot

plot <- plot_grid(total.length.plot, ten_kb.length.plot, ncol = 1)
plot


