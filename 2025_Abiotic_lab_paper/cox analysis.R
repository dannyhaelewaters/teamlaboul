# Load necessary libraries
library(ggplot2)
library(survival)
library(broom)
library(gt)

# Load the dataset
data <- read.csv("average_thalli_per_group_per_day.csv")

# Convert the Date column to Date type
data$Date <- as.Date(data$Date)

# Create the ggplot
ggplot(data, aes(x = Date, y = Average_Thalli, color = Group)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    title = "Average Number of Thalli Over Time",
    x = "Date",
    y = "Average Number of Thalli",
    color = "Group"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

#Infection curves (survival curve)
inftime_coxanalysis_olda <- read.table("inftime_coxanalysis_olda.txt", header = TRUE, sep = "\t")
infectioncurve <- survfit(Surv(timestop, inf) ~ Group, data = inftime_coxanalysis_olda)
ggsurvplot(infectioncurve, risk.table = TRUE, risk.table.height = 0.3, fun = "event",
           break.x.by = 4, censor.shape="+", censor.size = 6, risk.table.title = "Number of uninfected/live ladybirds", ylab = "Parasite prevalence",
           xlim = c(16,61), legend.title = "",
           legend.labs = c("Control", "Low humidity", "High humidity", "Low temperature", "High temperature"), 
           palette = "Set2")

#Cox proportional hazard test (to get numbers)
res.cox.full <- coxph(Surv(timestart, timestop, inf) ~ Group + sex, data = inftime_coxanalysis_olda)
summary(res.cox.full)
res.cox.melanic <- coxph(Surv(timestart, timestop, inf) ~ Group, data = inftimemelvsnonmel)
summary(res.cox.melanic) #Melanic vs. non-melanic

# Fit the Cox proportional hazards model
cox_model <- coxph(Surv(timestart, timestop, inf) ~ Group + sex, data = inftime_coxanalysis_olda)

# Extract the results with confidence intervals
cox_results <- tidy(cox_model, exponentiate = TRUE, conf.int = TRUE)

# Select and rename relevant columns
cox_results <- cox_results %>%
  select(term, estimate, conf.low, conf.high, p.value) %>%
  rename(
    Hazard_Ratio = estimate,
    Lower_95_CI = conf.low,
    Upper_95_CI = conf.high,
    P_value = p.value
  )

# Create a formatted table
cox_results_table <- cox_results %>%
  gt() %>%
  tab_header(
    title = "Cox Proportional Hazards Model Results",
    subtitle = "Hazard Ratios with 95% Confidence Intervals"
  ) %>%
  fmt_number(
    columns = c(Hazard_Ratio, Lower_95_CI, Upper_95_CI, P_value),
    decimals = 3
  ) %>%
  cols_label(
    term = "Variable",
    Hazard_Ratio = "Hazard Ratio",
    Lower_95_CI = "Lower 95% CI",
    Upper_95_CI = "Upper 95% CI",
    P_value = "P-value"
  )

# Print the table
print(cox_results_table)
