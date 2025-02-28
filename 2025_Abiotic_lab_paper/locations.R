# Load required libraries
library(ggplot2)
library(tidyr)
library(dplyr)

# Read the dataset
data <- readxl::read_excel("Infection location.xlsx")

# Split multiple locations into separate rows
data_clean <- data %>%
  separate_rows(First, sep = ",") %>%      # Split locations
  mutate(First = trimws(First))            # Trim whitespace

# Count infections per location, grouped by Group
location_counts <- data_clean %>%
  group_by(First, Group) %>%
  summarise(Count = n(), .groups = 'drop')

# Plot the stacked bar graph
ggplot(location_counts, aes(x = First, y = Count, fill = Group)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Infection Locations by Group",
       x = "Infection Location",
       y = "Number of Infections",
       fill = "Group") +
  theme_minimal()

# Plot the grouped bar graph
ggplot(location_counts, aes(x = First, y = Count, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Infection Locations by Group",
       x = "Infection Location",
       y = "Number of Infections",
       fill = "Group") +
  theme_minimal()

# Create a contingency table: Locations (rows) by Groups (columns)
contingency_table <- table(data_clean$First, data_clean$Group)

# Print the contingency table
print("Contingency Table:")
print(contingency_table)

# Perform the chi-squared test
chi_test <- chisq.test(contingency_table)

# Print the test results
print("Chi-Squared Test Results:")
print(chi_test)

# Visualize results with a mosaic plot
mosaicplot(contingency_table, main = "Mosaic Plot of Infection Locations by Group",
           xlab = "Infection Location", ylab = "Group", color = TRUE)

# Convert the contingency table to a data frame
contingency_df <- as.data.frame(as.table(contingency_table))

# Plot a heatmap using ggplot2
ggplot(contingency_df, aes(x = Var2, y = Var1, fill = Freq)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue") +
  labs(title = "Heatmap of Infection Locations by Group",
       x = "Group", y = "Infection Location", fill = "Count") +
  theme_minimal()

# Print the expected counts from the chi-squared test
print("Expected Counts:")
print(chi_test$expected)

# Check if any expected count is less than 5
if(any(chi_test$expected < 5)) {
  print("Warning: Some expected counts are less than 5. Chi-squared assumptions may be violated.")
} else {
  print("All expected counts are greater than or equal to 5. Assumptions are met.")
}

# Perform Fisher's Exact Test
fisher_test <- fisher.test(contingency_table)

# Print Fisher's Exact Test results
print("Fisher's Exact Test Results:")
print(fisher_test)

# Convert the contingency table to a data frame
contingency_df <- as.data.frame(as.table(contingency_table))

# Heatmap of observed counts
library(ggplot2)
ggplot(contingency_df, aes(x = Var2, y = Var1, fill = Freq)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue") +
  labs(title = "Heatmap of Observed Counts",
       x = "Group", y = "Infection Location", fill = "Count") +
  theme_minimal()

# Mosaic plot for contingency table
mosaicplot(contingency_table, main = "Mosaic Plot of Infection Locations by Group",
           xlab = "Group", ylab = "Infection Location", color = TRUE)

# Combine categories into "head"
data_combined <- data_clean %>%
  mutate(First = ifelse(First %in% c("head", "eye", "antenna", "labrum"), "head", First)) %>%
  mutate(First = ifelse(First %in% c("elytron", "elytron_tip"), "elytron", First))
# Create a new contingency table
contingency_table_combined <- table(data_combined$First, data_combined$Group)

# Print the updated contingency table
print("Updated Contingency Table:")
print(contingency_table_combined)

# Perform Fisher's Exact Test on the updated table
fisher_test_combined <- fisher.test(contingency_table_combined)

# Print Fisher's Exact Test results
print("Fisher's Exact Test Results (After Combining):")
print(fisher_test_combined)

# Visualize the updated results with a mosaic plot
mosaicplot(contingency_table_combined, main = "Mosaic Plot of Infection Locations (Combined)",
           xlab = "Group", ylab = "Infection Location", color = TRUE)

# Heatmap of observed counts
contingency_df_combined <- as.data.frame(as.table(contingency_table_combined))
ggplot(contingency_df_combined, aes(x = Var2, y = Var1, fill = Freq)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue") +
  labs(title = "Locations of first mature thalli on ladybird body",
       x = "Group", y = "Infection Location", fill = "Count") +
  theme_minimal()

# Count the frequencies of infection locations
location_counts <- table(data_combined$First)

# Print observed counts
print("Observed Counts of Locations:")
print(location_counts)

# Perform chi-squared goodness-of-fit test
# Null hypothesis: All locations are equally likely
chi_test_locations <- chisq.test(location_counts)

# Print test results
print("Chi-Squared Goodness-of-Fit Test Results:")
print(chi_test_locations)

# Load libraries
library(ggplot2)

# Count frequencies of infection locations
location_counts <- data_combined %>%
  count(First, name = "Count")

# Bar plot of infection location counts
ggplot(location_counts, aes(x = reorder(First, -Count), y = Count, fill = First)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  labs(title = "Counts of Infection Locations",
       x = "Infection Location",
       y = "Number of Infections") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Count frequencies of infection locations by group
location_group_counts <- data_combined %>%
  group_by(First, Group) %>%
  summarise(Count = n(), .groups = "drop")

# Stacked bar plot of infection locations by group
ggplot(location_group_counts, aes(x = reorder(First, -Count), y = Count, fill = Group)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Counts of Infection Locations by Group",
       x = "Infection Location",
       y = "Number of Infections",
       fill = "Group") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Define surface area proportions for each location
surface_areas <- data.frame(
  First = c("elytron", "pronotum", "head", "leg", "sternum"), # Adjust as needed
  SurfaceArea = c(0.60, 0.09, 0.06, 0.12, 0.13) # Example proportions, sum must be 1
)

# Merge infection data with surface areas
data_with_area <- data_combined %>%
  left_join(surface_areas, by = "First")

# Calculate adjusted counts
adjusted_counts <- data_with_area %>%
  group_by(First) %>%
  summarise(AdjustedCount = n() / unique(SurfaceArea), .groups = "drop")

# Create a contingency table with adjusted counts
adjusted_contingency_table <- table(data_with_area$First, data_with_area$Group) / surface_areas$SurfaceArea

# Perform chi-squared test
adjusted_chi_test <- chisq.test(adjusted_contingency_table)

# Print the results
print("Chi-Squared Test Results with Surface Area Adjustment:")
print(adjusted_chi_test)

print(adjusted_chi_test$expected)

# Bar plot for adjusted counts
ggplot(adjusted_counts, aes(x = reorder(First, -AdjustedCount), y = AdjustedCount)) +
  geom_bar(stat = "identity") +
  labs(title = "Adjusted Infection Counts by Location",
       x = "Infection Location",
       y = "Adjusted Count") +
  theme_minimal()

adjusted_fisher_test <- fisher.test(adjusted_contingency_table)

# Print Fisher's Exact Test results
print("Fisher's Exact Test Results with Surface Area Adjustment:")
print(adjusted_fisher_test)

fisher_test <- fisher.test(adjusted_contingency_table, simulate.p.value = TRUE, B = 10000)
print("Fisher's Exact Test with Simulation Results:")
print(fisher_test)

scaled_contingency_table <- round(adjusted_contingency_table * 1000)
fisher_test <- fisher.test(scaled_contingency_table, simulate.p.value = TRUE, B = 10000)
print(fisher_test)

residuals <- adjusted_chi_test$residuals
print(residuals)


# Install required package if not already installed
if (!requireNamespace("reshape2", quietly = TRUE)) install.packages("reshape2")

# Load required libraries
library(ggplot2)
library(reshape2)

# Convert the contingency table to a data frame
adjusted_contingency_df <- as.data.frame(as.table(adjusted_contingency_table))

# Rename columns for clarity
colnames(adjusted_contingency_df) <- c("Body_Part", "Group", "Expected_Count")

# Create a heatmap
ggplot(adjusted_contingency_df, aes(x = Group, y = Body_Part, fill = Expected_Count)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "blue", name = "Expected Count") +
  labs(title = "Heatmap of Expected Counts (Adjusted for Surface Area)",
       x = "Group", y = "Body Part") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

###############

# Summarize expected counts by body part
total_expected_by_location <- adjusted_contingency_df %>%
  group_by(Body_Part) %>%
  summarise(Total_Expected = mean(Expected_Count, na.rm = TRUE))

# Create a bar plot of total expected counts by location
ggplot(total_expected_by_location, aes(x = reorder(Body_Part, -Total_Expected), y = Total_Expected, fill = Body_Part)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  labs(title = "Total Expected Counts by Location (All Groups Combined)",
       x = "Body Part", y = "Total Expected Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(total_expected_by_location, aes(x = reorder(Body_Part, -Total_Expected), y = Total_Expected, fill = Body_Part)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Counts of Infection Locations by Group",
       x = "Infection Location",
       y = "Number of Infections",
       fill = "Group") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(total_expected_by_location)
print(adjusted_counts)

--
# Load necessary library
library(readxl)

# Read the data
file_path <- "fractions locations.xlsx"
data <- read_excel(file_path)

# Extract actual and expected counts
actual_counts <- data$actual_counts
expected_counts <- data$expected

# Perform the chi-square test
chi_test <- chisq.test(x = actual_counts, p = expected_counts / sum(expected_counts))

# Print results
print(chi_test)

#fisher's
fisher_test <- fisher.test(rbind(actual_counts, expected_counts))
print(fisher_test)

#monte carlo
chi_test_mc <- chisq.test(x = actual_counts, p = expected_counts / sum(expected_counts), simulate.p.value = TRUE, B = 10000)
print(chi_test_mc)

# Visualization
barplot(rbind(actual_counts, expected_counts),
        beside = TRUE,
        col = c("red", "blue"),
        legend = c("Actual", "Expected"),
        names.arg = c("Elytron", "Pronotum", "Head", "Leg", "Sternum"),
        main = "Observed vs. Expected Fungal Infections",
        ylab = "Counts")

# Load necessary library
library(readxl)

# Read the data
file_path <- "fractions locations.xlsx"
data <- read_excel(file_path)

# Check column names
print(colnames(data))  # See if there's an extra column

# Remove the first column if it's an index
if (colnames(data)[1] == "...1") {
  data <- data[, -1]  # Drop the first column
}

# Extract numeric counts
actual_counts <- as.numeric(data$actual_counts)
expected_counts <- as.numeric(data$expected)

# Ensure the matrix has the correct dimensions (2 x 5)
counts_matrix <- rbind(actual_counts, expected_counts)

# Define correct body part names
body_parts <- c("Elytron", "Pronotum", "Head", "Leg", "Sternum")

# Bar plot
barplot(counts_matrix, 
        beside = TRUE, 
        col = c("red", "blue"), 
        legend.text = c("Actual", "Expected"), 
        names.arg = body_parts, 
        main = "Observed vs. Expected Fungal Infections", 
        ylab = "Counts")

# Perform chi-square test with Monte Carlo simulation
chi_test_mc <- chisq.test(x = actual_counts, p = expected_counts / sum(expected_counts), simulate.p.value = TRUE, B = 10000)

# Extract standardized residuals
residuals <- chi_test_mc$stdres

# Create a data frame
residuals_df <- data.frame(
  BodyPart = c("Elytron", "Pronotum", "Head", "Leg", "Sternum"),
  StandardizedResiduals = residuals
)

# Print residuals
print(residuals_df)

# Highlight significant differences
residuals_df$Significance <- ifelse(abs(residuals_df$StandardizedResiduals) > 2, "Significant", "Not Significant")

# Print table with significance labels
print(residuals_df)

# Initialize empty vector for p-values
p_values <- numeric(length(actual_counts))

# Perform binomial test for each body part
for (i in 1:length(actual_counts)) {
  p_values[i] <- binom.test(actual_counts[i], sum(actual_counts), expected_counts[i] / sum(expected_counts))$p.value
}

# Adjust for multiple comparisons using Bonferroni correction
p_adj <- p.adjust(p_values, method = "bonferroni")

# Create a results table
binom_results <- data.frame(
  BodyPart = c("Elytron", "Pronotum", "Head", "Leg", "Sternum"),
  P_Value = p_values,
  Adjusted_P_Value = p_adj
)

# Print results
print(binom_results)

# Determine if each body part has more or fewer infections than expected
binom_results$Direction <- ifelse(actual_counts > expected_counts, "More Infected", "Less Infected")

# Print the updated results table
print(binom_results)

# Load necessary library
library(ggplot2)

# Create a data frame for plotting
plot_data <- data.frame(
  BodyPart = rep(c("Elytron", "Pronotum", "Head", "Leg", "Sternum"), each = 2),
  CountType = rep(c("Actual", "Expected"), times = 5),
  Count = c(actual_counts, expected_counts),
  Significance = rep(ifelse(binom_results$Adjusted_P_Value < 0.05, "Significant", "Not Significant"), each = 2)
)

# Bar plot using ggplot2
ggplot(plot_data, aes(x = BodyPart, y = Count, fill = interaction(CountType, Significance))) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Actual.Significant" = "red", "Actual.Not Significant" = "gray",
                               "Expected.Significant" = "blue", "Expected.Not Significant" = "lightblue"),
                    labels = c("Actual (Significant)", "Actual (Not Significant)", 
                               "Expected (Significant)", "Expected (Not Significant)")) +
  labs(title = "Observed vs. Expected Fungal Infections",
       x = "Body Part", y = "Counts", fill = "Legend") +
  theme_minimal()


# Load necessary libraries
library(dplyr)
library(knitr)
library(gt)

# Format the table using kable (for Markdown/HTML)
kable(results_table, digits = 4, col.names = c("Body Part", "P-Value", "Adjusted P-Value", "Direction"))

# Format the table using gt (for a styled output)
gt_table <- results_table %>%
  gt() %>%
  fmt_number(columns = c("P_Value", "Adjusted_P_Value"), decimals = 4) %>%
  tab_header(
    title = "Post-Hoc Binomial Test Results",
    subtitle = "Comparison of Observed vs. Expected Infections"
  ) %>%
  data_color(
    columns = "Adjusted_P_Value",
    colors = scales::col_numeric(
      palette = c("red", "orange", "green"),
      domain = c(0, 1)
    )
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = vars(Direction)
    )
  ) %>%
  cols_label(
    BodyPart = "Body Part",
    P_Value = "P-Value",
    Adjusted_P_Value = "Adjusted P-Value",
    Direction = "More/Less Infected"
  )

# Print gt table
print(gt_table)

# Load necessary libraries
library(dplyr)
library(writexl)

# Create results table
results_table <- binom_results %>%
  mutate(Direction = ifelse(actual_counts > expected_counts, "More Infected", "Less Infected")) %>%
  select(BodyPart, P_Value, Adjusted_P_Value, Direction)

# Rename columns for clarity
colnames(results_table) <- c("Body Part", "P-Value", "Adjusted P-Value", "More/Less Infected")

# Save to an Excel file
write_xlsx(results_table, "Binomial_Test_Results.xlsx")

# Print message
print("Results have been saved as 'Binomial_Test_Results.xlsx' in your working directory.")
