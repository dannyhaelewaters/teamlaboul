# Load required libraries
library(tidyverse)
library(plotly)
library(lubridate)
library(ggplot2)

# Load the dataset
data <- read.csv("Combine_Hesp_growth.csv", header = FALSE, stringsAsFactors = FALSE)

# Cleaning the dataset
# Remove unnecessary rows and assign column names
data_cleaned <- data[-c(1, 2), ]  # Removing the first two rows
colnames(data_cleaned) <- c('Group', 'Sex', 'Temperature', 'Humidity', 'Comments', as.character(data[1, 6:ncol(data)]))
data_cleaned <- data_cleaned[-1, ]  # Removing the header row introduced into the data

# Check column names and identify date columns
date_columns <- grep("^\\d{1,2}/\\d{1,2}/\\d{2,4}$", colnames(data_cleaned), value = TRUE)

# Convert wide format to long format
data_long <- data_cleaned %>%
  pivot_longer(
    cols = all_of(date_columns),  # Use identified date columns
    names_to = "Date",
    values_to = "Fungal_Growth"
  ) %>%
  mutate(
    Date = as.Date(Date, format = "%d/%m/%y"),  # Convert Date column to Date type
    Fungal_Growth = as.numeric(Fungal_Growth),  # Ensure Fungal_Growth is numeric
    Temperature = as.numeric(Temperature),     # Ensure Temperature is numeric
    Humidity = as.numeric(Humidity)            # Ensure Humidity is numeric
  ) %>%
  drop_na(Date)  # Drop rows where Date conversion failed

# Save cleaned and reshaped data to a new CSV file
write.csv(data_long, "cleaned_fungal_growth_data.csv", row.names = FALSE)

# View the cleaned and reshaped data
head(data_long)

# Scatterplot of Fungal Growth vs. Temperature
ggplot(data_long, aes(x = Temperature, y = Fungal_Growth, color = Humidity)) +
  geom_point() +
  labs(title = "Fungal Growth vs. Temperature and Humidity",
       x = "Temperature (°C)", y = "Fungal Growth") +
  theme_minimal()

# Scatterplot of Fungal Growth vs. Humidity
ggplot(data_long, aes(x = Humidity, y = Fungal_Growth, color = Temperature)) +
  geom_point() +
  labs(title = "Fungal Growth vs. Humidity and Temperature",
       x = "Humidity (%)", y = "Fungal Growth") +
  theme_minimal()

# Fit a quadratic regression model
model <- lm(Fungal_Growth ~ poly(Temperature, 2) + poly(Humidity, 2), data = data_long)

# Summary of the model
summary(model)

# Filter data until 6/6/2023 (the point of complete saturation)
data_filtered <- data_long %>%
  filter(Date <= as.Date("2023-06-06"))

# Fit a quadratic regression model on the filtered data
model_filtered <- lm(Fungal_Growth ~ poly(Temperature, 2) + poly(Humidity, 2), data = data_filtered)

# Summary of the filtered model
summary(model_filtered)

# Create a grid of Temperature and Humidity values
temperature_seq <- seq(min(data_filtered$Temperature, na.rm = TRUE),
                       max(data_filtered$Temperature, na.rm = TRUE), length.out = 100)
humidity_seq <- seq(min(data_filtered$Humidity, na.rm = TRUE),
                    max(data_filtered$Humidity, na.rm = TRUE), length.out = 100)

grid_filtered <- expand.grid(Temperature = temperature_seq, Humidity = humidity_seq)
grid_filtered$Fungal_Growth <- predict(model_filtered, newdata = grid_filtered)

# Heatmap of Fungal Growth for filtered data
ggplot(grid_filtered, aes(x = Temperature, y = Humidity, fill = Fungal_Growth)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(title = "Optimum Curve of Fungal Growth (Until 6/6/2023)",
       x = "Temperature (°C)", y = "Humidity (%)", fill = "Fungal Growth") +
  theme_minimal()

optimum_filtered <- grid_filtered[which.max(grid_filtered$Fungal_Growth), ]
print(optimum_filtered)

# Define the optimum point (update these values as necessary)
optimum_temp <- 23.17172  # Replace with the actual optimum temperature
optimum_humidity <- 89.39394  # Replace with the actual optimum humidity
optimum_growth <- 19.21372

# Create the heatmap
ggplot(grid_filtered, aes(x = Temperature, y = Humidity, fill = Fungal_Growth)) +
  geom_tile() +
  scale_fill_viridis_c() +
  scale_y_continuous(
    limits = c(30, 90),  # Set the range of the y-axis
    breaks = seq(30, 90, by = 10) # Set specific break points 
  ) +
  labs(
    title = "Optimum Curve of Fungal Growth",
    x = "Temperature (°C)", 
    y = "Humidity (%)", 
    fill = "Fungal Growth"
  ) +
  # Add the optimum point
  geom_point(aes(x = optimum_temp, y = optimum_humidity), color = "red", size = 3) +
  # Customize the axis font
  theme_minimal() +
  theme(
    axis.title = element_text(size = 14, face = "bold"),  # Bold and larger font for axis titles
    axis.text = element_text(size = 12, face = "bold")    # Bold and larger font for axis text
  )

# Plot the response surface
ggplot(grid_filtered, aes(x = Temperature, y = Humidity, fill = Fungal_Growth)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(title = "Optimal Fungal Growth Conditions",
       x = "Temperature", y = "Humidity", fill = "Growth")

# Define the optimum point
optimum_temp <- 23.17172
optimum_humidity <- 89.39394
optimum_growth <- 19.21372

# Create the 3D plot
plot_ly(
  data = grid_filtered,
  x = ~Temperature,
  y = ~Humidity,
  z = ~Fungal_Growth,
  type = "scatter3d",
  mode = "markers",
  marker = list(
    size = 5,
    color = ~Fungal_Growth,
    colorscale = "Viridis",
    showscale = TRUE
  )
) %>%
  add_markers(
    x = optimum_temp, 
    y = optimum_humidity, 
    z = optimum_growth,
    marker = list(size = 10, color = "red", symbol = "circle"),
    name = "Optimum Point"
  ) %>%
  layout(
    title = "3D Plot of Fungal Growth with Optimum Point",
    scene = list(
      xaxis = list(title = "Temperature (°C)", titlefont = list(size = 14, color = "black")),
      yaxis = list(title = "Humidity (%)", titlefont = list(size = 14, color = "black")),
      zaxis = list(title = "Fungal Growth", titlefont = list(size = 14, color = "black"))
    )
  )

# Add a week column to the dataset
data_long <- data_long %>%
  mutate(Week = floor_date(Date, unit = "week"))  # Group dates into weeks

# Aggregate fungal growth by week
weekly_growth <- data_long %>%
  group_by(Week) %>%
  summarize(Average_Growth = mean(Fungal_Growth, na.rm = TRUE),
            SD_Growth = sd(Fungal_Growth, na.rm = TRUE),
            .groups = "drop")

# Line plot of weekly fungal growth
ggplot(weekly_growth, aes(x = Week, y = Average_Growth)) +
  geom_line(color = "blue", size = 1) +
  geom_point(color = "red") +
  geom_ribbon(aes(ymin = Average_Growth - SD_Growth, ymax = Average_Growth + SD_Growth),
              fill = "blue", alpha = 0.2) +
  labs(title = "Weekly Fungal Growth Over Time",
       x = "Week",
       y = "Average Fungal Growth ± SD") +
  theme_minimal()

# Add a Week column to the dataset
data_long <- data_long %>%
  mutate(Week = floor_date(Date, unit = "week"))  # Group dates into weeks

# Aggregate fungal growth by Week and Group (Treatment)
weekly_growth_by_treatment <- data_long %>%
  group_by(Week, Group) %>%
  summarize(Average_Growth = mean(Fungal_Growth, na.rm = TRUE),
            SD_Growth = sd(Fungal_Growth, na.rm = TRUE),
            .groups = "drop")

# Line plot of weekly fungal growth by treatment
ggplot(weekly_growth_by_treatment, aes(x = Week, y = Average_Growth, color = Group)) +
  geom_line(size = 1) +
  geom_point() +
  geom_ribbon(aes(ymin = Average_Growth - SD_Growth, ymax = Average_Growth + SD_Growth, fill = Group),
              alpha = 0.2, color = NA) +
  labs(title = "Weekly Fungal Growth Over Time by Treatment",
       x = "Week",
       y = "Average Fungal Growth ± SD",
       color = "Treatment",
       fill = "Treatment") +
  theme_minimal()

# Line plot of daily fungal growth by treatment
ggplot(data_long, aes(x = Date, y = Fungal_Growth, color = Group)) +
  geom_line(size = 1, alpha = 0.8) +
  geom_point(alpha = 0.6) +
  labs(title = "Daily Fungal Growth Over Time by Treatment",
       x = "Date",
       y = "Fungal Growth",
       color = "Treatment") +
  theme_minimal()

ggplot(data_long, aes(x = Date, y = Fungal_Growth, color = Group)) +
  geom_line(size = 1, alpha = 0.8) +
  geom_point(alpha = 0.6) +
  geom_smooth(se = FALSE, linetype = "dashed") +
  labs(title = "Daily Fungal Growth Over Time by Treatment",
       x = "Date",
       y = "Fungal Growth",
       color = "Treatment") +
  theme_minimal()

# Line plot of daily fungal growth by treatment
ggplot(data_long, aes(x = Date, y = Fungal_Growth, color = Group)) +
  geom_line(size = 1, alpha = 0.8) +  # Line plot only
  labs(title = "Daily Fungal Growth Over Time by Treatment",
       x = "Date",
       y = "Fungal Growth",
       color = "Treatment") +
  theme_minimal()

# Smooth growth curves for daily fungal growth by treatment
ggplot(data_long, aes(x = Date, y = Fungal_Growth, color = Group)) +
  geom_smooth(se = TRUE, size = 1.2) +  # Smooth curves without standard error ribbons
  labs(title = "Daily Fungal Growth Over Time by Treatment",
       x = "Date",
       y = "Fungal Growth",
       color = "Treatment") +
  theme_minimal()

# Perform ANOVA
anova_model <- aov(Fungal_Growth ~ Group, data = data_long)

# Summary of the ANOVA results
summary(anova_model)

# Tukey HSD test for pairwise comparisons
tukey_results <- TukeyHSD(anova_model)

# View Tukey results
print(tukey_results)


