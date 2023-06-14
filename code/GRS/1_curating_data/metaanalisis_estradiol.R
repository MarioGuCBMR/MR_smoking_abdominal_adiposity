# METAANALISIS TO COMBINE FEMALE & MALE DATA
# install.packages("meta")
library(meta)
library(metafor)

# Set the random seed for reproducibility
set.seed(123)

# Generate hypothetical effect sizes for 5 male and 5 female studies
male_effect_sizes <- rnorm(5, 0.5, 0.2)
female_effect_sizes <- rnorm(5, 0.8, 0.3)

# Calculate hypothetical standard errors and sample sizes
male_std_errors <- sqrt(1/200 + 1/150 + 1/100) * sd(male_effect_sizes)
female_std_errors <- sqrt(1/300 + 1/200 + 1/100) * sd(female_effect_sizes)
male_sample_sizes <- rep(c(200, 150, 100), each = 5)
female_sample_sizes <- rep(c(300, 200, 100), each = 5)

# Combine the male and female data into a single data frame
data <- data.frame(
  study = 1:10,
  sex = rep(c("male", "female"), each = 5),
  logOR = c(male_effect_sizes, female_effect_sizes),
  se = c(male_std_errors, female_std_errors),
  n = c(male_sample_sizes, female_sample_sizes)
)

# Print the data frame
print(data)

# > head(data)
# study    sex     logOR         se   n
# 1     1   male 0.3879049 0.02387584 200
# 2     2   male 0.4539645 0.04726117 200
# 3     3   male 0.8117417 0.02387584 200
# 4     4   male 0.5141017 0.04726117 200
# 5     5   male 0.5258575 0.02387584 200
# 6     6 female 1.3145195 0.04726117 150

# Specify the meta-analysis model
model <- metareg(logOR, se, data = with(data, data.frame(logOR, se, sex, study)), byvar = "sex", studlab = "study")

# Print the meta-analysis results
summary(model)
