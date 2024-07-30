library("lmerTest")
library("readr")

data_dir <- "/Users/chloehampson/Desktop/projects/abide-analysis/dset/group/habenula/" # Make sure to leave the slash at the end
clusters <- c("1", "2", "3", "4")
roi <- "RSFC"

# Level-1 Predictors
group_var <- "Group"
categorical_vars  <- c("Sex")
numerical_vars <- c("Age")
phen_vars <- c("Phen1", "Phen2", "Phen3", "Phen4", "Phen5")

for (cluster in clusters) {
  data_path <- paste0(data_dir, "cluster-", cluster, "_data.csv")
  data <- read.table(file = data_path, sep = ',', header = TRUE)
  
  for (phen_var in phen_vars) {
    all_columns <- c(roi, categorical_vars, numerical_vars, phen_var, group_var, "Site")
    sub_data <- data[, all_columns]
    sub_data <- na.omit(sub_data)
    
    # Print column names to ensure 'Group' is included
    print(colnames(sub_data))
    
    # Convert categorical variables to factors
    for (var in categorical_vars) {
      sub_data[[var]] <- factor(sub_data[[var]])
    }
    
    # Relevel 'Group' to make ASD the reference category
    sub_data[[group_var]] <- relevel(factor(sub_data[[group_var]]), ref = "ASD")

    # Check the levels of 'Group' to ensure releveling worked
    print(levels(sub_data[[group_var]]))
    
    # Scale continuous predictors
    for (var in numerical_vars) {
      sub_data[[var]] <- scale(sub_data[[var]], center = TRUE, scale = TRUE)
    }
    
    sub_data[[phen_var]] <- scale(sub_data[[phen_var]], center = TRUE, scale = TRUE)
    
    # Initialize equation string
    fixed_effects <- paste(c(numerical_vars), collapse = " + ")
    equation_lme <- paste(roi, "~", fixed_effects, "+", group_var, "*", phen_var)
    
    # Conditionally add random effect for Site, except for Phen4
    if (phen_var != "Phen4" && length(unique(sub_data$Site)) > 1) {
      equation_lme <- paste(equation_lme, "+ (1|Site)")
    }
    
    # Print the equation to check
    message(equation_lme)
    
    # Run model
    if (phen_var == "Phen4") {
      model <- lm(as.formula(equation_lme), data = sub_data)
    } else {
      model <- lmer(as.formula(equation_lme), data = sub_data)
    }
    x <- summary(model)
    print(x)
    
    # Write results of the model to csv file
    out_file <- paste0(data_dir, "cluster-", cluster, "_", phen_var, "_table.csv")
    model_table <- as.data.frame(coef(summary(model)))
    write.csv(model_table, file = out_file, row.names = TRUE)
  }
}

