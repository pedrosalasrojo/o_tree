# Authors: Domenico Moramarco, Paolo Brunori, Pedro Salas-Rojo
# Date: June 2024 
# Purpose: Opportunity Tree function + Scores function

source("C:/Users/user/Desktop/POWER_CODE/Function_OTree.R")

data <- read.csv("C:/Users/user/Desktop/POWER_CODE/data.csv")

# Define baseline set of types, only with sex and fathers education
data <- data %>%
  mutate(types_np = as.numeric(factor(interaction(fedu, sex))))

# Define complete model to be tested. It includes ethnicity and mothers education
circum <- c("sex", "eth", "fedu", "medu")
model <- income ~ sex + eth + fedu + medu

# Estimate O-tree
otree <- o_tree(data, model,
                types = "types_np", outcome = "income", 
                circum = circum,
                alp = 0.01, mu = 0, minbucket = 50,
                weights = "weights", merge = "ctree", try.new = TRUE)

# Estimate IOp with original type partition and "new" type partitions
otree <- otree %>%
  ungroup() %>%
  group_by(types_np) %>%
  mutate(y_tilde_1 = weighted.mean(income, weights)) %>%
  ungroup() %>%
  group_by(m.types) %>%
  mutate(y_tilde_2 = weighted.mean(income, weights)) %>%
  ungroup()%>%
  group_by(o.types) %>%
  mutate(y_tilde_3 = weighted.mean(income, weights)) %>%
  ungroup()

# Higher IOp in the type-partition of the O-tree
gini1 <- gini.wtd(otree$y_tilde_1, otree$weights)
gini2 <- gini.wtd(otree$y_tilde_2, otree$weights)
gini3 <- gini.wtd(otree$y_tilde_3, otree$weights)

# Fewer types in the type-partition of the O-tree
types1 <- length(unique(otree$types_np))
types2 <- length(unique(otree$m.types))
types3 <- length(unique(otree$o.types))

# Higher reward score in the O-tree
rew1 <- reward_score(otree, "types_np", "income")$score
rew2 <- reward_score(otree, "m.types", "income")$score
rew3 <- reward_score(otree, "o.types", "income")$score

# Higher compensation score in the O-tree
com1 <- compen_score(otree, types = "types_np", outcome = "income", model = model, circum = circum, eps = 1200)$score
com2 <- compen_score(otree, types = "m.types", outcome = "income", model = model, circum = circum, eps = 1200)$score
com3 <- compen_score(otree, types = "o.types", outcome = "income", model = model, circum = circum, eps = 1200)$score
