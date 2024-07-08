# Opportunity tree and Reward/Compensation Scores

Set of R-functions to estimate the Opportunity Tree (O-tree) and the Reward and Compensation Scores as defined in Moramarco et al. (2024).

## Files

function_otree.R: this script contains the function producing the O-tree and the reward and compensation scores. Some comments clarify the arguments to be plugged.

to_run.R: this script runs the functions and computes and example exercise of an ad-hoc type partition. The script computes an O-tree and then estimates the associated (gini) Inequality of Opportunity, and the reward and compensation scores.

data.csv: example data used to run the example.

## Description of the example exercise

The example exercise is aimed at showing how to compute an O-tree from a given type partition. This type partition can be ad-hoc or produced by a splitting algorithm, such as Conditional Inference Trees (see Moramarco et al., 2024 for details).

Let's imagine that we estimate IOp on a sample using only two circumstances, sex (sex) and fathers education (fedu). The script first creates the type partition "types_np" as the interaction of both circumstances.

```
# Define baseline set of types, only with sex and fathers education
data <- data %>%
  mutate(types_np = as.numeric(factor(interaction(fedu, sex))))
```

We want to see two things. First, if the compensation and reward scores of this type partition are good enough to measure IOp. Second, if another type partition, produced with the O-tree, produces better (higher) scores. We let the type partition in the O-tree to explore other circumstances, like the ethnicity (eth) and the mothers education (medu).

```
# Define complete model to be tested in the O-tree. It now includes ethnicity and mothers education
circum <- c("sex", "eth", "fedu", "medu")
model <- income ~ sex + eth + fedu + medu
```

First, we estimate the O-tree. To run the function "o_tree" we need to define the data, the model (see above), the variable defining the original type partition (types), the name of the outcome (outcome), the complete vector of circumstance names that we want to check (circum), the significance level in the statistical tests (alp), the difference-in-means to be tested (mu), the minbucket parameter in the C-tree (minbucket), the name of the vector of weights (weights, if you have no weights, substitute by a vector of 1), the method used to merge types with the same statistical mean (p-values are obtained from a c-tree or a difference-in-means t-test; "ctree" or "t.test"), and try.new set as TRUE (default) if, once the types have been merged, you want to produce new types. 

```
# Estimate O-tree
otree <- o_tree(data, model,
                types = "types_np", outcome = "income", 
                circum = circum,
                alp = 0.01, mu = 0, minbucket = 50,
                weights = "weights", merge = "ctree", try.new = TRUE)
```

The o_tree function delivers a dataframe with three type partitions:

The *first type partition* (types_np) delivers a Gini of 0.112, obtained after assigning the average outcome from 24 types. The reward score ascends to 0.805 and the compensation score is essentially zero. 
```
gini1 <- gini.wtd(otree$y_tilde_1, otree$weights)
types1 <- length(unique(otree$types_np))
rew1 <- reward_score(otree, "types_np", "income")$score
com1 <- compen_score(otree, types = "types_np", outcome = "income", model = model, circum = circum, eps = 1200)$score
```
The *second type partition* (m.types) is obtained from the first step of the O-tree, and is the type partition obtained from merging types with the (statistically) same mean. The number of types descends from 24 to 5, but the Gini remains at 0.111 points. It seems clear that the first type partition overfitted the data. The more sparse type partition rises the reward score to 0.999, but maintains the compensation score in zero. It seems there are other circumstances that may be playing a role to estimate IOp.
```
gini2 <- gini.wtd(otree$y_tilde_2, otree$weights)
types2 <- length(unique(otree$m.types))
rew2 <- reward_score(otree, "m.types", "income")$score
com2 <- compen_score(otree, types = "m.types", outcome = "income", model = model, circum = circum, eps = 1200)$score
```
The *third and final type partition* (o.types) results from the O-tree, and is obtained after expanding the second type partition as described in Moramarco et al. (2024). We are left with 11 types, that now rise IOp to 0.139 Gini poitns. The reward score of the new type partition ascends to 0.956, while the compensation score rises to 0.047.
This way, the partition delivered by the O-tree produces higher IOp estimates using fewer types than the original type partition (from 24 to 11), while rising both scores. The resulting type partition is more satisfactory in terms of measurement and normative principles.
```
gini3 <- gini.wtd(otree$y_tilde_3, otree$weights)
types3 <- length(unique(otree$o.types))
rew3 <- reward_score(otree, "o.types", "income")$score
com3 <- compen_score(otree, types = "o.types", outcome = "income", model = model, circum = circum, eps = 1200)$score
```

_Disclaimer_: This is a work in progress version. Do not use or cite without checking.

*References*:
Moramarco, D., Brunori, P. and Salas-Rojo, P. (2024) Biases in inequality of opportunity estimates: measures and solutions (mimeo)

