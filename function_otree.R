
# Authors: Domenico Moramarco, Paolo Brunori, Pedro Salas-Rojo
# Date: June 2024 
# Purpose: Opportunity Tree function + Scores function

# The function o_tree is used to fit an opportunity tree as in Moramarco et al. (2024).
# The function takes as input a data frame, a model formula, the name of the type partition variable
# (e.g. "types"), the name of the outcome variable (e.g. "income"), the name of the circumstance variables
# (e.g. "C("sex", "educ_father", ...), the name of the weights variable (e.g. "weights"), the significance level 
# (default = 0.01), the value of the difference in means to be tested (default = 0), the value of minbucket (default = 100),
# and the method used to merge (ctree or t-test). An option to try new splits (default = TRUE) is also included. 

# The function returns the original dataset with two new columns:
# "m.types" showing the types after merging
# "o.types" showing the types after merging and trying a further split.

# The Reward and Compensation Scores deliver a "results" dataframe with the main
# Results (types, pvalues, power, and so on) to compute the score

# Check libraries, install if necessary, and open.
rm(list=ls())

packages <- c("dplyr", "pwr", "data.table", "partykit", "tm", "tidyverse", "dineq", "strucchange") 

for (package in packages) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package)
  }
  library(package, character.only = TRUE)
}

# Functions to calculate types with o_tree, reward score and compensation score ----

o_tree <- function(data, model, types, outcome, circum, weights,
                   alp = 0.05, mu = 0, minbucket = 100,
                   merge = c("ctree", "t.test"), try.new = TRUE) {
  
  set.seed(1)
  
  # Set basic information
  minc=1-alp
  data$types <- data[[types]]
  data$outcome <- data[[outcome]]
  data$weights <- data[[weights]]
  if (all(is.na(data$weights))) {
    data$weights <- rep(1, length(data$weights))
  }
  data[circum] <- lapply(data[circum], as.factor)
  
  # Show original types  
  print(paste0("Original types: "))
  print(paste0(sort(unique(data$types))))
  
  # Check if there are missing values in types, outcome, circumstances or weights.
  if (any(is.na(data$types) | is.na(data$outcome) | is.na(data$weights))) {
    stop("Function stopped: there are missing values in the data.")
  }
  
  # Check if some types can be merged
  # For terminal nodes in "types" run a double loop to check
  # pairwise combinations. Generate vals_b so each pair
  # is only tested once. For each pair, launch a t-test and get the p-value.
  # Store results. Assign "TRUE" if p-value is higher
  # than alp, so the null hypothesis of equal means is not rejected.
  
  data$m.types <- data$types
  results_merge <- as.data.frame(0)
  
  while (nrow(results_merge)>0){
    
    results_merge <- NA
    vals <- sort(unique(data$m.types))
    
    for(i in vals){
      vals_b <- vals[vals>i]
      for(j in vals_b){
        
        if(merge == "t.test"){
          
          t_test <- t.test(data$outcome[data$m.types==i],
                           data$outcome[data$m.types==j],
                           mu = mu,
                           alternative = c("two.sided"),
                           paired = FALSE,
                           var.equal = FALSE)
          
          pval <- round(t_test[["p.value"]], 10)
          
        } else if(merge == "ctree"){
          
          ctcheck <- data %>%
            dplyr::filter(m.types == i | m.types == j) %>%
            dplyr::select(outcome, m.types)
          
          ctree <- partykit::ctree(outcome ~ m.types,
                                   data = ctcheck, 
                                   control = ctree_control(mincriterion = minc, 
                                                           minbucket = minbucket,
                                                           maxdepth = 1))
          if(length(ctree)==1){
            pval <- round(sctest(ctree)[2], 10)
          } else if (length(ctree) > 1){
            pval <- round(sctest(ctree)[["1"]][2], 10)
          }
          
        }
        
        check <- ifelse(pval > alp, "TRUE", "FALSE")
        res <- cbind(i, j, pval, check)
        results_merge <- as.data.frame(na.omit(rbind(results_merge, res)))
        
      } 
    }
    
    colnames(results_merge) <- c("Node 1", "Node 2", "Pval", "Check")
    results_merge <- results_merge %>%
      filter(Check == "TRUE") %>%
      arrange(Pval)
    
    results_merge <- results_merge[rev(order(results_merge$Pval)),]
    
    data$m.types=data$m.types
    
    t1<-as.numeric(results_merge$`Node 1`[1])
    t2<-as.numeric(results_merge$`Node 2`[1])
    
    data$m.types[data$m.types== t1 | data$m.types==t2]<- (max(data$m.types) + 1 )
    
  }
  
  # All merges have been done.  
  print(paste0("New types after merging: "))
  print(paste0(sort(unique(data$m.types))))
  
  if(try.new == TRUE){
    # Fit and store subtrees.
    newnodes<-sort(unique(data$m.types))
    subtrees<-NA
    data$o.types <- data$m.types
    numt<-1
    
    for (c in newnodes) {
      
      attempt <- data[data$m.types==c,]
      
      attempt[circum] <- lapply(attempt[circum], as.factor)
      
      sub.tree <- partykit::ctree(model,
                                  data = attempt, 
                                  control = ctree_control(mincriterion = minc, 
                                                          minbucket = minbucket))
      
      if (length(unique(predict(sub.tree, type="node")))>1){
        
        attempt$o.types <-predict(sub.tree, type="node", 
                                  newdata = data[data$m.types==c,])*100 + numt*1000
        
        data <- data %>%
          dplyr::filter(m.types != c)
        data <- rbind(data, attempt)
        numt<-numt + 1
      }
    }
    
    # All merges have been done.  
    print(paste0("New types after splitting: "))
    print(paste0(sort(unique(data$o.types))))
  }
  
  data <- data %>%
    dplyr::select(-c(outcome))
  
  return(data = data)
  
}

reward_score <- function(data, types, outcome,
                         alp = 0.01, mu = 0, eps = 200){
  
  data$types <- data[[types]]
  data$outcome <- data[[outcome]]
  
  # Check if there are missing values in types or outcome
  if (any(is.na(data$types) | is.na(data$outcome))) {
    stop("Function stopped: there are missing values in the data.")
  }
  
  results <- NA
  vals <- sort(unique(data$types))
  
  for(i in vals){
    
    vals_b <- vals[vals>i]
    
    for(j in vals_b){
      
      t_test <- t.test(data$outcome[data$types==i],
                       data$outcome[data$types==j],
                       mu = mu,
                       alternative = c("two.sided"),
                       paired = FALSE,
                       var.equal = FALSE)
      
      pval <- round(t_test[["p.value"]], 10)
      
      cohen_d <- eps/sqrt((var(data$outcome[data$types==i])/length(data$outcome[data$types==i])) + 
                            (var(data$outcome[data$types==j])/length(data$outcome[data$types==j])))
      
      f_a <- pt(q=qt(p=1-alp/2, 
                     df = length(data$outcome[data$types==j])+length(data$outcome[data$types==i])-2), 
                df = length(data$outcome[data$types==j])+length(data$outcome[data$types==i])-2, 
                ncp=cohen_d)
      
      f_b <- pt(q=qt(p=alp/2, 
                     df = length(data$outcome[data$types==j])+length(data$outcome[data$types==i])-2), 
                df = length(data$outcome[data$types==j])+length(data$outcome[data$types==i])-2, 
                ncp=cohen_d)
      
      pow_test <- round(1 - f_a - f_b, 4)
      
      pop <- (sum(data$types == i) + sum(data$types == j))
      res <- cbind(i, j, pval, (1-pval), pow_test, pop)
      results <- as.data.frame(na.omit(rbind(results, res)))
      
    } 
  }
  
  colnames(results) <- c("Node 1", "Node 2", "Pval", "minPval", "Power", "Pop")
  results$share <- results$Pop/sum(results$Pop)
  results$reward <- results$share * results$minPval
  reward_score <- sum(results$reward)
  print(paste0("Reward Score for this type partition: ", round(mean(reward_score),6)))
  
  return(list(`results` = results, `score` = reward_score))
  
}

compen_score <- function(data, types, outcome, model, circum, 
                         alp = 0.01, mu = 0, eps = 200){
  
  data$types <- data[[types]]
  data$outcome <- data[[outcome]]
  
  # Check if there are missing values in types, outcome, circumstances or weights.
  if (any(is.na(data$types) | is.na(data$outcome))) {
    stop("Function stopped: there are missing values in the data.")
  }
  
  data[circum] <- lapply(data[circum], as.factor)
  
  results <- NA
  vals <- sort(unique(data$types))
  
  for(i in vals){
    
    nextsp <- data %>%
      dplyr::filter(types == i) 
    
    ctree <- partykit::ctree(model,
                             data = nextsp, 
                             control = ctree_control(mincriterion = 0, 
                                                     minbucket = 1,
                                                     maxdepth = 1))
    
    nextsp$types = predict(ctree, type="node", newdata = nextsp)
    
    if(length(ctree)==1){
      
      pval <- NA
      pow_test <- NA
      
    } else if (length(ctree) > 1){
      
      pval <- round(min(sctest(ctree)[["1"]][2,]), 10)

      cohen_d <- eps/sqrt((var(nextsp$outcome[nextsp$types==2])/length(nextsp$outcome[nextsp$types==2])) + 
                            (var(nextsp$outcome[nextsp$types==3])/length(nextsp$outcome[nextsp$types==3])))
      
      f_a <- pt(q=qt(p=1-alp/2, 
                     df = length(nextsp$outcome[nextsp$types==3])+length(nextsp$outcome[nextsp$types==2])-2), 
                df = length(nextsp$outcome[nextsp$types==3])+length(nextsp$outcome[nextsp$types==2])-2, 
                ncp=cohen_d)
      
      f_b <- pt(q=qt(p=alp/2, 
                     df = length(nextsp$outcome[nextsp$types==3])+length(nextsp$outcome[nextsp$types==2])-2), 
                df = length(nextsp$outcome[nextsp$types==3])+length(nextsp$outcome[nextsp$types==2])-2, 
                ncp=cohen_d)
      
      pow_test <- round(1 - f_a - f_b, 4)
      
    }
    
    pop <- nrow(nextsp)
    res <- cbind(i, pval, (1-pval), pow_test, pop)
    results <- as.data.frame(na.omit(rbind(results, res)))
    
  }
  
  colnames(results) <- c("Node", "Pval", "minPval", "Power", "Pop")
  results$share <- results$Pop/sum(results$Pop)
  results$compen <- results$share * results$Power * results$Pval
  compen_score <- sum(results$compen)
  print(paste0("Compensation Score for this type partition: ", round(mean(compen_score),6)))
  
  return(list(`results` = results, `score` = compen_score))
  
}


