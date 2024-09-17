library(HydeNet)
library(rjags)
library(MatchIt)
library(ggplot2)
library(plyr)
library(dplyr)
library(gtable)
library(grid)
library(GGally)
library(ggridges)
library("reshape2")
library(gridExtra)


#Mechanical functions (hidden)######
#x is numeric
check_integer = function(x){
  check = all.equal(x, as.integer(x))
  return(isTRUE(check))
}
set_p = function(p,model){
  p2 = log(p/(1-p))
  b0 = p2-model
  return(b0)
}
tot_bind <- function(datalist) {
  require(plyr)
  temp = rbind.fill(datalist)
  rownames(temp) <- unlist(lapply(datalist, row.names))
  return(temp)
}

beta_sum = function(run, a=0.05){
  
  #call needed data
  
  
  ate = run[[9]][2]
  calc_ate = as.data.frame(run[[2]])
  lower_int = as.data.frame(run[[3]])
  upper_int = as.data.frame(run[[4]])
  p_val = as.data.frame(run[[5]])
  n=length(calc_ate[[1]])
  
  #parse data
  names= unlist(lapply(colnames(lower_int), function(y) rep(y, nrow(lower_int))))
  upper_int= as.numeric(unlist(upper_int))
  lower_int = as.numeric(unlist(lower_int))
  calc_ate = as.numeric(unlist(calc_ate))
  p_val = as.numeric(unlist(p_val))
  data = dplyr::tibble(names,lower_int, upper_int, calc_ate, p_val)
  data_ci = data %>% dplyr::filter(lower_int <= ate & upper_int >= ate)
  data_under_ci = data %>% dplyr::filter(upper_int < ate)
  data_over_ci = data %>% dplyr::filter(lower_int > ate)
  data_ci$in_ci = 1
  data_under_ci$in_ci = 0
  data_over_ci$in_ci = 0
  data_ci$over_ci = 0
  data_ci$under_ci = 0
  data_under_ci$over_ci = 0
  data_over_ci$over_ci = 1
  data_under_ci$under_ci = 1
  data_over_ci$under_ci = 0
  
  data_under_a = data %>% dplyr::filter(p_val <= a)
  
  data = rbind(data_ci, data_over_ci, data_under_ci)
  #max((data_ci %>% dplyr::filter(names == "logistic_or"))[4])
  
  #create summary functions
  max_ci_beta = sapply(unique(data$names), function(x) max((data_ci %>% dplyr::filter(names == x))[4]))
  min_ci_beta = sapply(unique(data$names), function(x) min((data_ci %>% dplyr::filter(names == x))[4]))
  max_beta = sapply(unique(data$names), function(x) max((data %>% dplyr::filter(names == x))[4]))
  min_beta = sapply(unique(data$names), function(x) min((data %>% dplyr::filter(names == x))[4]))
  mean_beta = sapply(as.list(unique(data$names)), function(x) mean((data %>% dplyr::filter(names == x))[[4]]))
  names(mean_beta) = colnames(run[[2]])
  sd_beta = sapply(as.list(unique(data$names)), function(x) sd((data %>% dplyr::filter(names == x))[[4]]))
  if(is.null(dim(run[[1]]))==T){
    prop_over_ci = sum(data$over_ci)/length(data$over_ci)
    prop_under_ci = sum(data$under_ci/length(data$under_ci))
    beta_bias_mcse =sapply(as.list(unique(data$names)), function(x) sqrt(sum(((data %>% dplyr::filter(names == x))[[4]]-mean_beta)^2) * 1/(n*(n-1))) )
    
  }else{
    
    prop_over_ci = sapply(unique(data$names), function(x) length((data_over_ci %>% dplyr::filter(names == x))[[4]])/length(run[[1]][,1]))
    prop_under_ci = sapply(unique(data$names), function(x) length((data_under_ci %>% dplyr::filter(names == x))[[4]])/length(run[[1]][,1]))
    beta_bias_mcse =sapply(as.list(unique(data$names)), function(x) sqrt(sum(((data %>% dplyr::filter(names == x))[[4]]-mean_beta[[x]])^2) * 1/(n*(n-1))) )
  }
  
  beta_se =  sapply(as.list(unique(data$names)), function(x) sd((data %>% dplyr::filter(names == x))[[4]])/sqrt(n))
  beta_bias = mean_beta-ate
  
  reject_per = sapply(unique(data$names), function(x) length((data_under_a %>% dplyr::filter(names == x))[[4]])/n)
  coverage = 1-(prop_over_ci+prop_under_ci)
  
  #MCMC std error of functions
  
  coverage_mcse = sqrt((coverage*(1-coverage))/n)
  reject_per_mcse = sqrt((reject_per*(1-reject_per))/n)
  
  #amalgamte into model
  beta_in_ci = as.data.frame(rbind(max_ci_beta, min_ci_beta, max_beta, min_beta, mean_beta, sd_beta, prop_over_ci, prop_under_ci,
                                   beta_se, beta_bias, reject_per, coverage,
                                   beta_bias_mcse, reject_per_mcse, coverage_mcse)
  )
  beta_in_ci[sapply(beta_in_ci, is.infinite)]=NA
  if(is.null(dim(run[[1]]))==T){
    colnames(beta_in_ci) = "regression"
    
  }else{
    beta_in_ci = beta_in_ci[,c(colnames(run[[1]]))] #reorder columns to original run order for indexing
    
    
  }
  return(beta_in_ci)
  
}
#####
#make_model now as flexible argument lengths, DAGs must be functions, will be internal cleaning function

#hidden also
make_model = function(dag, ...){
  arg_list = list(...)
  if(class(dag)== "HydeNetwork"){
    dag_1 = dag
  }else if(length(arg_list)==0){
    dag_1 = dag()
    
  }else{
    dag_1 = do.call(dag, arg_list)
  }
  
  writeNetworkModel(dag_1, pretty = TRUE)
  comp_dag = compileJagsModel(dag_1)
  return(comp_dag)
}


#hidden?
misdiagnosis = function(df, variable, under_rate=0, over_rate=0){
  index_1 = which(df[variable] == 1)
  index_0 = which(df[variable] == 0)

  over = round(length(index_0)*over_rate)
  under = round(length(index_1)*under_rate)

  if(under != 0){
    df[index_1, variable][1:under] = 0
  }
  if(over != 0){
    df[index_0, variable][1:over] = 1
  }

  return(df)
}

#prep for read out######

get_ps = function(exposure, covariates, df){
  if(class(df[,exposure]) == "numeric" ){
    ps_mod = lm(as.formula(paste(exposure, paste(covariates, collapse=" + "), sep=" ~ ")), data = df)
    ps = as.numeric(plogis(fitted.values(ps_mod)))
    num_mod = lm(as.formula(paste(exposure, 1, sep=" ~ ")), data = df)
    num = as.numeric(plogis(fitted.values(num_mod)))
  }else if(class(df[,exposure]) == "integer"){
    ps_mod <- glm(as.formula(paste(exposure, paste(covariates, collapse=" + "), sep=" ~ ")), data = df, family="binomial")
    ps = fitted(ps_mod)
    num_mod = glm(as.formula(paste(exposure, 1, sep=" ~ ")), data = df, family = "binomial")
    num = fitted(num_mod)
  }else{
    print("exposure must be numeric or integer class")
  }
  
  df_out = data.frame(ps = ps,
                      weights = num/ps)
  df$ps = df_out$ps
  df$weights = df_out$weights
  return(df)
}

#####

#Read outs perhaps should all be hidden#####
#bi outcome maybe cut risk_ratio

odds_ratio = function(exposure, outcome, covariates=NULL, df){
  vars = c(exposure, covariates)
  cont_glm = stats::glm(as.formula(paste(outcome, paste(vars, collapse=" + "), sep=" ~ ")), data = df, family = "quasibinomial")
  exp_coef = as.numeric(cont_glm$coefficients[2])
  exp_or = exp(exp_coef)
  or_confint = stats::confint.default(cont_glm, parm = exposure, trace = F)
  
  or_upper = or_confint[1,2]
  or_lower = or_confint[1,1]
  
  int_diff = or_upper - or_lower
  
  or_df = data.frame("odds_ratio" = exp_or,
                     beta = exp_coef,
                     lower_int = or_lower,
                     upper_int = or_upper,
                     confint_diff = abs(int_diff),
                     p_val = coef(summary(cont_glm))[2,4],
                     n = nrow(df))
  
  
  row.names(or_df) = "logistic_regression"
  return(or_df)
}


#cont outcome
lm_beta = function(exposure, outcome, covariates=NULL, df){
  vars = c(exposure, covariates)
  lm1 = stats::lm(as.formula(paste(outcome, paste(vars, collapse=" + "), sep=" ~ ")), data = df)
  confint = confint(lm1, parm = exposure, trace = F)
  upper_int = confint[1,2]
  lower_int = confint[1,1]
  beta = as.numeric(lm1$coefficients[2])
  regression_df = data.frame("odds_ratio" = NA,
                             beta = beta,
                             lower_int =lower_int,
                             upper_int = upper_int,
                             confint_diff = abs(upper_int-lower_int),
                             p_val = coef(summary(lm1))[2,4],
                             n = nrow(df))
  
  rownames(regression_df) = "linear_regression"
  return(regression_df)
}
#any
ps_weight = function(exposure, outcome, covariates, df, weights){
  vars = c(exposure, covariates)
  if(class(df[,outcome]) == "numeric" ){
    cont_glm = stats::lm(as.formula(paste(outcome, paste(vars, collapse=" + "), sep=" ~ ")), data = df, weights = weights)
  }else if(class(df[,outcome]) == "integer"){
    cont_glm = stats::glm(as.formula(paste(outcome, paste(vars, collapse=" + "), sep=" ~ ")), data = df, weights = weights, family = "quasibinomial")
  }else{
    warning("exposure must be numeric or integer")
  }
  
  exp_coef = as.numeric(cont_glm$coefficients[2])
  exp_or = exp(exp_coef)
  or_confint = stats::confint.default(cont_glm, parm = exposure, trace = F)
  or_upper = or_confint[1,2]
  or_lower = or_confint[1,1]
  int_diff = or_upper - or_lower
  if(class(df[,outcome]) == "numeric" ){
    exp_or = NA
  }
  
  
  or_df = data.frame("odds_ratio" = exp_or,
                     beta = exp_coef,
                     lower_int = or_lower,
                     upper_int = or_upper,
                     confint_diff = abs(int_diff),
                     p_val = coef(summary(cont_glm))[2,4],
                     n = nrow(df))
  row.names(or_df) = "propensity_score_weighting"
  return(or_df)
}
#####


#Data Parsing hidden#####
bi_strat = function(x, df, column){
  index_1 = which(df[column] == x)
  index_0 = which(df[column] != x)
  df[index_0, column] = 0
  df[index_1, column] = 1
  return(list(df[index_1,], df[index_0,]))
}
dichotomize = function(column, df, div){
  x = as.numeric(quantile(df[,column])[div])
  index_0 = which(df[column] < x)
  index_1 = which(df[column] >= x)
  df[index_0, column] = 0
  df[index_1, column] = 1
  df[,column] = as.integer(df[,column])
  df = rbind(df[index_1,], df[index_0,])
  rownames(df) = NULL
  return(df)
}

#####

create_data = function(dag, n, positivity = F, ...){
  reclassify = as.integer
  jag_dag = make_model(dag, ...)
  sim_df = bindSim(HydeSim(jag_dag, variable.names = colnames(jag_dag$dag), n.iter = n, bind = FALSE))
  sim_df = sim_df[c(-length(sim_df), -(length(sim_df)-1))]
  
  relabel = lapply(sim_df, check_integer) # JAGS labels integers as numeric, have to reclassify them
  relabel = relabel[relabel != FALSE]
  relabel = names(relabel)
  sim_df[relabel] = lapply(sim_df[relabel], reclassify)
  
  sep_check = length(which(duplicated(t(sim_df))==TRUE))
  if(sep_check != 0){
    stop("complete separation occured (one variable completly predicts another). A beta value may be too large, or a sample size may be too small")
  }
  if(positivity==T){
    binary_cols = names(which(lapply(sim_df, is.integer)==TRUE))
    cont_cols = names(which(lapply(sim_df, is.integer)==FALSE))
    bi_sim = sim_df[binary_cols]
    cont_sim = sim_df[cont_cols]
    pos_check = lapply(bi_sim, function(x) sum(x)/length(x))
    all_1s = names(which(pos_check==1))
    bi_sim[1, all_1s]=0
    all_0s = names(which(pos_check==0))
    bi_sim[1, all_0s]=1
    names_order = colnames(sim_df)
    sim_df = cbind(bi_sim, cont_sim)
    sim_df = sim_df[,names_order]
  }
  
  return(sim_df)
}

apply_methods = function(exposure, outcome, covariates=NULL, sb=NULL, df, ratio=1, match_methods = NULL){
  re = function(df, name){
    rownames(df) = name
    return(df)
  }
  div=4
  #create empty data frame
  tot_df = data.frame("odds_ratio" = as.numeric(),
                      beta = as.numeric(),
                      lower_int = as.numeric(),
                      upper_int = as.numeric(),
                      confint_diff = as.numeric(),
                      p_val = as.numeric(),
                      n = as.numeric())
  if(is.null(match_methods)==F & is.null(covariates)==T){
    warning("Matching methods cannot be used if no covariates are given and thus were not run.")
  }
  
  #set selection bias here
  if(length(names(which(unlist(lapply(df[,sb],class))=="numeric")))!=0){
    cont_sb = names(which(unlist(lapply(df[,sb],class))=="numeric"))
    
    warning(paste0("the following nodes are continuous and cannot be in the sb argument: ",
                   paste0(cont_sb, collapse = ", ")), call. = FALSE)
  }
  if(is.null(sb)==F){
    for(i in sb){
      sb_df = bi_strat(1, df, i)[[1]]
      df = sb_df
    }
    
  }
  #rewriting for less if/else conditions
  
  #risk_ratio_quals
  # if(class(df[,outcome])=="integer" & length(find.package("logisticRR", quiet=TRUE))==1){
  #   tot_df = tot_bind(list(tot_df, risk_ratio(exposure, outcome, covariates, df = df)))
  # }
  
  #odds ratio quals
  if(class(df[,outcome])=="integer"){
    tot_df = tot_bind(list(tot_df, odds_ratio(exposure, outcome, covariates, df = df)))
  }
  
  #lm beta quals
  if(class(df[,outcome])=="numeric"){
    tot_df = tot_bind(list(tot_df, lm_beta(exposure, outcome, covariates, df)))
  }
  
  #ps weighting quals
  if(is.null(covariates)==F){
    ps_df = get_ps(exposure, covariates, df)
    tot_df =  tot_bind(list(tot_df,
                            ps_weight(exposure, outcome, covariates, ps_df, "weights")))
    
  }
  
  #matching quals
  if(length(find.package("MatchIt", quiet=TRUE))!=1 & is.null(match_methods)==F){
    warning("You've entered matching methods without the MatchIt package installed. Matching Analyses have been skipped")
  }
  
  if(is.null(match_methods)==F & length(find.package("MatchIt", quiet=TRUE))==1){
    di_df = df
    if(class(df[,exposure]) == "numeric"){
      warning("Matching methods require a binary exposure. The exposure has been dichotomized so the top 25% quartile are 1s and the bottom 75% are 0s")
      di_df = dichotomize(exposure, df, div)
    }
    
    match_df_list = lapply(match_methods, function(x) matchit_matching(exposure, covariates, di_df, d = x, ratio))
    if(class(df[,outcome]) == "numeric"){
      match_dfs = tot_bind(lapply(c(1:length(match_df_list)), function(x) re(lm_beta(exposure, outcome, df = match_df_list[[x]]), match_methods[[x]])
      ))
    }
    if(class(df[,outcome]) == "integer"){
      match_dfs = tot_bind(lapply(c(1:length(match_df_list)), function(x) re(odds_ratio(exposure, outcome, df = match_df_list[[x]]), match_methods[[x]])
      ))
    }
    
    tot_df = tot_bind(list(tot_df,match_dfs))
  }
  
  tot_df = as.data.frame(apply(tot_df, 2, unlist))
  return(tot_df)
}

varied_runs = function(runs, dag, exposure, outcome, covariates=NULL, sb=NULL, n=10000, positivity = F, misdiagnosis_v = outcome, under_r = 0, over_r = 0, ratio=1, match_methods = NULL, ...){
  randomize = function(variable, rmodel){
    if(is.null(variable) == TRUE){
      variable = rmodel
    } else {
      variable = rep(variable, runs)
    }
  }
  
  n = randomize(n, as.integer(runif(runs, 1000, 100000)))
  under_r = randomize(under_r, runif(runs, 0, 1))
  over_r = randomize(over_r, runif(runs, 0, 1))
  
  value_df = data.frame(n = n,
                        under_r = under_r,
                        over_r = over_r)
  #capture the set ate
  temp_dag = make_model(dag, ...)
  ate_set = temp_dag$dag[exposure, outcome]
  b_value = capture.output(temp_dag$jags)
  outcome_eq = grep(paste0("   ", outcome), b_value, value = T)
  outcome_eq = gsub(" ", "", outcome_eq, fixed = T)
  outcome_eq = unlist(strsplit(outcome_eq, ")", fixed = T))
  outcome_eq = unlist(strsplit(outcome_eq, "(", fixed = T))
  outcome_eq = unlist(strsplit(outcome_eq, "+", fixed = T))
  outcome_val = grep(exposure, outcome_eq, value = T)
  outcome_val = gsub(exposure, 1, outcome_val)
  outcome_val = gsub("[[:alpha:]]", "0", outcome_val)
  ate = sum(unlist(lapply(outcome_val, function(x) eval(parse(text = x)))))
  
  
  #temp_dag = lapply(c(1:runs), function(x) make_model(dag, value_df[x,1], value_df[x,3], value_df[x,4]))
  
  #FIX DIMENSION PROBLEM
  temp_df = lapply(c(1:runs), function(x) create_data(dag, value_df[x,1], positivity = positivity, ...))
  temp_df = lapply(c(1:runs), function(x) misdiagnosis(temp_df[[x]], misdiagnosis_v, under_r[x], over_r[x]))
  temp_output = lapply(temp_df, apply_methods, exposure = exposure, outcome = outcome, covariates = covariates, sb = sb, ratio=ratio, match_methods=match_methods)
  
  one_dim = FALSE
  if(names(temp_output[[1]])[1]=="apply(tot_df, 2, unlist)"){
    temp_output = lapply(1:runs,function(x) t(temp_output[[x]]))
    one_dim = TRUE
  }
  
  if(class(temp_df[[1]][,outcome])=="integer"){
    out_p = unlist(lapply(c(1:runs), function(x) sum(temp_df[[x]][,outcome])/nrow(temp_df[[x]])))
    
  }else{
    out_p=rep(NA, runs)
    
  }
  if(class(temp_df[[1]][,exposure])=="integer"){
    exp_p = unlist(lapply(c(1:runs), function(x) sum(temp_df[[x]][,exposure])/nrow(temp_df[[x]])))
    #mde = unlist(lapply(c(1:runs), function(x) sum(temp_df[[x]][,outcome])/nrow(temp_df[[x]])))
    #unlist(lapply(colnames(x_val), function(x) 0.02*sqrt(1 / (run[[3]][,x]*(1 - run[[3]][,x])*run[[5]][,x] ) )))
  }else{
    exp_p=rep(NA, runs)
    #mde=rep(NA, runs)
  }
  temp_output =lapply(c(1:runs), function(x) cbind(temp_output[[x]], exp_prev = rep(exp_p[x], nrow(temp_output[[x]]))))
  temp_output =lapply(c(1:runs), function(x) cbind(temp_output[[x]], out_prev = rep(out_p[x], nrow(temp_output[[x]]))))
  ##FIX SET ATE HERE
  temp_output =lapply(c(1:runs), function(x) cbind(temp_output[[x]], set_ate = rep(ate, nrow(temp_output[[x]]))))
  temp_output =lapply(c(1:runs), function(x) cbind(temp_output[[x]], over_r = rep(over_r[x], nrow(temp_output[[x]]))))
  temp_output =lapply(c(1:runs), function(x) cbind(temp_output[[x]], under_r = rep(under_r[x], nrow(temp_output[[x]]))))
  partition = laply(temp_output, as.matrix)
  if(one_dim == TRUE){
    partition = list(odds_ratio = partition[, 1],
                     calculated_ate = partition[, 2],
                     lower_int = partition[, 3],
                     upper_int = partition[, 4],
                     p_values = partition[, 6],
                     exp_prevalence = partition[, 8],
                     out_prevalence = partition[, 9],
                     sample_population = partition[, 7],
                     set_ate = partition[, 10],
                     over_r = partition[,11],
                     under_r = partition[,12])
    # partition = as.matrix(partition)
    # class(partition) = c("ScenarioMatrix", "matrix")
    return(partition)
  }
  run_list = list(ratio = partition[,, 1],
                  calculated_ate = partition[,, 2],
                  lower_int = partition[,, 3],
                  upper_int = partition[,, 4],
                  p_values = partition[,, 6],
                  exp_prevalence = partition[,, 8],
                  out_prevalence = partition[,, 9],
                  sample_population = partition[,, 7],
                  set_ate = partition[,, 10],
                  over_r = partition[,,11],
                  under_r = partition[,,12])
  # class(run_list) = "ScenarioMatrix"
  
  
  return(run_list)
}

reparse_runs = function(run_list, method = NULL, list_names=as.character(1:length(run_list))){
  extract_method = function(run,method){
    m_list = lapply(c(1:length(run)),function(x) run[[x]][,colnames(run[[1]]) == method])
    names(m_list) = names(run)
    m_list = as.data.frame(m_list)
    return(m_list)
  }
  if(length(run_list)==1){
    stop("must compare at least 2 different SenarioMatrices in order to reparse")
    # rep_size = nrow(as.data.frame(run_list))
    # dim_holder = data.frame("ratio" = rep(NA, rep_size),
    #                         "calculated_ate" = rep(NA, rep_size),
    #                         "lower_int" = rep(NA, rep_size),
    #                         "upper_int" = rep(NA, rep_size),
    #                         "p_values" = rep(NA, rep_size),
    #                         "exp_prevalence" = rep(NA,rep_size),
    #                         "out_prevalence" = rep(NA,rep_size),
    #                         "sample_population" = rep(NA,rep_size),
    #                         "set_ate" = rep(NA,rep_size),
    #                         "over_r"= rep(NA,rep_size),
    #                         "under_r" = rep(NA,rep_size))
    # run_list = list(dim_holder, run_list)
  }
  names(run_list) = list_names
  index = unlist(lapply(1:length(run_list), function(x) class(run_list[[x]][[1]])[1]))
  names(index)=list_names
  
  list_list = run_list[which(index != "numeric")]
  mat_list = run_list[which(index == "numeric")]
  mat_names = names(mat_list)
  
  if(length(mat_list) !=0){
    mat_list = lapply(c(1:length(mat_list)), function(x) as.data.frame(mat_list[[x]]))
    names(mat_list) = mat_names
  }
  
  if(is.null(method) == F && length(list_list)>0){
    extraction = lapply(list_list, function(x) extract_method(x,method))
    list_names = names(extraction)
  }else if(is.null(method) == T && length(list_list)>0){
    method_names = lapply(1:length(list_list), function(x) dimnames(list_list[[x]][[1]])[[2]])
    names(method_names)=names(list_list)
    extraction = lapply(names(list_list), function(x)
      lapply( method_names[[x]], function(y) extract_method(list_list[[x]], y) )
    )
    names(extraction) = names(list_list)
    header_names = unlist(lapply(1:length(method_names), function(x) rep(names(method_names)[[x]], length(method_names[[x]]))))
    method_names_2 = unname(unlist(method_names))
    expand_names = unlist(lapply(1:length(method_names_2), function(x) paste0(header_names[[x]], "_", method_names_2[[x]])))
    extraction = unlist(extraction, recursive = F)
    names(extraction) = expand_names
    
  }else{
    extraction = list_list
  }
  
  
  
  tot_list = append(extraction, mat_list)
  
  tot_mat = laply(tot_list, as.matrix)
  
  partition = aperm(tot_mat, c(2,1,3))
  colnames(partition) = c(names(extraction),mat_names)
  out_df = list(ratio = partition[,, 1],
                calculated_ate = partition[,, 2],
                lower_int = partition[,, 3],
                upper_int = partition[,, 4],
                p_values = partition[,, 5],
                exp_prevalence = partition[,, 6],
                out_prevalence = partition[,, 7],
                sample_population = partition[,, 8],
                set_ate = partition[,, 9],
                over_r = partition[,,10],
                under_r = partition[,,11])
  return(out_df)
}

#visualization#####
ci_ridges = function(run, title =NULL, subtitle=NULL){
  ate_val = as.data.frame(run[[2]])
  drawn_ci = beta_sum(run)
  
  if(is.null(dim(run[[1]]))==T){
    colnames(ate_val) = "regression"
    
  }
  
  all_in_ci = unlist(lapply(c(1:ncol(ate_val)), function(x){
    if(drawn_ci[7,x]==0&drawn_ci[8,x]==0){
      return(TRUE)
    }else{
      return(FALSE)
    }
  } ))
  
  drawn_ci[drawn_ci==0] =NA
  drawn_ci[drawn_ci==1] =NA
  
  # drawn_ci$tag = c(1,1,0,0)
  # drawn_ci = drawn_ci[!(duplicated(drawn_ci[1:3]) | duplicated(drawn_ci[1:3], fromLast = TRUE)), ]
  # drawn_ci = drawn_ci[drawn_ci$tag==1,]
  # drawn_ci = drawn_ci[-length(drawn_ci)]
  upper = sapply(drawn_ci[7,], function(x) rep(x, length(ate_val[[1]])))
  lower = sapply(drawn_ci[8,], function(x) rep(x, length(ate_val[[1]])))
  #rep(drawn_ci[,1],length(comp_df[[1]]))
  
  
  names= unlist(lapply(colnames(ate_val), function(y) rep(y, nrow(ate_val))))
  value= as.numeric(unlist(ate_val))
  upper = 1- as.numeric(unlist(upper))
  lower = as.numeric(unlist(lower))
  
  data=dplyr::tibble(names,value, upper, lower)
  data$tag =0
  data[data$value>0,]$tag = 1
  data$names = as.factor(data$names)
  data$tag = as.character(data$tag)
  
  lines = drawn_ci[1:2, , drop = FALSE]
  nas_df = drawn_ci[7:8, , drop = FALSE]
  #na_dex
  inf_dex = which(is.na(nas_df[1,]))
  neg_inf_dex = which(is.na(nas_df[2,]))
  all_in_dex = which(all_in_ci==TRUE)
  
  lines[1, inf_dex]= Inf
  lines[2, neg_inf_dex]= -Inf
  
  
  
  factor_df = as.data.frame(unique(as.integer(data$names)))
  factor_df = t(factor_df)
  colnames(factor_df) = unique(as.character(data$names))
  rownames(factor_df) = "i"
  
  #FIX HERE, col name problem
  lines = rbind(lines, factor_df, all_in_ci)
  # rownames(lines[4, ]) = "all" #perhaps unneseccary??
  
  #lines = lines[,order(lines[nrow(lines),])]
  
  drawn_ci = rbind(drawn_ci, factor_df)
  #drawn_ci = drawn_ci[,order(drawn_ci[nrow(drawn_ci),])]
  
  ###WORKING LINES
  p <- ggplot2::ggplot(data, aes(x=value, y=names)) +
    ggridges::stat_density_ridges(scale = 0.95,
                                  quantile_lines = TRUE,
                                  quantile_fun = function(x, ...) quantile(x, probs =
                                                                             c(sort(c(mean(data[data$value == x,]$lower))), sort(c(mean(data[data$value == x,]$upper)))), na.rm = TRUE)
    ) +
    ggridges::theme_ridges(center = TRUE) +
    ggplot2::ylab("Anaylsis Performed") +
    ggplot2::xlab("Estimated Beta") +
    ggplot2::ggtitle(label = paste(title),
                     subtitle = paste(subtitle)) +
    ggplot2::theme(plot.title = element_text(size = 18, face = "bold"),
                   plot.subtitle = element_text(size = 12),
                   axis.text = element_text(size = 14, face = "bold"),
                   axis.title = element_text(size = 14, face = "bold")) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_discrete(expand = expansion(mult = c(0.01, .1))) +
    if(max(data$value) > 10 | min(data$value < -10)){
      ggplot2::xlim(limits[1],limits[2])
    }
  p
  
  d <- ggplot_build(p)$data[[1]]
  ribbon = function(upper, lower, i, all) {
    if(upper == Inf & lower == -Inf & all == 0){
      return()
    }
    q = ggplot2::geom_ribbon(
      data = transform(subset(d, x <= upper & x >= lower & ymin == i), names = group),
      ggplot2::aes(x, ymin = ymin, ymax = ymax, group = group),
      fill = "lightblue2")
    return(q)
  }
  # i = 3
  # p + ribbon(lines[1,i], lines[2,i], lines[3,i])
  # p + ribbon(-0.1, -0.2, 2)
  q = p + lapply(c(1:ncol(lines)), function(x) ribbon(lines[1,x],lines[2,x],lines[3,x], lines[4,x]))
  
  # q =p + geom_ribbon(
  #   data = transform(subset(d, x <= 0.0373545924 & x >= -0.0412792948 & ymin == 3), names = group),
  #   aes(x, ymin = ymin, ymax = ymax, group = group),
  #   fill = "lightblue2") +
  #   geom_ribbon(
  #     data = transform(subset(d, x <= 0.0348230349 & x >= -0.0357772047  & ymin == 2), names = group),
  #     aes(x, ymin = ymin, ymax = ymax, group = group),
  #     fill = "lightblue2") +
  #   geom_ribbon(
  #     data = transform(subset(d, ymin == 1), names = group),
  #     aes(x, ymin = ymin, ymax = ymax, group = group),
  #     fill = "lightblue2") +
  #   geom_segment( aes(y=1, yend=length(colnames(ate_val))+1, x=run[[8]][1], xend=run[[8]][1]), color="navy", linetype = "dashed", lwd = 1)
  
  r = q+ ggridges::stat_density_ridges(scale = 0.95,
                                       quantile_lines = TRUE,
                                       quantile_fun = function(x, ...) quantile(x, probs =
                                                                                  c(sort(c(mean(data[data$value == x,]$lower))), sort(c(mean(data[data$value == x,]$upper)))), na.rm = TRUE),
                                       fill = "lightblue2",
                                       alpha= 0.01) +
    ggplot2::geom_vline( xintercept = run[[9]][1], color="navy", linetype = "dashed", lwd = 1)
  
  
  
  
  ###
  
  
  
  # Construct the six grobs - three symbols and three labels
  L1 = grid::rectGrob(height = .5, width = .5, gp = gpar(fill = "lightblue2", col = NA))
  L2 = grid::rectGrob(height = .5, width = .5, gp = gpar(fill = "grey50", col = NA))
  T1 = grid::textGrob("Yes", x = .2, just = "left")
  T2 = grid::textGrob("No", x = .2, just = "left")
  
  
  # Construct a gtable - 2 columns X 4 rows
  leg = gtable::gtable(width = unit(c(1,1), "cm"), height = unit(c(1.8,1,1), "cm"))
  
  # Place the six grob into the table
  leg = gtable::gtable_add_grob(leg, L1, t=2, l=1)
  leg = gtable::gtable_add_grob(leg, L2, t=3, l=1)
  leg = gtable::gtable_add_grob(leg, T1, t=2, l=2)
  leg = gtable::gtable_add_grob(leg, T2, t=3, l=2)
  
  # Give it a title (if needed)
  leg = gtable::gtable_add_grob(leg, grid::textGrob(expression(bold("True B in\n95% CI?")), vjust = 2), t=1, l=1, r=2)
  #leg = gtable_add_grob(leg, textGrob(expression(bold("95% CI?"))), t=2, l=1, r=2)
  # Get the ggplot grob for plot1
  g = ggplot2::ggplotGrob(r)
  
  # Get the position of the panel,
  # add a column to the right of the panel,
  # put the legend into that column,
  # and then add another spacing column
  pos = g$layout[grepl("panel", g$layout$name), c('t', 'l')]
  g = gtable::gtable_add_cols(g, sum(leg$widths), pos$l)
  g = gtable::gtable_add_grob(g, leg, t = pos$t, l = pos$l + 1)
  g = gtable::gtable_add_cols(g, unit(6, "pt"), pos$l)
  
  # Draw it
  grid::grid.newpage()
  return(grid::grid.draw(g))
  
}

#####
#Ui Facing Functions
#####

dag_ui = function(dag_string){
  dag = HydeNetwork(eval(str2lang(dag_string)))

  return(plot(dag))
}

get_nodes = function(dag_string){
  if(class(tryCatch(HydeNetwork(eval(str2lang(dag_string))), error =  function(x) x=1)) =="HydeNetwork"){
    dag = HydeNetwork(eval(str2lang(dag_string)))
    node_test = dag[["nodes"]]
    return(node_test)
  }

}


get_arrows = function(dag_string){
  if(class(tryCatch(HydeNetwork(eval(str2lang(dag_string))), error =  function(x) x=1)) =="HydeNetwork"){
    dag = HydeNetwork(eval(str2lang(dag_string)))
    parents = unlist(dag[["parents"]])
    children = names(parents)
    children = gsub('[[:digit:]]+', '', children)
    arrow_list = unlist(lapply(1:length(children), function(x) paste0(parents[x], " -> ", children[x])))
    return(arrow_list)
  }
}

handler_df = function(x){
  test_df = data.frame(unlist(lapply(x(), function(handle) {
    handle()[["distribution"]]
  })), holder = 1)
  return(rownames(test_df[test_df[,1]=="binary",]))
}

get_dist = function(x, y){
  output = as.character(lapply(y(), function(handle) {
    handle()[["distribution"]]})[["exposure"]])

  # output = gsub("binary", "dbern(", output )
  # output = gsub("continuous", "dnorm(", output )
  output



}

get_sum_stats = function(x){
  as.character(lapply(handler(), function(handle) {
    handle()[["distribution"]]})[["exposure"]])
  as.double(lapply(handler(), function(handle) {
    handle()[[2]]})[["exposure"]])
}

run_code = function(out_code){

  eval(parse(text = out_code))
  try_plot = ci_ridges(run_1)



  return(try_plot)


}

