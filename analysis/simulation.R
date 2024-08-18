library(tidyverse)
library(mcmapper)
library(pROC)

set.seed(2024)

prevs <- seq(from = 0.01, to = 0.50, by = 0.01)
c_stats <- seq(from = 0.51, to = 0.99, by = 0.01)

param_grid <- expand.grid(prev = prevs, c_stat = c_stats,type=c("beta","logitnorm",'probitnorm'))

calculate_sol <- T

if(calculate_sol){

    tmp_sol <- apply(param_grid, 1, function(tmp_row) {
      c(tmp_row,as.numeric(mcmap(as.numeric(tmp_row[1:2]),tmp_row[3])$value))
    }) %>%
      t() %>%
      as.data.frame() %>%
      mutate(across(c(1:2,4:5),as.numeric))

    # sanity-check: unique values
    tmp_sol %>%
      group_by(type) %>%
      distinct(V4,V5) %>%
      tally() -> sanity_check

    sanity_check

  write_rds(tmp_sol,
            "results/algos_sol.rds")
}

df_algo_sol <- read_rds("results/algos_sol.rds") %>%
  rename(arg1=V4,
         arg2=V5)

sim_dir <- "simulation_results"
dir.create(sim_dir)

n_outer <- 10000

se <- function(x) {sqrt(var(x)/length(x))}

choose_ralgo <- function(type){
  ifelse(type=="logitnorm", rlogitnorm,
         ifelse(type=="probitnorm", rprobitnorm,
                ifelse(type=="beta",rbeta,NA)
         )
  )
}

alpha_coverage <- 0.90
alpha <- (1-alpha_coverage)/2
qcrit <- abs(qnorm(alpha))

# parallelization on mac
library(doParallel)
library(foreach)
registerDoParallel(cl <- makeCluster(7))

foreach(i=1:nrow(df_algo_sol),.combine = "+") %dopar% {

  library(mcmapper)
  set.seed(i)

  tmp_pi_check <- tmp_cstat_check <- c()
  tmp_row <- df_algo_sol[i,]
  tmp_true <- c(tmp_row$prev,tmp_row$c_stat)
  ralgo <- choose_ralgo(tmp_row$type)

  n_inner <- as.numeric(1/tmp_row$prev*10)

  for(j in 1:n_outer){

    tmp_pi <- ralgo(n_inner,tmp_row$arg1,tmp_row$arg2)
    tmp_pi_mean <- mean(tmp_pi)
    tmp_pi_se <- se(tmp_pi)
    tmp_pi_CI <- c(tmp_pi_mean-qcrit*tmp_pi_se,tmp_pi_mean + qcrit*tmp_pi_se)

    tmp_Y <- rbinom(n_inner,1,tmp_pi)

    if(length(unique(tmp_Y))!=2){
      next
    }

    tmp_roc <- pROC::roc(tmp_Y~tmp_pi, ci=T,quiet=T,ci.method = c("delong"),conf.level=alpha_coverage)
    tmp_cstat <- as.vector(tmp_roc$ci)
    tmp_cstat_mean <- tmp_cstat[2]
    tmp_cstat_CI <- c(tmp_cstat[1],tmp_cstat[3])

    tmp_pi_check <- c(tmp_pi_check,as.numeric((tmp_true[1]<=tmp_pi_CI[2]) & (tmp_true[1]>=tmp_pi_CI[1])))
    tmp_cstat_check <-c(tmp_cstat_check,as.numeric((tmp_true[2]<=tmp_cstat_CI[2]) & (tmp_true[2]>=tmp_cstat_CI[1])))
  }

  res <- c(unlist(tmp_row),coverage_prob_arg1=mean(tmp_pi_check,na.rm=T),coverage_prob_arg2=mean(tmp_cstat_check,na.rm=T))
  readr::write_rds(res,paste0(sim_dir,"/",i,".rds"))
  return(0)
}

stop(cl)

# process results
sim_dir <- "simulation_results"

sim_files <- list.files(sim_dir)

sim_results <- lapply(sim_files,function(tmp_sim){
  read_rds(paste0(sim_dir,"/",tmp_sim))
}) %>%
  do.call(rbind,.) %>%
  as.data.frame() %>%
  mutate(across(c(prev:c_stat,arg1:coverage_prob_arg2),as.numeric)) %>%
  mutate(coverage_prob_arg1=round(coverage_prob_arg1,2),
         coverage_prob_arg2 = round(coverage_prob_arg2,2)) %>%
  arrange(type) %>%
  filter(coverage_prob_arg1!=0) %>%
  mutate(type = case_when(type=="beta" ~ "Beta",
                          type == "logitnorm" ~ "Logit-Normal",
                          type == "probitnorm" ~ "Probit-Normal",
                          TRUE ~ NA)) %>%
  pivot_longer(cols=6:7,names_to="parameter",values_to="Coverage probability") %>%
  mutate(parameter = ifelse(parameter=="coverage_prob_arg1",
                            "Coverage for prevalence",
                            "Coverage for c-statistic"),
         parameter = factor(parameter,levels=c("Coverage for prevalence","Coverage for c-statistic")))

write_rds(sim_results,"results/simulation_results.rds")

sim_results <- read_rds("results/simulation_results.rds")

ggplot(data=sim_results,aes(x=prev,y=c_stat,fill=`Coverage probability`))+
  geom_tile() +
  facet_grid(parameter~type)+
  theme_classic() +
  xlab("Prevalence") +
  ylab("C-statistic") +
  theme(legend.title = element_text("Coverage probability"),
        text = element_text(size=20)) -> fig_sim

ggsave("figures/fig_sim.jpeg",plot=fig_sim,device = "jpg")
