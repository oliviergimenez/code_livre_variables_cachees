
get_init_moveHMM <- function(dta, nbStates = 2, id_point = NULL, zero_mass = FALSE){
  if(is.null(id_point)){
    id_point = 'id_point'
  }
  K <- nbStates
  fitted_data <- na.omit(dta)
  init <- fitted_data %>% 
    dplyr::select(angle, step) %>% 
    kmeans( centers = K, nstart = 500)
  cluster_kmeans <- tibble( cluster = init$cluster, !!id_point := as.numeric(names(init$cluster)))
  par_init <- fitted_data %>% 
    inner_join(cluster_kmeans) %>% 
    group_by(cluster) %>% 
    summarise(m_step = mean(step),
              sd_step = sd(step), 
              m_angle = mean(angle),
              zeros_prop = mean(step == 0))  %>% 
    arrange(m_step)
  
  mu0 <- par_init$m_step # step mean (three parameters: one for each state)
  sigma0 <- par_init$sd_step # step SD
  mass0 <- par_init$zeros_prop
  stepPar0 <- c(mu0,sigma0)
  if(zero_mass){
    stepPar0 <- c(mu0,sigma0, mass0)
  }
  angleMean0 <- par_init$m_angle # angle mean
  kappa0 <- rep(1,K) # angle concentration
  anglePar0 <- c(angleMean0,kappa0)
  ## call to fitting function
  return(list(stepPar0 = stepPar0,
           anglePar0 = anglePar0))
  
}

get_init_depmix <- function(dta, nbStates = 2){
  K <- nbStates
  fitted_data <- na.omit(dta)
  init <- fitted_data %>% 
    dplyr::select(v_p, v_r) %>% 
    na.omit() %>% 
    kmeans( centers = K, nstart = 500)
  fitted_data %>% 
    mutate(cluster = init$cluster) %>% 
    group_by(cluster) %>% 
    summarise(m_vp = mean(v_p),
              sd_vp = sd(v_p), 
              m_vr = mean(v_r),
              sd_vr = sd(v_r))  %>% 
    arrange(m_vp) %>% 
    dplyr::select(-cluster) %>% 
    as.matrix() %>% 
    t() %>% 
    as.numeric()
}
