#'
#'
#'
#'

vbgSim <- function(t = NULL, m_init = NULL, m_inf = NULL, k = NULL, cpi = NULL, m_approx = 0.99){
  if(is.null(m_init) | is.null(m_init)) stop("Error: both 'm_init' and 'm_inf' must be provided.")
  if(is.null(k) & is.null(cpi)) stop("Error: either 'k' or 'cpi' must be provided.")
  # if no time vector is provided, the function estimates all the vbg parameters 
  # based on the provided parameters
  if(is.null(t)){
    if(is.null(k)){
      # objective function
      f <- function(k){
        t0 = log(1 - (m_init/m_inf)^(1/3))/k
        
        mt <- m_inf * (1 - exp(-k * (cpi - t0)))^3
        
        mt - (m_approx * m_inf)
      }
      k = uniroot(f, interval = c(1e-6,10))$root
      
      t0 = log(1 - (m_init/m_inf)^(1/3)) / k
      
      return(list(m_init = m_init,
                  m_inf = m_inf,
                  cpi = cpi,
                  k = k,
                  t0 = t0))
    }
    if(is.null(cpi)){
      t0 = log(1 - (m_init/m_inf)^(1/3)) / k
      cpi = t0 - log(1 - m_approx^(1/3)) / k
      
      return(list(m_init = m_init,
                  m_inf = m_inf,
                  cpi = cpi,
                  k = k,
                  t0 = t0))
    }
  } else{
  # if a vector of times is provided estimate mass at each time
  if(is.null(k)){
    # objective function
    f <- function(k){
      t0 = log(1 - (m_init/m_inf)^(1/3))/k
      
      mt <- m_inf * (1 - exp(-k * (cpi - t0)))^3
      
      mt - (m_approx * m_inf)
    }
    k = uniroot(f, interval = c(1e-6,10))$root
    
    t0 = log(1 - (m_init/m_inf)^(1/3)) / k
  } else{ 
    t0 = log(1 - (m_init/m_inf)^(1/3)) / k
    }
   mtVec =  m_inf*(1 - exp(-k * (t + t0)))^3 
   return(mtVec)
  }
}
