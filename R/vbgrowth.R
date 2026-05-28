#'
#'
#'
#'

vbgSim <- function(t = NULL, m_init = NULL, m_inf = NULL, k = NULL, cpi = NULL, m_approx = 0.99){
  if(is.null(m_init) | is.null(m_init)) stop("Error: both 'm_init' and 'm_inf' must be provided.")
  if(is.null(k) & is.null(cpi)) stop("Error: either 'k' or 'cpi' must be provided.")
  if (m_init <= 0 || m_inf <= 0)stop("Error: 'm_init' and 'm_inf' must be positive.")
  if (m_init >= m_inf) stop("Error: 'm_init' must be less than 'm_inf'.")
  if (is.null(k) && is.null(cpi)) stop("Error: either 'k' or 'cpi' must be provided.")
  if (m_approx <= 0 || m_approx >= 1) stop("Error: 'm_approx' must be between 0 and 1.")
  
  # Positive-mass VBGF reparameterized by initial mass
  pred_mass <- function(t, k) {
    m_inf * (
      1 - (1 - (m_init / m_inf)^(1 / 3)) * exp(-k * t)
    )^3
  }
  
  # Estimate k from cpi if needed
  if (is.null(k)) {
    f <- function(k) {
      pred_mass(cpi, k) - m_approx * m_inf
    }
    
    k <- uniroot(f, interval = c(1e-8, 10))$root
  }
  
  # Estimate cpi from k if needed
  if (is.null(cpi)) {
    cpi <- -log(
      (1 - m_approx^(1 / 3)) /
        (1 - (m_init / m_inf)^(1 / 3))
    ) / k
  }
  
  # Return parameters if no time vector supplied
  if (is.null(t)) {
    return(list(
      m_init = m_init,
      m_inf = m_inf,
      cpi = cpi,
      k = k,
      m_approx = m_approx
    ))
  }
  
  # Return predicted masses
  pred_mass(t, k)
}
#   # if no time vector is provided, the function estimates all the vbg parameters 
#   # based on the provided parameters
#   if(is.null(t)){
#     if(is.null(k)){
#       # objective function
#       f <- function(k){
#         t0 = log(1 - (m_init/m_inf)^(1/3))/k
#         
#         mt <- m_inf * (1 - exp(-k * (cpi - t0)))^3
#         
#         mt - (m_approx * m_inf)
#       }
#       k = uniroot(f, interval = c(1e-6,10))$root
#       
#       t0 = log(1 - (m_init/m_inf)^(1/3)) / k
#       
#       return(list(m_init = m_init,
#                   m_inf = m_inf,
#                   cpi = cpi,
#                   k = k,
#                   t0 = t0))
#     }
#     if(is.null(cpi)){
#       t0 = log(1 - (m_init/m_inf)^(1/3)) / k
#       cpi = t0 - log(1 - m_approx^(1/3)) / k
#       
#       return(list(m_init = m_init,
#                   m_inf = m_inf,
#                   cpi = cpi,
#                   k = k,
#                   t0 = t0))
#     }
#   } else{
#   # if a vector of times is provided estimate mass at each time
#   if(is.null(k)){
#     # objective function
#     f <- function(k){
#       t0 = log(1 - (m_init/m_inf)^(1/3))/k
#       
#       mt <- m_inf * (1 - exp(-k * (cpi - t0)))^3
#       
#       mt - (m_approx * m_inf)
#     }
#     k = uniroot(f, interval = c(1e-6,10))$root
#     
#     t0 = log(1 - (m_init/m_inf)^(1/3)) / k
#   } else{ 
#     t0 = log(1 - (m_init/m_inf)^(1/3)) / k
#     }
#    mtVec =  m_inf*(1 - exp(-k * (t + t0)))^3 
#    mtVec
#    return(mtVec)
#   }
# }
