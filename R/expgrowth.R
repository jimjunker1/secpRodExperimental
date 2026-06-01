#'
#'
#'
#'

expSim <- function(t = NULL, m_init = NULL, Ce = NULL, cpi = NULL, m_cpi = NULL) {
  
  if (is.null(m_init)) {
    stop("Error: 'm_init' must be provided.")
  }
  
  if (m_init <= 0) {
    stop("Error: 'm_init' must be positive.")
  }
  
  if (is.null(Ce) && (is.null(cpi) || is.null(m_cpi))) {
    stop("Error: provide either 'Ce', or both 'cpi' and 'm_cpi' to estimate Ce.")
  }
  
  if (!is.null(Ce) && Ce <= 0) {
    stop("Error: 'Ce' must be positive.")
  }
  
  if (!is.null(m_cpi)) {
    if (m_cpi <= 0) {
      stop("Error: 'm_cpi' must be positive.")
    }
    
    if (m_cpi <= m_init) {
      stop("Error: 'm_cpi' must be greater than 'm_init' for positive exponential growth.")
    }
  }
  
  if (!is.null(cpi) && cpi <= 0) {
    stop("Error: 'cpi' must be positive.")
  }
  
  # Exponential mass growth
  pred_mass <- function(t, Ce) {
    m_init * exp(Ce * t)
  }
  
  # Estimate Ce from cpi and target mass if needed
  if (is.null(Ce)) {
    Ce <- log(m_cpi / m_init) / cpi
  }
  
  # Estimate cpi from Ce and target mass if needed
  if (is.null(cpi) && !is.null(m_cpi)) {
    cpi <- log(m_cpi / m_init) / Ce
  }
  
  # Return parameters if no time vector supplied
  if (is.null(t)) {
    return(list(
      m_init = m_init,
      m_cpi = m_cpi,
      cpi = cpi,
      Ce = Ce
    ))
  }
  
  # Return predicted masses
  pred_mass(t, Ce)
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
