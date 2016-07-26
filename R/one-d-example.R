# Hamiltonian example in one dimension
# H(q, p) = U(q) + K(p), U(q)=q^2/2, K(p)=p^2/2
# Corresponds to Gaussian distribution for q with mean zero and variance one

# eulers method
euler <- function(steps, epsilon, q.0=0, p.0=1) {
    # p_i(t+epsilon) = p_i(t) - epsilon partial U/partial q_i(q(t))
    # q_i(t+epsilon) = q_i(t) + epsilon p_i(t)
    
    # initialize
    p <- rep(p.0, steps + 1)
    q <- rep(q.0, steps + 1)
    
    for (i in 2:(steps+1)) {
        p[i] <- p[i-1] - epsilon * q[i-1]
        q[i] <- q[i-1] + epsilon * p[i-1]
    }
    
    return(cbind(q, p))
}

euler.values <- euler(20, 0.3)

# plot values
library(ggplot2)
theme_set(theme_classic())
ggplot(data.frame(euler.values), aes(x=q, y=p)) +
    geom_point() +
    expand_limits(x=c(-2, 2), y=c(-2, 2)) +
    annotate("path", 
             x=cos(seq(0, 2*pi, len=100)), 
             y=sin(seq(0, 2*pi, len=100)))
