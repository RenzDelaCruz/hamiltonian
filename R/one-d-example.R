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

# modified eulers method
euler.mod <- function(steps, epsilon, q.0=0, p.0=1) {
    # p_i(t+epsilon) = p_i(t) - epsilon partial U/partial q_i(q(t))
    # q_i(t+epsilon) = q_i(t) + epsilon p_i(t+epsilon)
    
    # initialize
    p <- rep(p.0, steps + 1)
    q <- rep(q.0, steps + 1)
    
    for (i in 2:(steps+1)) {
        p[i] <- p[i-1] - epsilon * q[i-1]
        q[i] <- q[i-1] + epsilon * p[i]
    }
    
    return(cbind(q, p))
}

# leapfrog method
leapfrog <- function(steps, epsilon, q.0=0, p.0=1) {
    # p_i(t+epsilon/2) = p_i(t) - epsilon/2 partial U/partial q_i(q(t))
    # q_i(t+epsilon) = q_i(t) + epsilon p_i(t+epsilon/2)
    # p_i(t+epsilon) = p_i(t+epsilon/2) - epsilon/2 partial U/partial q_i(q(t+epsilon))
    
    # initialize
    p <- rep(p.0, steps + 1)
    q <- rep(q.0, steps + 1)
    
    for (i in 2:(steps+1)) {
        p.half <- p[i-1] - epsilon / 2 *q[i-1]
        q[i] <- q[i-1] + epsilon * p.half
        p[i] <- p.half - epsilon / 2 * q[i]
    }
    
    return(cbind(q, p))
}

euler.values.1 <- euler(20, 0.3)
euler.values.2 <- euler(60, 0.1)
euler.mod.values.1 <- euler.mod(20, 0.3)
euler.mod.values.2 <- euler.mod(60, 0.1)
leapfrog.values.1 <- leapfrog(20, 0.3)
leapfrog.values.2 <- leapfrog(60, 0.1)

# plot values
data <- rbind(data.frame(euler.values.1, method="Euler's, epsilon=0.3"),
              data.frame(euler.values.2, method="Euler's, epsilon=0.1"),
              data.frame(euler.mod.values.1, method="Mod. Euler's, epsilon=0.3"),
              data.frame(euler.mod.values.2, method="Mod. Euler's, epsilon=0.1"),
              data.frame(leapfrog.values.1, method="Leap-frog, epsilon=0.3"),
              data.frame(leapfrog.values.2, method="Leap-frog, epsilon=0.1"))

library(ggplot2)
theme_set(theme_classic())
ggplot(data, aes(x=q, y=p)) +
    geom_point() +
    expand_limits(x=c(-2, 2), y=c(-2, 2)) +
    annotate("path", 
             x=cos(seq(0, 2*pi, len=100)), 
             y=sin(seq(0, 2*pi, len=100))) +
    facet_wrap(~method, ncol=2)
