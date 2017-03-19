bsc <- function(S, T, t, K, r, s, q)
{
    d_plus <- (log(S/K) + (r - q + 0.5 * s^2) * (T - t)) / (s * sqrt(T-t))
    d_minus <- d_plus - s * sqrt(T-t)
    S * exp(-q*(T-t)) * pnorm(d_plus) - K * exp(-r*(T-t)) * pnorm(d_minus)
}

sigmas <- seq(0.05, 0.5, by=0.01)
fsig <- bsc(50, 0.5, 0.0, 45, 0.06, sigmas, 0.02) - 7

print("")
print("Using brute-force:")
bsc(50, 0.5, 0.0, 45, 0.06, 0.25, 0.02) - 7

fsig <- function(sigma)
    bsc(50, 0.5, 0.0, 45, 0.06, sigma, 0.02) - 7

# epsilon is tolerance, which is the error threshold.
bisection <- function(f, a, b, epsilon=0.001)
{
    while (b - a > epsilon)
    {
        mid <- (b + a) / 2
        if (sign(f(mid)) == sign(f(a)))
            a <- mid
        else
            b <- mid
    }
    (a + b) / 2
}

print("")
print("Using bisection:")
bsc(50, 0.5, 0.0, 45, 0.06, bisection(fsig, 0.1, 0.3), 0.02) - 7

# derivative of bsc with respect to σ
vega <- function(S,T,t,K,r,s,q)
{
    d_plus <- (log(S/K) + (r - q + 0.5 * s^2) * (T - t)) / (s * sqrt(T-t))
    S * exp(-q*(T-t)) * (1 / sqrt(2*pi)) * exp(-d_plus^2/2) * sqrt(T-t)
}

dfsig <- function(sigma)
    vega(50, 0.5, 0.0, 45, 0.06, sigma, 0.02)

u <- 100
x <- 1
i <- 0
epsilon <- 1e-6
while (u / x > epsilon)
{
    u <- fsig(x) / dfsig(x)
    x <- x - u
    i <- i + 1
}

print("")
print("Using newton:")
fsig(x) 