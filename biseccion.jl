# Bisection algorithm for one variable real continuous functions
# Author: Arturo Erdely
# Last update: 2023-10-04

"""
    biseccion(f, a, b; δ = abs((a + b)/2) / 1_000_000, m = 10_000)

Bisection method: If `f` is a one variable real function over the closed interval [a,b]
such that `f(a)f(b)<0` and with only one root in such interval, with a maximum absolute
error `δ` a value `c` is to be found such that `f(c)≈0` but if a maximum numer of
iterations `m` is reached the algorithm stops with a warning.

## Example
```
f(x) = (x - 3) * (x - 1) * (x + 1)
biseccion(f, -1.9, 0)
```
"""
function biseccion(f, a, b; δ = abs((a + b)/2) / 1_000_000, m = 10_000)
    iter = 1
    z = (a + b) / 2
    while iter ≤ m
        if f(z) == 0.0 || (b - a)/2 ≤ δ
            break # stopping criteria
        elseif f(a) * f(z) > 0
            a = z
        else
            b = z
        end
        z = (a + b) / 2
        iter += 1
    end
    iter -= 1
    if iter == m
        @warn "The maximum number of $m iterations has been performed before reaching stopping criteria."
    end
    return (raíz = z, dif = f(z), numiter = iter, maxiter = m, tol = δ)
end
