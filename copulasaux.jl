### Funciones auxiliares para `copulas.jl`
### Autor: Arturo Erdely
### Última actualización: 2026-04-03


"""
    Fn(x::Real, xobs::Vector{<:Real})

Función de distribución empírica evaluada en el punto `x`
con base en una muestra aleatoria observada `xobs`.
## Ejemplo
```
Fn(0.0, randn(1_000)) # valor cercano a 0.5
```
"""
function Fn(x::Real, xobs::Vector{<:Real})
    return count(xobs .≤ x) / length(xobs)
end

"""
    Fn(x::Vector{<:Real}, xobs::Vector{<:Real})

Función de distribución empírica evaluada en cada punto del vector `x`
con base en una muestra aleatoria observada `xobs`.
## Ejemplo
```
Fn([-1.96, 0.0, 1.96], randn(10_000)) # valores cercanos a [0.025, 0.5, 0.975]
```
"""
function Fn(x::Vector{<:Real}, xobs::Vector{<:Real})
    m = length(x)
    n = length(xobs)
    v = zeros(m)
    for k ∈ 1:m
        v[k] = count(xobs .≤ x[k]) / n
    end
    return v
end


"""
    poligonal(x::Vector{<:Real}, xobs::Vector{<:Real}; mínimo = minimum(xobs), máximo = maximum(xobs))

Aproximación poligonal de una función de distribución continua en un vector dado de valores `x`
con base en una muestra aleatoria observada dada por el vector `xobs`. El parámetro opcional
`mínimo` puede escogerse como un valor menor que el mínimo de la muestra, y `máximo` como
un valor mayor que el máximo de la muestra, si así se requiere.
greater value than the sample maximum, if required.
## Ejemplo
```
xobs = randn(10_000)
x = [-1.96, 0, 1.96]
poligonal(x, xobs) # valores cercanos a [0.025, 0.5, 0.975]
```
"""
function poligonal(x::Vector{<:Real}, xobs::Vector{<:Real}; mínimo = minimum(xobs), máximo = maximum(xobs))
    n = length(xobs)
    xord = sort(xobs)
    xp = (xord[1:(n-1)] .+ xord[2:n]) ./ 2
    xp = vcat(mínimo, xp, máximo)
    function g(z, xp)
        a = 0.0
        if z > xp[n+1]
            a = 1.0
        elseif z > xp[1]
            k = findmax(z .≤ xp)[2]
            a = (1/n) * ((k-1) - (xp[k] - z)/(xp[k] - xp[k-1]))
        end
        return a
    end
    m = length(x)
    polivalores = zeros(m)
    for i ∈ 1:m
        polivalores[i] = g(x[i], xp)
    end
    return polivalores 
end

@info "Fn  poligonal"