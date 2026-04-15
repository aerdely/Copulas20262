### Aproximación poligonal univariada 
### Autor: Arturo Erdely
### Fecha: 2026-04-03
## Paquetes instalados previamente

using Plots, LaTeXStrings, Distributions

## Función de distribución empírica univariada

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

methods(Fn) 

@doc Fn 


## Poligonal
"""
    poligonal(x::Vector{<:Real}, xobs::Vector{<:Real}; mínimo = minimum(xobs), máximo = maximum(xobs))

Aproximación poligonal de una función de distribución continua en un vector dado de valores `x`
con base en una muestra aleatoria observada dada por el vector `xobs`. El parámetro opcional
`mínimo` puede escogerse como un valor menor que el mínimo de la muestra, y `máximo` como
un valor mayor que el máximo de la muestra, si así se requiere.

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

@doc poligonal


## Ejemplo: m.a ~ Gamma(2, 3)

# Simular muestra y calcular distribución teórica y empírica
begin
    n = 1000 # probar n ∈ {5, 100, 1_000}
    Xmodel = Gamma(2.0, 3.0)
    xobs = rand(Xmodel, n)
    x = collect(range(0.0, quantile(Xmodel, 0.995), length = 1_000))
    yFX = cdf.(Xmodel, x) # Distribución teórica
    yFn = Fn(x, xobs) # Distribución empírica
end;

# Aproximación poligonal + gráfica comparativa
begin
    yPoligonal = poligonal(x, xobs) # probar parámetros opcionales: mínimo = 0.0, máximo = ...
    plot(x, yFX, label = "Teórica", legend = :right, lw = 4.0)
    xaxis!("x"); yaxis!("F(x)"); title!("Función de distribución")
    scatter!(x, yFn, label = "Empírica", lw = 3.0, color = :black)
    plot!(x, yPoligonal, label = "Poligonal", lw = 2.0, color = :red)
end

