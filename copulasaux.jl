### Funciones auxiliares para `copulas.jl`
### Autor: Arturo Erdely
### Última actualización: 2026-05-04


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


# repetidos (1 method)
"""
    repetidos(x::Vector)

Returns a type `Int` vector with the positions in vector `x` of repeated values.

## Examples
```
repetidos([1, 2, 3, 4, 1, 2, 5, 2])

repetidos(['A', 'B', 'C', 'A', 'B', 'A'])

repetidos(rand(10_000)) # no repeated values
```
"""
function repetidos(x::Vector)
    if allunique(x)
        return Int[]
    else
        u = unique(x)
        n = length(u)
        xpos = Int[]
        for k ∈ 1:n
            ix = findall(x .== u[k])
            if length(ix) > 1
                append!(xpos, ix[2:end])
            end
        end
        sort!(xpos)
    end
    return xpos 
end


# resumen (4 methods)
"""
    resumen(v::Vector{<:Real}, desplegar = false)

Calculates a tuple of summary statistics for the values in vector `v`.
If `desplegar = true` it also displays a summary table. The calculated
statistics are: sample size, total number of repeated values, mean,
minimum, first quartile, median, third quartile and maximum.

*Dependencies*:
- function `repetidos` 
- package `Statistics` (Julia standard library)

## Example
```
v = rand(rand(50), 55);
r = resumen(v, true)
r.median
```
"""
function resumen(v::Vector{<:Real}, desplegar = false)
    # dependencies: `repetidos` `Statistics`
    n = length(v)
    rep = length(repetidos(v))
    μ = mean(v)
    m = minimum(v)
    Q1 = quantile(v, 0.25)
    Q2 = median(v)
    Q3 = quantile(v, 0.75)
    M = maximum(v)
    r = (samplesize = n, repeated = rep, mean = μ, min = m,
         quartile_1 = Q1, median = Q2, quartile_3 = Q3, max = M)
    if desplegar
        display(hcat(collect(keys(r)), collect(values(r))))
    end
    return r
end


"""
    resumen(matriz::Matrix{<:Real}, desplegar = false)

Calculates a tuple of summary statistics for each column in `matriz`.
If `desplegar = true` it also displays a summary table. The calculated
statistics are: sample size, total number of repeated values, mean,
minimum, first quartile, median, third quartile and maximum.

*Dependencies*:
- function `resumen(v::Vector{<:Real},...)`
- package `Statistics` (Julia standard library)

## Example
```
M = [[.1,.3,.2,.1,.2,.4] [.6,.4,.2,.5,.3,.1] [.4,.5,.4,.6,.7,.8]]
r = resumen(M, true)
r.repeated
r.table
```
"""
function resumen(matriz::Matrix{<:Real}, desplegar = false)
    # dependencies: `Statistics` `resumen(v::Vector{<:Real},...)`
    nfilas, ncols = size(matriz)
    v = collect(1:ncols)
    n = fill(nfilas, ncols); rep = fill(0, ncols)
    μ = fill(0.0, ncols); m = fill(0.0, ncols)
    Q1 = fill(0.0, ncols); Q2 = fill(0.0, ncols)
    Q3 = fill(0.0, ncols); M = fill(0.0, ncols)
    D = Matrix{Any}(undef, nfilas, 0)
    for j ∈ 1:ncols
        r = resumen(matriz[:, j])
        n[j] = r.samplesize; rep[j] = r.repeated
        μ[j] = r.mean; m[j] = r.min
        Q1[j] = r.quartile_1; Q2[j] = r.median
        Q3[j] = r.quartile_3; M[j] = r.max
        if j == 1
            D = vcat([:variable j], hcat(collect(keys(r)), collect(values(r))))
        else
            D = hcat(D, pushfirst!(collect(values(r)), j))
        end
    end
    rr = (variable = v, samplesize = n, repeated = rep, mean = μ, min = m,
         quartile_1 = Q1, median = Q2, quartile_3 = Q3, max = M, table = D)
    if desplegar
        display(D)
    end
    return rr
end


# desempate! (4 methods)   desempate (4 methods)
"""
    desempate!(v::Vector{<:AbstractFloat}, ε = minimum(abs.(v)) / 1_000_000.0)

Adds (in place) a uniform random noise ±`ε` to the repeated values in `v` (if any). 
If `ε` is not specified then it defaults to `ε = minimum(abs.(v)) / 1_000_000`.

*Dependencies*: 
- function `repetidos`

## Example
```
a = [.1, .2, .1, .1, .3, .2]
desempate!(a)
println(a)
```
"""
function desempate!(v::Vector{<:AbstractFloat}, ε = minimum(abs.(v)) / 1_000_000.0)
    # dependencies: `repetidos`
    ipos = repetidos(v)
    npos = length(ipos)
    v[ipos] .+= ε .* (2 .* rand(npos) .- 1)
    return v
end


"""
    desempate!(matriz::Matrix{<:AbstractFloat}, ε = minimum(abs.(matriz)) / 1_000_000.0)

Adds (in place) a uniform random noise ±`ε` to the repeated values in each column of `matriz`
(if any). If `ε` is not specified then it defaults to `ε = minimum(abs.(matriz)) / 1_000_000`.

*Dependencies*: 
- function `repetidos`
- function `desempate!(v::Vector{<:AbstractFloat},...)`

## Example
```
A = [[.1,.2,.1,.3] [.1,.2,.3,.4] [.1,.1,.2,.1]]
println(A)
desempate!(A)
println(A)
```
"""
function desempate!(matriz::Matrix{<:AbstractFloat}, ε = minimum(abs.(matriz)) / 1_000_000.0)
    # dependencies: `repetidos` `desempate!(v::Vector{Float64},...)`
    ncol = size(matriz)[2]
    for j ∈ 1:ncol
        matriz[:, j] = desempate!(matriz[:, j], ε)
    end
    return(matriz)
end


"""
    desempate(v::Vector{<:AbstractFloat}, ε = minimum(abs.(v)) / 1_000_000.0)

Adds a uniform random noise ±`ε` to the repeated values in a copy of `v` (if any).
If `ε` is not specified then it defaults to `ε = minimum(abs.(v)) / 1_000_000`.
Returns a tuple with such copy but with no repeated values, a vector with the 
positions where repeated values were found, and the total number of repeated values.

*Dependencies*: 
- function `repetidos`

## Example
```
a = [.1, .2, .1, .1, .3, .2]
b = desempate(a)
a ≠ b[1] 
```
"""
function desempate(v::Vector{<:AbstractFloat}, ε = minimum(abs.(v)) / 1_000_000.0)
    # dependencies: `repetidos`
    vv = copy(v)
    ipos = repetidos(vv)
    npos = length(ipos)
    vv[ipos] .+= ε .* (2 .* rand(npos) .- 1)
    return (vv, ipos, npos)
end

"""
    desempate(matriz::Matrix{<:AbstractFloat}, ε = minimum(abs.(matriz)) / 1_000_000.0)

Adds a uniform random noise ±`ε` to the repeated values in each column of
a copy of `matriz` (if any). If `ε` is not specified then it defaults to
`ε = minimum(abs.(v)) / 1_000_000`. Returns a tuple with such copy but with no repeated values,
a vector with the positions where repeated values were found in each column,
and the total number of repeated values in each column.

*Dependencies*: 
- function `repetidos`
- function `desempate(v::Vector{<:AbstractFloat},...)`
        
## Example
```
A = [[.1,.2,.1,.3] [.1,.2,.3,.4] [.1,.1,.2,.1]]
R = desempate(A)
R[1]
A .≠ R[1]
R[2]
R[3]
```
"""
function desempate(matriz::Matrix{<:AbstractFloat}, ε = minimum(abs.(matriz)) / 1_000_000.0)
    # dependencies: `repetidos` `desempate(v::Vector{<:AbstractFloat},...)`
    d = size(matriz)
    M = zeros(d)
    ipos = Vector{Vector{Int}}(undef, d[2])
    npos = Vector{Int}(undef, d[2])
    for j ∈ 1:d[2]
        des = desempate(matriz[:, j], ε)
        M[:, j] = des[1]
        ipos[j] = des[2]
        npos[j] = des[3]
    end
    return (M, ipos, npos)
end


@info "Fn  poligonal  repetidos  resumen  desempate  desempate!"
