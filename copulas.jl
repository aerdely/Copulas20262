### Biblioteca de funciones para cópulas 
### Autor: Arturo Erdely
### Última actualización: 2026-04-19

using Distributions, HCubature # instalados previamente

include("copulasaux.jl") # en el mismo directorio que este archivo


"""
    Cn(u::Real, v::Real, cobs::Matrix{<:Real})

Calcula la cópula empírica bivariada en el punto `(u,v)` con base en las observaciones
(o pseudo-observaciones) de una cópula subyacente, dadas por la matriz `cobs`
de dos columnas, donde cada fila de `cobs` representa un par de observaciones.

> Advertencia: no debe haber valores repetidos en una misma columna de `cobs`.
"""
function Cn(u::Real, v::Real, cobs::Matrix{<:Real})
    n = size(cobs, 1)
    return (1/n) * sum( (cobs[:,1] .≤ u) .* (cobs[:,2] .≤ v) )
end


"""
    CopBer(u::Real, v::Real, cobs::Matrix{<:Real}, m::Int = size(cobs, 1))

Calcula la cópula de Bernstein bivariada en el punto `(u,v)` con base en las observaciones
(o pseudo-observaciones) de una cópula subyacente, dadas por la matriz `cobs`
de dos columnas, donde cada fila de `cobs` representa un par de observaciones. 
El parámetro `m` es el orden de la cópula de Bernstein (igual al tamaño de la muestra, por default).

> Advertencia: no debe haber valores repetidos en una misma columna de `cobs`.

Utiliza: 
- función `Cn`
- paquete `Distributions`
"""
function CopBer(u::Real, v::Real, cobs::Matrix{<:Real}, m::Int = size(cobs, 1))
    valor = 0.0
    uu = range(0, 1, m+1)
    vv = range(0, 1, m+1)
    copemp = [Cn(u,v,cobs) for u in uu, v in vv]
    for i ∈ 1:m+1
        for j ∈ 1:m+1
            valor += copemp[i,j] * pdf(Binomial(m, u), i-1) * pdf(Binomial(m, v), j-1)
        end
    end
    return valor
end


"""
    CopCB(u::Real, v::Real, cobs::Matrix{<:Real}, m::Int = size(cobs, 1))

Calcula la cópula `checkerboard` bivariada en el punto `(u,v)` con base en las
observaciones (o pseudo-observaciones) de una cópula subyacente, dadas por la matriz `cobs`
de dos columnas, donde cada fila de `cobs` representa un par de observaciones. 
El parámetro `m` define la rejilla de valores {0,1/m,...,(m-1)/m,m}^2 en donde se
evalúa la cópula empírica, para luego extenderla sobre todo [0,1]^2
(igual al tamaño de la muestra, por default).

> Advertencia: no debe haber valores repetidos en una misma columna de `cobs`.

Utiliza: 
- función `Cn`
"""
function CopCB(u::Real, v::Real, cobs::Matrix{<:Real}, m::Int = size(cobs, 1))
    mu = m*u 
    fmu = floor(mu)
    mv = m*v 
    fmv = floor(mv)
    uinf = fmu / m
    usup = uinf + 1/m 
    vinf = fmv / m 
    vsup = vinf + 1/m 
    λu = mu - fmu
    λv = mv - fmv
    Cuv = (1-λu)*(1-λv)*Cn(uinf,vinf,cobs) + (1-λu)*λv*Cn(uinf,vsup,cobs) + λu*(1-λv)*Cn(usup,vinf,cobs) + λu*λv*Cn(usup,vsup,cobs)
    return Cuv 
end


"""
    CopLB(u::Real, v::Real, cobs::Matrix{<:Real}, m::Int = size(cobs, 1))

Calcula la cópula `Linear B-Spline` bivariada en el punto `(u,v)` con base en las
observaciones (o pseudo-observaciones) de una cópula subyacente, dadas por la matriz `cobs`
de dos columnas, donde cada fila de `cobs` representa un par de observaciones. 
El parámetro `m` define la rejilla de valores {0,1/m,...,(m-1)/m,m}^2 en donde se
evalúa la cópula empírica, para luego extenderla sobre todo [0,1]^2
(igual al tamaño de la muestra, por default).

> Advertencia: no debe haber valores repetidos en una misma columna de `cobs`.

Utiliza: 
- función `Cn`
"""
function CopLB(u::Real, v::Real, cobs::Matrix{<:Real}, m::Int = size(cobs, 1))
    valor = 0.0
    uu = range(0, 1, m+1)
    vv = range(0, 1, m+1)
    copemp = [Cn(u,v,cobs) for u in uu, v in vv]
    function N(i, u)
        if i == m
            if (m-1)/m ≤ u ≤ 1 
                return m*u - m + 1
            else
                return 0.0
            end
        else
            if (i-1)/m ≤ u < i/m 
                return m*u - i + 1
            elseif i/m ≤ u < (i+1)/m 
                return i + 1 - m*u 
            else
                return 0.0 
            end
        end
    end
    for i ∈ 1:m+1
        for j ∈ 1:m+1
            valor += copemp[i,j] * N(i-1, u) * N(j-1, v)
        end
    end
    return valor
end


"""
    Schweizer(C::Function)

Calcula la medida de dependencia de Schweizer-Wolff para la cópula dada por la función `C(u,v)`.

Utiliza:
- paquete `HCubature`
"""
function Schweizer(C::Function)
    g(x) = abs(C(x[1], x[2]) - x[1]*x[2])
    integral, error = hcubature(g, [0.0, 0.0], [1.0, 1.0])
    return (sw = 12 * integral, error = error)
end


"""
    Schweizer(cobs::Matrix{<:Real}, m::Int = size(cobs, 1))

Calcula la medida de dependencia de Schweizer-Wolff empírica para las observaciones
(o pseudo-observaciones) de una cópula subyacente, dadas por la matriz `cobs` de dos columnas, 
donde cada fila de `cobs` representa un par de observaciones. El parámetro `m` define la rejilla
de valores {0,1/m,...,(m-1)/m,m}^2 en donde se evalúa la cópula empírica (igual al tamaño de la
muestra, por default).

Utiliza:
- función `Cn`
"""
function Schweizer(cobs::Matrix{<:Real}, m::Int = size(cobs, 1))
    uu = range(1/m, 1, m)
    vv = range(1/m, 1, m)
    s = 0.0
    for i ∈ 1:m
        for j ∈ 1:m
            s += abs(Cn(uu[i], vv[j], cobs) - uu[i]*vv[j])
        end
    end
    return 12 * s / m^2
end


"""
    Spearman(C::Function)

Calcula la medida de concordancia de Spearman para la cópula dada por la función `C(u,v)`

Utiliza:
- paquete `HCubature`
"""
function Spearman(C::Function)
    g(x) = C(x[1], x[2]) - x[1]*x[2]
    integral, error = hcubature(g, [0.0, 0.0], [1.0, 1.0])
    return (sp = 12 * integral, error = error)
end


"""
    Spearman(cobs::Matrix{<:Real}, m::Int = size(cobs, 1))

Calcula la medida de concordancia de Spearman empírica para las observaciones (o pseudo-observaciones)
de una cópula subyacente, dadas por la matriz `cobs` de dos columnas, donde cada fila de `cobs`
representa un par de observaciones. El parámetro `m` define la rejilla de valores
{0,1/m,...,(m-1)/m,m}^2 en donde se evalúa la cópula empírica (igual al tamaño de la muestra, por default).

Utiliza:
- función `Cn`
"""
function Spearman(cobs::Matrix{<:Real}, m::Int = size(cobs, 1))
    uu = range(1/m, 1, m)
    vv = range(1/m, 1, m)
    s = 0.0
    for i ∈ 1:m
        for j ∈ 1:m
            s += Cn(uu[i], vv[j], cobs) - uu[i]*vv[j]
        end
    end
    return 12 * s / m^2
end


@info "Cn  CopBer  CopCB  CopLB  Schweizer  Spearman"
