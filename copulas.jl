### Biblioteca de funciones para cópulas 
### Autor: Arturo Erdely
### Última actualización: 2026-04-03

using Distributions # instalado previamente

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
El parámetro `m` es el orden de la cópula de Bernstein (igual al tamaño de la muestra por default).

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

@info "Cn  CopBer"