### Biblioteca de funciones para subcópulas 
### Autor: Arturo Erdely
### Última actualización: 2026-05-20

include("02probestim.jl")


"""
    subcopem(xobs::Vector{<:Real}, yobs::Vector{<:Real})

Calcula la subcópula empírica bivariada a partir de pares de observaciones `(x,y)` proporcionados
por los vectores `xobs` y `yobs`. Devuelve una tupla etiquetada con los siguientes elementos:
- `mat` = Matriz de valores de la subcópula empírica.
- `u` = Vector de valores de la función de masa acumulada para `xobs`.
- `v` = Vector de valores de la función de masa acumulada para `yobs`.

> Dependencias: función `masaprob` del archivo `02probestim.jl`

## Ejemplo:
```
xobs = [3, 3, 1, 4, 1, 5, 1, 4, 6, 2, 4, 4];
yobs = [2, 1, 1, 0, 1, 2, 0, 0, 0, 1, 0, 0];
Sn = subcopem(xobs, yobs);
Sn.mat
Sn.u
Sn.v 
```
"""
function subcopem(xobs::Vector{<:Real}, yobs::Vector{<:Real})
    n = length(xobs)
    if n ≠ length(yobs)
        error("Ambos vectores deben ser de igual longitud")
        return nothing
    end
    X = masaprob(xobs)
    Y = masaprob(yobs)
    r = X.valores
    m1 = length(r)
    s = Y.valores
    m2 = length(s)
    px = X.probs
    py = Y.probs
    qx = pushfirst!(cumsum(px), 0.0)
    qy = pushfirst!(cumsum(py), 0.0)
    Sn = zeros(m1+1, m2+1)
    for i ∈ 2:(m1+1), j ∈ 2:(m2+1)
        Sn[i,j] = sum((xobs .≤ r[i-1]) .* (yobs .≤ s[j-1])) / n 
    end
    return (mat = Sn, u = qx, v = qy)
end


"""
    dsub(subcop)

Calcula la medida de dependencia para subcópulas propuesta por Erdely (2017)
```math
d(S) = 4 * ( sup{S - Π} - sup{Π - S} )
```
Toma como entrada `subcop` una subcópula empírica representada por la tupla
devuelta por `subcopem` y devuelve una tupla etiquetada con los siguientes elementos:
- `W` = Mínima dependencia inferior.
- `inf` = Extremo de dependencia inferior.
- `d` = Medida de dependencia total.
- `sup` = Extremo de dependencia superior.
- `M` = Máxima dependencia superior.

## Ejemplo:
```
xobs = [3, 3, 1, 4, 1, 5, 1, 4, 6, 4, 4];
yobs = [2, 1, 1, 0, 1, 2, 0, 0, 1, 0, 0];
Sn = subcopem(xobs, yobs); 
dsub(Sn)
```
"""
function dsub(subcop)
    Sn = subcop.mat
    uu = subcop.u
    nu = length(uu)
    vv = subcop.v
    nv = length(vv)
    Π = subcop.u * transpose(subcop.v)
    W = similar(Π)
    M = similar(Π)
    for i ∈ 1:nu, j ∈ 1:nv
        W[i,j] = max(uu[i] + vv[j] - 1, 0)
        M[i,j] = min(uu[i], vv[j])
    end
    dW = minimum(W - Π)
    dM = maximum(M - Π)
    D = Sn - Π
    dsup = maximum(D)
    dinf = minimum(D)
    dS = dsup + dinf
    return (W = 4*dW, inf = 4*dinf, d = 4*(dsup + dinf), sup = 4*dsup, M =4*dM)
end

@info "subcopem  dsub"
