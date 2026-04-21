### (pseudo)observaciones de una cópula
### Autor: Arturo Erdely
### Fecha: 2026-04-18

## Paquetes instalados previamente

using Plots, LaTeXStrings

## Scripts requeridos

# varinfo()
include("copulas.jl") # en el mismo directorio que este archivo
# varinfo()
## Cópula paramétrica subfamilia Clayton (Nelsen 4.2.1, θ > 0)

C(u, v, θ) = (u^(-θ) + v^(-θ) - 1)^(-1/θ)

# Simular cópula paramétrica
function simulaCopula(θ, n)
    u = rand(Uniform(0, 1), n)
    ψi(t,u,θ) = (1 - (u^(-θ))*(1 - t^(-θ/(1+θ))))^(-1/θ)
    t = rand(Uniform(0, 1), n)
    v = zeros(n)
    for j ∈ 1:n
        v[j] = ψi(t[j], u[j], θ)
    end
    return u, v
end

## Simular vector aleatorio (X,Y) con cópula Clayton y marginales dadas
begin
    n = 3_000
    θ = 3.5
    U, V = simulaCopula(θ, n)
    X = quantile.(Gamma(2,1), U)
    Y = quantile.(Normal(0,1), V)
    scatter(X, Y, title = "cópula Clayton($θ) + Gamma(2,1) + Normal(0,1)",
            xlabel = "X", ylabel = "Y", label = false, mc = :black,
            ms = 1, size = (600, 600)
    )
end


## Observaciones versus pseudo-observaciones de la cópula

begin
    scatter(U, V, title = "observaciones de la cópula", xlabel = "U", ylabel = "V",
            label = false, mc = :black, ms = 1, size = (600, 600))
end

begin
    pU = poligonal(X, X)
    pV = poligonal(Y, Y)
    scatter(pU, pV, title = "pseudo-observaciones de la cópula", xlabel = "U", ylabel = "V",
            label = false, mc = :black, ms = 1, size = (600, 600))
end


## Ejemplo de evaluación de las cópulas

u, v = 0.358, 0.817  # verificar también en la frontera de [0,1]^2
C(u, v, θ) # cópula teórica
Cn(u, v, [pU pV]) # cópula empírica
CopBer(u, v, [pU pV], 30) # cópula Bernstein
CopCB(u, v, [pU pV], 30) # cópula checkerboard
CopLB(u, v, [pU pV], 30) # cópula Linear B-Spline


## Comparar curvas de nivel de las tres cópulas

# malla para curvas de nivel
begin
    m = 50
    uu = range(0, 1, m+1)
    vv = range(0, 1, m+1)
end;

#  Cópula teórica 
begin
    Copula = [C(u,v,θ) for u in uu, v in vv] 
    nivCopula = contour(uu, vv, Copula, xlabel = L"u", ylabel = L"v",
                        xticks = [0, 0.5, 1], yticks = [0, 0.5, 1], 
                        fill = true, size = (550, 500), title = "Cópula teórica"
    )
end

# Cópula empírica
begin
    Copula = [Cn(u,v, [pU pV]) for u in uu, v in vv] 
    nivCopula = contour(uu, vv, Copula, xlabel = L"u", ylabel = L"v",
                        xticks = [0, 0.5, 1], yticks = [0, 0.5, 1], 
                        fill = true, size = (550, 500), title = "Cópula empírica"
    )
end

# Cópula Bernstein (tarda aprox 30 segundos)
@time begin
    Copula = [CopBer(u,v, [pU pV], 30) for u in uu, v in vv] 
    nivCopula = contour(uu, vv, Copula, xlabel = L"u", ylabel = L"v",
                        xticks = [0, 0.5, 1], yticks = [0, 0.5, 1], 
                        fill = true, size = (550, 500), title = "Cópula Bernstein"
    )
end

# Cópula Linear B-Spline (tarda aprox 30 segundos)
@time begin
    Copula = [CopLB(u,v, [pU pV], 30) for u in uu, v in vv] 
    nivCopula = contour(uu, vv, Copula, xlabel = L"u", ylabel = L"v",
                        xticks = [0, 0.5, 1], yticks = [0, 0.5, 1], 
                        fill = true, size = (550, 500), title = "Cópula Linear B-Spline"
    )
end

# Cópula checkerboard (tarda menos de 1 segundo)
@time begin
    Copula = [CopCB(u,v, [pU pV], 30) for u in uu, v in vv] 
    nivCopula = contour(uu, vv, Copula, xlabel = L"u", ylabel = L"v",
                        xticks = [0, 0.5, 1], yticks = [0, 0.5, 1], 
                        fill = true, size = (550, 500), title = "Cópula checkerboard"
    )
end
