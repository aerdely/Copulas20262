### Regresión vía cópulas
### Autor: Arturo Erdely
### Fecha: 2026-05-10

#=
    Cópula paramétrica subfamilia Clayton (Nelsen 4.2.1, θ > 0)
=#

## Paquetes instalados previamente

using Plots, LaTeXStrings, Distributions


# Simular cópula paramétrica
ψi(t,u,θ) = (1 - (u^(-θ))*(1 - t^(-θ/(1+θ))))^(-1/θ) # ψ⁻¹(t|u) para la cópula Clayton
function simulaCopula(θ, n)
    u = rand(Uniform(0, 1), n)
    t = rand(Uniform(0, 1), n)
    v = zeros(n)
    for j ∈ 1:n
        v[j] = ψi(t[j], u[j], θ)
    end
    return u, v
end

begin
    θ = 3
    U, V = simulaCopula(θ, 1000)
    scatter(U, V, xlabel = L"U", ylabel = L"V", size = (500,470), 
            ms = 2.0, mcolor = :blue, label = "")
end


## Simular (X,Y) con cópula arquimediana y marginales dadas 

function simular_XY(n, θ, distX, distY)
    U, V = simulaCopula(θ, n)
    X = quantile(distX, U)
    Y = quantile(distY, V)
    return X, Y
end

begin
    θ = 5
    distX = Beta(0.8,0.5)
    distY = Gamma(2,3)
    X, Y = simular_XY(1000, θ, distX, distY)
    scatter(X, Y, xlabel = L"X", ylabel = L"Y", size = (500,470), 
            ms = 2.0, mcolor = :green, label = "")
end

# comprobando marginales

histogram(X, bins = 30, xlabel = L"X", ylabel = "densidad",
          title = "Histograma de X ~ Beta(0.8,0.5)", size = (500,400), label = "empírica",
          mcolor = :orange, linecolor = :black, normalize = true)
plot!(x -> pdf(distX, x), 0:0.01:1, color = :red, lw = 2.0, label = "teórica")

histogram(Y, bins = 30, xlabel = L"Y", ylabel = "denidad",
          title = "Histograma de Y ~ Gamma(2,3)", size = (500,400), label = "empírica",
          mcolor = :purple, linecolor = :black, normalize = true)
plot!(x -> pdf(distY, x), 0:0.1:20, color = :red, lw = 2.0, label = "teórica")


## Curvas de regresión
qsup = 0.975 # cuantil superior
qinf = 0.025 # cuantil inferior

# Para las (pseudo) observaciones de la copula
scatter(U, V, xlabel = L"U", ylabel = L"V", size = (500,470), 
            ms = 2.0, mcolor = :blue, label = ""
)
uu = collect(range(0, 1, length = 1_000));
plot!(uu, [ψi(qsup, u, θ) for u ∈ uu], color = :green, lw = 2.0, label = "cuantil $qsup")
plot!(uu, [ψi(0.5, u, θ) for u ∈ uu], color = :red, lw = 2.0, label = "mediana")
plot!(uu, [ψi(qinf, u, θ) for u ∈ uu], color = :orange, lw = 2.0, label = "cuantil $qinf")


# Para las variables originales
scatter(X, Y, xlabel = L"X", ylabel = L"Y", size = (500,470), 
            ms = 2.0, mcolor = :blue, label = ""
)
xx = collect(range(minimum(X), maximum(X), length = 1_000));
plot!(xx, [quantile(distY, ψi(qsup, cdf(distX, x), θ)) for x ∈ xx], 
      color = :green, lw = 2.0, label = "cuantil $qsup")
plot!(xx, [quantile(distY, ψi(0.5, cdf(distX, x), θ)) for x ∈ xx], 
      color = :red, lw = 2.0, label = "mediana")    
plot!(xx, [quantile(distY, ψi(qinf, cdf(distX, x), θ)) for x ∈ xx], 
      color = :orange, lw = 2.0, label = "cuantil $qinf")
