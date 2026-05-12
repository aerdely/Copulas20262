### Medidas de dependencia y concordancia
### Autor: Arturo Erdely
### Fecha: 2026-05-07

## Paquetes instalados previamente

using Plots, LaTeXStrings
## Scripts requeridos

begin
    include("copulas.jl") # en el mismo directorio que este archivo
    include("biseccion.jl") # en el mismo directorio que este archivo
end


#=
    Cópula arquimediana con parámetro θ ≥ 1, definida mediante
    el generador φ(t) = (1-t)^θ, t ∈ [0,1] (familia 4.2.2 de Nelsen)
    
    C(u,v,θ) = max{1 - [(1-u)^θ + (1-v)^θ]^(1/θ), 0}, u,v ∈ [0,1]

    C(u,v,1) = W(u,v) = max{u + v - 1, 0} (cota inferior Fréchet-Hoeffding)
    C(u,v,∞) = M(u,v) = min{u, v} (cota superior Fréchet-Hoeffding)

    C(u,v,θ) ≠ Π(u,v)=uv para todo valor admisible θ ≥ 1
=#

C(u,v,θ) = max(1 - ((1-u)^θ + (1-v)^θ)^(1/θ), 0)

## Simular cópula paramétrica

# curva frontera de región cero 
g(u,θ) = 1 - (1 - (1-u)^θ)^(1/θ) 

# ψ⁻¹(t|u) como función de t, para u y θ fijos
function ψi(t, u, θ)
    if t > (1-u)^(θ-1) 
        return 1 - (1-u)*(1/(t^(θ/(θ-1))) - 1)^(1/θ)
    else
        return g(u, θ)
    end
end

# Función para simular n observaciones de la cópula arquimediana con parámetro θ
function simular_copula(n, θ)
    U = rand(n)
    V = similar(U)
    for i ∈ 1:n
        u = U[i]
        w = rand()
        V[i] = ψi(w, u, θ)
    end
    return U, V
end

begin
    θ = 1.5
    U, V = simular_copula(10_000, θ)
    scatter(U, V, xlabel = L"U", ylabel = L"V", size = (500,470), 
            ms = 1.0, mcolor = :black, label = "")
end
# agregar curva frontera de región cero
plot!(sort(U), g.(sort(U), θ), color = :red, lw = 1.0, label = L"g_{θ}(u)")

# dependencia y concordancia teóricas y empíricas
σ = Schweizer((u,v) -> C(u,v,θ)) # teórica
ρ = Spearman((u,v) -> C(u,v,θ)) # teórica
σn = Schweizer([U V], 100) # empírica
ρn = Spearman([U V], 100) # empírica

# secciones diagonales
begin
    tt = collect(range(0, 1, length = 100))
    δn = [Cn(t,t,[U V]) for t ∈ tt]
    λn = [Cn(t,1-t,[U V]) for t ∈ tt]
    δM = tt 
    δΠ = tt .^ 2
    δW = @. (2*tt - 1)*(tt > 1/2) # = máx{2u-1, 0}
    λM = min.(tt, 1 .- tt) # o bien: λM = @. min(tt, 1 - tt)
    λΠ = @. tt*(1 - tt)
    λW = zeros(length(tt))
end;
begin # diagonal principal
    plot(tt, δM, color = :green, lw = 2, label = "M", title = "diagonal principal", size = (500,500))
    xaxis!(L"t"); yaxis!(L"δ(t) = C(t,t)")
    plot!(tt, δΠ, color = :gray, lw = 2, label = "Π")
    plot!(tt, δW, color = :red, lw = 2, label = "W")
    plot!(tt, δn, color = :black, lw = 2, label = L"δ_n")
end
begin # diagonal secundaria
    plot(tt, λM, color = :green, lw = 2, label = "M", title = "diagonal secundaria", size = (500,500))
    xaxis!(L"t"); yaxis!(L"λ(t) = C(t,1-t)")
    plot!(tt, λΠ, color = :gray, lw = 2, label = "Π")
    plot!(tt, λW, color = :red, lw = 2, label = "W")
    plot!(tt, λn, color = :black, lw = 2, label = L"λ_n")
end


## Medidas de dependencia y concordancia en función de θ
begin
    θs = collect(1:20)
    σs = [Schweizer((u,v) -> C(u,v,θ))[1] for θ in θs]
    ρs = [Spearman((u,v) -> C(u,v,θ))[1] for θ in θs]
    plot(θs, σs, label = "Schweizer-Wolff", xlabel = L"θ", ylabel = "medida", lw = 2, size = (500,470))
    plot!(θs, ρs, label = "Spearman", lw = 2)
end

begin
    θs = collect(1:0.2:5)
    σs = [Schweizer((u,v) -> C(u,v,θ))[1] for θ in θs]
    ρs = [Spearman((u,v) -> C(u,v,θ))[1] for θ in θs]
    plot(θs, σs, label = "Schweizer-Wolff", xlabel = L"θ", ylabel = "medida", lw = 2, size = (500,470))
    plot!(θs, ρs, label = "Spearman", lw = 2)
end
plot!(θs, abs.(ρs), label = "abs(Spearman)", lw = 1, ls = :dash)

# parámetro correspondiente a una concordancia de Spearman dada
function θdeRho(C::Function, ρ::Float64, intervalo::Tuple{Float64, Float64})
    f(θ) = Spearman((u,v) -> C(u,v,θ))[1] - ρ
    return biseccion(f, intervalo[1], intervalo[2], δ = 0.0001, m = 1000)
end

θestim = θdeRho(C, ρn, (1.0, 10.0))
