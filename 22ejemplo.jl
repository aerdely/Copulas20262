### Medidas de dependencia y concordancia
### Autor: Arturo Erdely
### Fecha: 2026-05-04

## Paquetes instalados previamente

using Plots, LaTeXStrings
## Scripts requeridos

begin
    include("copulas.jl") # en el mismo directorio que este archivo
    include("biseccion.jl") # en el mismo directorio que este archivo
end


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


## Medidas de dependencia y concordancia teóricas 

@doc Schweizer 

@doc Spearman

begin
    θ = 3.5
    println("Cópula Clayton(θ = $θ)")
    σ = Schweizer((u,v) -> C(u,v,θ))
    println("Schweizer-Wolff = ", σ[1])
    ρ = Spearman((u,v) -> C(u,v,θ))
    println("Spearman = ", ρ[1])
end

# parámetro correspondiente a una concordancia de Spearman dada
function θdeRho(C::Function, ρ::Float64, intervalo::Tuple{Float64, Float64})
    f(θ) = Spearman((u,v) -> C(u,v,θ))[1] - ρ
    return biseccion(f, intervalo[1], intervalo[2])
end

p = θdeRho(C, 0.8203, (0.01, 10.0))

p = θdeRho(C, 0.5, (0.01, 10.0))
Spearman((u,v) -> C(u,v,p.raíz))
Schweizer((u,v) -> C(u,v,p.raíz))


## Dependencia y concordancia empíricas

begin
    n = 3_000
    θ = 3.5
    U, V = simulaCopula(θ, n)
    scatter(U, V, title = "cópula Clayton(θ = $θ)",
            xlabel = "U", ylabel = "V", label = false, mc = :black,
            ms = 1, size = (600, 600)
    )
end

σ = Schweizer((u,v) -> C(u,v,θ)) # valor teórico
ρ = Spearman((u,v) -> C(u,v,θ)) # valor teórico

CopEst(u,v) = CopCB(u,v, hcat(U,V), 30)
@time Schweizer((u,v) -> CopEst(u,v)) # (46 segundos aprox)
@time Spearman((u,v) -> CopEst(u,v)) # (46 segundos aprox)

@time Schweizer(hcat(U,V), 100) # (menos de 1 segundo)

begin
    n = 10_000
    θ = 3.5
    U, V = simulaCopula(θ, n)
    scatter(U, V, title = "cópula Clayton(θ = $θ)",
            xlabel = "U", ylabel = "V", label = false, mc = :black,
            ms = 1, size = (600, 600)
    )
end

@time Schweizer(hcat(U,V), 100) 
@time Spearman(hcat(U,V), 100)


## Clayton rotada 90°

begin
    scatter(U, 1 .- V, title = "cópula Clayton(θ = $θ) rotada 90°",
            xlabel = "U", ylabel = "V", label = false, mc = :black,
            ms = 1, size = (600, 600)
    )
end

coprot = hcat(U, 1 .- V);
@time Schweizer(coprot, 100) 
@time Spearman(coprot, 100) 
