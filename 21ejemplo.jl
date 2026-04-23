### Estimación paramétrica de una cópula
### Autor: Arturo Erdely
### Fecha: 2026-04-22

## Paquetes instalados previamente

using Plots, LaTeXStrings
## Scripts requeridos

include("copulas.jl") # en el mismo directorio que este archivo


## Cópula paramétrica subfamilia Clayton (Nelsen 4.2.1, θ > 0)

begin
    C(u, v, θ) = (u^(-θ) + v^(-θ) - 1)^(-1/θ)
    c(u, v, θ) = (1+θ) * (u^(-θ) + v^(-θ) - 1)^(-1/θ-2) * (u^(-θ-1)) * (v^(-θ-1))
    logc(u, v, θ) = log(1+θ) - (1+θ)*(log(u) + log(v)) - (1/θ + 2) * log(u^(-θ) + v^(-θ)-1)
end;

# función de log-verosimilitud
function LogVer(θ::Float64, cobs::Matrix{Float64})
    valor = 0.0
    for i ∈ 1:n 
        valor += logc(cobs[i, 1], cobs[i, 2], θ)
    end
    return valor 
end

# función a minimizar (método de mínimos cuadrados ordinarios)
function MCO(θ::Float64, cobs::Matrix{Float64})
    valor = 0.0
    for i ∈ 1:n
        valor += (C(cobs[i, 1], cobs[i, 2], θ) - Cn(cobs[i, 1], cobs[i, 2], cobs))^2
    end
    return valor
end

# Approximate Bayesian Computation (ABC)
function métodoABC(θmin, θmax, m, k, cobs::Matrix{Float64})
    n = size(cobs, 1)
    CnObs = zeros(n)
    for i ∈ 1:n
        CnObs[i] = Cn(cobs[i, 1], cobs[i, 2], cobs)
    end
    θsim = rand(Uniform(θmin, θmax), m)
    ε = zeros(m)
    for j ∈ 1:m
        csim = simulaCopula(θsim[j], n)
        for i ∈ 1:n
            ε[j] += (CnObs[i] - Cn(cobs[i, 1], cobs[i, 2], csim))^2
        end
    end
    jopt = sortperm(ε)[1:k]
    return θsim[jopt], ε[jopt]
end


## Simular cópula paramétrica

function simulaCopula(θ, n)
    u = rand(Uniform(0, 1), n)
    ψi(t,u,θ) = (1 - (u^(-θ))*(1 - t^(-θ/(1+θ))))^(-1/θ)
    t = rand(Uniform(0, 1), n)
    v = zeros(n)
    for j ∈ 1:n
        v[j] = ψi(t[j], u[j], θ)
    end
    return [u v] # matriz de dos columnas
end

begin
    n = 3_000
    θ = 2.5
    cobs = simulaCopula(θ, n)
end


## Estimaciones paramétricas de la cópula

# Por máxima verosimilitud
@time begin
    th = range(0.01, 5.0, length = 1_000)
    Lth = [LogVer(θ, cobs) for θ in th]
    k = findmax(Lth)[2]
    plot(th, Lth, xlabel = L"θ", ylabel = L"\log\,L(θ\,\,|D)",
         label = "", size = (500, 500), lw = 3
    )
    θmv = th[k]
    vline!([θ], color = :black, label = "θ = $(round(θ, digits = 3))")
    vline!([θmv], color = :red, label = "θ* = $(round(θmv, digits = 3))")
end

# Por mínimos cuadrados ordinarios (30 segundos aprox.)
@time begin
    th = range(0.01, 5.0, length = 1_000)
    hth = [MCO(θ, cobs) for θ in th]
    k = findmin(hth)[2]
    plot(th, hth, xlabel = L"θ", ylabel = L"MCO(θ\,\,|D)",
         label = "", size = (500, 500), lw = 3
    )
    θh = th[k]
    vline!([θ], color = :black, label = "θ = $(round(θ, digits = 3))")
    vline!([θh], color = :red, label = "θ* = $(round(θh, digits = 3))")
end

# Por método ABC (30 segundos aprox.)
@time begin
    θmin = 0.01
    θmax = 5.0
    m = 1_000
    k = 50
    θopt, εopt = métodoABC(θmin, θmax, m, k, cobs)
end;

extrema(εopt)
mean(θopt)
quantile(θopt, [0.025, 0.5, 0.975])

begin
    histogram(θopt, normalize = true, xlabel = L"θ", ylabel = "Densidad a posteriori",
             label = "", size = (500, 500), color = :lightyellow
    )
    vline!([θ], color = :violet, lw = 2, label = "θ = $(round(θ, digits = 3))")
    vline!([quantile(θopt, 0.975)], lw = 3, color = :blue, label = "Cuantil 97.5% = $(round(quantile(θopt, 0.975), digits = 3))")   
    vline!([median(θopt)], lw = 3.5, color = :red, label = "Mediana = $(round(median(θopt), digits = 3))")
    vline!([quantile(θopt, 0.025)], lw = 3, color = :green, label = "Cuantil 2.5% = $(round(quantile(θopt, 0.025), digits = 3))")
end
