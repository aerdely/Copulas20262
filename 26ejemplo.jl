### Subcópula de una regresión binaria
### Autor: Arturo Erdely
### Fecha: 2026-05-23

#=
    (X,Y) con X absolutamente continua, Y ~ Bernoulli(θ)

    Y|X=x ~ Bernoulli(θ(x))

    θ(x) = θf(x|Y=1)/f(x) = 1 / (1 + (1-θ)f(x|Y=0)/(θf(x|Y=1))) 

    X|Y=0 ~ Beta(α₀, β₀)
    X|Y=1 ~ Beta(α₁, β₁)

=#

using Plots, LaTeXStrings, Distributions, Optim, SpecialFunctions

include("biseccion.jl")
include("subcopulas.jl")

begin
    α₀ = 4  # 4, 5, 20, 20
    β₀ = 7  # 7, 4, 20, 20
    α₁ = 5  # 5, 4, 0.5, 20
    β₁ = 4  # 4, 7, 0.2, 20
    X₀ = Beta(α₀, β₀) # X|Y=0
    X₁ = Beta(α₁, β₁) # X|Y=1
    k = beta(α₁, β₁) / beta(α₀, β₀)
    regθ(x) = 1 / (1 + k*(1-θ)*(x^(α₀-α₁))*((1-x)^(β₀-β₁))/θ)
    xx = collect(0.0:0.005:1.0)
    plot(xx, pdf.(X₀, xx), label=L"f_{\!X|Y}\,(x\,|\,0)", title="Densidades condicionales",
         xlabel=L"x", ylabel=L"f_{\!X|Y}\,(x\,|\,y)", lw =2
    )
    plot!(xx, pdf.(X₁, xx), label=L"f_{\!X|Y}\,(x\,|\,1)", lw = 2)
end

begin 
    θ = 0.3
    n = 10_000
    Y = rand(Bernoulli(θ), n)
    X = zeros(n)
    for i in 1:n
        if Y[i] == 0
            X[i] = rand(X₀)
        else
            X[i] = rand(X₁)
        end
    end  
    FX(x) = θ*cdf(X₁, x) + (1-θ)*cdf(X₀, x)
    function FXinv(u)
        g(x) = FX(x) - u
        biseccion(g, 0.0, 1.0).raíz
    end
    scatter(X, Y, title="Regresión binaria", 
            xlabel=L"x", ylabel=L"\theta_x = \mathbb{P}(Y=1)",
            legend = false
    )
    plot!(xx, regθ.(xx), lw=3)
end

begin
    W(u) = (u - θ) * (u > θ)
    Π(u) = (1 - θ) * u 
    M(u) = min(u, 1-θ)
    S(u) = (1 - θ) * cdf(X₀, FXinv(u))
    uu = collect(range(0, 1, length=1_000))
    plot(uu, M.(uu), lw = 2, color = :green, label = L"M(u,1-θ)",
         xlabel = L"u", ylabel = L"S_{\!X,Y\,}(u,1-θ)",
         xticks = [0, θ, 1-θ, 1], yticks = [0, 1-θ, 1],
         size = (400,400)
    )
    hline!([1], lw = 0.1, color = :lightgray, label = "")
    plot!(uu, Π.(uu), lw = 2, color = :gray, label = L"\Pi(u,1-θ)")
    plot!(uu, W.(uu), lw = 2, color = :red, label = L"W(u,1-θ)")
    plot!(uu, S.(uu), lw = 3, color = :black, label = L"S_{\!X,Y\,}(u,1-θ)")
end

begin # dependencia teórica
    dsup(u) = Π(u[1]) - S(u[1])
    dinf(u) = S(u[1]) - Π(u[1])
    ds = (-1) * Optim.minimum(optimize(dsup, [0.5]))
    di = Optim.minimum(optimize(dinf, [0.5]))
    d = ds + di
    d = (W = -4*θ*(1-θ), inf = 4*di, S = 4*d, sup = 4*ds, M = 4*θ*(1-θ))
end;

# subcópula empírica
Sn = subcopem(X, Y);
dn = dsub(Sn)
d
plot!(Sn.u, Sn.mat[:, 2], lw = 1.2, color = :violet, label = L"S_n(u,1-θ)")
