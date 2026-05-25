### Regresión binaria
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

begin # Densidades condicionales X|Y=y, y ∈ {0,1}
    α₀ = 4  # 4, 5, 20
    β₀ = 7  # 7, 4, 20
    α₁ = 5  # 5, 4, 0.5
    β₁ = 4  # 4, 7, 0.2
    X₀ = Beta(α₀, β₀)  # X|Y=0
    X₁ = Beta(α₁, β₁)  # X|Y=1
    k = beta(α₁, β₁) / beta(α₀, β₀)
    regθ(x) = 1 / (1 + k*(1-θ)*(x^(α₀-α₁))*((1-x)^(β₀-β₁))/θ)
    xx = collect(0.0:0.005:1.0)
    plot(xx, pdf.(X₀, xx), label=L"f_{\!X|Y}\,(x\,|\,0)", title="Densidades condicionales",
         xlabel=L"x", ylabel=L"f_{\!X|Y}\,(x\,|\,y)", lw =2
    )
    plot!(xx, pdf.(X₁, xx), label=L"f_{\!X|Y}\,(x\,|\,1)", lw = 2)
end


## Regresión binaria versus regresión logística

begin # ejecutar varias veces con los mismos parámetros
    θ = 0.1
    n = 10_000
    Y = rand(Bernoulli(θ), n)
    X = zeros(n)
    for i in 1:n # simular X|Y=y, y ∈ {0,1}
        if Y[i] == 0
            X[i] = rand(X₀)
        else
            X[i] = rand(X₁)
        end
    end    
    reglog1(x, a, b) = 1 / (1 + exp(-a - b*x)) # regresión logística
    reglog2(x, a, b) = 1 / (1 + exp(-a - b*log(x) + b*log(1-x))) # regresión logística con logit transformada
    function logL(a, b) # log-verosimilitud de la regresión logística
        sum(log.(reglog.(X, a, b)) .* Y + log.(1 .- reglog.(X, a, b)).*(1 .- Y))
    end
    reglog = reglog1
    res = optimize(x -> -logL(x[1], x[2]), [0.0, 0.0]) # maximización de la log-verosimilitud
    a₁, b₁ = res.minimizer
    reglog = reglog2
    res = optimize(x -> -logL(x[1], x[2]), [0.0, 0.0]) # maximización de la log-verosimilitud
    a₂, b₂ = res.minimizer
    # Graficar y comparar
    scatter(X, Y, title="Regresión binaria", 
            xlabel=L"x", ylabel=L"\theta_x = \mathbb{P}(Y=1)",
            label = "", legend = :left
    )
    plot!(xx, regθ.(xx), label=L"θ_x", lw=3)
    plot!(xx, reglog1.(xx, a₁, b₁), label="RegLog 1", lw=2)
    plot!(xx, reglog2.(xx, a₂, b₂), label="RegLog 2", lw=2)
end
