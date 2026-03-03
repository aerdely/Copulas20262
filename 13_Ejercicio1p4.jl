### Ejercicio 1.4 de la tarea
### Autor: Dr. Arturo Erdely
### Versión: 2026-02-07


## Cargar paquetes necesarios:

begin
    using Distributions, Plots, LaTeXStrings # previamente instalados
    include("02probestim.jl")
end


## a) Comprobación empírica de las marginales X,Y 

@doc Exponential 

# Utilizarmos la reparametrización λ = 1/θ

function simXY(n, λ = [1,1,1])
    # simulador de (X,Y) = (min{Z1,Z12} , min{Z2,Z12})
    Z1 = Exponential(1 / λ[1])
    Z2 = Exponential(1 / λ[2])
    Z12 = Exponential(1 / λ[3])
    simZ1 = rand(Z1, n)
    simZ2 = rand(Z2, n)
    simZ12 = rand(Z12, n)
    simX = minimum([simZ1 simZ12], dims = 2)[:, 1] # [:, 1] para que sea vector columna, no matriz
    simY = minimum([simZ2 simZ12], dims = 2)[:, 1] # [:, 1] para que sea vector columna, no matriz
    return (X = simX, Y = simY)
end

function marginales(n, λ = [1,1,1])
    sim = simXY(n, λ)
    X = Exponential(1 / (λ[1] + λ[3])) # distribución marginal teórica de X
    Y = Exponential(1 / (λ[2] + λ[3])) # distribución marginal teórica de Y
    # densidades marginales (teóricas y empíricas)
    fX(x) = pdf(X, x)
    fY(y) = pdf(Y, y)
    dempX = densprob(sim.X)
    dempY = densprob(sim.Y)
    xx = collect(range(0, dempX.max, length = 1_000))
    plot(xx, dempX.fdp.(xx), lw = 2, label = "empírica", xlabel = L"X", ylabel = "densidad")
    G01 = plot!(xx, fX.(xx), lw = 2, label = "teórica")
    yy = collect(range(0, dempY.max, length = 1_000))
    plot(yy, dempY.fdp.(yy), lw = 2, label = "empírica", xlabel = L"Y", ylabel = "densidad")
    G02 = plot!(yy, fY.(yy), lw = 2, label = "teórica")
    G03 = plot(G01, G02, layout = (2,1), size = (400, 600)) # juntando G01 y G02
    # Funciones de distribución (teóricas y empíricas)
    FX(x) = cdf(X, x)
    FY(y) = cdf(Y, y)
    empX = distprob(sim.X)
    empY = distprob(sim.Y)
    plot(xx, empX.fda.(xx), lw = 4, label = "empírica", xlabel = L"X", ylabel = "función de distribución")
    G04 = plot!(xx, FX.(xx), lw = 2, label = "teórica")
    plot(yy, empY.fda.(yy), lw = 4, label = "empírica", xlabel = L"Y", ylabel = "función de distribución")
    G05 = plot!(yy, FY.(yy), lw = 2, label = "teórica")
    G06 = plot(G04, G05, layout = (2,1), size = (400, 600)) # juntando G04 y G05
    return (n = n, λ = λ, sim = sim, fX = fX, FX = FX, fY = fY, FY = FY, dempX, empX, dempY, empY, gfX = G01, gfY = G02, gf = G03, gFX = G04, gFY = G05, gF = G06)
end

# ---> especificar parámetros λ <--- #
marg = marginales(10_000, [2,5,10]);
keys(marg)
marg.gf
marg.gF


## b) Función de distribución conjunta de (X,Y)

function conjunta(marg) # teórica y empírica
    # marg: resultado de la función `marginales`
    function FXY(x,y)
        if x > 0 && y > 0 
            marg.FX(x) + marg.FY(y) - 1 + (1 - marg.FX(x))*(1 - marg.FY(y))*min(exp(marg.λ[3]*x), exp(marg.λ[3]*y))
        else
            return 0.0
        end
    end
    function empFXY(x,y)
        iX = findall(marg.sim.X .≤ x)
        iY = findall(marg.sim.Y .≤ y)
        iXY = iX ∩ iY 
        return length(iXY) / marg.n
    end
    return (teórica = FXY, empírica = empFXY)
end

begin # teórica
    FXY = conjunta(marg); # keys(FXY)
    xvec = collect(range(0, marg.empX.max, length = 100));
    yvec = collect(range(0, marg.empY.max, length = 100));
    teórica = transpose([FXY.teórica(x,y) for x ∈ xvec, y ∈ yvec]);
    contour(xvec, yvec, teórica, fill = true, xlabel = L"x", ylabel = L"y", title = L"F_{X,Y}\,(x,y)")
    plot!([0],[0], color = :black, label = "teórica")
end
begin # agregar recta y = x
    vmax = min(xvec[end], yvec[end])
    plot!([0, vmax], [0, vmax], lw = 1, color = :red, label = L"y = x")
end
begin # empírica (tarda unos segundos)
    empírica = transpose([FXY.empírica(x,y) for x ∈ xvec, y ∈ yvec]);
    contour!(xvec, yvec, empírica, fill = false, color = :cyan3, lw = 1, label = "teórica")
    plot!([0],[0], color = :cyan3, label = "empírica")
end


## c) Distribución condicional de Y|X=x

function FYdX(y, x, λ = [1,1,1]) # P(Y ≤ y | X = x, λ)
    if 0 < y < x 
        return 1- exp(-λ[2]*y) 
    elseif 0 < x ≤ y 
        return 1 - exp(-λ[2]*y - λ[3]*(y - x))*λ[1]/(λ[1] + λ[3])
    else
        return 0.0
    end
end

begin
    x = 0.1
    yvec = collect(range(0, quantile(marg.sim.Y, 0.99), length = 1_000))
    yval = [FYdX(y, x, marg.λ) for y ∈ yvec]
    scatter(yvec, yval, ms = 2, mc = :black, xlabel = L"y", title = L"F_{Y|X}\,(y\,\,|x)", label = "")
    vline!([x], color = :pink, label = "x = $x")
end


## d) No existe densidad conjunta de (X,Y)

# P(Y = X) > 0
marg.λ[3] / sum(marg.λ) # teórica 
mean(marg.sim.X .== marg.sim.Y) # estimación empírica


## e) Gráfico de dispersión (scatter plot) de (X,Y)

function gdisper(n, λ = [1,1,1])
    sim = simXY(n, λ)
    PXY = λ[3] / sum(λ)
    scatter(sim.X, sim.Y, ms = 1, mc = :black, xlabel = L"X", ylabel = L"Y",
    label = "")
    iguales = findall(sim.X .== sim.Y)
    scatter!(sim.X[iguales], sim.Y[iguales], ms = 1, mc = :blue, msw = 0, label = "P(Y = X) = $(round(PXY, digits = 3))")
end

gdisper(10_000, [2,5,10]) # ejecutar varias veces, y con distintos parámetros
