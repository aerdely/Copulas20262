### Simulación y estimación no paramétrica de una cópula
### Autor: Arturo Erdely
### Fecha: 2026-04-02
## Paquetes instalados previamente

using Plots, LaTeXStrings, Distributions

## Cópula paramétrica subfamilia Clayton (Nelsen 4.2.1, θ > 0)

C(u, v, θ) = (u^(-θ) + v^(-θ) - 1)^(-1/θ)
∂uC(u, v, θ) = (u^(-θ) + v^(-θ) - 1)^(-1/θ-1) * u^(-θ-1)

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


## Estimaciones no paramétricas de la cópula 

begin # simulación de la cópula paramétrica
    n = 3_000
    θ = 2.5
    U, V = simulaCopula(θ, n)
end;

# Cópula empírica
Cn(u,v) = (1/n) * sum( (U .≤ u) .* (V .≤ v) )

# Cópula Bernstein
function CopBer(u, v, m)
    valor = 0.0
    uu = range(0, 1, m+1)
    vv = range(0, 1, m+1)
    copemp = [Cn(u,v) for u in uu, v in vv]
    for i ∈ 1:m+1
        for j ∈ 1:m+1
            valor += copemp[i,j] * pdf(Binomial(m, u), i-1) * pdf(Binomial(m, v), j-1)
        end
    end
    return valor
end

# Ejemplo de evaluación de las cópulas
u, v = 0.4, 0.7
C(u, v, θ) # valor teórico
Cn(u, v) # valor empírico
CopBer(u, v, 100) # valor Bernstein


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
    Copula = [Cn(u,v) for u in uu, v in vv] 
    nivCopula = contour(uu, vv, Copula, xlabel = L"u", ylabel = L"v",
                        xticks = [0, 0.5, 1], yticks = [0, 0.5, 1], 
                        fill = true, size = (550, 500), title = "Cópula empírica"
    )
end

# Cópula Bernstein (tarda aprox 12 segundos)
@time begin
    Copula = [CopBer(u,v,50) for u in uu, v in vv] 
    nivCopula = contour(uu, vv, Copula, xlabel = L"u", ylabel = L"v",
                        xticks = [0, 0.5, 1], yticks = [0, 0.5, 1], 
                        fill = true, size = (550, 500), title = "Cópula Bernstein"
    )
end



## Funciones auxiliares: encontrar solución (raíz) de f(c)=0 

# Método de bisección para resolver f(x)=0

include("biseccion.jl") 

@doc biseccion

begin # ejemplo
    f(x) = (x - 3) * (x - 1) * (x + 1)
    println("¿signos opuestos?", (f(-1.9), f(0)))
    b = biseccion(f, -1.9, 0)
end 
b.raíz
f(b.raíz) # valor de f en la raíz encontrada
b.numiter

# Minimización de f(x)^2 para resolver f(x)=0

using Optim 
begin    
    ff(x) = f(x[1])^2
    res = optimize(ff, [-1.9])
end
res.minimizer
res.minimum
res.iterations

## Simulación de la cópula Bernstein

function ∂γ(u, m, i)
    if i == m
        return m*(u^(m-1))
    elseif i == 0
        return -m*((1-u)^(m-1))
    else
        return pdf(Binomial(m-1, u), i-1)*m*(1 - u*(m-i)/(i*(1-u)))
    end    
end

function ∂uCopBer(u, v, m)
    valor = 0.0
    uu = range(0, 1, m+1)
    vv = range(0, 1, m+1)
    copemp = [Cn(u,v) for u in uu[2:end], v in vv[2:end]]
    for i ∈ 1:m
        for j ∈ 1:m
            valor += copemp[i,j] * ∂γ(u, m, i) * pdf(Binomial(m, v), j-1)
        end
    end
    return valor
end

∂uC(0.3, 0.6, 2.5) # valor teórico
∂uCopBer(0.3, 0.6, 100) # valor Bernstein

function ψiBer(t, u, m)
    g(v) = ∂uCopBer(u, v, m) - t
    return biseccion(g, 0.0, 1.0).raíz
    # gg(v) = g(v[1])^2
    # res = optimize(gg, [0.5])
    # return res.minimizer[1]
end

m = 100
u = 0.3
v = 0.6 
t = ∂uCopBer(u, v, m)
ψiBer(t, u, m) # valor de ψ⁻¹(t|u) para t, u, m dados

begin
    θ = 5
    n = 500
    U, V = simulaCopula(θ, n)
    scatter(U, V, xlabel = L"U", ylabel = L"V", size = (500,470), 
            ms = 1.5, mcolor = :black, legend = false, 
            title = "Simulación cópula paramétrica")
end

@time begin  # tarda aprox 1 minuto 
    m = 50
    Ub = rand(Uniform(0, 1), n)
    Vb = similar(Ub)
    for i ∈ eachindex(Ub)
        u = Ub[i]
        w = rand()
        Vb[i] = ψiBer(w, u, m)
    end
    scatter(Ub, Vb, xlabel = L"U", ylabel = L"V", size = (500,470), 
            ms = 1.5, mcolor = :black, legend = false, 
            title = "Simulación cópula Bernstein")
end

#=
  Conclusiones:

- La cópula Bernstein es una aproximación suave de la cópula teórica,
  pero solo aproxima la parte absolutamente continua.
- La simulación vía la cópula Bernstein es computacionalmente costosa.
- Hay que controlar más a detalle la precisión de las aproximaciones numéricas.
- No reproduce otras características que se verán más adelante, como tail-dependence.

=#
