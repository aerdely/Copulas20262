### Simular cópula Marshall-Olkin
### Autor: Arturo Erdely
### Fecha: 2026-03-21

#=
    Cópula Marshall-Olkin con parámetros 0<α<1 y 0<β<1, dada por
    C(u,v) = mín{u^(1-α)v, uv^(1-β)}
=#

## Paquetes instalados previamente

using Plots, LaTeXStrings, Distributions


## Funciones auxiliares para simular

# curva de parte singular 
h(u,α,β) = u^(α/β) 

# ψ(v|u) = ∂C(u,v)/∂u como función de v, para u, α, β fijos
function ∂C∂u(v, u, α, β) 
    if 0 < v < h(u, α, β)
        return (1-α)*u^(-α)*v
    elseif h(u, α, β) < v < 1
        return v^(1-β)
    else
        return NaN
    end
end

# ψ⁻¹(t|u) como función de t, para u, α, β fijos
function ψi(t, u, α, β)
    if 0 < t < (1-α)*u^(α/β - α) 
        return t*(u^α) / (1-α)
    elseif (1-α)*u^(α/β - α) ≤ t ≤ u^(α/β - α)
        return u^(α/β)
    elseif u^(α/β - α) < t < 1
        return t^(1/(1-β))
    else
        return NaN
    end
end

# Gráfica de ψ(v|u) para u, α, β fijos
begin
    α, β = 0.3, 0.8   # 0 < α < 1, 0 < β < 1
    u = 0.2   # 0 < u < 1
    vv = collect(range(0, 1, length = 1_000))
    ψ = [∂C∂u(v, u, α, β) for v ∈ vv]
    scatter(vv, ψ, xlabel = L"v", ylabel = L"\partial C_{\theta}(u,v) / \partial u",
         title = "u = $u, (α,β) = $((α,β))", size = (500,400), label = L"v\mapsto \partial C(u,v) / \partial u",
         ms = 2.0, mcolor = :black)
    hu = h(u,α,β)
    vg1 = (1-α)*u^(α/β - α)
    vg2 = u^(α/β - α)
    vline!([hu], ls = :dash, color = :red, label = "u^(α/β) = $(round(hu, digits=4))")
    hline!([vg1], ls = :dash, color = :gray, label = "(1-α)*u^(α/β - α) = $(round(vg1, digits=4))")
    hline!([vg2], ls = :dash, color = :gray, label = "u^(α/β - α) = $(round(vg2, digits=4))")
end

# Agregar gráfica de ψ⁻¹(t|u) para u, α, β fijos (para comprobar)
begin
    tt = collect(range(0, 1, length = 1_000))
    plot!([ψi(t, u, α, β) for t ∈ tt], tt, color = :cyan,
               label = L"\psi^{(-1)}_{θ}(t|u)")
end


## Simulación

# Función para simular n observaciones de la cópula con parámetros α y β
function simular_copula(n, α, β)
    U = rand(n)
    V = similar(U)
    for i ∈ 1:n
        u = U[i]
        w = rand()
        V[i] = ψi(w, u, α, β)
    end
    return U, V
end

#= 
   Función para calcular las probabilidades teóricas y empíricas
   de estar por encima, por debajo o en la parte singular de la cópula
=#
function probabilidad(ε = 0.0001)
    teórica_superior = α*(1-β) / (α + β - α*β)
    teórica_inferior = (1-α)*β / (α + β - α*β)
    teórica_singular = α*β / (α + β - α*β)
    conteo_superior = 0
    conteo_inferior = 0
    conteo_singular = 0
    for i ∈ 1:n
        if V[i] > h(U[i], α, β) + ε
            conteo_superior += 1
        elseif V[i] < h(U[i], α, β) - ε
            conteo_inferior += 1
        else
            conteo_singular += 1
        end
    end
    empírica_superior = conteo_superior / n
    empírica_inferior = conteo_inferior / n
    empírica_singular = conteo_singular / n
    println("Probabilidad estar por encima de la parte singular:")
    println("--------------------------------------------------")
    println("Teórica: ", round(teórica_superior, digits=4), " | Empírica: ", round(empírica_superior, digits=4))
    println()
    println("Probabilidad estar por debajo de la parte singular:")
    println("--------------------------------------------------")
    println("Teórica: ", round(teórica_inferior, digits=4), " | Empírica: ", round(empírica_inferior, digits=4))
    println()
    println("Probabilidad estar en la parte singular: ")    
    println("---------------------------------------")
    println("Teórica: ", round(teórica_singular, digits=4), " | Empírica: ", round(empírica_singular, digits=4))
    println()
    return (teórica = (superior = teórica_superior, inferior = teórica_inferior, singular = teórica_singular), empírica = (superior = empírica_superior, inferior = empírica_inferior, singular = empírica_singular))
end

begin
    n = 10_000
    α, β = 1/3, 1/4   # 0 < α < 1, 0 < β < 1
    U, V = simular_copula(n, α, β)
    scatter(U, V, xlabel = L"U", ylabel = L"V", size = (500,470), 
            ms = 1.5, mcolor = :black, label = "")
end
# agregar curva parte singular
plot!(sort(U), h.(sort(U), α, β), color = :red, lw = 1.0, label = L"u^{α/β}")

pp = probabilidad();
