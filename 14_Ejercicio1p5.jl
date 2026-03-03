### Ejercicio 1.5 de la tarea
### Autor: Dr. Arturo Erdely
### Versión: 2026-02-16


## Cargar paquetes necesarios:

using Plots, LaTeXStrings # previamente instalados


## Familia de cópulas Marshall-Olkin

begin
    K(u,v,α,β) = min((u^(1-α))*v, u*(v^(1-β)))
    g(u,α,β) = u^(α/β) # parte singular
    uu = collect(range(0, 1, length = 100))
    vv = collect(range(0, 1, length = 100))
end;


## Curvas de nivel de K(u,v,α,β)

begin 
    α = 0.5
    β = 0.3
    gg = [g(u,α,β) for u ∈ uu]
    ww = transpose([K(u,v,α,β) for u ∈ uu, v ∈ vv])
    contour(uu, vv, ww, fill = true, xlabel = L"u", ylabel = ylabel = L"v", 
            title = "K(u,v)    α = $α    β = $β", size = (430,400))
end
# identificar parte singular:
plot!(uu, gg, lw = 2, color = :green, label = L"v = u^{\alpha/\beta}")
