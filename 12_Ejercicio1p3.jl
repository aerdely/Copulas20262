### Ejercicio 1.3 de la tarea
### Autor: Dr. Arturo Erdely
### Versión: 2026-02-07


## Cargar paquetes necesarios:

using Plots, LaTeXStrings # previamente instalados


## a) Función de densidad conjunta de (X,Y)

begin 
    f(x,y) = exp(-y) * (0 < x < y)
    xvec = collect(range(0, 3.5, length = 100))
    yvec = collect(range(0, 3.5, length = 100))
    zmat = transpose([f(x,y) for x in xvec, y in yvec])
    contour(xvec, yvec, zmat, xlabel = L"x", ylabel = L"y", fill = true, 
            size = (450,450), title = "Conjuntos de nivel de f(x,y)")
    graf_a = plot!([0, 3.5], [0, 3.5], lw = 2, lc = :red, 
                   label = L"y = x", legend = :bottomright)
end


## b) Función de distribución conjunta de (X,Y)

leyenda = true
begin 
    function F(x,y)
        if 0 < x ≤ y
            return 1 - exp(-x) - x*exp(-y)
        elseif 0 < y < x 
            return 1 - exp(-y) - y*exp(-y)
        else
            return 0.0
        end
    end
    Zmat = transpose([F(x,y) for x in xvec, y in yvec])
    contour(xvec, yvec, Zmat, xlabel = L"x", ylabel = L"y", fill = true, 
            size = (450,450), title = L"F_{X,Y}\,(x,y)", legend = leyenda)
    graf_b = plot!([0, 3.5], [0, 3.5], lw = 2, lc = :red, 
                   label = L"y = x", legend = :bottomright)
end

begin # Si X,Y fuesen independientes
    FX(x) = F(x, Inf)
    FY(y) = F(Inf, y)
    Findep(x,y) = FX(x)*FY(y)
    Zindep = transpose([Findep(x,y) for x in xvec, y in yvec])
    contour(xvec, yvec, Zindep, xlabel = L"x", ylabel = L"y", fill = true, 
            size = (450,450), title = L"F_X(x)F_Y(y)", legend = leyenda)
    graf_bextra = plot!([0, 3.5], [0, 3.5], lw = 2, lc = :red, 
                   label = L"y = x", legend = :bottomright)
end

# Comparando F(x,y) versus X⟂Y
plot(graf_b, graf_bextra, layout = (1,2), size = (900, 450))


## c) Scatter plot de simulaciones de (X,Y)

function simXY(n)
    uu = rand(n)
    xx = -log.(1 .- uu)
    vv = rand(n)
    yy = xx .- log.(1 .- vv)
    return (x = xx, y = yy)
end

begin
    n = 3_000
    sim = simXY(n);
    vmax = maximum(sim.x ∪ sim.y)
    graf_c = scatter(sim.x, sim.y, mc = :lightgreen, ms = 2,
                     xlims = (0, vmax), ylims = (0, vmax),
                     size = (450,450), label = "")
    graf_cbis = plot(graf_c, [0, vmax], [0, vmax], lw = 2, lc = :red, 
                     xlabel = L"x", ylabel = L"y", title = "Simulaciones de (X,Y)", 
                     label = L"y = x", legend = :bottomright)
end


## d) Agregar gráfica c) a la de a) 

plot(graf_a)
scatter!(sim.x, sim.y, mc = :green, ms = 1, markerstrokewidth = 0,
         xlims = (0, 3.5), ylims = (0, 3.5), label = "")
