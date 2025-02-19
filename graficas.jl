using Plots

# Matriz de prueba (reemplaza con tu matriz real)
A = rand(ComplexF64, 80, 80)

# Tomamos la magnitud de la matriz
A_mag = abs.(A)

# Graficamos la superficie
p = surface(A_mag, title="Análisis de contingencias n-1", xlabel="Línea m-n", ylabel="Línea p-q", 
            zlabel="Aumento de corriente en magnitud", color=:viridis)

# Guardamos la gráfica en un archivo PNG
savefig(p, "grafica_contingencias_n1.png")