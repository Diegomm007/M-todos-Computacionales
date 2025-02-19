using LinearAlgebra
using DataFrames
using CSV
using Plots
using SparseArrays
using DelimitedFiles, Printf
using PrettyTables
# Calcular la matriz de admitancia nodal

function calcular_ybus(lines,nodes)
    
"""
    Entradas: lines (DataFrames)
              nodes (DataFrames)

    Salida:   Ybus : matriz
    
"""
num_nodes = nrow(nodes)
num_lines = nrow(lines)
Ybus = zeros(num_nodes,num_nodes)*1im
for k = 1:num_lines
    n1 = lines.FROM[k]
    n2 = lines.TO[k]
    yL = 1/(lines.R[k]+lines.X[k]*1im)
    Bs = lines.B[k]*1im/2
    Ybus[n1, n1] += yL + Bs # diagonal
    Ybus[n1, n2] -= yL # fuera
    Ybus[n2, n1] -= yL # fuera
    Ybus[n2, n2] += yL + Bs # diagonal
end
return Ybus
end

## Función principal

lines = DataFrame(CSV.File("lines.csv"))
nodes = DataFrame(CSV.File("nodes.csv"))

Ybus = calcular_ybus(lines, nodes)
#Ybus = sparse(Ybus)
## GRÁFICA
# Función para formatear los números con ancho fijo
function formatear(z)
    return @sprintf("%10f + %10fim", real(z), imag(z))
end

# Convertir cada número a una cadena alineada
matriz_str = join([join(formatear.(fila), ", ") for fila in eachrow(Ybus)], "\n")

# Guardar en CSV
open("Y_bus.csv", "w") do file
    write(file, matriz_str)
end

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

## Función para hacer TAP = 1
function modificar_taps(lines)
        
"""
Entradas: lines (DataFrames)

Salida:   linesmod (DataFrames) con Tap = 1

"""

# Vamos a reemplazar los valores "0" de la columna TAP de lines.csv  por "1"
num_lines = nrow(lines)
for k in 1:num_lines
    if lines[k, :TAP] == 0
        lines[k, :TAP] = 1
    end
end
CSV.write("linesmodf.csv", lines)
end


## Función para Ybus con taps

lines_modf = DataFrame(CSV.File("linesmodf.csv"))

function calcular_ybus_taps(linesmodf,nodes)

"""
    Entradas: linesmodf (DataFrames)
              nodes (DataFrames)

    Salida:   Ybus : matriz con modelo de taps
    
"""
## Se añade el modelo pi del transformador con el tap modificado

num_nodes = nrow(nodes)
num_lines_mod = nrow(linesmodf)
Ybus_taps = zeros(num_nodes,num_nodes)*1im
for k = 1:num_lines_mod
    n1 = linesmodf.FROM[k]
    n2 = linesmodf.TO[k]
    yL = 1/(linesmodf.R[k]+linesmodf.X[k]*1im)
    t = linesmodf.TAP[k]
    yL_t = yL*t
    Bs = linesmodf.B[k]*1im/2
    Bs_1 = (t*t-t)*yL
    Bs_2 = (1-t)*yL
    Ybus_taps[n1, n1] += yL_t + Bs + Bs_1 # diagonal
    Ybus_taps[n1, n2] -= yL_t # fuera
    Ybus_taps[n2, n1] -= yL_t # fuera
    Ybus_taps[n2, n2] += yL_t + Bs + Bs_2 # diagonal
end
return Ybus_taps

end
Ybus_taps = calcular_ybus_taps(lines_modf,nodes)

matriz_str1 = join([join(formatear.(fila), ", ") for fila in eachrow(Ybus_taps)], "\n")

# Guardar en CSV
open("Y_bus_taps.csv", "w") do file
    write(file, matriz_str1)
end

#- - - - - - - - - - - - - - - - - - -- - - - - - - - - - 

##FLUJO DC

# 1. Declaramos una función para eliminar el nodo slack


function eliminar_slack(nodes)
    
"""
    Entradas: nodes (DataFrames)

    Salida:   nodesmod : DataFrame sin el nodo slack
    
"""
nodes = filter(row -> row.TYPE != 3, nodes) #Se recorre cada fila con row  y se filtra para dejar las filas conde TYPE no sea 3
CSV.write("nodesmodf.csv", nodes)   
end

## Ejemplo para elminar el nodo slack de nodes.csv

eliminar_slack(nodes)  # Se crea nodesmodf.csv

# 2. Declaramos una función para calcular la matriz B

function calcular_B(nodes, lines)
    
"""
    Entradas: nodes (DataFrames)
              lines (DataFrames)

    Salida:   B : matriz
    
"""
num_nodes = nrow(nodes)
num_lines = nrow(lines)
B = zeros(num_nodes,num_nodes)
for k = 1:num_lines
    n1 = lines.FROM[k]
    n2 = lines.TO[k]
    yL = 1/(lines.X[k])
    B[n1, n1] += yL # diagonal
    B[n1, n2] -= yL # fuera
    B[n2, n1] -= yL # fuera
    B[n2, n2] += yL # diagonal
end
for k = 1:num_nodes
    j = nodes.TYPE[k]

    if j == 3
        slack = nodes.NUMBER[k] 
        B = B[setdiff(1:end, slack), setdiff(1:end, slack)]
    end
end
return B
end
# Ejemplo para calcular la matriz B

B = calcular_B(nodes,lines) 

# Convertir la matriz a un DataFrame con los nombres personalizados
df = DataFrame(B, :auto)

# Guardar el DataFrame como un archivo CSV
CSV.write("B.csv", df; writeheader=false)

## 3. Calcular la matriz de angulos

function calcular_theta(nodes,nodesmodf,B)
    
"""
    Entradas: nodes (DataFrames)
              nodesmodf (DataFrames)
              B (matriz de admitancia DC)

    Salida:   Theta : vector de angulos 
    
""" 
num_nodes = nrow(nodesmodf)
P = zeros(num_nodes,1)
for k = 1:num_nodes
    P_neta = nodesmodf.PGEN[k]-nodesmodf.PLOAD[k]
    P[k,1] = P_neta
    end
Thetha = inv(B)*P
#Agregamos el angulo del nodo slack
num_nodes_slack = nrow(nodes)
for k = 1:num_nodes_slack
    j = nodes.TYPE[k]

    if j == 3
        slack = nodes.NUMBER[k] 
        Thetha = vcat(Thetha[1:slack-1, :], [0.0], Thetha[slack:end, :]) 
    end
end
return Thetha
end

# Ejemplo para calcular la matriz de angulos
nodes_modf = DataFrame(CSV.File("nodesmodf.csv"))
Theta = calcular_theta(nodes,nodes_modf,B)

# Convertir la matriz a un DataFrame con los nombres personalizados
df = DataFrame(Theta, :auto)

# Guardar el DataFrame como un archivo CSV
CSV.write("Theta.csv", df; writeheader=false)




## 4. Calcular los flujos de Potencia
function calcular_flujos(lines,Theta)
    
    """
        Entradas: lines (DataFrames)
                  Theta (vector de angulos)
    
        Salida:   P_flujos : vector de flujos de potencia 
        
    """
num_lines = nrow(lines)
P = zeros(num_lines,1)
for k = 1:num_lines
    n1 = lines.FROM[k]
    n2 = lines.TO[k]
    angulo = Theta[n1] - Theta[n2]
    P_k = angulo/(lines.X[k])
    P[k, 1] = P_k
end
return P
end
    
# Ejemplo para calcular los flujos de potencia
P_flujos = calcular_flujos(lines,Theta)

# Convertir la matriz a un DataFrame con los nombres personalizados
df = DataFrame(P_flujos, :auto)

# Guardar el DataFrame como un archivo CSV
CSV.write("P_flujos.csv", df; writeheader=false)


#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

## Análisis de contingencias n-1

function contingencias(ybus,lines,nodes)

    """
        Entradas: ybus (matriz de admitancia)
                  lines (DataFrames)
                  nodes (DataFrames)
    
        Salida:   contingencias : DataFrame con las contingencias
        
    """
zbus = inv(ybus)
num_lines = nrow(lines)
contingencias = zeros(num_lines,num_lines)*1im
for k = 1:num_lines
    m = lines.FROM[k]
    n = lines.TO[k]
    # Se define el factor de distribución de la línea
    for j = 1:num_lines
        if j != k
        p = lines.FROM[j]
        q = lines.TO[j]
        z_pm = zbus[p,m]
        z_pn = zbus[p,n]
        z_qm = zbus[q,m]
        z_qn = zbus[q,n]
        z_mm = zbus[m,m]
        z_nn = zbus[n,n]
        z_mn = zbus[m,n]
        L_pq_mn = (-z_mn/z_qm)*((z_pm-z_pn-z_qm+z_qn)/(z_mm+z_nn+2*z_mn-z_mn))
    # Se calcula la corriente previa a la salida de la linea
        I_pq_old = (nodes.VPU[p]-nodes.VPU[q])/(zbus[p,q])
        I_mn_old = (nodes.VPU[m]-nodes.VPU[n])/(zbus[m,n])
    # Se calcula la corriente después de la salida de la línea
        I_pq_new = I_pq_old + L_pq_mn*I_mn_old
        contingencias[k,j] = I_pq_new
        end
    end
end
return contingencias
end

# Ejemplo para calcular las contingencias

cont = contingencias(Ybus,lines,nodes)
display(cont)

# Convertir la matriz a un DataFrame con los nombres personalizados
df = DataFrame(cont, :auto)

# Guardar el DataFrame como un archivo CSV
CSV.write("Cont.csv", df; writeheader=false)




# -----------------------------------------
A = cont
using Graphs, GraphPlot, Cairo, Compose

g = SimpleGraph(80)  # Crear un grafo con 80 nodos (líneas de transmisión)
     
     # Agregar conexiones donde el aumento de corriente es significativo
for i in 1:80, j in 1:80
     if abs(A[i, j]) > 0.1  # Filtro para conexiones significativas
         add_edge!(g, i, j)
    end
end
     
# Dibujar el grafo
p = gplot(g, nodelabel=1:80, nodesize=5, title="")
draw(PNG("grafo.png", 800, 800), p)

A_mag = abs.(A)

# Graficamos la superficie
p = surface(A_mag, title="Análisis de contingencias n-1", xlabel="Línea m-n", ylabel="Línea p-q", 
            zlabel="Aumento de corriente en magnitud", color=:viridis)

# Guardamos la gráfica en un archivo PNG
savefig(p, "contingencias.png")
