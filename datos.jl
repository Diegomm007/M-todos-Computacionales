using CSV, DataFrames

df = CSV.read("Cont.csv", DataFrame)  # Carga el CSV en un DataFrame


function dataframe_a_md(df)
    encabezado = "| " * join(names(df), " | ") * " |"
    separador = "|-" * join(repeat(["-"], size(df,2)), "-|")
    filas = [ "| " * join(row, " | ") * " |" for row in eachrow(df) ]
    return join([encabezado, separador, filas...], "\n")
end

tabla_md = dataframe_a_md(df)  # Convertir CSV a tabla Markdown



function insertar_csv_en_md(md_file, csv_texto, marcador)
    texto = read(md_file, String)  # Leer todo el Markdown como String
    nuevo_texto = replace(texto, marcador => marcador * "\n\n" * csv_texto * "\n\n")  # Insertar CSV después del marcador
    write(md_file, nuevo_texto)  # Guardar los cambios
end

insertar_csv_en_md("README.md", tabla_md, "<!-- INSERTAR_CSV_AQUI -->")
