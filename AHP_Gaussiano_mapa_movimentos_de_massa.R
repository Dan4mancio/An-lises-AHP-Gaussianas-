# Pacotes necessários
library(terra)

# Caminho da pasta com os .tif
pasta_rasters <- "C:/Usuario/Seu_User/Pasta"

# Lista dos arquivos
arquivos <- list.files(pasta_rasters, pattern = "\\.tif$", full.names = TRUE)
nomes_variaveis <- tools::file_path_sans_ext(basename(arquivos))

# Carrega todos os rasters
rasters <- lapply(arquivos, rast)
names(rasters) <- nomes_variaveis

# Raster de referência (índice i)
ref <- rasters[[i]]
# para as análises de movimentos de massa provavelmente seu raster referência será o de declividade, recomendo colocar um número na frente
# para ele ter i= 1 se você for fazer diversas análises.

# Define quais são categóricos (ajuste conforme seus dados)
categoricos <- c("Geomorfologia", "Litologia/Geologia",
                 "Uso_e_Cobertura", "Vias_de_acesso")

# Função que alinha (reprojeta + reamostra) cada raster ao ref
alinhar_raster <- function(r, nome, ref) {
  # Se já estiver exatamente igual, retorna sem modificar
  if (compareGeom(r, ref, stopOnError = FALSE, res = TRUE, ext = TRUE, crs = TRUE)) {
    cat(nome, "já está alinhado.\n")
    return(r)
  }
  
  # Escolhe método: "near" para categóricos, "bilinear" para contínuos
  metodo <- if (nome %in% categoricos) "near" else "bilinear"
  cat(nome, "projetado/reamostrado com método", metodo, "\n")
  
  # project() já faz tudo: muda CRS, extensão e resolução para ficar igual ao ref
  r_align <- project(r, ref, method = metodo)
  return(r_align)
}

# Aplica o alinhamento para todos
rasters_align <- list()
for (i in seq_along(rasters)) {
  rasters_align[[i]] <- alinhar_raster(rasters[[i]], names(rasters)[i], ref)
}
names(rasters_align) <- names(rasters)

# -------------------------------------------------------------------
# Cálculo dos pesos gaussianos a partir dos rasters JÁ ALINHADOS
# -------------------------------------------------------------------
fator_gaussiano <- function(r) {
  vals <- values(r, mat = FALSE)
  vals <- vals[!is.na(vals)]
  if (length(vals) == 0) return(0)
  norm_vals <- vals / sum(vals)        # normaliza pela soma
  mu <- mean(norm_vals)
  if (mu == 0) return(0)
  sigma <- sd(norm_vals)
  gf <- sigma / mu
  return(gf)
}

fatores <- sapply(rasters_align, fator_gaussiano)
pesos_gauss <- fatores / sum(fatores)

print(data.frame(variavel = nomes_variaveis, peso_gauss = pesos_gauss))

# -------------------------------------------------------------------
# Combinação linear ponderada
# -------------------------------------------------------------------
# Opcional: normalizar cada raster para [0,1] antes de somar (evita que variáveis com
# magnitudes muito diferentes dominem o resultado). Descomente se desejar:
# for (i in seq_along(rasters_align)) {
#   r <- rasters_align[[i]]
#   r <- (r - global(r, "min", na.rm = TRUE)[[1]]) /
#        (global(r, "max", na.rm = TRUE)[[1]] - global(r, "min", na.rm = TRUE)[[1]])
#   rasters_align[[i]] <- r
# }

# Inicializa com o primeiro raster ponderado
r_gauss <- rasters_align[[1]] * pesos_gauss[1]

# Soma os demais
for (i in 2:length(rasters_align)) {
  r_gauss <- r_gauss + rasters_align[[i]] * pesos_gauss[i]
}

# Normaliza o resultado final para [0,1]
min_g <- global(r_gauss, "min", na.rm = TRUE)[[1]]
max_g <- global(r_gauss, "max", na.rm = TRUE)[[1]]
r_gauss <- (r_gauss - min_g) / (max_g - min_g)

# Salva o raster final
writeRaster(r_gauss, "Movimentos _de_Massa_gaussiano.tif", overwrite = TRUE)

cat("Processo concluído sem erros de CRS!\n")

#---------------------------------------------------
# Para visualizar o mapa
#---------------------------------------------------

library(terra)

# Carregar o raster gerado
r <- rast("Movimentos _de_Massa_gaussiano.tif")

# Plotar
plot(r, 
     main = "Movimentos _de_Massa_gaussiano.tif",
     col = hcl.colors(100, "viridis"),  # paleta viridis
     range = c(0,1),                    # já normalizado entre 0 e 1
     axes = TRUE)
