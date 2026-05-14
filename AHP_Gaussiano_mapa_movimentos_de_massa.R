# =========================================================
# SUSCETIBILIDADE A MOVIMENTOS DE MASSA
# MÉTODO GAUSSIANO
# =========================================================

# ---------------------------------------------------------
# PACOTES
# ---------------------------------------------------------

library(terra)

# ---------------------------------------------------------
# CAMINHO DA PASTA DOS RASTERS
# ---------------------------------------------------------

pasta_rasters <- "C:/Usuario/Seu_User/Pasta"

# ---------------------------------------------------------
# LISTA DOS ARQUIVOS .TIF
# ---------------------------------------------------------

arquivos <- list.files(
  pasta_rasters,
  pattern = "\\.tif$",
  full.names = TRUE
)

# ---------------------------------------------------------
# NOMES DAS VARIÁVEIS
# ---------------------------------------------------------

nomes_variaveis <- tools::file_path_sans_ext(
  basename(arquivos)
)

# ---------------------------------------------------------
# CARREGAR RASTERS
# MANTENDO APENAS 1 BANDA
# ---------------------------------------------------------

rasters <- lapply(arquivos, function(x){

  r <- rast(x)

  # Mantém somente a primeira banda
  r[[1]]

})

names(rasters) <- nomes_variaveis

# ---------------------------------------------------------
# VERIFICAÇÃO DAS BANDAS
# ---------------------------------------------------------

cat("\n====== VERIFICAÇÃO DAS BANDAS ======\n")

for(i in seq_along(rasters)){

  cat(
    names(rasters)[i],
    "->",
    nlyr(rasters[[i]]),
    "banda(s)\n"
  )

}

# ---------------------------------------------------------
# RASTER DE REFERÊNCIA
# ---------------------------------------------------------
# Recomenda-se utilizar o raster de declividade
# como referência para movimentos de massa.
# Exemplo:
# 01_Declividade.tif

ref <- rasters[[1]]

# ---------------------------------------------------------
# VARIÁVEIS CATEGÓRICAS
# ---------------------------------------------------------

categoricos <- c(
  "Geomorfologia",
  "Litologia_Geologia",
  "Uso_e_Cobertura",
  "Vias_de_acesso"
)

# ---------------------------------------------------------
# FUNÇÃO DE ALINHAMENTO
# ---------------------------------------------------------

alinhar_raster <- function(r, nome, ref){

  # Verifica se já está alinhado
  if(compareGeom(
    r,
    ref,
    stopOnError = FALSE,
    res = TRUE,
    ext = TRUE,
    crs = TRUE
  )){

    cat(nome, "já está alinhado.\n")

    return(r)

  }

  # Método de reamostragem
  metodo <- if(nome %in% categoricos){
    "near"
  } else {
    "bilinear"
  }

  cat(
    nome,
    "projetado/reamostrado com método",
    metodo,
    "\n"
  )

  # Reprojeção + reamostragem
  r_align <- project(
    r,
    ref,
    method = metodo
  )

  return(r_align)

}

# ---------------------------------------------------------
# ALINHAMENTO DOS RASTERS
# ---------------------------------------------------------

cat("\n====== ALINHAMENTO DOS RASTERS ======\n")

rasters_align <- list()

for(i in seq_along(rasters)){

  rasters_align[[i]] <- alinhar_raster(
    rasters[[i]],
    names(rasters)[i],
    ref
  )

}

names(rasters_align) <- names(rasters)

# ---------------------------------------------------------
# NORMALIZAÇÃO [0-1]
# ---------------------------------------------------------
# Evita que variáveis com magnitudes maiores
# dominem o resultado final.

normalizar <- function(r){

  min_r <- global(
    r,
    "min",
    na.rm = TRUE
  )[1,1]

  max_r <- global(
    r,
    "max",
    na.rm = TRUE
  )[1,1]

  r_norm <- (r - min_r) /
    (max_r - min_r)

  return(r_norm)

}

cat("\n====== NORMALIZAÇÃO ======\n")

for(i in seq_along(rasters_align)){

  cat(
    "Normalizando:",
    names(rasters_align)[i],
    "\n"
  )

  rasters_align[[i]] <- normalizar(
    rasters_align[[i]]
  )

}

# ---------------------------------------------------------
# FUNÇÃO GAUSSIANA
# ---------------------------------------------------------

fator_gaussiano <- function(r){

  vals <- values(r, mat = FALSE)

  vals <- vals[!is.na(vals)]

  if(length(vals) == 0){

    return(0)

  }

  norm_vals <- vals / sum(vals)

  mu <- mean(norm_vals)

  if(mu == 0){

    return(0)

  }

  sigma <- sd(norm_vals)

  gf <- sigma / mu

  return(gf)

}

# ---------------------------------------------------------
# CÁLCULO DOS PESOS
# ---------------------------------------------------------

cat("\n====== CÁLCULO DOS PESOS ======\n")

fatores <- sapply(
  rasters_align,
  fator_gaussiano
)

pesos_gauss <- fatores / sum(fatores)

# ---------------------------------------------------------
# TABELA DE PESOS
# ---------------------------------------------------------

tabela_pesos <- data.frame(
  Variavel = nomes_variaveis,
  Peso_Gaussiano = round(pesos_gauss, 4)
)

print(tabela_pesos)

# ---------------------------------------------------------
# COMBINAÇÃO LINEAR PONDERADA
# ---------------------------------------------------------

cat("\n====== COMBINAÇÃO LINEAR ======\n")

# Raster vazio
r_gauss <- rast(
  rasters_align[[1]]
)

values(r_gauss) <- 0

# Soma ponderada
for(i in seq_along(rasters_align)){

  r_temp <- rasters_align[[i]][[1]]

  r_gauss <- r_gauss +
    (r_temp * pesos_gauss[i])

}

# ---------------------------------------------------------
# NORMALIZAÇÃO FINAL
# ---------------------------------------------------------

min_g <- global(
  r_gauss,
  "min",
  na.rm = TRUE
)[1,1]

max_g <- global(
  r_gauss,
  "max",
  na.rm = TRUE
)[1,1]

r_gauss <- (
  r_gauss - min_g
) / (
  max_g - min_g
)

# ---------------------------------------------------------
# VERIFICA RESULTADO FINAL
# ---------------------------------------------------------

cat("\n====== RESULTADO FINAL ======\n")

print(r_gauss)

cat(
  "\nNúmero de bandas:",
  nlyr(r_gauss),
  "\n"
)

# ---------------------------------------------------------
# EXPORTAR RASTER FINAL
# ---------------------------------------------------------

writeRaster(
  r_gauss,
  "Movimentos_de_Massa_gaussiano.tif",
  overwrite = TRUE
)

# ---------------------------------------------------------
# CLASSIFICAÇÃO OPCIONAL
# ---------------------------------------------------------

classes <- classify(
  r_gauss,
  matrix(
    c(
      0.0, 0.2, 1,
      0.2, 0.4, 2,
      0.4, 0.6, 3,
      0.6, 0.8, 4,
      0.8, 1.0, 5
    ),
    ncol = 3,
    byrow = TRUE
  )
)

# ---------------------------------------------------------
# EXPORTAR MAPA CLASSIFICADO
# ---------------------------------------------------------

writeRaster(
  classes,
  "Classes_Movimentos_de_Massa.tif",
  overwrite = TRUE
)

# ---------------------------------------------------------
# VISUALIZAÇÃO DO MAPA CONTÍNUO
# ---------------------------------------------------------

plot(
  r_gauss,
  main = "Suscetibilidade a Movimentos de Massa",
  col = hcl.colors(100, "viridis"),
  range = c(0,1),
  axes = TRUE
)

# ---------------------------------------------------------
# VISUALIZAÇÃO DO MAPA CLASSIFICADO
# ---------------------------------------------------------

plot(
  classes,
  main = "Classes de Suscetibilidade",
  col = hcl.colors(5, "inferno"),
  axes = TRUE
)

# ---------------------------------------------------------
# FINALIZAÇÃO
# ---------------------------------------------------------

cat("\n=====================================\n")

cat("PROCESSAMENTO CONCLUÍDO COM SUCESSO\n")

cat("\nRaster contínuo salvo como:\n")

cat("Movimentos_de_Massa_gaussiano.tif\n")

cat("\nRaster classificado salvo como:\n")

cat("Classes_Movimentos_de_Massa.tif\n")

cat("=====================================\n")
