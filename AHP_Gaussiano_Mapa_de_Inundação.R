# =========================================================
# SUSCETIBILIDADE À INUNDAÇÃO - MÉTODO GAUSSIANO
# Código corrigido
# =========================================================

# ---------------------------------------------------------
# PACOTES
# ---------------------------------------------------------

library(terra)

# ---------------------------------------------------------
# PASTA DOS RASTERS
# ---------------------------------------------------------

pasta_rasters <- "C:/Usuario/Seu_Usuraio/Pasta"

# ---------------------------------------------------------
# LISTA DOS ARQUIVOS
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
# CARREGA RASTERS
# MANTENDO APENAS 1 BANDA
# ---------------------------------------------------------

rasters <- lapply(arquivos, function(x){

  r <- rast(x)

  # Mantém somente a primeira banda
  r[[1]]

})

names(rasters) <- nomes_variaveis

# ---------------------------------------------------------
# VERIFICA NÚMERO DE BANDAS
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

ref <- rasters[[1]]

# ---------------------------------------------------------
# VARIÁVEIS CATEGÓRICAS
# ---------------------------------------------------------

categoricos <- c(
  "Geomorfologia", 
  "Litologia/Geologia",
  "Uso_e_Cobertura",
  "Drenagens"
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

    cat(nome, "já alinhado.\n")

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
    "-> reprojetado com método",
    metodo,
    "\n"
  )

  r_align <- project(
    r,
    ref,
    method = metodo
  )

  return(r_align)

}

# ---------------------------------------------------------
# ALINHA TODOS OS RASTERS
# ---------------------------------------------------------

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
# VARIÁVEIS QUE DEVEM SER INVERTIDAS
# (MENOR VALOR = MAIOR INUNDAÇÃO)
# ---------------------------------------------------------

variaveis_invertidas <- c(
  "Slope",
  "Altitude",
  "HAND",
  "Distancia_drenagem"
)

# ---------------------------------------------------------
# FUNÇÃO DE NORMALIZAÇÃO
# ---------------------------------------------------------

normalizar <- function(r, inverter = FALSE){

  min_r <- global(r, "min", na.rm = TRUE)[1,1]
  max_r <- global(r, "max", na.rm = TRUE)[1,1]

  r_norm <- (r - min_r) / (max_r - min_r)

  # Inverte lógica
  if(inverter){

    r_norm <- 1 - r_norm

  }

  return(r_norm)

}

# ---------------------------------------------------------
# NORMALIZA TODOS OS RASTERS
# ---------------------------------------------------------

cat("\n====== NORMALIZAÇÃO ======\n")

for(i in seq_along(rasters_align)){

  nome <- names(rasters_align)[i]

  inverter <- nome %in% variaveis_invertidas

  if(inverter){

    cat(nome, "-> invertido\n")

  } else {

    cat(nome, "-> normal\n")

  }

  rasters_align[[i]] <- normalizar(
    rasters_align[[i]],
    inverter = inverter
  )

}

# ---------------------------------------------------------
# FUNÇÃO FUZZY GAUSSIANA
# ---------------------------------------------------------

fuzzy_gauss <- function(r){

  vals <- values(r, mat = FALSE)

  vals <- vals[!is.na(vals)]

  media <- mean(vals)

  desvio <- sd(vals)

  g <- exp(
    -((r - media)^2) /
      (2 * desvio^2)
  )

  return(g)

}

# ---------------------------------------------------------
# APLICA GAUSSIANA
# ---------------------------------------------------------

cat("\n====== APLICANDO FUNÇÃO GAUSSIANA ======\n")

rasters_gauss <- list()

for(i in seq_along(rasters_align)){

  cat(
    "Processando:",
    names(rasters_align)[i],
    "\n"
  )

  rasters_gauss[[i]] <- fuzzy_gauss(
    rasters_align[[i]]
  )

}

names(rasters_gauss) <- names(rasters_align)

# ---------------------------------------------------------
# CÁLCULO DOS PESOS
# ---------------------------------------------------------

fator_gaussiano <- function(r){

  vals <- values(r, mat = FALSE)

  vals <- vals[!is.na(vals)]

  mu <- mean(vals)

  sigma <- sd(vals)

  if(mu == 0){

    return(0)

  }

  return(sigma / mu)

}

# Calcula fatores
fatores <- sapply(
  rasters_gauss,
  fator_gaussiano
)

# Normaliza pesos
pesos_gauss <- fatores / sum(fatores)

# ---------------------------------------------------------
# TABELA DE PESOS
# ---------------------------------------------------------

tabela_pesos <- data.frame(
  Variavel = names(rasters_gauss),
  Peso = round(pesos_gauss, 4)
)

print(tabela_pesos)

# ---------------------------------------------------------
# COMBINAÇÃO LINEAR PONDERADA
# ---------------------------------------------------------

cat("\n====== COMBINAÇÃO LINEAR ======\n")

# Raster vazio
r_gauss <- rast(rasters_gauss[[1]])

values(r_gauss) <- 0

# Soma ponderada
for(i in seq_along(rasters_gauss)){

  r_temp <- rasters_gauss[[i]]

  # Garante apenas 1 banda
  r_temp <- r_temp[[1]]

  r_gauss <- r_gauss +
    (r_temp * pesos_gauss[i])

}

# ---------------------------------------------------------
# NORMALIZA RESULTADO FINAL
# ---------------------------------------------------------

min_g <- global(r_gauss, "min", na.rm = TRUE)[1,1]

max_g <- global(r_gauss, "max", na.rm = TRUE)[1,1]

r_gauss <- (r_gauss - min_g) /
  (max_g - min_g)

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
# EXPORTA RASTER FINAL
# ---------------------------------------------------------

writeRaster(
  r_gauss,
  "C:/Usuario/Seu_Usuario/Desktop/Perigo_de_inundacao_gaussiana.tif",
  overwrite = TRUE
)

# ---------------------------------------------------------
# PLOT FINAL
# ---------------------------------------------------------

plot(
  r_gauss,
  main = "Perigo de  Inundação (Gaussiano)"
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

# Plot classificado
plot(
  classes,
  main = "Classes de Suscetibilidade"
)

# ---------------------------------------------------------
# EXPORTA MAPA CLASSIFICADO
# ---------------------------------------------------------

writeRaster(
  classes,
  "C:/Users/santo/Desktop/classes_inundacao_gaussiana.tif",
  overwrite = TRUE
)

# ---------------------------------------------------------
# FINAL
# ---------------------------------------------------------

cat("\n=====================================\n")
cat("PROCESSAMENTO CONCLUÍDO COM SUCESSO\n")
cat("=====================================\n")
