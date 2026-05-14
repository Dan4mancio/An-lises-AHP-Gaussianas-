# =========================================================
# COMBINAÇÃO DE MAPAS DE SUSCETIBILIDADE
# MOVIMENTO DE MASSA + INUNDAÇÃO
# MÉTODO GAUSSIANO CORRIGIDO
# =========================================================

# ---------------------------------------------------------
# PACOTES
# ---------------------------------------------------------

library(terra)

# ---------------------------------------------------------
# CAMINHOS DOS RASTERS
# ---------------------------------------------------------

caminho_mapa1 <- "C:/Usuário/Seu_Usuario/Desktop/Mapas/ Seu_mapa_de_Perigo_a_movimento_de_massa.tif"

caminho_mapa2 <- "C:/Usuário/Seu_Usuario/Desktop/Mapas/ Seu_mapa_de_Perigo_a_inundação.tif"

# ---------------------------------------------------------
# NOMES DOS MAPAS
# ---------------------------------------------------------

nome_mapa1 <- "Movimento de Massa"

nome_mapa2 <- "Inundacao"

# ---------------------------------------------------------
# CARREGAR RASTERS
# MANTENDO APENAS 1 BANDA
# ---------------------------------------------------------

r1 <- rast(caminho_mapa1)[[1]]

r2 <- rast(caminho_mapa2)[[1]]

# Lista nomeada
rasters <- list(r1, r2)

names(rasters) <- c(
  nome_mapa1,
  nome_mapa2
)

# ---------------------------------------------------------
# VERIFICA BANDAS
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
# RASTERS CATEGÓRICOS
# (AQUI NÃO EXISTEM)
# ---------------------------------------------------------

categoricos <- c()

# ---------------------------------------------------------
# FUNÇÃO DE ALINHAMENTO
# ---------------------------------------------------------

alinhar_raster <- function(r, nome, ref){

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

  r_align <- project(
    r,
    ref,
    method = metodo
  )

  return(r_align)

}

# ---------------------------------------------------------
# ALINHAMENTO
# ---------------------------------------------------------

cat("\n====== ALINHAMENTO ======\n")

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
# MUITO IMPORTANTE
# ---------------------------------------------------------

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
# CÁLCULO DOS PESOS GAUSSIANOS
# ---------------------------------------------------------

fator_gaussiano <- function(r){

  vals <- values(r, mat = FALSE)

  vals <- vals[!is.na(vals)]

  if(length(vals) == 0){

    return(0)

  }

  mu <- mean(vals)

  sigma <- sd(vals)

  if(mu == 0){

    return(0)

  }

  gf <- sigma / mu

  return(gf)

}

# Calcula fatores
fatores <- sapply(
  rasters_gauss,
  fator_gaussiano
)

# Pesos normalizados
pesos_gauss <- fatores / sum(fatores)

# ---------------------------------------------------------
# TABELA DE PESOS
# ---------------------------------------------------------

cat("\n====== PESOS CALCULADOS ======\n")

tabela_pesos <- data.frame(
  Variavel = names(rasters_gauss),
  Fator = round(fatores, 4),
  Peso = round(pesos_gauss, 4)
)

print(tabela_pesos)

# ---------------------------------------------------------
# COMBINAÇÃO LINEAR PONDERADA
# ---------------------------------------------------------

cat("\n====== COMBINAÇÃO LINEAR ======\n")

# Raster vazio
r_combinado <- rast(
  rasters_gauss[[1]]
)

values(r_combinado) <- 0

# Soma ponderada
for(i in seq_along(rasters_gauss)){

  r_temp <- rasters_gauss[[i]][[1]]

  r_combinado <- r_combinado +
    (r_temp * pesos_gauss[i])

}

# ---------------------------------------------------------
# NORMALIZA RESULTADO FINAL
# ---------------------------------------------------------

min_final <- global(
  r_combinado,
  "min",
  na.rm = TRUE
)[1,1]

max_final <- global(
  r_combinado,
  "max",
  na.rm = TRUE
)[1,1]

r_combinado <- (
  r_combinado - min_final
) / (
  max_final - min_final
)

# ---------------------------------------------------------
# VERIFICA RESULTADO FINAL
# ---------------------------------------------------------

cat("\n====== RESULTADO FINAL ======\n")

print(r_combinado)

cat(
  "\nNúmero de bandas:",
  nlyr(r_combinado),
  "\n"
)

# ---------------------------------------------------------
# EXPORTAR RASTER FINAL
# ---------------------------------------------------------

writeRaster(
  r_combinado,
  "C:/Usuario/Seu_Usuario/Desktop/Mapas/Mapa_combinado_Gaussiano_corrigido.tif",
  overwrite = TRUE
)

# ---------------------------------------------------------
# CLASSIFICAÇÃO EM 5 CLASSES
# ---------------------------------------------------------

classes <- classify(
  r_combinado,
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
  "C:/Usuario/Seu_Usuario/Desktop/Mapas/Classes_combinadas.tif",
  overwrite = TRUE
)

# ---------------------------------------------------------
# PLOT MAPA CONTÍNUO
# ---------------------------------------------------------

plot(
  r_combinado,
  main = "Suscetibilidade Combinada"
)

# ---------------------------------------------------------
# PLOT MAPA CLASSIFICADO
# ---------------------------------------------------------

plot(
  classes,
  main = "Classes de Suscetibilidade"
)

# ---------------------------------------------------------
# FINALIZAÇÃO
# ---------------------------------------------------------

cat("\n=====================================\n")

cat("PROCESSAMENTO CONCLUÍDO COM SUCESSO\n")

cat("Raster contínuo salvo em:\n")

cat(
  "combinacao_ahp_gaussiana_corrigida.tif\n"
)

cat("\nRaster classificado salvo em:\n")

cat(
  "classes_combinadas.tif\n"
)

cat("=====================================\n")
