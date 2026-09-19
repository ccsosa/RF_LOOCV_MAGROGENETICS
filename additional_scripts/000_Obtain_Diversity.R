procesar_especie <- function(file_path, sheet_name, out_dir) {
  tab <- readxl::read_xlsx(file_path, sheet = sheet_name)
  colnames(tab) <- gsub("[.-]", "_", colnames(tab))
  
  # Identificar dinámicamente columnas de loci (excluyendo metadatos conocidos)
  meta_cols <- c("specie", "Population", "Locality", "Latitude", "Longitude", "Individual_Collection_number", "site")
  mic <- setdiff(colnames(tab), meta_cols)
  
  # 1. Objeto genind y estadísticas por sitio (si existen coordenadas)
  if (all(c("Longitude", "Latitude") %in% colnames(tab))) {
    tab$site <- paste0(tab$Longitude, ";", tab$Latitude)
    
    tabAA <- adegenet::df2genind(tab[, mic], sep = "/", loc.names = mic, 
                                 pop = factor(tab$site), NA.char = "0")
    
    hfstat <- hierfstat::genind2hierfstat(tabAA)
    bs <- hierfstat::basic.stats(hfstat)
    
    het_table <- data.frame(
      site = names(colMeans(bs$Ho, na.rm = TRUE)),
      Ho   = colMeans(bs$Ho, na.rm = TRUE),
      He   = colMeans(bs$Hs, na.rm = TRUE)
    )
    
    # Agregar metadatos agregados por sitio
    het_table$Longitude <- as.numeric(sapply(strsplit(het_table$site, ";"), `[`, 1))
    het_table$Latitude  <- as.numeric(sapply(strsplit(het_table$site, ";"), `[`, 2))
    
    write.csv(het_table, file.path(out_dir, paste0(sheet_name, "_site.csv")), row.names = FALSE)
  }
  
  # 2. Convertir a inbreedR y calcular sMLH individual
  geno_split <- lapply(mic, function(l) {
    alleles <- do.call(rbind, strsplit(as.character(tab[[l]]), "/"))
    colnames(alleles) <- paste0(l, c("_a1", "_a2"))
    as.data.frame(alleles, stringsAsFactors = FALSE)
  })
  
  genotypes <- do.call(cbind, geno_split)
  genotypes[genotypes == "0"] <- NA
  genotypes_conv <- inbreedR::convert_raw(genotypes)
  # table(table(tab$site))  # si la mayoría de sitios tiene freq=1, confirmado
  tab$MLH  <- inbreedR::MLH(genotypes_conv)
  tab$sMLH <- inbreedR::sMLH(genotypes_conv)
  
  write.csv(tab, file.path(out_dir, paste0(sheet_name, "_individual.csv")), row.names = FALSE)
}

dir <- "D:/PROGRAMAS/Dropbox/TESIS_JORGE/Datos Genéticos"
file <- "Hoja_Maestra_Datos_Geneticos__CS.xlsx"
# Ejemplo de uso:
file_path <- paste0(dir, "/", file)
procesar_especie(file_path, " acutus_SSR", dir)
procesar_especie(file_path, "intermedius_SSR", dir)
procesar_especie(file_path, "moreletti_SSR", dir)
