# Mendelian segregation calculations
#  AA Aa aa
#AA
#Aa
#aa
AA <- c(1, 1/2, 0,
        0, 1/2, 1,
        0,   0, 0) %>% 
  matrix(nrow = 3, byrow = TRUE)

aa <- c(0,   0, 0,
        1, 1/2, 0,
        0, 1/2, 1) %>% 
  matrix(nrow = 3, byrow = TRUE)

Aa <- c(1/2, 1/4,   0,
        1/2, 1/2, 1/2,
        0,   1/4, 1/2) %>% 
  matrix(nrow = 3, byrow = TRUE)

# Self
S  <- c(1, 1/4, 0,
        0, 1/2, 0,
        0, 1/4, 1) %>% 
  matrix(nrow = 3, byrow = TRUE)

#       AA Aa aa
f1 <- c( 0, 1, 0)

add.names <- function(x){ 
  names(x) <- c("AA","Aa", "aa")
  x}

# Crossing against an heterozygous donor is just like an extra round backcross <<<<<

# The BC1S4s --------------------------------------


bc1  <- AA %*% f1
bc2 <-  AA %*% bc1
bc3 <-  AA %*% bc2

# for AA x aa founder
bc1s4 <- ( S %*% S %*% S %*% S %*% bc1)[,1]

bc1s4

# for AA x Aa founder
bc2s4 <- ( S %*% S %*% S %*% S %*% bc2)[,1]
bc2s4

# The BC2S3s --------------------------------------
bc2s3 <- (S %*% S %*% S %*% bc2)[,1]
bc3s3 <- (S %*% S %*% S %*% bc3)[,1]

bc2s3
bc3s3

#-------------------------------------------------

add.names(bc1s4)
add.names(bc2s4)
add.names(bc2s3)
add.names(bc3s3)

x <- results[["diplo_chr1"]][["genotype_counts_df"]]
