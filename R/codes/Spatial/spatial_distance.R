dist1 <- function(pt1, pt2) {
  ((2 * asin(sqrt((sin((pt1[1] - pt2[, 1]) / 2))^2 + cos(pt1[1]) * cos(pt2[, 1]) *
                    (sin((pt1[2] - pt2[, 2]) / 2))^2))) * 1.852 * (180 * 60) / pi)
}

#pt1 c'est une latitude et longitude
# pt2 ca peut être la même chose ou un vecteur de points
# Defining a matrix with distances between each individual
# Row 1 corresponds to the dstance between child 1 and all other children

dist_matrix = matrix(NA, nrow = nrow(Vdist_Data), ncol = nrow(Vdist_Data))

for(i in seq(nrow(Vdist_Data))) {
  x1 <- unlist(Vdist_Data[i, c("Lat_Rad", "Long_Rad")])
  x2 <- Vdist_Data[seq(nrow(Vdist_Data)),
                   c("Lat_Rad", "Long_Rad")]
  
  dist_matrix[i,] <- dist1(x1, x2)
}

#Vdist_data c'est juste les lat et long de chaque individu