oceania <- read.csv("Oceania-towers.csv")

library(dplyr)
library(sf)
library(spatstat)
library(leaflet)

sydney_bbox <- c(lon_min = 150.5, lon_max = 151.5, lat_min = -34.2, lat_max = -33.5)
sydney_towers <- oceania %>% 
  filter(
    LON >= sydney_bbox["lon_min"] & LON <= sydney_bbox["lon_max"],
    LAT >= sydney_bbox["lat_min"] & LAT <= sydney_bbox["lat_max"]
  ) %>% 
  filter(!is.na(LON), !is.na(LAT))  # Eliminar NA

sydney_sf <- st_as_sf(sydney_towers, coords = c("LON", "LAT"), crs = 4326)
sydney_utm <- st_transform(sydney_sf, 32756)  # UTM para Sydney (EPSG:32756)
coords <- st_coordinates(sydney_utm)
sydney_ppp <- as.ppp(coords, W = owin(
  xrange = range(coords[, "X"]),
  yrange = range(coords[, "Y"])
))

plot(sydney_ppp, main = "Distribución de torres celulares en Sydney", pch = 20, cex = 0.3, cols = "darkblue")

quadrat_test <- quadrat.test(sydney_ppp, nx = 6, ny = 6)
plot(quadrat_test, main = "Prueba Quadrat para Aleatoriedad Espacial")
print(quadrat_test)

set.seed(123)
K_env <- envelope(sydney_ppp, Kest, nsim = 20, correction = "Ripley")
plot(K_env, main = "Función K con bandas de confianza de 20 simulaciones")

set.seed(123)
sydney_sample <- sydney_towers %>% 
  slice_sample(n = 1000)

leaflet(sydney_sample) %>%
  addTiles() %>%
  addCircleMarkers(
    lng = ~LON, 
    lat = ~LAT,
    radius = 2,
    color = "red",
    popup = ~paste("Operador:", Network, "<br>MCC:", MCC, "| MNC:", MNC)
  ) %>%
  addControl("Torres celulares en Sydney (muestra de 1000)", position = "topright")

# 1. Kernel Smoothing
h_opt <- bw.diggle(sydney_ppp)
lambda_kernel <- density(sydney_ppp, sigma = h_opt)
plot(lambda_kernel, main = "Intensidad por kernel")

# 2. Modelo Log-Lineal (distancia al centro)
centro <- centroid.owin(sydney_ppp$window)

dist_center <- function(x, y) {
  sqrt((x - centro$x)^2 + (y - centro$y)^2)
}

# Ajustar el modelo de Poisson inhomogéneo con covariable espacial
model <- ppm(sydney_ppp ~ dist_center)

summary(model)
plot(predict(model, type = "trend"), main = "Intensidad estimada por Log-Lineal dependiente de la distancia al centro")

grid <- as.data.frame(expand.grid(
  x = seq(sydney_ppp$window$xrange[1], sydney_ppp$window$xrange[2], length.out = 100),
  y = seq(sydney_ppp$window$yrange[1], sydney_ppp$window$yrange[2], length.out = 100)
)
grid$dist_center <- dist_center(grid$x, grid$y)
grid$lambda <- predict(model, locations = grid, type = "trend")

library(ggplot2)
ggplot(grid, aes(x, y, fill = lambda)) +
  geom_tile() +
  scale_fill_viridis_c(name = "λ(x,y)") +
  labs(title = "Intensidad predicha (log-lineal)") +
  coord_equal()