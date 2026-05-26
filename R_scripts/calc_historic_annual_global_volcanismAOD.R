library(ncdf4)
library(dplyr)

nc_file <- "/Users/tylerbagwell/Downloads/eVolv2k_v3_EVA_AOD_-500_1900_1.nc"
nc_file <- "/Users/tylerbagwell/Downloads/SAOD_merged_EVA_eVolv2k_v3p1_CMIP6_v3_stitched_at_1900.nc"

nc <- nc_open(nc_file)

aod550 <- ncvar_get(nc, "aod550")  # [lat, time]
lat <- ncvar_get(nc, "lat")
time <- ncvar_get(nc, "time")

nc_close(nc)

# Area weights for latitude bands
w_lat <- cos(lat * pi / 180)
w_lat <- w_lat / sum(w_lat, na.rm = TRUE)

# Global weighted mean at each time step
aod_global_time <- as.numeric(crossprod(w_lat, aod550))

# Convert fractional years to integer years
year <- floor(time)

# Annual global mean
df_aod_annual <- data.frame(
  year = year,
  aod550_global = aod_global_time
) %>%
  group_by(year) %>%
  summarise(
    aod550_global = mean(aod550_global, na.rm = TRUE),
    .groups = "drop"
  )

summary(df_aod_annual)


### plot
head(df_aod_annual)
tail(df_aod_annual)

plot(
  df_aod_annual$year,
  df_aod_annual$aod550_global,
  type = "l",
  col = "black",
  lty = 1,
  lwd = 1.0,
  xlab = "Year",
  ylab = "Global annual mean AOD550",
  xlim = c(0,2000)
)

lines(
  df$year,
  df$V,
  col = "red",
  lwd = 1.5,
  lty = 2
)

legend(
  "topleft",
  legend = c("lat-weighted SAOD_merged_EVA_eVolv2k_v3p1_CMIP6_v3_stitched_at_1900", "Barboza 2019"),
  col = c("black", "red"),
  lty = c(1, 2),
  bty = "n"
)




?rgamma
as.integer(rgamma(10, shape = 1.08, rate = 19.5))
