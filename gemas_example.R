# Code to reproduce the results of Sipilä et al. (2023).

# packages
library("SpatialBSS")
library("StatDA")
library("ggmap")
library("sf")
library("sp")
library("gstat")
library("tidyr")
library("xtable")
library("ggplot2")
library("kernelshap")
devtools::install_github("mikasip/NonlinearBSS")
library("NonlinearBSS")

# functions
contras_mat <- function(x) {
    c_mat <- matrix(0, nrow = ncol(x), ncol = ncol(x) - 1)
    for (i in 1:ncol(c_mat)) {
        c_mat[1:i, i] <- 1 / i
        c_mat[i + 1, i] <- (-1)
        c_mat[, i] <- c_mat[, i] * sqrt(i / (i + 1))
    }
    return(c_mat)
}

clr_func <- function(x) {
    x_clr <- apply(x, 1, function(row) row / prod(row)^(1 / length(row)))
    x_clr <- log(x_clr)
    return(t(x_clr))
}

plot_map <- function(coords, variable, map, quant = TRUE) {
    x <- data.frame(x = coords[, 1], y = coords[, 2], variable = variable)

    if (quant) {
        q <- quantile(variable, c(0, 0.05, 0.25, 0.75, 0.95, 1))
        leg.round <- 2
        leg.wid <- 4
        leg.just <- "right"
        leg <- rep(NA, length(q) - 1)
        leg[1] <- paste(
            "  ", roundpretty(q[1], leg.round),
            "-", format(roundpretty(q[1 + 1], leg.round),
                width = leg.wid, justify = leg.just
            )
        )
        for (i in 2:(length(q) - 1)) {
            leg[i] <- paste(
                ">", roundpretty(q[i], leg.round),
                "-", format(roundpretty(q[i + 1], leg.round),
                    width = leg.wid, justify = leg.just
                )
            )
        }

        x$variable <- cut(variable,
            breaks = q,
            include.lowest = TRUE
        )

        shapes <- rev(c(17, 17, 20, 15, 15))
        shape_size <- rev(c(1.7, 1.0, 0.3, 1.2, 1.8))
        shape_size <- 1.6 * rev(c(1.3, 0.6, 0.2, 0.9, 1.5))
        color <- c("royalblue3", "royalblue1", "black", "brown1", "brown3")

        g <- ggmap::ggmap(map,
            base_layer = ggplot(aes(
                x = x, y = y, color = variable, shape = variable,
                size = variable
            ), data = x)
        ) +
            geom_point(alpha = 1) +
            theme_bw() +
            xlab(expression("Longitude" * ~ degree * W)) +
            ylab(expression("Latitude" * ~ degree * N)) +
            theme(legend.position = "bottom") +
            scale_shape_manual(
                values = shapes,
                labels = leg
            ) +
            scale_size_manual(
                values = shape_size,
                labels = leg
            ) +
            scale_color_manual(
                values = color,
                labels = leg
            ) +
            guides(shape = guide_legend(nrow = 2, byrow = TRUE)) +
            theme(legend.title = element_blank()) +
            theme(legend.text = element_text(size = 11))
    } else {
        g <- ggmap::ggmap(map,
            base_layer = ggplot(aes(x = x, y = y, color = variable), data = x)
        ) +
            geom_point(alpha = 1) +
            theme_bw() +
            xlab(expression("Longitude" * ~ degree * W)) +
            ylab(expression("Latitude" * ~ degree * N)) +
            theme(legend.position = "bottom") +
            # scale_colour_gradientn(colours =rainbow(7))	+
            # scale_colour_gradientn(low = "grey", high = "brown") +
            scale_colour_gradientn(colors = c("turquoise", "yellow", "red", "black")) +
            theme(legend.title = element_blank()) +
            theme(legend.text = element_text(size = 11))
    }
    plot(g)
    return(g)
}

# data
library("robCompositions")
data("gemas")
head(gemas)
gemas <- gemas[-which(gemas$Ycoord == min(gemas$Ycoord)), ]
coords <- as.matrix(gemas[, c("Xcoord", "Ycoord")])
row.names(coords) <- NULL

coords_sf <- sf::st_as_sf(data.frame(coords),
    coords = c(1, 2)
)
st_crs(coords_sf) <- "EPSG:3035"
coords_sf <- st_transform(coords_sf, "+proj=longlat +datum=WGS84")
coords_ll <- st_coordinates(coords_sf, "+proj=longlat +datum=WGS84")
head(coords_ll)

aux_grid_start <- c(min(coords[, 1]), min(coords[, 2]))
size <- 100000
aux_grid_x <- seq(aux_grid_start[1] - size * 5, max(coords[, 1]) + size * 5, by = size)
aux_grid_y <- seq(aux_grid_start[2] - size * 5, max(coords[, 2]) + size * 5, by = size)
grid_coords <- data.frame(matrix(NA, ncol = 4))
colnames(grid_coords) <- c("x", "y", "order", "grid_ind")
ind <- 1
grid_ind <- 1
for (i in 1:(length(aux_grid_x) - 1)) {
    for (j in 1:(length(aux_grid_y) - 1)) {
        grid_coords[ind, ] <- c(aux_grid_x[i], aux_grid_y[j], 1, grid_ind)
        grid_coords[ind + 1, ] <- c(aux_grid_x[i], aux_grid_y[j + 1], 2, grid_ind)
        grid_coords[ind + 2, ] <- c(aux_grid_x[i + 1], aux_grid_y[j + 1], 3, grid_ind)
        grid_coords[ind + 3, ] <- c(aux_grid_x[i + 1], aux_grid_y[j], 4, grid_ind)
        grid_ind <- grid_ind + 1
        ind <- ind + 4
    }
}
grid_coords_sf <- st_as_sf(data.frame(grid_coords),
    coords = c(1, 2)
)
st_crs(grid_coords_sf) <- "EPSG:3035"
grid_coords_sf <- st_transform(grid_coords_sf, "+proj=longlat +datum=WGS84")
grid_coords_ll <- st_coordinates(grid_coords_sf, "+proj=longlat +datum=WGS84")
head(grid_coords_ll)
grid_coords$lon <- grid_coords_ll[, 1]
grid_coords$lat <- grid_coords_ll[, 2]

# register_stadiamaps("YOUR_API_KEY", write = TRUE)
europe_map <- get_stadiamap(c(left = min(coords_ll[, 1]) - 0.5, bottom = min(coords_ll[, 2]) - 0.5, right = max(coords_ll[, 1]) + 0.5, top = max(coords_ll[, 2]) + 0.5),
    maptype = "stamen_terrain_background", zoom = 4
)

g <- ggmap(europe_map,
    base_layer = ggplot(
        data = data.frame(coords_ll),
        aes(y = Y, x = X)
    )
) + geom_point(color = "darkred") +
    geom_polygon(aes(lon, lat, group = grid_ind), alpha = 0, linewidth = 0.1, colour = "black", data = grid_coords)
plot(g)

vars <- c("Al", "Ba", "Ca", "Cr", "Fe", "K", "Mg", "Mn", "Na", "Nb", "P", "Si", "Sr", "Ti", "V", "Y", "Zn", "Zr")

colnames(gemas)
gemas_dat <- as.matrix(gemas[, vars])
gemas_df <- data.frame(gemas_dat)
gemas_dat <- as.matrix(gemas_df)
row.names(gemas_dat) <- NULL

clr_dat <- clr_func(gemas_dat)
c_mat <- contras_mat(gemas_dat)

data_ilr <- log(gemas_dat) %*% c_mat

p <- 17
seed <- 11092023
resiVAE <- iVAE_spatial(data_ilr, coords, c(100000, 100000), c(1, 1), p, epochs = 1000, seed = seed, batch_size = 64)
resSBSS <- sbss(as.matrix(data_ilr), as.matrix(coords), "ring", c(0, 167000))
resSNSS <- snss_sjd(as.matrix(data_ilr), as.matrix(coords), n_block = 2, "ring", c(0, 167000))

select_points_sparsely <- function(coords, n) {
    distance_matrix <- as.matrix(dist(coords))
    selected_inds <- numeric(n)
    selected_inds[1] <- 1
    for (i in 2:(n)) {
        dists_unselected_to_selected <- as.matrix(distance_matrix[-selected_inds, selected_inds])
        row_min_dists <- apply(dists_unselected_to_selected, 1, min)
        max_dist_ind <- which(row_min_dists == max(row_min_dists))
        selected_inds[i] <- as.numeric(names(max_dist_ind)[1])
    }
    return(list(coords = coords[selected_inds, ], inds = selected_inds))
}
# Interpretations of ICs

X <- cbind(clr_dat, resiVAE$aux_data)
X <- as.data.frame(X)
bg_points <- select_points_sparsely(coords, 200)

g_bg <- ggmap(europe_map,
    base_layer = ggplot(
        data = data.frame(coords_ll[bg_points$inds, ]),
        aes(y = Y, x = X)
    )
) + geom_point(color = "darkred") + xlab("Longitude") + ylab("Latitude")
plot(g_bg)

bg_X <- X[bg_points$inds, ]
explainer <- kernelshap(resiVAE, X,
    bg_X = bg_X, orig_dim = 18, pred_fun = function(object, X, orig_dim) {
        ilr_dat <- as.matrix(X[, 1:18]) %*% c_mat
        predict(object, newdata = as.matrix(ilr_dat), aux_data = as.matrix(X[, (orig_dim + 1):ncol(X)]))
    }, feature_names = colnames(X)[1:18]
)
orig_dim <- 17

mashap <- data.frame(matrix(NA, ncol = 17, nrow = 18))
rownames(mashap) <- sapply(colnames(gemas_dat), function(name) paste0("clr(", name, ")"))
colnames(mashap) <- sapply(1:17, function(i) paste0("IC", i))

i <- 1
for (l in explainer$S) {
    mashap[, i] <- apply(abs(l), 2, mean)
    i <- i + 1
}
mashap_scaled <- sweep(mashap, 2, colSums(mashap), "/")

X <- as.data.frame(resiVAE$IC)
bg_X <- X[bg_points$inds, ]
explainer2 <- kernelshap(resiVAE, X, bg_X = bg_X, pred_fun = function(object, X) {
    pred <- predict(object, newdata = as.matrix(X), IC_to_data = TRUE)
    return(as.matrix(pred) %*% t(c_mat))
})

explainer2$baseline
mashap2 <- data.frame(matrix(NA, ncol = 17, nrow = 18))
rownames(mashap2) <- sapply(colnames(gemas_dat), function(name) paste0("clr(", name, ")"))
colnames(mashap2) <- sapply(1:17, function(i) paste0("IC", i))
i <- 1
for (l in explainer2$S) {
    mashap2[i, ] <- apply(abs(l), 2, mean)
    i <- i + 1
}
mashap_scaled2 <- sweep(mashap2, 1, rowSums(mashap2), "/")
mashap_scaled2
importance_vals <- colMeans(mashap_scaled2)
threshold <- -sort(-importance_vals)[13]
important_inds <- which(importance_vals >= threshold)
xtable(mashap_scaled2, digits = 3)
xtable(data.frame(t(colMeans(mashap_scaled2))), digits = 3)
importance_order <- order(importance_vals, decreasing = TRUE)
sorted_mashap_scaled2 <- mashap_scaled2[, importance_order]
xtable(sorted_mashap_scaled2, digits = 3)
xtable(data.frame(t(colMeans(sorted_mashap_scaled2))), digits = 3)

sorted_mashap_scaled <- mashap_scaled[, importance_order]
xtable(sorted_mashap_scaled, digits = 3)

g1 <- plot_map(coords_ll, resiVAE$IC[, importance_order[1]], map = europe_map, quant = FALSE)
g2 <- plot_map(coords_ll, resiVAE$IC[, importance_order[2]], map = europe_map, quant = FALSE)

X <- as.data.frame(resSBSS$s)
X <- X
bg_X <- X[bg_points$inds, ]
explainer2_sbss <- kernelshap(resSBSS, X, bg_X = bg_X, pred_fun = function(object, X) {
    sweep(tcrossprod(as.matrix(X), object$w_inv), 2, object$x_mu, "+") %*% t(c_mat)
})

explainer2_sbss$baseline
mashap2_sbss <- data.frame(matrix(NA, ncol = 17, nrow = 18))
rownames(mashap2_sbss) <- sapply(colnames(gemas_dat), function(name) paste0("clr(", name, ")"))
colnames(mashap2_sbss) <- sapply(1:17, function(i) paste0("IC", i))
i <- 1
for (l in explainer2_sbss$S) {
    mashap2_sbss[i, ] <- apply(abs(l), 2, mean)
    i <- i + 1
}
mashap_scaled2_sbss <- sweep(mashap2_sbss, 1, rowSums(mashap2_sbss), "/")
mashap_scaled2_sbss
importance_vals_sbss <- colMeans(mashap_scaled2_sbss)
xtable(mashap_scaled2_sbss, digits = 3)
xtable(data.frame(t(colMeans(mashap_scaled2_sbss))), digits = 3)
importance_order_sbss <- order(importance_vals_sbss, decreasing = TRUE)
sorted_mashap_scaled2_sbss <- mashap_scaled2_sbss[, importance_order_sbss]
xtable(sorted_mashap_scaled2_sbss, digits = 3)
xtable(data.frame(t(colMeans(sorted_mashap_scaled2_sbss))), digits = 3)

g1_sbss <- plot_map(coords_ll, resSBSS$s[, importance_order_sbss[1]], map = europe_map, quant = FALSE)
g2_sbss <- plot_map(coords_ll, resSBSS$s[, importance_order_sbss[2]], map = europe_map, quant = FALSE)

snss_ics <- tcrossprod(data_ilr, resSNSS$w)
X <- as.data.frame(snss_ics)
X <- X
bg_X <- X[bg_points$inds, ]
explainer2_snss <- kernelshap(resSNSS, X, bg_X = bg_X, pred_fun = function(object, X) {
    tcrossprod(as.matrix(X), solve(resSNSS$w)) %*% t(c_mat)
})

explainer2_snss$baseline
mashap2_snss <- data.frame(matrix(NA, ncol = 17, nrow = 18))
rownames(mashap2_snss) <- sapply(colnames(gemas_dat), function(name) paste0("clr(", name, ")"))
colnames(mashap2_snss) <- sapply(1:17, function(i) paste0("IC", i))
i <- 1
for (l in explainer2_snss$S) {
    mashap2_snss[i, ] <- apply(abs(l), 2, mean)
    i <- i + 1
}
mashap_scaled2_snss <- sweep(mashap2_snss, 1, rowSums(mashap2_snss), "/")
mashap_scaled2_snss
importance_vals_snss <- colMeans(mashap_scaled2_snss)
xtable(mashap_scaled2_snss, digits = 3)
xtable(data.frame(t(colMeans(mashap_scaled2_snss))), digits = 3)
importance_order_snss <- order(importance_vals_snss, decreasing = TRUE)
sorted_mashap_scaled2_snss <- mashap_scaled2_snss[, importance_order_snss]
xtable(sorted_mashap_scaled2_snss, digits = 3)
xtable(data.frame(t(colMeans(sorted_mashap_scaled2_snss))), digits = 3)

g1_snss <- plot_map(coords_ll, snss_ics[, importance_order_snss[1]], map = europe_map, quant = FALSE)
g2_snss <- plot_map(coords_ll, snss_ics[, importance_order_snss[2]], map = europe_map, quant = FALSE)


# Prediction part

# Method to predict using ordinary kriging
kriging_predict <- function(train_data, train_coords, test_coords, universal = FALSE) {
    vg_df <- data.frame(cbind(train_coords, train_data))
    vg_df_test <- data.frame(cbind(test_coords))
    names(vg_df) <- c("x", "y", "var")
    names(vg_df_test) <- c("x", "y")
    coordinates(vg_df) <- ~ x + y
    coordinates(vg_df_test) <- ~ x + y
    vg_var <- variogram(var ~ 1, data = vg_df)
    vg_s <- fit.variogram(object = vg_var, model = vgm(0.1, c("Sph", "Mat", "Exp")), fit.kappa = seq(0.1, 2, by = 0.1), fit.method = 6)
    pred1 <- NULL
    if (universal) {
        pred1 <- krige(var ~ x + y,
            vg_df, vg_df_test,
            model = vg_s,
            nmax = 50,
        )
    } else {
        pred1 <- krige(var ~ 1,
            vg_df, vg_df_test,
            model = vg_s,
            nmax = 50,
        )
    }
    return(pred1$var1.pred)
}

# Method to predict using cokriging
cokriging_predict <- function(train_data, train_coords, test_coords) {
    vg_df <- data.frame(cbind(train_coords, train_data))
    vg_df_test <- data.frame(cbind(test_coords))
    var_names <- sapply(1:ncol(train_data), FUN = function(i) paste0("var", i))
    names(vg_df) <- c("x", "y", var_names)
    names(vg_df_test) <- c("x", "y")
    coordinates(vg_df) <- ~ x + y
    coordinates(vg_df_test) <- ~ x + y
    g <- gstat(NULL, id = "var1", form = var1 ~ 1, data = vg_df, nmax = 50)
    g <- gstat(g, id = "var2", form = var2 ~ 1, data = vg_df, nmax = 50)
    g <- gstat(g, id = "var3", form = var3 ~ 1, data = vg_df, nmax = 50)
    g <- gstat(g, id = "var4", form = var4 ~ 1, data = vg_df, nmax = 50)
    g <- gstat(g, id = "var5", form = var5 ~ 1, data = vg_df, nmax = 50)
    g <- gstat(g, id = "var6", form = var6 ~ 1, data = vg_df, nmax = 50)
    g <- gstat(g, id = "var7", form = var7 ~ 1, data = vg_df, nmax = 50)
    g <- gstat(g, id = "var8", form = var8 ~ 1, data = vg_df, nmax = 50)
    g <- gstat(g, id = "var9", form = var9 ~ 1, data = vg_df, nmax = 50)
    g <- gstat(g, id = "var10", form = var10 ~ 1, data = vg_df, nmax = 50)
    g <- gstat(g, id = "var11", form = var11 ~ 1, data = vg_df, nmax = 50)
    g <- gstat(g, id = "var12", form = var12 ~ 1, data = vg_df, nmax = 50)
    g <- gstat(g, id = "var13", form = var13 ~ 1, data = vg_df, nmax = 50)
    g <- gstat(g, id = "var14", form = var14 ~ 1, data = vg_df, nmax = 50)
    g <- gstat(g, id = "var15", form = var15 ~ 1, data = vg_df, nmax = 50)
    g <- gstat(g, id = "var16", form = var16 ~ 1, data = vg_df, nmax = 50)
    g <- gstat(g, id = "var17", form = var17 ~ 1, data = vg_df, nmax = 50)
    vg1 <- fit.variogram(object = variogram(var1 ~ 1, data = vg_df), model = vgm(0.1, c("Mat")), fit.method = 6)
    vg2 <- fit.variogram(object = variogram(var2 ~ 1, data = vg_df), model = vgm(0.1, c("Mat")), fit.method = 6)
    vg3 <- fit.variogram(object = variogram(var3 ~ 1, data = vg_df), model = vgm(0.1, c("Mat")), fit.method = 6)
    vg4 <- fit.variogram(object = variogram(var4 ~ 1, data = vg_df), model = vgm(0.1, c("Mat")), fit.method = 6)
    vg5 <- fit.variogram(object = variogram(var5 ~ 1, data = vg_df), model = vgm(0.1, c("Mat")), fit.method = 6)
    vg6 <- fit.variogram(object = variogram(var6 ~ 1, data = vg_df), model = vgm(0.1, c("Mat")), fit.method = 6)
    vg7 <- fit.variogram(object = variogram(var7 ~ 1, data = vg_df), model = vgm(0.1, c("Mat")), fit.method = 6)
    vg8 <- fit.variogram(object = variogram(var8 ~ 1, data = vg_df), model = vgm(0.1, c("Mat")), fit.method = 6)
    vg9 <- fit.variogram(object = variogram(var9 ~ 1, data = vg_df), model = vgm(0.1, c("Mat")), fit.method = 6)
    vg10 <- fit.variogram(object = variogram(var10 ~ 1, data = vg_df), model = vgm(0.1, c("Mat")), fit.method = 6)
    vg11 <- fit.variogram(object = variogram(var11 ~ 1, data = vg_df), model = vgm(0.1, c("Mat")), fit.method = 6)
    vg12 <- fit.variogram(object = variogram(var12 ~ 1, data = vg_df), model = vgm(0.1, c("Mat")), fit.method = 6)
    vg13 <- fit.variogram(object = variogram(var13 ~ 1, data = vg_df), model = vgm(0.1, c("Mat")), fit.method = 6)
    vg14 <- fit.variogram(object = variogram(var14 ~ 1, data = vg_df), model = vgm(0.1, c("Mat")), fit.method = 6)
    vg15 <- fit.variogram(object = variogram(var15 ~ 1, data = vg_df), model = vgm(0.1, c("Mat")), fit.method = 6)
    vg16 <- fit.variogram(object = variogram(var16 ~ 1, data = vg_df), model = vgm(0.1, c("Mat")), fit.method = 6)
    vg17 <- fit.variogram(object = variogram(var17 ~ 1, data = vg_df), model = vgm(0.1, c("Mat")), fit.method = 6)
    models <- list(vg1, vg2, vg3, vg4, vg5, vg6, vg7, vg8, vg9, vg10, vg11, vg12, vg13, vg14, vg15, vg16, vg17)
    ranges <- numeric(17)
    for (i in 1:17) {
        ranges[i] <- models[[i]]$range[1]
    }
    common_range <- mean(ranges)
    g <- gstat(g, id = "var1", model = vgm(vg1$psill[1], "Mat", common_range, nugget = 0), fill.all = T)
    v.cross <- variogram(g)
    g <- fit.lmc(v.cross, g, fit.method = 6, correct.diagonal = 1.01)
    pred_all <- predict(g, vg_df_test)
    pred_var_names <- sapply(1:17, FUN = function(i) paste0("var", i, ".pred"))
    return(as.matrix(as.data.frame(pred_all[, pred_var_names])[, 1:17]))
}

######################
## CROSSVALIDATION ###
######################

set.seed(01122023)
K <- 10
p <- ncol(data_ilr)
n <- nrow(clr_dat)
suffle_ind <- sample(1:n)
epochs <- 1000
test_data_inds <- split(suffle_ind, cut(seq_along(suffle_ind), K, labels = FALSE))
preds <- list()
res_error <- data.frame(matrix(ncol = 4, nrow = 7 * K))
names(res_error) <- c("Method", "Mean squared error", "Mean absolute error", "Fold")
ind <- 1
for (fold in 1:K) {
    print(paste0("Fold ", fold))
    test_indices <- test_data_inds[[fold]]
    n_test <- length(test_indices)
    clr_train_dat <- clr_dat[-test_indices, ]
    clr_test_dat <- clr_dat[test_indices, ]
    coords_test <- coords[test_indices, ]
    coords_train <- coords[-test_indices, ]
    data_ilr_test <- data_ilr[test_indices, ]
    data_ilr_train <- data_ilr[-test_indices, ]
    seed <- 01122023 + fold * 10
    set.seed(seed)
    co_kriging_pred <- cokriging_predict(data_ilr_train, coords_train, coords_test)
    co_kriging_pred <- as.matrix(co_kriging_pred) %*% t(c_mat)
    preds <- append(preds, list(list(method = "Cokriging", pred = co_kriging_pred, fold = fold)))
    mse_cokrig_obs <- mean(as.matrix(co_kriging_pred - clr_test_dat)^2)
    mae_cokrig_obs <- mean(as.matrix(abs(co_kriging_pred - clr_test_dat)))

    resiVAE_train <- iVAE_spatial(data_ilr_train, coords_train, c(100000, 100000), c(1, 1), p, epochs = epochs, seed = seed, error_dist_sigma = 0.01, batch_size = 64)
    simple_kriging_pred_ICs <- data.frame(matrix(NA, nrow = nrow(data_ilr_test), ncol = p))
    for (i in seq_len(p)) {
        simple_kriging_pred_ICs[, i] <- kriging_predict(resiVAE_train$IC[, i], coords_train, coords_test)
    }
    pred_obs_ic_sk <- predict(resiVAE_train, as.matrix(simple_kriging_pred_ICs), IC_to_data = TRUE)
    pred_obs_ic_sk <- as.matrix(pred_obs_ic_sk) %*% t(c_mat)
    preds <- append(preds, list(list(method = "iVAEOrdKriging", pred = pred_obs_ic_sk, fold = fold)))
    mse_krig_ic <- mean(as.matrix(pred_obs_ic_sk - clr_test_dat)^2)
    mae_krig_ic <- mean(as.matrix(abs(pred_obs_ic_sk - clr_test_dat)))

    uni_kriging_pred_ICs <- data.frame(matrix(NA, nrow = nrow(data_ilr_test), ncol = p))
    for (i in seq_len(p)) {
        uni_kriging_pred_ICs[, i] <- kriging_predict(resiVAE_train$IC[, i], coords_train, coords_test, universal = TRUE)
    }
    pred_obs_ic_sk_uni <- predict(resiVAE_train, as.matrix(uni_kriging_pred_ICs), IC_to_data = TRUE)
    pred_obs_ic_sk_uni <- as.matrix(pred_obs_ic_sk_uni) %*% t(c_mat)
    preds <- append(preds, list(list(method = "iVAEUniKriging", pred = pred_obs_ic_sk_uni, fold = fold)))
    mse_unikrig_ic <- mean(as.matrix(pred_obs_ic_sk_uni - clr_test_dat)^2)
    mae_unikrig_ic <- mean(as.matrix(abs(pred_obs_ic_sk_uni - clr_test_dat)))

    res_error[ind, ] <- c("iVAEOrdKriging", mse_krig_ic, mae_krig_ic, fold)
    ind <- ind + 1
    res_error[ind, ] <- c("iVAEUniKriging", mse_unikrig_ic, mae_unikrig_ic, fold)
    ind <- ind + 1
    res_error[ind, ] <- c("Cokriging", mse_cokrig_obs, mae_cokrig_obs, fold)
    ind <- ind + 1

    resSBSS_train <- sbss(as.matrix(data_ilr_train), as.matrix(coords_train), "ring", c(0, 167000), maxiter = 1000)

    krig_pred_ICs_sbss <- matrix(NA, nrow = n_test, ncol = p)
    for (i in 1:p) {
        IC_i <- as.matrix(resSBSS_train$s[, i], ncol = 1)
        krig_pred_ICs_sbss[, i] <- kriging_predict(IC_i, coords_train, coords_test)
    }
    krig_sbss_pred_obs <- sweep(tcrossprod(krig_pred_ICs_sbss, resSBSS_train$w_inv), 2, resSBSS_train$x_mu, "+") %*% t(c_mat)
    preds <- append(preds, list(list(method = "SBSSOrdKriging", pred = krig_sbss_pred_obs, fold = fold)))
    mse_krig_ic_sbss <- mean(as.matrix(krig_sbss_pred_obs - clr_test_dat)^2)
    mae_krig_ic_sbss <- mean(as.matrix(abs(krig_sbss_pred_obs - clr_test_dat)))

    res_error[ind, ] <- c("SBSSOrdKriging", mse_krig_ic_sbss, mae_krig_ic_sbss, fold)
    ind <- ind + 1

    unikrig_pred_ICs_sbss <- matrix(NA, nrow = n_test, ncol = p)
    for (i in 1:p) {
        IC_i <- as.matrix(resSBSS_train$s[, i], ncol = 1)
        unikrig_pred_ICs_sbss[, i] <- kriging_predict(IC_i, coords_train, coords_test, universal = TRUE)
    }
    unikrig_sbss_pred_obs <- sweep(tcrossprod(unikrig_pred_ICs_sbss, resSBSS_train$w_inv), 2, resSBSS_train$x_mu, "+") %*% t(c_mat)
    preds <- append(preds, list(list(method = "SBSSUniKriging", pred = unikrig_sbss_pred_obs, fold = fold)))
    unimse_krig_ic_sbss <- mean(as.matrix(unikrig_sbss_pred_obs - clr_test_dat)^2)
    unimae_krig_ic_sbss <- mean(as.matrix(abs(unikrig_sbss_pred_obs - clr_test_dat)))
    res_error[ind, ] <- c("SBSSUniKriging", unimse_krig_ic_sbss, unimae_krig_ic_sbss, fold)
    ind <- ind + 1

    resSNSS_train <- snss_sjd(as.matrix(data_ilr_train), as.matrix(coords_train), n_block = 2, "ring", c(0, 167000), maxiter = 1000)
    snss_ics <- tcrossprod(data_ilr_train, resSNSS_train$w)
    krig_pred_ICs_snss <- matrix(NA, nrow = n_test, ncol = p)
    for (i in 1:p) {
        IC_i <- as.matrix(snss_ics[, i], ncol = 1)
        krig_pred_ICs_snss[, i] <- kriging_predict(IC_i, coords_train, coords_test)
    }
    krig_snss_pred_obs <- tcrossprod(krig_pred_ICs_snss, solve(resSNSS_train$w)) %*% t(c_mat)
    preds <- append(preds, list(list(method = "SNSSOrdKriging", pred = krig_snss_pred_obs, fold = fold)))
    mse_krig_ic_snss <- mean(as.matrix(krig_snss_pred_obs - clr_test_dat)^2)
    mae_krig_ic_snss <- mean(as.matrix(abs(krig_snss_pred_obs - clr_test_dat)))
    res_error[ind, ] <- c("SNSSOrdKriging", mse_krig_ic_snss, mae_krig_ic_snss, fold)
    ind <- ind + 1

    unikrig_pred_ICs_snss <- matrix(NA, nrow = n_test, ncol = p)
    for (i in 1:p) {
        IC_i <- as.matrix(snss_ics[, i], ncol = 1)
        unikrig_pred_ICs_snss[, i] <- kriging_predict(IC_i, coords_train, coords_test, universal = TRUE)
    }
    unikrig_snss_pred_obs <- tcrossprod(unikrig_pred_ICs_snss, solve(resSNSS_train$w)) %*% t(c_mat)
    preds <- append(preds, list(list(method = "SNSSUniKriging", pred = unikrig_snss_pred_obs, fold = fold)))
    mse_unikrig_ic_snss <- mean(as.matrix(unikrig_snss_pred_obs - clr_test_dat)^2)
    mae_unikrig_ic_snss <- mean(as.matrix(abs(unikrig_snss_pred_obs - clr_test_dat)))
    res_error[ind, ] <- c("SNSSUniKriging", mse_unikrig_ic_snss, mae_unikrig_ic_snss, fold)
    ind <- ind + 1
}
save(preds, file = "res_preds.RData")
save(res_error, file = "res_error.RData")

names(res_error) <- c("Method", "MSE", "MAE", "Fold")
res_error$RMSE <- sqrt(as.numeric(res_error$MSE))
ag1 <- aggregate(as.numeric(MSE) ~ Method, data = res_error, mean)
ag2 <- aggregate(as.numeric(MAE) ~ Method, data = res_error, mean)
ag3 <- aggregate(as.numeric(RMSE) ~ Method, data = res_error, mean)
ag <- cbind(ag1, ag2[, 2], ag3[, 2])
names(ag) <- c("Method", "MSE", "MAE", "RMSE")

xtable(ag, digits = 4)
