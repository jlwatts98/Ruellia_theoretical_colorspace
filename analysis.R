##### Analysis for Ruellia theoretical colorspace Paper #####


##### Load in and Process Data #####
# load libraries
library(readxl)
library(tidyr)
library(ggplot2)
library(patchwork)
library(corrplot)
library(psych)
library(ggborderline)
library(ape)
library(ggbiplot)
library(ggtree)
library(pavo)
library(readr)



# xl files
antho_db = read_excel("color data/color data/antho data/antho_db.xlsx") |>
    dplyr::filter(tissue == "corolla") |>
    dplyr::rename(
        Pelargonidin = `Mass of Antho per Tissue [mg/g]_Pelargonidin`,
        Cyanidin = `Mass of Antho per Tissue [mg/g]_Cyanidin`,
        Peonidin = `Mass of Antho per Tissue [mg/g]_Peonidin`,
        Delphinidin = `Mass of Antho per Tissue [mg/g]_Delphinidin`,
        Petunidin = `Mass of Antho per Tissue [mg/g]_Petunidin`,
        Malvidin = `Mass of Antho per Tissue [mg/g]_Malvidin`
    )

# make duplicates unique
antho_db$hplc_id1 = make.unique(antho_db$hplc_id)
duplicated(antho_db$hplc_id1)

## remove specimens with all zeros. These specimens do not have any antho cyanins
antho_db = antho_db |>
    dplyr::filter(
        Pelargonidin > 0 |
            Cyanidin > 0 |
            Peonidin > 0 |
            Petunidin > 0 |
            Malvidin > 0 |
            Delphinidin > 0)

spec_db <- read_excel("color data/color data/spec data/spec_db.xlsx")
load("color data/color data/spec data/03-mean.specs.RData")

# log normalize and scale anthocyanin concentrations
# log the distribution + 1, to keep 0s as is [log(1) = 0]
antho_db$Pelargonidin_log = log1p(antho_db$Pelargonidin)
antho_db$Cyanidin_log = log1p(antho_db$Cyanidin)
antho_db$Peonidin_log = log1p(antho_db$Peonidin)
antho_db$Petunidin_log = log1p(antho_db$Petunidin)
antho_db$Delphinidin_log = log1p(antho_db$Delphinidin)
antho_db$Malvidin_log = log1p(antho_db$Malvidin)

# scale so mean is 0 and sd is 1 and predictor distributions are directly comparable
antho_db$Pelargonidin_standard = scale(antho_db$Pelargonidin_log)
antho_db$Cyanidin_standard = scale(antho_db$Cyanidin_log)
antho_db$Peonidin_standard = scale(antho_db$Peonidin_log)
antho_db$Petunidin_standard = scale(antho_db$Petunidin_log)
antho_db$Delphinidin_standard = scale(antho_db$Delphinidin_log)
antho_db$Malvidin_standard = scale(antho_db$Malvidin_log)


##### Make a table with Collection number, species, antho id, spec id, and phylo id #####
table_df = antho_db |>
    dplyr::select(`Sample Collection #`,
                  species,
                  `North South`,
                  `East West`,
                  hplc_id,
                  spec_ID,
                  phylo_ID
                  )
table_df$species = paste0("R. ", table_df$species)
table_df$`North South` = round(as.numeric(ifelse(table_df$`North South` == "?", 
                                                 "",table_df$`North South` )), 5)
table_df$`East West` = round(as.numeric(ifelse(table_df$`East West` == "?", 
                                                 "",table_df$`East West` )), 5)

table_df$`North South` = ifelse(is.na(table_df$`North South`), "",table_df$`North South`)
table_df$`East West` = ifelse(is.na(table_df$`East West`), "",table_df$`East West`)

table_df$hplc_id = ifelse(is.na(table_df$hplc_id), "no", "yes")

table_df$spec_ID = ifelse(is.na(table_df$spec_ID), "no", "yes")

table_df$phylo_ID = ifelse(is.na(table_df$phylo_ID), "no", "yes")

# save table
write.csv(table_df, "objects/supp_table_1.csv")

# merge and filter data to only include species with anthos and corolla
antho_spec = merge(spec_db, antho_db, by = "spec_ID") 

# look at the structure of the dataframe
str(antho_spec)

# filter out the useless variables and rename them
antho_spec = antho_spec[, c(1:23, 30, 37, 44, 51, 58, 65, 73:78)]

# transpose mean.specs and rename as mean_specs
mean_specs = as.data.frame(t(mean.specs))

# assign first row as a column
colnames(mean_specs) <- mean_specs[1,]

# delete first row
mean_specs <- mean_specs[-1, ] 

# assign rownames as first column of the dataframe
mean_specs <- tibble::rownames_to_column(mean_specs, "spec_ID")

# merge the two data frames by spec_ID
antho_spec = merge(antho_spec, mean_specs, by = "spec_ID")

# remove yellow and white
antho_spec = antho_spec |>
    dplyr::filter(color != "yellow", color != "white")

# view the structure
str(antho_spec)

# add binary predictors
anthos1 = as.data.frame(antho_spec[, 24:29])

# make binary data
anthos_bin = matrix(nrow = nrow(antho_spec), ncol = 6)

for ( i in 1:6 ) {
    anthos_bin[, i] = ifelse(anthos1[, i] > 0, 1, anthos1[, i])
}

anthos_bin = as.data.frame(anthos_bin) |>
    dplyr::mutate_all(function(x) as.factor(x))
names(anthos_bin) <- names(antho_spec)[24:29]

# add suffix to each column
original_cols <- colnames(anthos_bin)
colnames(anthos_bin) <- paste(original_cols,"01",sep="_")

# add to antho_spec dataframe
antho_spec = cbind(antho_spec, anthos_bin)

# add binary branch data
antho_spec$del_branch = ifelse(antho_spec$Delphinidin_01 == 1 | 
                                   antho_spec$Petunidin_01 == 1 |
                                   antho_spec$Malvidin_01 == 1, 1, 0)
antho_spec$cya_branch = ifelse(antho_spec$Cyanidin_01 == 1 | 
                                   antho_spec$Peonidin_01 == 1, 1, 0)
antho_spec$pel_branch = ifelse(antho_spec$Pelargonidin_01 == 1, 1, 0)


# Spectra

# first, break everything up and stack the data
specs = antho_spec[, c(36:436)]

specs = stack(specs)
colnames(specs)[1] = "reflectance"
colnames(specs)[2] = "wavelength"
specs$wavelength = as.numeric(specs$wavelength) + 299

specs$sample = rep(antho_spec$spec_ID, 401)
specs$colorcode = rep(antho_spec$hex_color, 401)
specs$color = rep(antho_spec$color, 401)

# add anthocyanin concentrations to the specs data.frame
specs$Pelargonidin_standard = rep(antho_spec$Pelargonidin_standard, 401)
specs$Cyanidin_standard = rep(antho_spec$Cyanidin_standard, 401)
specs$Peonidin_standard = rep(antho_spec$Peonidin_standard, 401)
specs$Petunidin_standard = rep(antho_spec$Petunidin_standard, 401)
specs$Delphinidin_standard = rep(antho_spec$Delphinidin_standard, 401)
specs$Malvidin_standard = rep(antho_spec$Malvidin_standard, 401)

# add binary data
specs$Pelargonidin_01 = rep(antho_spec$Pelargonidin_01, 401)
specs$Cyanidin_01 = rep(antho_spec$Cyanidin_01, 401)
specs$Peonidin_01 = rep(antho_spec$Peonidin_01, 401)
specs$Petunidin_01 = rep(antho_spec$Petunidin_01, 401)
specs$Delphinidin_01 = rep(antho_spec$Delphinidin_01, 401)
specs$Malvidin_01 = rep(antho_spec$Malvidin_01, 401)

# add binary branch data
specs$del_branch = rep(antho_spec$del_branch, 401)
specs$pel_branch = rep(antho_spec$pel_branch, 401)
specs$cya_branch = rep(antho_spec$cya_branch, 401)

# look at structure
str(specs)

##### Model Predicting Reflectance Spectra from Anthocyanins - Original Submission #####

# regression with binary indicators to account for zero-inflation
lms = list()

for ( i in 300:700 ) {
    lms[[i]] = lm(reflectance ~ Pelargonidin + Cyanidin + 
                      Peonidin + Petunidin +
                      Delphinidin + Malvidin +
                      Pelargonidin_01 + Cyanidin_01 + 
                      Peonidin_01 + Petunidin_01 +
                      Delphinidin_01 + Malvidin_01, 
                      data = specs |> 
                      dplyr::filter(wavelength == i))
}

lms = lms[300:700]

summary(lms[[150]])

specs = specs[order(specs$wavelength),]

preds = list()
resids = list()
coefficients = list()


for ( i in 1:401 ) {
    preds[i] = lms[[i]][5]
    
}


for ( i in 1:401 ) {
    resids[i] = lms[[i]][2]
    
}

for ( i in 1:401 ) {
    coefficients[i] = lms[[i]][1]
    
}


preds = unlist(preds)
resids = unlist(resids)

# extract model coefficients and add to dataframe
coefficients = unlist(coefficients)

intercepts = rep(coefficients[seq.int(1L, length(coefficients), 13L)], each = 61)
Pelargonidin = rep(coefficients[seq.int(2L, length(coefficients), 13L)], each = 61)
Pelargonidin_01 = rep(coefficients[seq.int(8L, length(coefficients), 13L)], each = 61)
Cyanidin = rep(coefficients[seq.int(3L, length(coefficients), 13L)], each = 61)
Cyanidin_01 = rep(coefficients[seq.int(9L, length(coefficients), 13L)], each = 61)
Peonidin = rep(coefficients[seq.int(4L, length(coefficients), 13L)], each = 61)
Peonidin_01 = rep(coefficients[seq.int(10L, length(coefficients), 13L)], each = 61)
Petunidin = rep(coefficients[seq.int(5L, length(coefficients), 13L)], each = 61)
Petunidin_01 = rep(coefficients[seq.int(11L, length(coefficients), 13L)], each = 61)
Delphinidin = rep(coefficients[seq.int(6L, length(coefficients), 13L)], each = 61)
Delphinidin_01 = rep(coefficients[seq.int(12L, length(coefficients), 13L)], each = 61)
Malvidin = rep(coefficients[seq.int(7L, length(coefficients), 13L)], each = 61)
Malvidin_01 = rep(coefficients[seq.int(13L, length(coefficients), 13L)], each = 61)

# add wavelength and sample
wls = rep(300:700, each = 61)
samples = unique(specs$sample)
sample = rep(samples, 401)


model_preds = data.frame(preds, resids, wls, sample, intercepts,
                         Pelargonidin, Cyanidin, Peonidin, Petunidin, Delphinidin,
                         Malvidin,
                         Pelargonidin_01, Cyanidin_01, Peonidin_01, Petunidin_01, Delphinidin_01,
                         Malvidin_01) |>
    dplyr::rename(
        reflect_pred = `preds`,
        resid = resids,
        wavelength = wls
    )

model_coefficients = model_preds[seq.int(1L, nrow(model_preds), 61L), c(5:17)] |>
    dplyr::mutate_all(~(scale(.) %>% as.vector))
model_coefficients = stack(model_coefficients) |>
    dplyr::rename(
        variable = ind,
        standard_coef = values
    )
model_coefficients$wavelength = rep(300:700, 13)

coefficients_plot = ggplot(data = model_coefficients, mapping =
                               aes(x = wavelength, y = standard_coef, color = variable)) +
    geom_line() +
    labs(y = "standardized coefficent value") +
    theme_bw()
coefficients_plot

# set negative predictions to 0 and larger than 1 to 1
model_preds$reflect_pred = ifelse(
    model_preds$reflect_pred < 0, 0, 
    ifelse(model_preds$reflect_pred > 1, 1, 
           model_preds$reflect_pred))

# add colorcodes
model_preds_wide = model_preds[, c(1,3,4)] %>%
    pivot_wider(names_from = sample,
                values_from = reflect_pred)

model_preds_rspec = as.rspec(model_preds_wide)

pred_colorcodes = spec2rgb(model_preds_rspec)
pred_colorcodes

model_preds$colorcode = rep(pred_colorcodes, 401)


# plot
pred_specs_plot = ggplot(data = model_preds, mapping = aes(
    x = wavelength, y = reflect_pred, group = sample
)) +
    geom_line(aes(color = sample), linewidth = 1) +
    theme_bw() +
    scale_color_manual(values = c(pred_colorcodes)) +
    theme(legend.position = "none")
pred_specs_plot

# poly regression with binary indicators (of three branches) to account for zero-inflation
lms2 = list()

for ( i in 300:700 ) {
    lms2[[i]] = lm(reflectance ~ poly(Pelargonidin,2) + poly(Cyanidin,2) + 
                      poly(Peonidin,2) + poly(Petunidin,2) +
                      poly(Delphinidin,2) + poly(Malvidin,2) +
                      pel_branch +
                      del_branch + cya_branch, 
                      data = specs |> 
                      dplyr::filter(wavelength == i))
}

lms2 = lms2[300:700]

summary(lms2[[350]])

specs = specs[order(specs$wavelength),]

preds2 = list()
resids2 = list()
coefficients2 = list()


for ( i in 1:401 ) {
    preds2[i] = lms2[[i]][5]
    
}


for ( i in 1:401 ) {
    resids2[i] = lms2[[i]][2]
    
}

for ( i in 1:401 ) {
    coefficients2[i] = lms2[[i]][1]
    
}


preds2 = unlist(preds2)
resids2 = unlist(resids2)

# extract model coefficients2 and add to dataframe
coefficients2 = unlist(coefficients2)

intercepts = rep(coefficients2[seq.int(1L, length(coefficients2), 16L)], each = 61)
Pelargonidin = rep(coefficients2[seq.int(2L, length(coefficients2), 16L)], each = 61)
Pelargonidin2 = rep(coefficients2[seq.int(3L, length(coefficients2), 16L)], each = 61)
Cyanidin = rep(coefficients2[seq.int(4L, length(coefficients2), 16L)], each = 61)
Cyanidin2 = rep(coefficients2[seq.int(5L, length(coefficients2), 16L)], each = 61)
Peonidin = rep(coefficients2[seq.int(6L, length(coefficients2), 16L)], each = 61)
Peonidin2 = rep(coefficients2[seq.int(7L, length(coefficients2), 16L)], each = 61)
Petunidin = rep(coefficients2[seq.int(8L, length(coefficients2), 16L)], each = 61)
Petunidin2 = rep(coefficients2[seq.int(9L, length(coefficients2), 16L)], each = 61)
Delphinidin = rep(coefficients2[seq.int(10L, length(coefficients2), 16L)], each = 61)
Delphinidin2 = rep(coefficients2[seq.int(11L, length(coefficients2), 16L)], each = 61)
Malvidin = rep(coefficients2[seq.int(12L, length(coefficients2), 16L)], each = 61)
Malvidin2 = rep(coefficients2[seq.int(13L, length(coefficients2), 16L)], each = 61)
pel_branch = rep(coefficients2[seq.int(14L, length(coefficients2), 16L)], each = 61)
del_branch = rep(coefficients2[seq.int(15L, length(coefficients2), 16L)], each = 61)
cya_branch = rep(coefficients2[seq.int(16L, length(coefficients2), 16L)], each = 61)

# add wavelength and sample
wls = rep(300:700, each = 61)
samples = unique(specs$sample)
sample = rep(samples, 401)


model_preds2 = data.frame(preds2, resids2, wls, sample, intercepts,
                         Pelargonidin, Cyanidin, Peonidin, Petunidin, Delphinidin,
                         Malvidin,
                         Pelargonidin2, Cyanidin2, Peonidin2, Petunidin2, Delphinidin2,
                         Malvidin2,
                         pel_branch,
                         del_branch,
                         cya_branch
                         ) |>
    dplyr::rename(
        reflect_pred = `preds2`,
        resid = resids2,
        wavelength = wls
    )

model_coefficients2 = model_preds2[seq.int(1L, nrow(model_preds2), 61L), c(5:20)] |>
    dplyr::mutate_all(~(scale(.) %>% as.vector))
model_coefficients2 = stack(model_coefficients2) |>
    dplyr::rename(
        variable = ind,
        standard_coef = values
    )
model_coefficients2$wavelength = rep(300:700, 16)

coefficients_plot2 = ggplot(data = model_coefficients2, mapping =
                               aes(x = wavelength, y = standard_coef, color = variable)) +
    geom_line(linewidth = 1) +
    labs(y = "standardized coefficent value") +
    theme_bw()
coefficients_plot2

# set negative predictions to 0 and larger than 1 to 1
model_preds2$reflect_pred = ifelse(
    model_preds2$reflect_pred < 0, 0, 
    ifelse(model_preds2$reflect_pred > 1, 1, 
           model_preds2$reflect_pred))

# add colorcodes
model_preds2_wide = model_preds2[, c(1,3,4)] %>%
    pivot_wider(names_from = sample,
                values_from = reflect_pred)

model_preds2_rspec = as.rspec(model_preds2_wide)

pred_colorcodes = spec2rgb(model_preds2_rspec)
pred_colorcodes

model_preds2$colorcode = rep(pred_colorcodes, 401)


# plot
pred_specs_plot2 = ggplot(data = model_preds2, mapping = aes(
    x = wavelength, y = reflect_pred, group = sample
)) +
    geom_line(aes(color = sample), linewidth = 1) +
    theme_bw() +
    ylab("reflectance") +
    xlim(400, 700) +
    scale_color_manual(values = c(pred_colorcodes)) +
    theme(legend.position = "none")
pred_specs_plot2

# compare to real distribution

cc = specs$colorcode

specs_plot = ggplot(data = specs, mapping = aes(
    x = wavelength, y = reflectance, group = sample
)) +
    geom_line(aes(color = sample),linewidth = 1) +
    theme_bw() +
    xlim(400, 700) +
    scale_color_manual(values = c(cc)) +
    theme(legend.position = "none")
specs_plot

pred_specs_plot2

ggsave(filename = "objects/raw_spectra.jpg", plot = specs_plot,
       device = "jpg", width = 7, height = 5, dpi = 600)
ggsave(filename = "objects/pred_spectra.jpg", plot = pred_specs_plot2,
       device = "jpg", width = 7, height = 5, dpi = 600)


# model two is more representative of the real distribution
# linear relationships
# squared relationships
# branch present or absent

# use this model to predict a theoretical colorspace

##### model three, zero inflated lognormal distributions - Revised Submission #####

# run model (a model containing only the ABP biosynthesis branches was also included, but had lower adjusted R2 values all around)

library(dplyr)
library(glmnet)

# set seed
set.seed(123)

# Initialize list to store models
lasso_models = list()

# Loop over wavelengths
for (i in 400:700) {
    
    # Subset data
    df = specs |> dplyr::filter(wavelength == i)
    
    # Predictor matrix
    X = model.matrix(
        reflectance ~ pel_branch + del_branch + cya_branch +
            poly(Pelargonidin_standard, 2) + poly(Cyanidin_standard, 2) + poly(Peonidin_standard, 2) +
            poly(Petunidin_standard, 2) + poly(Delphinidin_standard, 2) + poly(Malvidin_standard, 2),
        data = df
    )[,-1]
    
    y = df$reflectance
    
    # Fit Lasso with 10-fold CV
    lasso_mod = cv.glmnet(X, y, alpha = 0.5, nfolds = 30, lambda = c(0, exp(seq(log(1e-4), log(10), length = 100))))
    
    # Store the model object
    lasso_models[[i]] = lasso_mod
}


lasso_models = lasso_models[400:700]



library(dplyr)

# Initialize results dataframe
lasso_metrics = data.frame(
    wavelength = numeric(),
    R2 = numeric(),
    RMSE = numeric()
)

# Loop over all stored models
for (i in 1:301) {
    
    # Subset data for this wavelength
    df = specs |> filter(wavelength == i + 399)

    # Predictor matrix
    X <- model.matrix(
            reflectance ~ pel_branch + del_branch + cya_branch +
                poly(Pelargonidin_standard, 2) + poly(Cyanidin_standard, 2) + poly(Peonidin_standard, 2) +
                poly(Petunidin_standard, 2) + poly(Delphinidin_standard, 2) + poly(Malvidin_standard, 2),
            data = df
    )[,-1]
    
    y_true <- df$reflectance
    
    # Get stored model
    lasso_mod = lasso_models[[i]]
    
    # get lam
    lam = lasso_mod$lambda.min / 10
    
    # Predict on training data
    y_pred = predict(lasso_mod, X, s = lam)
    
    # Compute metrics
    R2 = cor(y_true, y_pred)^2
    RMSE = sqrt(mean((y_true - y_pred)^2))
    
    # Append to results
    lasso_metrics = rbind(lasso_metrics, data.frame(
        wavelength = as.numeric(i + 399),
        R2 = R2,
        RMSE = RMSE
    ))
}

library(dplyr)
library(glmnet)

# Corrected CV function for one wavelength
cv_lasso_wavelength <- function(wl, df, k = 5, repeats = 10) {
    
    # subset to this wavelength
    dat = df |> filter(wavelength == as.numeric(wl))

    
    X = model.matrix(
        reflectance ~ pel_branch + del_branch + cya_branch +
            poly(Pelargonidin_standard, 2) + poly(Cyanidin_standard, 2) + poly(Peonidin_standard, 2) +
            poly(Petunidin_standard, 2) + poly(Delphinidin_standard, 2) + poly(Malvidin_standard, 2),
        data = dat
    )[,-1]
    
    y = dat$reflectance
    
    # store all metrics across folds and repeats
    R2_all = c()
    RMSE_all = c()
    
    for (r in 1:repeats) {
        # randomly assign folds
        folds = sample(rep(1:k, length.out = nrow(dat)))
        
        y_pred_all = numeric(length(y))
        
        for (f in 1:k) {
            train_idx = which(folds != f)
            test_idx  = which(folds == f)
            
            X_train = X[train_idx, ]
            y_train = y[train_idx]
            
            X_test = X[test_idx, ]
            
            # REFIT Lasso on training fold
            lasso_mod_fold = cv.glmnet(X_train, y_train, alpha = 0.5, nfolds = 5)
            
            # lam
            lam = lasso_mod$lambda.min / 10
            
            # Predict on held-out fold
            y_pred_all[test_idx] = predict(lasso_mod_fold, X_test, s = lam)
        }
        
        # compute metrics for this repeat
        R2_all = c(R2_all, cor(y, y_pred_all)^2)
        RMSE_all = c(RMSE_all, sqrt(mean((y - y_pred_all)^2)))
    }
    
    # return mean ± SD
    data.frame(
        wavelength = as.numeric(wl),
        R2_mean   = mean(R2_all),
        R2_sd     = sd(R2_all),
        RMSE_mean = mean(RMSE_all),
        RMSE_sd   = sd(RMSE_all)
    )
}

# Apply across all wavelengths
cv_results <- lapply(400:700, cv_lasso_wavelength, df = specs, k = 5, repeats = 10)

# Combine into one dataframe
cv_results_df <- bind_rows(cv_results)
head(cv_results_df)


# smooth
library(zoo)

cv_results_df$R2_smooth = rollmean(cv_results_df$R2_mean, k = 15, fill = 'extend', align = "center")
cv_results_df$R2_sd_smooth = rollmean(cv_results_df$R2_sd, k = 15, fill = 'extend', align = "center")
cv_results_df$RMSE_smooth = rollmean(cv_results_df$RMSE_mean, k = 15, fill = 'extend', align = "center")
cv_results_df$RMSE_sd_smooth = rollmean(cv_results_df$RMSE_sd, k = 15, fill = 'extend', align = "center")


# calculate confidence intervals
cv_results_df$r2conf_low = cv_results_df$R2_smooth - (1.96 * cv_results_df$R2_sd_smooth)
cv_results_df$r2conf_high = cv_results_df$R2_smooth + (1.96 * cv_results_df$R2_sd_smooth)


# calculate confidence intervals
cv_results_df$rmseconf_low = cv_results_df$RMSE_smooth - (1.96 * cv_results_df$RMSE_sd_smooth)
cv_results_df$rmseconf_high = cv_results_df$RMSE_smooth + (1.96 * cv_results_df$RMSE_sd_smooth)

# save
write.csv(cv_results_df, "objects/LASSO_CV_results.csv")

# plot
ggplot(cv_results_df, aes(x = wavelength, y = R2_smooth)) +
    geom_line(size = 1) +
    theme_classic() + 
    geom_ribbon(aes(ymin = r2conf_low, ymax = r2conf_high), 
                alpha = 0.2, color = NA, show.legend = F) + 
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "Wavelength (nm)", y = expression("Cross Validated Out-of-Sample R"^2)) +
    theme(axis.line = element_line(size = 1),
          axis.text = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 15),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 15))

#save
ggsave("objects/R2_plot.jpg", height = 6, width = 6, units = "in", dpi = 600)


ggplot(cv_results_df, aes(x = wavelength, y = RMSE_smooth)) +
    geom_line(size = 1) +
    theme_classic() + 
    geom_ribbon(aes(ymin = rmseconf_low, ymax = rmseconf_high), 
                alpha = 0.2, color = NA, show.legend = F) + 
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "Wavelength (nm)", y = "Cross Validated Out-of-Sample RMSE") +
    theme(axis.line = element_line(size = 1),
          axis.text = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 15),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 15))

#save
ggsave("objects/RMSE_plot.jpg", height = 6, width = 6, units = "in", dpi = 600)

# wavelengths corresponding to each model
wavelengths = 400:700

# extract coefficients for all wavelengths
lasso_coefs = lapply(seq_along(lasso_models), function(i) {
    mod = lasso_models[[i]]
    
    # extract coefficients at lambda.min
    coefs = coef(mod, s = "lambda.min")
    
    # convert sparse matrix to data frame
    coefs_df = as.data.frame(as.matrix(coefs))
    coefs_df$variable = rownames(coefs_df)
    coefs_df$wavelength = wavelengths[i]
    colnames(coefs_df)[1] = "Estimate"
    
    coefs_df
}) |> dplyr::bind_rows()

lasso_coefs


ggplot(lasso_coefs |> dplyr::filter(variable %in% c(
    "Pelargonidin_standard",
    "Cyanidin_standard",
    "Peonidin_standard",
    "Petunidin_standard",
    "Delphinidin_standard",
    "Malvidin_standard"
)), aes(x = wavelength, 
        y = Estimate,   # or Estimate
        color = variable,
        group = variable)) +  # group ensures one line per variable
    geom_line(size = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1) +
    labs(y = "LM Coefficient (Anthocyanin Concentration)", x = "Wavelength (nm)", color = "Anthocyanin") +
    scale_color_manual(values = c(
        "Pelargonidin_standard" = "#E41A1C",
        "Cyanidin_standard"   = "#377EB8",
        "Peonidin_standard"     = "#4DAF4A",
        "Petunidin_standard"    = "#984EA3",
        "Delphinidin_standard"  = "#F8766D",
        "Malvidin_standard"   = "#FFD700"
    ),    labels = c(
        "Pelargonidin_standard" = "Pelargonidin",
        "Cyanidin_standard"   = "Cyanidin",
        "Peonidin_standard"     = "Peonidin",
        "Petunidin_standard"    = "Petunidin",
        "Delphinidin_standard"  = "Delphinidin",
        "Malvidin_standard"   = "Malvidin"
    )) +
    theme_classic() +
    theme(axis.line = element_line(size = 1),
          axis.text = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 15),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 15))

ggsave("objects/antho_conc_coefficients.jpg", 
       width = 7.75, height = 6, units = "in", dpi = 600)

ggplot(lasso_coefs |> dplyr::filter(variable %in% c(
    "Pelargonidin_011",
    "Cyanidin_011",
    "Peonidin_011",
    "Petunidin_011",
    "Delphinidin_011",
    "Malvidin_011"
)), aes(x = wavelength, 
        y = Estimate,   # or Estimate
        color = variable,
        group = variable)) +  # group ensures one line per variable
    geom_line(size = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1) +
    labs(y = "LM Coefficient (Anthocyanin Presence)", x = "Wavelength (nm)", color = "Anthocyanin") +
    scale_color_manual(values = c(
        "Pelargonidin_011" = "#E41A1C",
        "Cyanidin_011"   = "#377EB8",
        "Peonidin_011"     = "#4DAF4A",
        "Petunidin_011"    = "#984EA3",
        "Delphinidin_011"  = "#F8766D",
        "Malvidin_011"   = "#FFD700"
    ),    labels = c(
        "Pelargonidin_011" = "Pelargonidin",
        "Cyanidin_011"   = "Cyanidin",
        "Peonidin_011"     = "Peonidin",
        "Petunidin_011"    = "Petunidin",
        "Delphinidin_011"  = "Delphinidin",
        "Malvidin_011"   = "Malvidin"
    )) +
    theme_classic() +
    theme(axis.line = element_line(size = 1),
          axis.text = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 15),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 15))

ggsave("objects/antho_pres_coefficients.jpg", 
       width = 7.75, height = 6, units = "in", dpi = 600)

##### Simpler model from which to extract coefficients for visualization - Revised Submission #####
library(glmmTMB)
lms3 = list()

for ( i in 400:700 ) {
    lms3[[i]] = glmmTMB(reflectance ~ 
                   # presence/absence of each pigment to deal with 0 inflation
                   Pelargonidin_01 + Cyanidin_01 + Peonidin_01 +
                   Petunidin_01 + Delphinidin_01 + Malvidin_01 +
                   # logged and scaled continuous predictors
                   Pelargonidin_standard + Cyanidin_standard + Peonidin_standard +
                   Petunidin_standard + Delphinidin_standard + Malvidin_standard,
               # gaussian as dependent variable is continuous
               family = gaussian,
               data = specs |> dplyr::filter(wavelength == i)
                   )
    
}

lms3 = lms3[400:700]

# predict
preds3 = lapply(1:301, function(i) {
    mod = lms3[[i]]
    d   = specs |> dplyr::filter(wavelength == i + 399)
    
    d$fit = predict(mod, newdata = d)
    d$wavelength = i + 399
    d
})

# compile into df
preds3_df = dplyr::bind_rows(preds3)

# calculate residuals
preds3_df$residual = preds3_df$reflectance - preds3_df$fit

summ_squared_resids = preds3_df |>
    dplyr::group_by(wavelength) |>
    dplyr::summarise(SSR = sum(residual^2))

ggplot(summ_squared_resids, aes(x = wavelength, y = SSR)) +
    geom_line(size = 1) +
    labs(y = "Sum of Squared Residuals", x = "Wavelength (nm)") +
    theme_classic()

library(performance)

R2_list = lapply(lms3, function(mod) {
    r2(mod)  # returns a data.frame with R2_marginal and R2_conditional
})

R2 = list()
adj_R2 = list()

for ( i in 1:301 ) {
    R2[[i]] = R2_list[[i]][[1]]
    adj_R2[[i]] = R2_list[[i]][[2]]
}

R2 = unlist(R2)
adj_R2 = unlist(adj_R2)

R2_df = data.frame(R2, adj_R2)
R2_df$wavelength = 400:700

ggplot(R2_df, aes(x = wavelength, y = adj_R2)) +
    geom_line(size = 1) +
    theme_classic()

coefficients3_df = lapply(1:301, function(i) {
    coefs = summary(lms3[[i]])$coefficients[["cond"]] |> as.data.frame()
    coefs$variable = rownames(coefs)   # keep term names
    coefs$wavelength = i + 399                    # model index
    coefs
}) |> dplyr::bind_rows()

# calculate confidence intervals
coefficients3_df$conf_low = coefficients3_df$Estimate - (1.96 * coefficients3_df$`Std. Error`)
coefficients3_df$conf_high = coefficients3_df$Estimate + (1.96 * coefficients3_df$`Std. Error`)

# create a significance variable
coefficients3_df$sig = ifelse(coefficients3_df$`Pr(>|z|)` < 0.05, "significant", "not significant")

ggplot(coefficients3_df |> dplyr::filter(variable %in% c(
    "Pelargonidin_standard",
    "Cyanidin_standard",
    "Peonidin_standard",
    "Petunidin_standard",
    "Delphinidin_standard",
    "Malvidin_standard"
)), aes(x = wavelength, 
                             y = Estimate,   # or Estimate
                             color = variable,
                             group = variable)) +  # group ensures one line per variable
    geom_line(aes(size = sig)) +
    geom_ribbon(aes(ymin = conf_low, ymax = conf_high, fill = variable), 
                alpha = 0.1, color = NA, show.legend = F) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1) +
    scale_size_manual(values = c("significant" = 1.5, "not significant" = 0.8), guide = "none") +
    labs(y = "LM Coefficient (Anthocyanin Concentration)", x = "Wavelength (nm)", color = "Anthocyanin") +
    scale_color_manual(values = c(
        "Pelargonidin_standard" = "#E41A1C",
        "Cyanidin_standard"   = "#377EB8",
        "Peonidin_standard"     = "#4DAF4A",
        "Petunidin_standard"    = "#984EA3",
        "Delphinidin_standard"  = "#F8766D",
        "Malvidin_standard"   = "#FFD700"
    ),    labels = c(
        "Pelargonidin_standard" = "Pelargonidin",
        "Cyanidin_standard"   = "Cyanidin",
        "Peonidin_standard"     = "Peonidin",
        "Petunidin_standard"    = "Petunidin",
        "Delphinidin_standard"  = "Delphinidin",
        "Malvidin_standard"   = "Malvidin"
    )) +
    scale_fill_manual(values = c(
        "Pelargonidin_standard" = "#E41A1C",
        "Cyanidin_standard"   = "#377EB8",
        "Peonidin_standard"     = "#4DAF4A",
        "Petunidin_standard"    = "#984EA3",
        "Delphinidin_standard"  = "#F8766D",
        "Malvidin_standard"   = "#FFD700"
    )) +
    theme_classic() +
    theme(axis.line = element_line(size = 1),
          axis.text = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 15),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 15))

ggsave("objects/antho_conc_coefficients.jpg", 
       width = 7.75, height = 6, units = "in", dpi = 600)

ggplot(coefficients3_df |> dplyr::filter(variable %in% c(
    "Pelargonidin_011",
    "Cyanidin_011",
    "Peonidin_011",
    "Petunidin_011",
    "Delphinidin_011",
    "Malvidin_011"
)), aes(x = wavelength, 
        y = Estimate,   # or Estimate
        color = variable,
        group = variable)) +  # group ensures one line per variable
    geom_line(aes(size = sig)) +
    geom_ribbon(aes(ymin = conf_low, ymax = conf_high, fill = variable), 
                alpha = 0.1, color = NA, show.legend = F) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1) +
    scale_size_manual(values = c("significant" = 1.5, "not significant" = 0.8), guide = "none") +
    labs(y = "LM Coefficient (Anthocyanin Presence)", x = "Wavelength (nm)", color = "Anthocyanin\nPresence") +
    scale_color_manual(values = c(
        "Pelargonidin_011" = "#E41A1C",
        "Cyanidin_011"   = "#377EB8",
        "Peonidin_011"     = "#4DAF4A",
        "Petunidin_011"    = "#984EA3",
        "Delphinidin_011"  = "#F8766D",
        "Malvidin_011"   = "#FFD700"
    ),    labels = c(
        "Pelargonidin_011" = "Pelargonidin",
        "Cyanidin_011"   = "Cyanidin",
        "Peonidin_011"     = "Peonidin",
        "Petunidin_011"    = "Petunidin",
        "Delphinidin_011"  = "Delphinidin",
        "Malvidin_011"   = "Malvidin"
    )) +
    scale_fill_manual(values = c(
        "Pelargonidin_011" = "#E41A1C",
        "Cyanidin_011"   = "#377EB8",
        "Peonidin_011"     = "#4DAF4A",
        "Petunidin_011"    = "#984EA3",
        "Delphinidin_011"  = "#F8766D",
        "Malvidin_011"   = "#FFD700"
    )) +
    theme_classic() +
    theme(axis.line = element_line(size = 1),
          axis.text = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 15),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 15))

ggsave("objects/antho_pres_coefficients.jpg", 
       width = 7.75, height = 6, units = "in", dpi = 600)

# conduct cross validation to calculate out-of-sample R2 and RMSE
library(caret)

# set seed for reproducibility
set.seed(123)

# define the CV function for one wavelength
cv_glmmTMB_wavelength <- function(wl, df, k = 5) {
    
    # subset to this wavelength
    dat <- df |> filter(wavelength == wl)
    
    # create folds for this subset
    folds <- createFolds(dat$reflectance, k = k)
    
    # loop through folds
    results <- lapply(folds, function(idx) {
        train <- dat[-idx, ]
        test  <- dat[idx, ]
        
        m <- glmmTMB(
            reflectance ~ 
                Pelargonidin_01 + Cyanidin_01 + Peonidin_01 +
                Petunidin_01 + Delphinidin_01 + Malvidin_01 +
                Pelargonidin_standard + Cyanidin_standard + Peonidin_standard +
                Petunidin_standard + Delphinidin_standard + Malvidin_standard,
            family = gaussian,
            data = train
        )
        
        preds = predict(m, newdata = test, type = "response")
        obs   = test$reflectance
        
        data.frame(
            wavelength = wl,
            R2   = cor(obs, preds)^2,
            RMSE = sqrt(mean((obs - preds)^2))
        )
    })
    
    bind_rows(results)
}

# apply across wavelengths 400–700
cv_results = lapply(400:700, cv_glmmTMB_wavelength, df = specs, k = 5)

# summarise per wavelength
cv_summary <- bind_rows(cv_results) %>%
    group_by(wavelength) %>%
    summarise(
        mean_R2 = mean(R2, na.rm = TRUE),
        sd_R2   = sd(R2, na.rm = TRUE),
        mean_RMSE = mean(RMSE, na.rm = TRUE),
        sd_RMSE   = sd(RMSE, na.rm = TRUE),
        .groups = "drop"
    )

cv_summary

# save
write.csv(cv_summary, "objects/CV_results.csv")

# bind with full model R2 df
model_diagnostics = merge(R2_df, cv_summary, by = "wavelength")

# smooth
library(zoo)

model_diagnostics$R2_smooth = rollmean(model_diagnostics$mean_R2, k = 15, fill = 'extend', align = "center")
model_diagnostics$R2_sd_smooth = rollmean(model_diagnostics$sd_R2, k = 15, fill = 'extend', align = "center")
model_diagnostics$RMSE_smooth = rollmean(model_diagnostics$mean_RMSE, k = 15, fill = 'extend', align = "center")
model_diagnostics$RMSE_sd_smooth = rollmean(model_diagnostics$sd_RMSE, k = 15, fill = 'extend', align = "center")


# calculate confidence intervals
model_diagnostics$r2conf_low = model_diagnostics$R2_smooth - (1.96 * model_diagnostics$R2_sd_smooth)
model_diagnostics$r2conf_high = model_diagnostics$R2_smooth + (1.96 * model_diagnostics$R2_sd_smooth)


# calculate confidence intervals
model_diagnostics$rmseconf_low = model_diagnostics$RMSE_smooth - (1.96 * model_diagnostics$RMSE_sd_smooth)
model_diagnostics$rmseconf_high = model_diagnostics$RMSE_smooth + (1.96 * model_diagnostics$RMSE_sd_smooth)


ggplot(model_diagnostics, aes(x = wavelength, y = R2_smooth)) +
    geom_line(size = 1) +
    theme_classic() + 
    geom_ribbon(aes(ymin = r2conf_low, ymax = r2conf_high), 
                alpha = 0.2, color = NA, show.legend = F) + 
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "Wavelength (nm)", y = expression("Cross Validated Out-of-Sample R"^2)) +
    theme(axis.line = element_line(size = 1),
          axis.text = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 15),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 15))

ggsave("objects/R2_plot.jpg", height = 6, width = 6, units = "in", dpi = 600)
    
ggplot(model_diagnostics, aes(x = wavelength, y = RMSE_smooth)) +
    geom_line(size = 1) +
    theme_classic() + 
    geom_ribbon(aes(ymin = rmseconf_low, ymax = rmseconf_high), 
                alpha = 0.2, color = NA, show.legend = F) + 
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "Wavelength (nm)", y = "Cross Validated Out-of-Sample RMSE") +
    theme(axis.line = element_line(size = 1),
          axis.text = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 15),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 15))
    
ggsave("objects/RMSE_plot.jpg", height = 6, width = 6, units = "in", dpi = 600)

##### Theoretical Morphospace - Original Submission #####
# Build a Theoretical Morphospace

# remove outliers
antho_db = antho_db |>
    dplyr::filter(Pelargonidin < 1.5 &
                  Cyanidin < .6 &
                  Peonidin < .03 &
                  Petunidin < .1 &
                  Malvidin < .6 &
                  Delphinidin < .1)

# remove specimens with all zeros. These specimens do not have any antho cyanins
antho_db = antho_db |>
    dplyr::filter(
        Pelargonidin > 0 |
            Cyanidin > 0 |
            Peonidin > 0 |
            Petunidin > 0 |
            Malvidin > 0 |
            Delphinidin > 0)

# nested for loop of 6 anthocyanins

# using quantiles
len = 3

pelarg_grid = quantile(antho_db$Pelargonidin, probs = c(0.05, 0.75, 0.95))
delph_grid = quantile(antho_db$Delphinidin, probs = c(0.05, 0.75, 0.95))
cya_grid = quantile(antho_db$Cyanidin, probs = c(0.05, 0.75, 0.95))
mal_grid = quantile(antho_db$Malvidin, probs = c(0.05, 0.75, 0.95))
peo_grid = quantile(antho_db$Peonidin, probs = c(0.05, 0.75, 0.95))
pet_grid = quantile(antho_db$Petunidin, probs = c(0.05, 0.75, 0.95))

# for loop
antho_concs = matrix(nrow = len^6, ncol = 6)

counter = 1
for ( i in delph_grid ) {
    for ( j in pelarg_grid ) {
        for ( k in cya_grid ) {
            for ( l in mal_grid ) {
                for ( m in peo_grid ) {
                    for ( n in pet_grid ) {
                        
                        
                        
                        # assign
                        delph <- i
                        pelarg <- j
                        cya <- k
                        mal <- l
                        peo <- m
                        pet <- n
                        
                        # fill in concentrations
                        antho_concs[counter, 1] = delph
                        antho_concs[counter, 2] = pelarg
                        antho_concs[counter, 3] = cya
                        antho_concs[counter, 4] = mal
                        antho_concs[counter, 5] = peo
                        antho_concs[counter, 6] = pet
                        
                        # update counter
                        counter <- counter + 1
                    }
                }
            }
        }
    }
}


# make into a dataframe
antho_concs = as.data.frame(antho_concs) |>
    dplyr::rename('Delphinidin' = V1,
                  'Pelargonidin' = V2,
                  'Cyanidin' = V3,
                  'Malvidin' = V4,
                  'Peonidin' = V5,
                  'Petunidin' = V6)
antho_concs$sim = c(1:nrow(antho_concs))

# add binary branch data
antho_concs$del_branch = ifelse(antho_concs$Delphinidin > 0 | 
                                   antho_concs$Petunidin > 0 |
                                   antho_concs$Malvidin > 0, 1, 0)
antho_concs$cya_branch = ifelse(antho_concs$Cyanidin > 0 | 
                                   antho_concs$Peonidin > 0, 1, 0)
antho_concs$pel_branch = ifelse(antho_concs$Pelargonidin > 0, 1, 0)

# remove all zeros
antho_concs = antho_concs |>
    dplyr::filter(
        Pelargonidin > 0 |
            Cyanidin > 0 |
            Peonidin > 0 |
            Petunidin > 0 |
            Malvidin > 0 |
            Delphinidin > 0)

# theoretical colorspace concentrations
# anthocyanin concentration distributions
pel_plot = ggplot(data = antho_db, mapping = aes(
    x = Pelargonidin
)) +
    geom_histogram(bins = 45) +
    xlim(-.01, (pelarg_grid[3] + (.25 * pelarg_grid[3]))) +
    geom_vline(data = as.data.frame(pelarg_grid), 
               aes(xintercept = pelarg_grid)) +
    theme_bw() 
pel_plot

cya_plot = ggplot(data = antho_db, mapping = aes(
    x = Cyanidin
)) +
    geom_histogram(bins = 45) +
    xlim(-.005, (cya_grid[3] + (.35 * cya_grid[3]))) +
    geom_vline(data = as.data.frame(cya_grid), 
               aes(xintercept = cya_grid)) +
    theme_bw() 
cya_plot

peo_plot = ggplot(data = antho_db, mapping = aes(
    x = Peonidin
)) +
    geom_histogram(bins = 45) +
    xlim(-.0005, (peo_grid[3] + (.25 * peo_grid[3]))) +
    geom_vline(data = as.data.frame(peo_grid), 
               aes(xintercept = peo_grid)) +
    theme_bw()
peo_plot

del_plot = ggplot(data = antho_db, mapping = aes(
    x = Delphinidin
)) +
    geom_histogram(bins = 45) +
    xlim(-.001, (delph_grid[3] + (.25 * delph_grid[3]))) +
    geom_vline(data = as.data.frame(delph_grid), 
               aes(xintercept = delph_grid)) +
    theme_bw() 
del_plot

pet_plot = ggplot(data = antho_db, mapping = aes(
    x = Petunidin
)) +
    geom_histogram(bins = 45) +
    xlim(-.001, (pet_grid[3] + (.25 * pet_grid[3]))) +
    geom_vline(data = as.data.frame(pet_grid), 
               aes(xintercept = pet_grid)) +
    theme_bw() 
pet_plot

mal_plot = ggplot(data = antho_db, mapping = aes(
    x = Malvidin
)) +
    geom_histogram() +
    xlim(-.005, (mal_grid[3] + (.25 * mal_grid[3]))) +
    geom_vline(data = as.data.frame(mal_grid), 
               aes(xintercept = mal_grid)) +
    theme_bw()
mal_plot

distributions = pel_plot + cya_plot + peo_plot +
    pet_plot + del_plot + mal_plot
distributions

ggsave(filename = "objects/concentrations.jpg", plot = distributions,
       device = "jpg", width = 7, height = 5, dpi = 600)

# use model and anthocyanin concentration trait space to predict colorspace

preds = matrix(nrow = 401, ncol = nrow(antho_concs))

for ( j in 1:nrow(antho_concs) ) {
    for ( i in 1:401 ) {
        preds[i,j] = predict(lms2[[i]], newdata = antho_concs[j,])
    }
}

# modify into dataframe and add wavelength, sim #, colorcode,  column
theo_preds = as.data.frame(preds)

stack_theo_preds = stack(theo_preds)

# set negative predictions to 0 - > 1 to 1
stack_theo_preds$values = ifelse(
    stack_theo_preds$values < 0, 0, 
    ifelse(stack_theo_preds$values > 1, 1, 
           stack_theo_preds$values))

theo_preds = unstack(stack_theo_preds)
theo_preds = tibble::rownames_to_column(theo_preds, "wavelength")
theo_preds$wavelength = as.numeric(theo_preds$wavelength) + 299

theo_morpho_rspec = as.rspec(theo_preds)

theo_colorcodes = spec2rgb(theo_morpho_rspec)

theo_specs = stack(theo_preds[,2:ncol(theo_preds)])
theo_specs$colorcode = rep(theo_colorcodes, each = 401)
theo_specs$wavelength = rep(300:700, length(theo_colorcodes))
theo_specs = theo_specs %>%
    dplyr::rename(
        simsamp = ind,
        reflectance = values
    )


# add anthocyanin concentrations
theo_specs$Pelargonidin = rep(antho_concs$Pelargonidin, each = 401)
theo_specs$Delphinidin = rep(antho_concs$Delphinidin, each = 401)
theo_specs$Cyanidin = rep(antho_concs$Cyanidin, each = 401)
theo_specs$Petunidin = rep(antho_concs$Petunidin, each = 401)
theo_specs$Peonidin = rep(antho_concs$Peonidin, each = 401)
theo_specs$Malvidin = rep(antho_concs$Malvidin, each = 401)

# add binary anthocyanin branch predictors
theo_specs$pel_branch = rep(antho_concs$pel_branch, each = 401)
theo_specs$del_branch = rep(antho_concs$del_branch, each = 401)
theo_specs$cya_branch = rep(antho_concs$cya_branch, each = 401)


# save
readr::write_csv(theo_specs, "objects/theo_specs.csv")

# predict colors from real distribution
# use model and anthocyanin concentration trait space to predict colorspace

preds = matrix(nrow = 401, ncol = nrow(antho_db))

# add binary branch data to antho_db
antho_db$del_branch = ifelse(antho_db$Delphinidin > 0 | 
                                    antho_db$Petunidin > 0 |
                                    antho_db$Malvidin > 0, 1, 0)
antho_db$cya_branch = ifelse(antho_db$Cyanidin > 0 | 
                                    antho_db$Peonidin > 0, 1, 0)
antho_db$pel_branch = ifelse(antho_db$Pelargonidin > 0, 1, 0)

for ( j in 1:nrow(antho_db) ) {
    for ( i in 1:401 ) {
        preds[i,j] = predict(lms2[[i]], newdata = antho_db[j,])
    }
}

# modify into dataframe and add wavelength, sim #, colorcode,  column
real_preds = as.data.frame(preds)

stack_real_preds = stack(real_preds)

# set negative predictions to 0 - > 1 to 1
stack_real_preds$values = ifelse(
    stack_real_preds$values < 0, 0, 
    ifelse(stack_real_preds$values > 1, 1, 
           stack_real_preds$values))

real_preds = unstack(stack_real_preds)
real_preds = tibble::rownames_to_column(real_preds, "wavelength")
real_preds$wavelength = as.numeric(real_preds$wavelength) + 299

real_morpho_rspec = as.rspec(real_preds)

real_colorcodes = spec2rgb(real_morpho_rspec)

real_specs = stack(real_preds[,2:ncol(real_preds)])
real_specs$colorcode = rep(real_colorcodes, each = 401)
real_specs$wavelength = rep(300:700, length(real_colorcodes))
real_specs = real_specs %>%
    dplyr::rename(
        simsamp = ind,
        reflectance = values
    )

real_specs$simsamp = rep(antho_db$hplc_id1, each = 401)


# add anthocyanin concentrations
real_specs$Pelargonidin = rep(antho_db$Pelargonidin, each = 401)
real_specs$Delphinidin = rep(antho_db$Delphinidin, each = 401)
real_specs$Cyanidin = rep(antho_db$Cyanidin, each = 401)
real_specs$Petunidin = rep(antho_db$Petunidin, each = 401)
real_specs$Peonidin = rep(antho_db$Peonidin, each = 401)
real_specs$Malvidin = rep(antho_db$Malvidin, each = 401)

# add binary anthocyanin branch predictors
real_specs$pel_branch = rep(antho_db$pel_branch, each = 401)
real_specs$del_branch = rep(antho_db$del_branch, each = 401)
real_specs$cya_branch = rep(antho_db$cya_branch, each = 401)


# save
readr::write_csv(real_specs, "objects/real_specs.csv")


##### Theoretical Colorspace V2 - Revised Submission #####
# remove specimens with all zeros. These specimens do not have any anthocyanins
antho_db = antho_db |>
    dplyr::filter(
        Pelargonidin > 0 |
            Cyanidin > 0 |
            Peonidin > 0 |
            Petunidin > 0 |
            Malvidin > 0 |
            Delphinidin > 0)


# define a low, medium, and high value, create a gridded theoretical anthocyanin space
# which includes all combinations of the low, medium, and high values
# low = min (as this is always 0 (standardized), because the actual distribution of anthocyanins is strongly zero inflated)
# medium = 0 because the anthocyanin concentrations are standardized to mean 0, sd = 1,
# high = 3, because this represents 99.7% of the data

# nested for loop of 6 anthocyanins

# using quantiles
len = 3

pelarg_grid = c(min(antho_db$Pelargonidin_standard), 0, 3)
delph_grid = c(min(antho_db$Delphinidin_standard, na.omit = T), 0, 3)
cya_grid = c(min(antho_db$Cyanidin_standard), 0, 3)
mal_grid = c(min(antho_db$Malvidin_standard), 0, 3)
peo_grid = c(min(antho_db$Peonidin_standard), 0, 3)
pet_grid = c(min(antho_db$Petunidin_standard), 0, 3)




# for loop
antho_concs = matrix(nrow = len^6, ncol = 6)

counter = 1
for ( i in delph_grid ) {
    for ( j in pelarg_grid ) {
        for ( k in cya_grid ) {
            for ( l in mal_grid ) {
                for ( m in peo_grid ) {
                    for ( n in pet_grid ) {
                        
                        
                        
                        # assign
                        delph <- i
                        pelarg <- j
                        cya <- k
                        mal <- l
                        peo <- m
                        pet <- n
                        
                        # fill in concentrations
                        antho_concs[counter, 1] = delph
                        antho_concs[counter, 2] = pelarg
                        antho_concs[counter, 3] = cya
                        antho_concs[counter, 4] = mal
                        antho_concs[counter, 5] = peo
                        antho_concs[counter, 6] = pet
                        
                        # update counter
                        counter <- counter + 1
                    }
                }
            }
        }
    }
}


# make into a dataframe
antho_concs = as.data.frame(antho_concs) |>
    dplyr::rename('Delphinidin_standard' = V1,
                  'Pelargonidin_standard' = V2,
                  'Cyanidin_standard' = V3,
                  'Malvidin_standard' = V4,
                  'Peonidin_standard' = V5,
                  'Petunidin_standard' = V6)

# remove row one - all = 0, no anthocyanins
antho_concs = antho_concs[-1,]

# add unique simulation number
antho_concs$sim = c(1:nrow(antho_concs))

# add presence absence for each anthocyanin
antho_concs$Pelargonidin_01 = ifelse(antho_concs$Pelargonidin_standard > -0.001, 1, 0)
antho_concs$Delphinidin_01 = ifelse(antho_concs$Delphinidin_standard > -0.001, 1, 0)
antho_concs$Cyanidin_01 = ifelse(antho_concs$Cyanidin_standard > -0.001, 1, 0)
antho_concs$Malvidin_01 = ifelse(antho_concs$Malvidin_standard > -0.001, 1, 0)
antho_concs$Peonidin_01 = ifelse(antho_concs$Peonidin_standard > -0.001, 1, 0)
antho_concs$Petunidin_01 = ifelse(antho_concs$Petunidin_standard > -0.001, 1, 0)

# add binary branch data
antho_concs$del_branch = ifelse(antho_concs$Delphinidin_01 > 0 | 
                                    antho_concs$Petunidin_01 > 0 |
                                    antho_concs$Malvidin_01 > 0, 1, 0)
antho_concs$cya_branch = ifelse(antho_concs$Cyanidin_01 > 0 | 
                                    antho_concs$Peonidin_01 > 0, 1, 0)
antho_concs$pel_branch = ifelse(antho_concs$Pelargonidin_01 > 0, 1, 0)

# predict, new ridge model, less overfit
# make into matrix with only necessary predictors
antho_concs_red = antho_concs |>
    dplyr::select(-sim, -Pelargonidin_01, -Delphinidin_01, -Cyanidin_01,
     -Malvidin_01, -Peonidin_01, -Petunidin_01)

antho_concs_red = as.matrix(antho_concs_red)

# make a model matrix
X_poly = model.matrix(
    ~ del_branch + cya_branch + pel_branch +
        poly(Pelargonidin_standard, 2) +
        poly(Cyanidin_standard, 2) +
        poly(Peonidin_standard, 2) +
        poly(Petunidin_standard, 2) +
        poly(Delphinidin_standard, 2) +
        poly(Malvidin_standard, 2),
    data = antho_concs
)

# glmnet already adds an intercept, so you can drop it if you like
X_poly <- X_poly[, -1]

# initiate empty matrix
preds = matrix(nrow = 301, ncol = 728)

for (i in 1:301) {
    mod = lasso_models[[i]]
    
    lam = lasso_mod$lambda.min / 10
    
    # predict all rows at once
    preds[i, ] = predict(mod, newx = X_poly, s = lam)
}

preds_df = as.data.frame(preds) |>
    dplyr::mutate(across(everything(), ~ pmin(pmax(., 0), 1)))

preds_df$wavelength = 400:700


write.csv(preds_df, file = 'objects/theoretical_colorspace.csv')

preds_df = read.csv('objects/theoretical_colorspace.csv')[, 2:729]

# put into long format and add associated anthocyanin data for each replicate
preds_stack = stack(preds_df) |>
    dplyr::rename(reflectance = values,
                  simsamp = ind) |>
    dplyr::mutate(wavelength = rep(400:700, 728))

preds_stack$Pelargonidin_standard = rep(antho_concs$Pelargonidin_standard, each = 301)
preds_stack$Cyanidin_standard = rep(antho_concs$Cyanidin_standard, each = 301)
preds_stack$Peonidin_standard = rep(antho_concs$Peonidin_standard, each = 301)
preds_stack$Petunidin_standard = rep(antho_concs$Petunidin_standard, each = 301)
preds_stack$Delphinidin_standard = rep(antho_concs$Delphinidin_standard, each = 301)
preds_stack$Malvidin_standard = rep(antho_concs$Malvidin_standard, each = 301)

preds_stack$Pelargonidin_01 = rep(antho_concs$Pelargonidin_01, each = 301)
preds_stack$Cyanidin_01 = rep(antho_concs$Cyanidin_01, each = 301)
preds_stack$Peonidin_01 = rep(antho_concs$Peonidin_01, each = 301)
preds_stack$Petunidin_01 = rep(antho_concs$Petunidin_01, each = 301)
preds_stack$Delphinidin_01 = rep(antho_concs$Delphinidin_01, each = 301)
preds_stack$Malvidin_01 = rep(antho_concs$Malvidin_01, each = 301)

preds_stack$del_branch = rep(antho_concs$del_branch, each = 301)
preds_stack$cya_branch = rep(antho_concs$cya_branch, each = 301)
preds_stack$pel_branch = rep(antho_concs$pel_branch, each = 301)

# calculate color codes

preds_df300 = as.data.frame(matrix(data = .1, nrow = 100, ncol = 728))
preds_rspec = rbind(preds_df300, preds_df)
preds_rspec = as.rspec(preds_rspec |> dplyr::mutate(wavelength = 300:700), whichwl = 'wavelength')
 
theo_colorcodes = spec2rgb(preds_rspec)
theo_colorcodes

preds_stack$colorcode = rep(theo_colorcodes, each = 301)

# save
write.csv(preds_stack, "objects/theoretical_colorspace_long.csv")

# predict the reflectance spectra of the observed flowers

# clean up predictors dataframe
antho_predictors = antho_db |>
    dplyr::select(Cyanidin_standard, Pelargonidin_standard, Delphinidin_standard,
                  Peonidin_standard, Malvidin_standard, Petunidin_standard)

antho_predictors$Pelargonidin_01 = ifelse(antho_predictors$Pelargonidin_standard > min(antho_predictors$Pelargonidin_standard), 1, 0)
antho_predictors$Delphinidin_01 = ifelse(antho_predictors$Delphinidin_standard > min(antho_predictors$Delphinidin_standard), 1, 0)
antho_predictors$Cyanidin_01 = ifelse(antho_predictors$Cyanidin_standard > min(antho_predictors$Cyanidin_standard), 1, 0)
antho_predictors$Malvidin_01 = ifelse(antho_predictors$Malvidin_standard > min(antho_predictors$Malvidin_standard), 1, 0)
antho_predictors$Peonidin_01 = ifelse(antho_predictors$Peonidin_standard > min(antho_predictors$Peonidin_standard), 1, 0)
antho_predictors$Petunidin_01 = ifelse(antho_predictors$Petunidin_standard > min(antho_predictors$Petunidin_standard), 1, 0)

# add binary branch data
antho_predictors$del_branch = ifelse(antho_predictors$Delphinidin_01 > 0 | 
                                    antho_predictors$Petunidin_01 > 0 |
                                    antho_predictors$Malvidin_01 > 0, 1, 0)
antho_predictors$cya_branch = ifelse(antho_predictors$Cyanidin_01 > 0 | 
                                    antho_predictors$Peonidin_01 > 0, 1, 0)
antho_predictors$pel_branch = ifelse(antho_predictors$Pelargonidin_01 > 0, 1, 0)

X_poly_real = model.matrix(
        ~ del_branch + cya_branch + pel_branch +
            poly(Pelargonidin_standard, 2) +
            poly(Cyanidin_standard, 2) +
            poly(Peonidin_standard, 2) +
            poly(Petunidin_standard, 2) +
            poly(Delphinidin_standard, 2) +
            poly(Malvidin_standard, 2),
        data = antho_predictors
    )

# glmnet already adds an intercept, so you can drop it if you like
X_poly_real <- X_poly_real[, -1]

# initiate empty matrix
real_preds = matrix(nrow = 301, ncol = 197)

for (i in 1:301) {
    mod = lasso_models[[i]]
    
    lam = lasso_mod$lambda.min / 10
    # predict all rows at once
    real_preds[i, ] = predict(mod, newx = X_poly_real, s = lam)
}


real_preds_df = as.data.frame(real_preds) |>
    dplyr::mutate(across(everything(), ~ pmin(pmax(., 0), 1)))

real_preds_df$wavelength = 400:700

write.csv(real_preds_df, file = 'objects/observed_colorspace.csv')

real_preds_df = read.csv('objects/observed_colorspace.csv')[, 2:198]

# put into long format and add associated anthocyanin data for each replicate
real_preds_stack = stack(real_preds_df) |>
    dplyr::rename(reflectance = values,
                  simsamp = ind) |>
    dplyr::mutate(wavelength = rep(400:700, 197))


real_preds_stack$simsamp = rep(antho_db$hplc_id1, each = 301)

real_preds_stack$Pelargonidin_standard = rep(antho_predictors$Pelargonidin_standard, each = 301)
real_preds_stack$Cyanidin_standard = rep(antho_predictors$Cyanidin_standard, each = 301)
real_preds_stack$Peonidin_standard = rep(antho_predictors$Peonidin_standard, each = 301)
real_preds_stack$Petunidin_standard = rep(antho_predictors$Petunidin_standard, each = 301)
real_preds_stack$Delphinidin_standard = rep(antho_predictors$Delphinidin_standard, each = 301)
real_preds_stack$Malvidin_standard = rep(antho_predictors$Malvidin_standard, each = 301)

real_preds_stack$Pelargonidin_01 = rep(antho_predictors$Pelargonidin_01, each = 301)
real_preds_stack$Cyanidin_01 = rep(antho_predictors$Cyanidin_01, each = 301)
real_preds_stack$Peonidin_01 = rep(antho_predictors$Peonidin_01, each = 301)
real_preds_stack$Petunidin_01 = rep(antho_predictors$Petunidin_01, each = 301)
real_preds_stack$Delphinidin_01 = rep(antho_predictors$Delphinidin_01, each = 301)
real_preds_stack$Malvidin_01 = rep(antho_predictors$Malvidin_01, each = 301)

real_preds_stack$del_branch = rep(antho_predictors$del_branch, each = 301)
real_preds_stack$cya_branch = rep(antho_predictors$cya_branch, each = 301)
real_preds_stack$pel_branch = rep(antho_predictors$pel_branch, each = 301)

# calculate color codes

real_preds_df300 = as.data.frame(matrix(data = .1, nrow = 100, ncol = 197))
real_preds_rspec = rbind(real_preds_df300, real_preds_df)
real_preds_rspec = as.rspec(real_preds_rspec |> dplyr::mutate(wavelength = 300:700), whichwl = 'wavelength')


real_colorcodes = spec2rgb(real_preds_rspec)
real_colorcodes

real_preds_stack$colorcode = rep(real_colorcodes, each = 301)

# save
write.csv(real_preds_stack, "objects/real_colorspace_long.csv")


##### Colorspace Comparisons and Trends - Original Submission #####

# read in specs
real_specs = readr::read_csv("objects/real_specs.csv")

theo_specs = readr::read_csv("objects/theo_specs.csv")


# just every 401st row
theo_anthos = theo_specs[c(seq(1,nrow(theo_specs), by = 401)),]
# just every 401st row
real_anthos = real_specs[c(seq(1,nrow(real_specs), by = 401)),]

anthos = rbind(theo_anthos, real_anthos) |>
    dplyr::mutate(
        realtheo = c(rep("theo", nrow(theo_anthos)),
                     rep("real", nrow(real_anthos)))
    ) |>
    dplyr::select(
        simsamp,
        colorcode,
        realtheo,
        Pelargonidin,
        Delphinidin,
        Cyanidin,
        Petunidin,
        Peonidin,
        Malvidin,
        pel_branch,
        del_branch,
        cya_branch
    )


# save
readr::write_csv(anthos, "objects/anthos.csv")


# calculate colorspace values from spectra

RYGB = function(
        spec_df
) {
    
    # calculate brightness
    brightness = spec_df |>
        dplyr::group_by(simsamp) |>
        dplyr::summarise(brightness = sum(reflectance), .groups = "drop")
    
    spec_df$brightness = rep(brightness$brightness, each = 301)
    
    # calculate R
    R = spec_df |>
        dplyr::group_by(simsamp) |>
        dplyr::filter(wavelength >= 625) |>
        dplyr::summarise(R = sum(reflectance) / mean(brightness))
    
    # calculate Y
    Y = spec_df |>
        dplyr::group_by(simsamp) |>
        dplyr::filter(wavelength <= 624 & wavelength >= 550) |>
        dplyr::summarise(Y = sum(reflectance) / mean(brightness))
    
    # calculate G
    G = spec_df |>
        dplyr::group_by(simsamp) |>
        dplyr::filter(wavelength <= 549 & wavelength >= 475) |>
        dplyr::summarise(G = sum(reflectance) / mean(brightness))
    
    # calculate B
    B = spec_df |>
        dplyr::group_by(simsamp) |>
        dplyr::filter(wavelength <= 474 & wavelength >= 400) |>
        dplyr::summarise(B = sum(reflectance) / mean(brightness))
    
    # final dataframe
    output = merge(brightness, R, by = "simsamp")
    output = merge(output, Y, by = "simsamp")
    output = merge(output, G, by = "simsamp")
    output = merge(output, B, by = "simsamp")
    
    return(output)
}


# from Stacey's paper: Smith, S. D. (2014). Quantifying color variation: 
# improved formulas for calculating hue with segment classification. 
# Applications in Plant Sciences, 2(3), 1300088.

calc_colorspace <- function(
        RYGB_out,
        format="degrees"
        ) {
            
            # define variables
            R = RYGB_out[, 3]
            Y = RYGB_out[, 4]
            G = RYGB_out[, 5]
            B = RYGB_out[, 6]
            
            MS=Y-B
            LM=R-G
            chroma=sqrt(LM^2+MS^2)
            unmod<-sign(MS)*acos(LM/chroma)
            hue<-unmod%%(2*pi)
            if (format=="degrees") {
                hue<-hue*(180/pi)
                output = data.frame(LM, MS, hue, chroma)
                output = cbind(output, RYGB_out)
                return(output)
            }
            else {
                output = data.frame(LM, MS, hue, chroma)
                output = cbind(output, RYGB_output)
                return(output)
            }
}

# run function to calculate color parameters
real_colors = RYGB(real_specs)

real_colors = calc_colorspace(real_colors)

theo_colors = RYGB(theo_specs)

theo_colors = calc_colorspace(theo_colors)

# append anthos data to the colors dataset
real_col_anth = merge(real_colors, real_anthos, by = "simsamp")
real_col_anth = real_col_anth[, c(1:10, 12, 14:22)]

theo_col_anth = merge(theo_colors, theo_anthos, by = "simsamp")
theo_col_anth = theo_col_anth[, c(1:10, 12, 14:22)]

col_anth = rbind(real_col_anth, theo_col_anth) |>
    dplyr::mutate(
        realtheo = c(rep("real", nrow(real_col_anth)),
                     rep("theo", nrow(theo_col_anth))
                     ))

# fix hue values (add 360 degrees to those less than 100 degrees 
# to avoid artificial separation between 0 and 360)
real_col_anth$hue = ifelse(real_col_anth$hue < 100, 
                           real_col_anth$hue + 360, 
                           real_col_anth$hue)
theo_col_anth$hue = ifelse(theo_col_anth$hue < 100, 
                           theo_col_anth$hue + 360, 
                           theo_col_anth$hue)
col_anth$hue = ifelse(col_anth$hue < 100, 
                           col_anth$hue + 360, 
                           col_anth$hue)

# add column describing which branches are present
col_anth$branches = ifelse(col_anth$del_branch == 1 &
                           col_anth$cya_branch == 0 &
                           col_anth$pel_branch == 0, 
                           "Delphinidin",
                    ifelse(col_anth$del_branch == 0 &
                           col_anth$cya_branch == 1 &
                           col_anth$pel_branch == 0, 
                           "Cyanidin",
                    ifelse(col_anth$del_branch == 0 &
                           col_anth$cya_branch == 0 &
                           col_anth$pel_branch == 1, 
                           "Pelargonidin",
                    ifelse(col_anth$del_branch == 1 &
                           col_anth$cya_branch == 1 &
                           col_anth$pel_branch == 0, 
                           "del_cya",
                    ifelse(col_anth$del_branch == 0 &
                           col_anth$cya_branch == 1 &
                           col_anth$pel_branch == 1, 
                           "cya_pel",
                    ifelse(col_anth$del_branch == 1 &
                           col_anth$cya_branch == 0 &
                           col_anth$pel_branch == 1, 
                           "del_pel",
                    ifelse(col_anth$del_branch == 1 &
                            col_anth$cya_branch == 1 &
                           col_anth$pel_branch == 1, 
                           "all", NA)))))))

unique(col_anth$branches)

# save
readr::write_csv(col_anth, "objects/col_anth.csv")

# plot color parameters in colorspace
cc = real_col_anth$colorcode

real_color_plot = ggplot(data = col_anth |> dplyr::filter(realtheo == "real"),
                         mapping = aes(x = MS, y = LM, color = simsamp)) +
    geom_point(aes(shape = branches), size = 2) +
    scale_shape_manual(values = c(8, 19, 15, 17, 10, 13, 14, 1),
                       breaks = c("all",
                                  "Pelargonidin",
                                  "Cyanidin",
                                  "Delphinidin",
                                  "cya_pel",
                                  "del_pel",
                                  "del_cya",
                                  "none"),
                       labels = c("all",
                                 "Pelargonidin",
                                 "Cyanidin",
                                 "Delphinidin",
                                 "C & P",
                                 "D & P",
                                 "D & C",
                                 "none"),
                       name = "ABP Branches") +
    theme_bw() +
    xlim(-.4, .3) +
    ylim(.1, .85) +
    scale_color_manual(values = c(cc), guide = "none") +
    theme(legend.position = "none")
real_color_plot  

ggsave(filename = "objects/real_color_plot.jpg", plot = real_color_plot,
       device = "jpg", width = 7, height = 5, dpi = 600)

# plot real spectra
real_specs_plot = ggplot(data = real_specs, mapping = aes(
    x = wavelength, y = reflectance, group = simsamp
)) +
    geom_line(aes(color = simsamp),linewidth = 1) +
    theme_bw() +
    scale_color_manual(values = c(cc)) +
    theme(legend.position = "none")
real_specs_plot

# do the same for theoretical spectra
cc = theo_col_anth$colorcode

theo_color_plot = ggplot(data = col_anth |> dplyr::filter(realtheo == "theo"),
                         mapping = aes(x = MS, y = LM, color = simsamp)) +
    geom_point(aes(shape = branches), size = 2) +
    scale_shape_manual(values = c(8, 19, 15, 17, 10, 13, 14, 1),
                       breaks = c("all",
                                  "Pelargonidin",
                                  "Cyanidin",
                                  "Delphinidin",
                                  "cya_pel",
                                  "del_pel",
                                  "del_cya",
                                  "none"),
                       labels = c("all",
                                  "Pelargonidin",
                                  "Cyanidin",
                                  "Delphinidin",
                                  "C & P",
                                  "D & P",
                                  "D & C",
                                  "none"),
                       name = "ABP Branches") +
    theme_bw() +
    xlim(-.4, .3) +
    ylim(.1, .85) +
    scale_color_manual(values = c(cc), guide = "none") +
    theme(legend.position = "right")
theo_color_plot  

ggsave(filename = "objects/theo_color_plot.jpg", plot = theo_color_plot,
       device = "jpg", width = 7, height = 5, dpi = 600)

duo_color = real_color_plot + theo_color_plot
duo_color

ggsave(filename = "objects/duo_color.jpg", plot = duo_color,
       device = "jpg", width = 7, height = 4, dpi = 600)


# plot random subset of theo spectra

# random subset
set.seed(1)
rand_subset = sample(theo_col_anth$simsamp, size = 272, replace = F)

cc = theo_col_anth |>
        dplyr::filter(simsamp %in% rand_subset)
cc = cc$colorcode

theo_specs_plot = ggplot(data = theo_specs |>
                             dplyr::filter(simsamp %in% rand_subset), 
                         mapping = aes(
    x = wavelength, y = reflectance, group = simsamp
)) +
    geom_line(aes(color = simsamp),linewidth = 1) +
    theme_bw() +
    scale_color_manual(values = c(cc)) +
    theme(legend.position = "none")
theo_specs_plot


# 3-d plot including brightness
library(scatterplot3d)
jpeg("objects/color_plot_real.jpeg", quality = 100, width = 7.5,
     height = 5, units = "in", res = 600)
scatterplot3d(real_col_anth[,c(3,2,6)], pch = 16, type="h", 
              color = real_col_anth$colorcode, angle =225, scale.y = .5)
dev.off()

jpeg("objects/color_plot_theo.jpeg", quality = 100, width = 7.5,
     height = 5, units = "in", res = 600)
scatterplot3d(theo_col_anth[,c(3,2,6)], pch = 16, type="h", 
              color = theo_col_anth$colorcode, angle = 225, scale.y = .5)
dev.off()

par(mfrow = c(1,2),
    mar = c(.1, .1,.1, .1))
scatterplot3d(theo_col_anth[,c(3,2,6)], pch = 16, type="h", 
              color = theo_col_anth$colorcode, angle = 135, scale.y = .5)
scatterplot3d(real_col_anth[,c(3,2,6)], pch = 16, type="h", 
              color = real_col_anth$colorcode, angle =135, scale.y = .5)


##### Colorspace Comparisons and Trends - Revised Submission #####

# read in data from theoretical colorspace section
rcol = read.csv('objects/real_colorspace_long.csv')[,-1]
tcol = read.csv('objects/theoretical_colorspace_long.csv')[,-1]

# rbind
colors = rbind(rcol, tcol)

# add column to indicate observed or theoretical
colors$realtheo = c(rep("real", nrow(rcol)), rep("theo", nrow(tcol)))

# extract only anthocyanin data
anthos = colors[c(seq(1,nrow(colors), by = 301)),] |>
    dplyr::select(-reflectance)

# segment color classification
colors_seg = RYGB(colors)

# get colorspace variables
color_vars = calc_colorspace(colors_seg)

# append antho data
colorspace = merge(color_vars, anthos, by = "simsamp")

# fix hue values (add 360 degrees to those less than 120 degrees 
# to avoid artificial separation between 0 and 360)
colorspace$hue = ifelse(colorspace$hue < 100, 
                           colorspace$hue + 360, 
                           colorspace$hue)

# add branches present
# add column describing which branches are present
colorspace$branches = ifelse(colorspace$del_branch == 1 &
                               colorspace$cya_branch == 0 &
                               colorspace$pel_branch == 0, 
                           "Delphinidin",
                           ifelse(colorspace$del_branch == 0 &
                                      colorspace$cya_branch == 1 &
                                      colorspace$pel_branch == 0, 
                                  "Cyanidin",
                                  ifelse(colorspace$del_branch == 0 &
                                             colorspace$cya_branch == 0 &
                                             colorspace$pel_branch == 1, 
                                         "Pelargonidin",
                                         ifelse(colorspace$del_branch == 1 &
                                                    colorspace$cya_branch == 1 &
                                                    colorspace$pel_branch == 0, 
                                                "del_cya",
                                                ifelse(colorspace$del_branch == 0 &
                                                           colorspace$cya_branch == 1 &
                                                           colorspace$pel_branch == 1, 
                                                       "cya_pel",
                                                       ifelse(colorspace$del_branch == 1 &
                                                                  colorspace$cya_branch == 0 &
                                                                  colorspace$pel_branch == 1, 
                                                              "del_pel",
                                                              ifelse(colorspace$del_branch == 1 &
                                                                         colorspace$cya_branch == 1 &
                                                                         colorspace$pel_branch == 1, 
                                                                     "all", NA)))))))

# plot
# plot color parameters in colorspace

# get colorcodes
theocc = colorspace |> dplyr::filter(realtheo == "theo") |> dplyr::pull(colorcode)
realcc = colorspace |> dplyr::filter(realtheo == "real") |> dplyr::pull(colorcode)


# define similar lims for both plots
xmin = min(colorspace$MS)
xmax = max(colorspace$MS)
ymin = min(colorspace$LM)
ymax = max(colorspace$LM)

theo_color_plot = ggplot(data = colorspace |> dplyr::filter(realtheo == "theo"),
                         mapping = aes(x = MS, y = LM, color = simsamp)) +
    geom_point(aes(shape = branches), size = 2) +
    scale_shape_manual(values = c(8, 19, 15, 17, 10, 13, 14),
                       breaks = c("all",
                                  "Pelargonidin",
                                  "Cyanidin",
                                  "Delphinidin",
                                  "cya_pel",
                                  "del_pel",
                                  "del_cya"),
                       labels = c("all",
                                  "Pelargonidin",
                                  "Cyanidin",
                                  "Delphinidin",
                                  "C & P",
                                  "D & P",
                                  "D & C"),
                       name = "ABP Branches") +
    lims(x = c(xmin, xmax), y = c(ymin, ymax)) +
    theme_classic() +
    scale_color_manual(values = c(theocc), guide = "none") +
    theme(legend.position = "right",
            axis.line = element_line(size = 1),
                axis.text = element_text(size = 12, color = "black"),
                axis.title = element_text(size = 15),
                legend.text = element_text(size = 12),
                legend.title = element_text(size = 15))
theo_color_plot 

ggsave(filename = "objects/theoretical_colorspace.jpg", plot = theo_color_plot,
       device = "jpg", width = 7.75, height = 6, dpi = 600)

real_color_plot = ggplot(data = colorspace |> dplyr::filter(realtheo == "real"),
                         mapping = aes(x = MS, y = LM, color = simsamp)) +
    geom_point(aes(shape = branches), size = 2) +
    scale_shape_manual(values = c(8, 19, 15, 17, 10, 13, 14),
                       breaks = c("all",
                                  "Pelargonidin",
                                  "Cyanidin",
                                  "Delphinidin",
                                  "cya_pel",
                                  "del_pel",
                                  "del_cya"),
                       labels = c("all",
                                  "Pelargonidin",
                                  "Cyanidin",
                                  "Delphinidin",
                                  "C & P",
                                  "D & P",
                                  "D & C"),
                       name = "ABP Branches") +
    theme_classic() +
    lims(x = c(xmin, xmax), y = c(ymin, ymax)) +
    scale_color_manual(values = c(realcc), guide = "none") +
    theme(legend.position = "right",
          axis.line = element_line(size = 1),
          axis.text = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 15),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 15))
real_color_plot

ggsave(filename = "objects/real_colorspace.jpg", plot = real_color_plot,
       device = "jpg", width = 7.75, height = 6, dpi = 600)



real_color_plot_hc = ggplot(data = colorspace |> dplyr::filter(realtheo == "real"),
                         mapping = aes(x = hue, y = chroma, color = simsamp)) +
    geom_point(aes(shape = branches), size = 2) +
    scale_shape_manual(values = c(8, 19, 15, 17, 10, 13, 14),
                       breaks = c("all",
                                  "Pelargonidin",
                                  "Cyanidin",
                                  "Delphinidin",
                                  "cya_pel",
                                  "del_pel",
                                  "del_cya"),
                       labels = c("all",
                                  "Pelargonidin",
                                  "Cyanidin",
                                  "Delphinidin",
                                  "C & P",
                                  "D & P",
                                  "D & C"),
                       name = "ABP Branches") +
    theme_classic() +
    scale_color_manual(values = c(realcc), guide = "none") +
    theme(legend.position = "right",
          axis.line = element_line(size = 1),
          axis.text = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 15),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 15))
real_color_plot_hc


theo_color_plot_hc = ggplot(data = colorspace |> dplyr::filter(realtheo == "theo"),
                            mapping = aes(x = hue, y = chroma, color = simsamp)) +
    geom_point(aes(shape = branches), size = 2) +
    scale_shape_manual(values = c(8, 19, 15, 17, 10, 13, 14),
                       breaks = c("all",
                                  "Pelargonidin",
                                  "Cyanidin",
                                  "Delphinidin",
                                  "cya_pel",
                                  "del_pel",
                                  "del_cya"),
                       labels = c("all",
                                  "Pelargonidin",
                                  "Cyanidin",
                                  "Delphinidin",
                                  "C & P",
                                  "D & P",
                                  "D & C"),
                       name = "ABP Branches") +
    theme_classic() +
    scale_color_manual(values = c(theocc), guide = "none") +
    theme(legend.position = "right",
          axis.line = element_line(size = 1),
          axis.text = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 15),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 15))
theo_color_plot_hc



real_specs_plot = ggplot(data = colors |> dplyr::filter(realtheo == "real"), mapping = aes(
    x = wavelength, y = reflectance, group = simsamp
)) +
    geom_line(aes(color = simsamp),linewidth = 1, alpha = 0.5) +
    theme_classic() +
    labs(x = "Wavelength (nm)", y = "Reflectance") +
    scale_color_manual(values = c(realcc)) +
    theme(legend.position = "none",
          axis.line = element_line(size = 1),
          axis.text = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 15),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 15)) 
    
real_specs_plot

ggsave(filename = "objects/real_specs.jpg", plot = real_specs_plot,
       device = "jpg", width = 7.75, height = 6, dpi = 600)

theo_specs_plot = ggplot(data = colors |> dplyr::filter(realtheo == "theo"), mapping = aes(
    x = wavelength, y = reflectance, group = simsamp
)) +
    geom_line(aes(color = simsamp),linewidth = 1, alpha = 0.5) +
    theme_classic() +
    labs(x = "Wavelength (nm)", y = "Reflectance") +
    scale_color_manual(values = c(theocc)) +
    theme(legend.position = "none",
          axis.line = element_line(size = 1),
          axis.text = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 15),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 15)) 

theo_specs_plot

ggsave(filename = "objects/theo_specs.jpg", plot = theo_specs_plot,
       device = "jpg", width = 7.75, height = 6, dpi = 600)


##### ABP Branches - Original Submission #####
# check that real and theoretical distributions are similar
col_anth$Pelargonidin_01 = ifelse(col_anth$Pelargonidin > 0, 1, 0)
col_anth$Delphinidin_01 = ifelse(col_anth$Delphinidin > 0, 1, 0)
col_anth$Cyanidin_01 = ifelse(col_anth$Cyanidin > 0, 1, 0)
col_anth$Petunidin_01 = ifelse(col_anth$Petunidin > 0, 1, 0)
col_anth$Peonidin_01 = ifelse(col_anth$Peonidin > 0, 1, 0)
col_anth$Malvidin_01 = ifelse(col_anth$Malvidin > 0, 1, 0)

col_anth$number_of_anthos = col_anth$Malvidin_01 + 
    col_anth$Peonidin_01 +
    col_anth$Petunidin_01 +
    col_anth$Cyanidin_01 +
    col_anth$Delphinidin_01 +
    col_anth$Pelargonidin_01


col_anth$number_of_branches = col_anth$del_branch +
    col_anth$cya_branch + col_anth$pel_branch

# plot
library(viridis)


num_anthos_plot = ggplot(data = col_anth,
                   mapping = aes(x = number_of_anthos, fill = realtheo)) +
    geom_bar(aes(y = ..prop..),
             position = "identity", alpha = .5) +
    scale_fill_manual(breaks = c("theo", "real"),
                        values = c(viridis::viridis(2)),
                      name = "",
                      labels = c("Theoretical", "Observed")
                      ) +
    scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6)) +
    xlab("Unique Anthocyanins") +
    ylab("Percent") +
    theme_bw() +
    theme(legend.position = "none")
num_anthos_plot

num_branches_plot = ggplot(data = col_anth,
                         mapping = aes(x = number_of_branches, fill = realtheo)) +
    geom_bar(aes(y = ..prop..),
             position = "identity",
             alpha = .5) +
    scale_fill_manual(breaks = c("theo", "real"),
                      values = c(viridis::viridis(2)),
                      name = "",
                      labels = c("Theoretical", "Observed")
    ) +
    xlab("Unique ABP Branches") +
    ylab("Percent") +
    theme_bw()
num_branches_plot

nums = num_anthos_plot + num_branches_plot
nums

ggsave(filename = "objects/number_plot.jpg", plot = nums,
       device = "jpg", width = 7, height = 4, dpi = 600)


# poisson regression model - number of anthos
# Fit a Poisson regression model
poisson_model <- glm(number_of_anthos ~ realtheo, family = poisson(link = "log"), data = col_anth)

# Display model summary
summary(poisson_model)

# poisson regression model - number of branches
# Fit a Poisson regression model
poisson_model2 <- glm(number_of_branches ~ realtheo, family = poisson(link = "log"), data = col_anth)

# Display model summary
summary(poisson_model2)

##### ABP Branches - Revised Submission ######

# run a chisquared test
# create a table
real_antho_pres = anthos |> 
    dplyr::filter(realtheo == "real") |>
    dplyr::select(9:14) |>
    dplyr::rename(Pelargonidin = Pelargonidin_01,
                  Cyanidin = Cyanidin_01,
                  Delphinidin = Delphinidin_01,
                  Peonidin = Peonidin_01,
                  Petunidin = Petunidin_01,
                  Malvidin = Malvidin_01)

real_antho_branch = anthos |> 
    dplyr::filter(realtheo == "real") |>
    dplyr::select(15:17) |>
    dplyr::rename(`Pelargonidin Branch` = pel_branch,
                  `Delphinidin Branch` = del_branch,
                  `Cyanidin Branch` = cya_branch)

# a function to do chi-squared tests for all combinations of anthocyanins
compute_chi_residuals = function(data) {
    traits = colnames(data)
    results_list = list()  # Store results for each pair
    
    for (i in 1:(length(traits) - 1)) {
        for (j in (i + 1):length(traits)) {
            trait1 = traits[i]
            trait2 = traits[j]
            
            cont_table = table(data[[trait1]], data[[trait2]])
            
            # Only proceed if there’s enough data
            if (min(dim(cont_table)) > 1) {
                chi_test = chisq.test(cont_table)
                residuals = chi_test$stdres  # Get standardized residuals
                
                pair_results = expand.grid(Level1 = rownames(residuals), Level2 = colnames(residuals))
                pair_results$Trait1 = trait1
                pair_results$Trait2 = trait2
                pair_results$P_Value = chi_test$p.value
                pair_results$Residual = as.vector(residuals)
                
                # Interpret direction
                pair_results$Direction = ifelse(pair_results$Residual > 1.96, "Positive",
                                                ifelse(pair_results$Residual < -1.96, "Negative", "Neutral"))
                
                results_list[[paste(trait1, trait2, sep = "_vs_")]] = pair_results
            }
        }
    }
    
    # Combine all results into one big dataframe
    results_df = do.call(rbind, results_list)
    rownames(results_df) = NULL  # Reset row names
    return(results_df)
}


chi_result = compute_chi_residuals(real_antho_pres) |>
    dplyr::filter(Level1 == 1, Level2 == 1)

chi_result


# Define the order of anthocyanins
ac_order <- c("Delphinidin","Petunidin","Malvidin", # del_branch
             "Cyanidin", "Peonidin", # cya branch
             "Pelargonidin" # Pel branch
                )

# Initialize empty matrix
res_matrix <- matrix(NA, nrow = 6, ncol = 6)
rownames(res_matrix) <- ac_order
colnames(res_matrix) <- ac_order

# Fill *both* halves to make it symmetric
for(i in 1:nrow(chi_result)) {
    r <- chi_result$Trait1[i]
    c <- chi_result$Trait2[i]
    val <- chi_result$Residual[i]
    res_matrix[r, c] <- val
    res_matrix[c, r] <- val   # mirror fill
}

# Keep only one triangle
res_matrix[upper.tri(res_matrix, diag = TRUE)] <- NA

library(tibble)
# Convert matrix to long format
res_df <- as.data.frame(res_matrix) %>%
    tibble::rownames_to_column(var = "Trait1") %>%
    pivot_longer(cols = -Trait1, names_to = "Trait2", values_to = "Residual") %>%
    # Keep only upper triangle
    filter(!is.na(Residual))

# Optional: order factors to keep the matrix shape
res_df <- res_df %>%
    mutate(Trait1 = factor(Trait1, levels = ac_order),
           Trait2 = factor(Trait2, levels = ac_order))

# Plot triangular heatmap
ggplot(res_df, aes(x = Trait2, y = Trait1, fill = Residual)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                         name = "Chi-Squared\nResidual",
                         limits = c(-8.7, 8.7),
                         breaks = c(-6, -4, -2, 0, 2, 4, 6)) +
    theme_classic() +
    coord_fixed() +
        labs(x = "", y = "") +
    theme(
          axis.line = element_line(size = 1),
          axis.text = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 15),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 15))

ggsave("objects/chi_squared.jpg", width = 7.75, height = 6, units = "in", dpi = 600)


# rerun on branches of the ABP Pathway


branch_result = compute_chi_residuals(real_antho_branch) |>
    dplyr::filter(Level1 == 1, Level2 == 1)
# Optional: order factors to keep the matrix shape

br_order = c("Delphinidin Branch", "Cyanidin Branch", "Pelargonidin Branch")

branch_result <- branch_result %>%
    mutate(Trait1 = factor(Trait1, levels = br_order),
           Trait2 = factor(Trait2, levels = br_order))

# Plot triangular heatmap
ggplot(branch_result, aes(x = Trait2, y = Trait1, fill = Residual)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                         name = "Chi-Squared\nResidual",
                         limits = c(-8.5, 8.5),
                         breaks = c(-6, -4, -2, 0, 2, 4, 6)) +
    theme_classic() +
    coord_fixed() +
    scale_x_discrete(labels = c("Cyanidin\nBranch", "Pelargonidin\nBranch")) +
    scale_y_discrete(labels = c("Delphinidin\nBranch", "Cyanidin\nBranch")) +
    labs(x = "", y = "") +
    theme(axis.line = element_line(size = 1),
          axis.text = element_text(size = 25, color = "black"),
          axis.title = element_text(size = 15),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 15))

ggsave("objects/chi_squared_branches.jpg", width = 7.75, height = 6, units = "in", dpi = 600)

# part two - random sampling to produce independent anthocyanin combinations
# compare anthocyanin numbers in null distribution to simulated distribution

# get df with only real anthos
r_anthos = anthos |>
    dplyr::filter(realtheo == "real") |>
    dplyr::select(Pelargonidin_01, Delphinidin_01, Cyanidin_01,
                  Peonidin_01, Petunidin_01, Malvidin_01)

# set seed
set.seed(123)
# simulate 10000 random distribution for each anthocyanin, derived from the actual distribution
pel_sim = sample(r_anthos$Pelargonidin_01, size = 10000, replace = T)
del_sim = sample(r_anthos$Delphinidin_01, size = 10000, replace = T)
cya_sim = sample(r_anthos$Cyanidin_01, size = 10000, replace = T)
peo_sim = sample(r_anthos$Peonidin_01, size = 10000, replace = T)
pet_sim = sample(r_anthos$Petunidin_01, size = 10000, replace = T)
mal_sim = sample(r_anthos$Malvidin_01, size = 10000, replace = T)

# create dataframe
sim_anthos = data.frame(pel_sim, del_sim, cya_sim, peo_sim, pet_sim, mal_sim) |>
    dplyr::rename(Pelargonidin_01 = pel_sim,
                  Delphinidin_01 = del_sim,
                  Cyanidin_01 = cya_sim,
                  Peonidin_01 = peo_sim,
                  Petunidin_01 = pet_sim,
                  Malvidin_01 = mal_sim)

rsim_anthos = rbind(r_anthos, sim_anthos)

# row sum to get number of anthocyanins
rsim_anthos$number_of_anthos = rowSums(rsim_anthos)

rsim_anthos$realsim = c(rep("real", 197), rep("sim", 10000))

# exclude rowsums of 0
rsim_anthos = rsim_anthos |>
    dplyr::filter(number_of_anthos != 0)

# add branches
# add binary branch data
rsim_anthos$del_branch = ifelse(rsim_anthos$Delphinidin_01 == 1 | 
                                   rsim_anthos$Petunidin_01 == 1 |
                                   rsim_anthos$Malvidin_01 == 1, 1, 0)
rsim_anthos$cya_branch = ifelse(rsim_anthos$Cyanidin_01 == 1 | 
                                   rsim_anthos$Peonidin_01 == 1, 1, 0)
rsim_anthos$pel_branch = ifelse(rsim_anthos$Pelargonidin_01 == 1, 1, 0)

# get number of unique ABP Branches
rsim_anthos$num_branches = rowSums(rsim_anthos[,9:11])

# plot
library(viridis)


num_anthos_plot = ggplot(data = rsim_anthos,
                         mapping = aes(x = number_of_anthos, fill = realsim)) +
    geom_bar(aes(y = ..prop..),
             position = "identity", alpha = .5) +
    scale_fill_manual(breaks = c("real", "sim"),
                      values = c(viridis::viridis(2)),
                      name = "",
                      labels = c("Observed", "Theoretical")
    ) +
    scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6)) +
    xlab("Unique Anthocyanins") +
    ylab("Proportion") +
    theme_classic() +
    theme(legend.position = "none",
          axis.line = element_line(size = 1),
          axis.text = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 15),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 15))
num_anthos_plot

ggsave(filename = "objects/num_anthos.jpg", plot = num_anthos_plot,
       device = "jpg", width = 6, height = 6, units = "in", dpi = 600)

num_branches_plot = ggplot(data = rsim_anthos,
                         mapping = aes(x = num_branches, fill = realsim)) +
    geom_bar(aes(y = ..prop..),
             position = "identity", alpha = .5) +
    scale_fill_manual(breaks = c("real", "sim"),
                      values = c(viridis::viridis(2)),
                      name = "",
                      labels = c("Observed", "Theoretical")
    ) +
    scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6)) +
    xlab("Unique ABP Branches") +
    ylab("Proportion") +
    theme_classic() +
    theme(legend.position = "none",
          axis.line = element_line(size = 1),
          axis.text = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 15),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 15))
num_branches_plot

ggsave(filename = "objects/num_branches.jpg", plot = num_branches_plot,
       device = "jpg", width = 6, height = 6, units = "in", dpi = 600)

real_an = rsim_anthos |>
    dplyr::filter(realsim == "real")

table(real_an$num_branches)

# statistics
# poisson regression model - number of anthos
# Fit a Poisson regression model
poisson_model <- glm(number_of_anthos ~ realsim, family = poisson(link = "log"), data = rsim_anthos)

# Display model summary
summary(poisson_model)

# poisson regression model - number of branches
# Fit a Poisson regression model
poisson_model2 <- glm(num_branches ~ realsim, family = poisson(link = "log"), data = rsim_anthos)

# Display model summary
summary(poisson_model2)

##### Hue and Brightness - Original Submission #####
hue_plot = ggplot(data = col_anth,
                  mapping = aes(x = realtheo,
                                y = hue,
                                fill = realtheo)) +
    geom_violin() +
    scale_fill_manual(breaks = c("theo", "real"),
                      values = c(viridis::viridis(2)),
                      name = "",
                      labels = c("Theoretical", "Observed")
    ) +
    ylab("Hue") +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank())
hue_plot

brightness_plot = ggplot(data = col_anth,
                         mapping = aes(x = realtheo,
                                       y = brightness,
                                       fill = realtheo)) +
    geom_violin() +
    scale_fill_manual(breaks = c("theo", "real"),
                      values = c(viridis::viridis(2)),
                      name = "",
                      labels = c("Theoretical", "Observed")
    ) +
    ylab("Brightness") +
    theme_bw() +
    theme(legend.position = "right",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank())
brightness_plot

hue_brightness = hue_plot + brightness_plot
hue_brightness

ggsave(filename = "objects/hue_brightness.jpg", plot = hue_brightness,
       device = "jpg", width = 7, height = 4, dpi = 600)

# compare variance using a levine's test
# Import required package
library(car)

# Using leveneTest()
levene = leveneTest(hue ~ realtheo, data = col_anth)

# print the result
levene

# Using leveneTest()
levene = leveneTest(brightness ~ realtheo, data = col_anth)

# print the result
levene

t.test(brightness ~ realtheo, data = col_anth)

(149.09 - 138.93) / 138.93

variance = col_anth |>
    dplyr::group_by(realtheo) |>
    dplyr::summarise(sd_hue = sd(hue),
                     sd_brightness = sd(brightness),
                     mean_hue = mean(hue),
                     mean_brightness = mean(brightness))
variance


(27.1 - 16.8) / 16.8

##### Visualize Model Accuracy in Colorspace - Original Submission #####
# compare training data to model prediction
specs$simsamp = specs$sample

# run function to calculate color parameters
colors = RYGB(specs)

colors = calc_colorspace(colors)
colors$hue = ifelse(colors$hue < 100, 
                           colors$hue + 360, 
                           colors$hue)
samples = colors$simsamp
colors$colorcode = specs$colorcode[1:61]

pred_colors = real_col_anth |>
    dplyr::mutate(simsamp = antho_db$spec_ID) |>
    dplyr::filter(simsamp %in% samples)
pred_colors = pred_colors[, c(1:11)]
samples = pred_colors$simsamp

colors = colors |>
    dplyr::filter(simsamp %in% samples)

real_pred = rbind(colors, pred_colors)
real_pred$realpred = c(rep("real", nrow(colors)),
                       rep("pred", nrow(pred_colors)))

cc = real_pred$colorcode

accuracy_plot = ggplot(data = real_pred,
                         mapping = aes(x = MS, y = LM, color = simsamp)) +
    geom_line(aes(group = simsamp), color = "black") +
    geom_point(aes(shape = realpred)) +
    theme_bw() +
    scale_color_manual(values = c(cc), guide = "none") +
    scale_shape_manual(values = c(19, 17),
                        labels = c("Predicted", "Measured")) +
    theme(legend.position = "none",
          legend.title = element_blank())
accuracy_plot

ggsave(filename = "objects/model_accuracy.jpg", plot = accuracy_plot,
       device = "jpg", width = 7, height = 5, dpi = 600)

# hue, brightness, and chroma
cc = colors$colorcode

hue_real_pred <- ggplot(data=real_pred, aes(x=realpred, y=hue, group=simsamp)) +
    geom_line(aes(color = simsamp),
              size=.75,
              alpha=1) +
    scale_color_manual(values = c(cc),
                       guide = "none") +
    scale_x_discrete(labels=c("pred" = "Predicted", "real" = "Measured"),
                     expand = c(0,0.1)) +
    theme_bw() +
    ylab("Hue") +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(size = 11,
                                     color = "black"))
          #axis.text.x = element_blank(),
          #axis.ticks.x = element_blank())
    
hue_real_pred

brightness_real_pred <- ggplot(data=real_pred, aes(x=realpred, y=brightness, group=simsamp)) +
    geom_line(aes(color = simsamp),
              size=.75,
              alpha=1) +
    scale_color_manual(values = c(cc),
                       guide = "none") +
    scale_x_discrete(labels=c("pred" = "Predicted", "real" = "Measured"),
                     expand = c(0,0.1)) +
    theme_bw() +
    ylab("Brightness") +
    theme(axis.title.x = element_blank(), 
          axis.text.x = element_text(size = 11,
                           color = "black"))

brightness_real_pred

chroma_real_pred <- ggplot(data=real_pred, aes(x=realpred, y=chroma, group=simsamp)) +
    geom_line(aes(color = simsamp),
              size=.75,
              alpha=1) +
    scale_color_manual(values = c(cc),
                       guide = "none") +
    scale_x_discrete(labels=c("pred" = "Predicted", "real" = "Measured"),
                     expand = c(0,0.1)) +
    theme_bw() +
    ylab("Chroma") +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(size = 11,
                                     color = "black"))

chroma_real_pred

model_accuracy = accuracy_plot + hue_real_pred + brightness_real_pred +
    chroma_real_pred + plot_layout(ncol=2) + plot_annotation(tag_levels = 'A')
model_accuracy

ggsave(filename = "objects/model_accuracy_4panel.jpg", plot = model_accuracy,
       device = "jpg", width = 7, height = 5, dpi = 600)

# table 1 - model diagnostics
# first pivot wider
# Reshape data to wide format
real = real_pred |>
    dplyr::filter(realpred == "real")
pred = real_pred |>
    dplyr::filter(realpred == "pred")
realpred_wide = merge(real, pred, by = "simsamp")

# linear models
hue_lm = lm(realpred_wide$hue.x ~ realpred_wide$hue.y)
summary(hue_lm)

chroma_lm = lm(realpred_wide$chroma.x ~ realpred_wide$chroma.y)
summary(chroma_lm)
plot(realpred_wide$chroma.x, realpred_wide$chroma.y)

brightness_lm = lm(realpred_wide$brightness.x ~ realpred_wide$brightness.y)
summary(brightness_lm)

LM_lm = lm(realpred_wide$LM.x ~ realpred_wide$LM.y)
summary(LM_lm)

MS_lm = lm(realpred_wide$MS.x ~ realpred_wide$MS.y)
summary(MS_lm)



##### Calculate Hypervolumes and Intersection - Original Submission #####
# library
library(hypervolume)
library(rgl)
library(alphahull)


# split the dataset on theoretical vs. real
col_anth_split = split(col_anth, col_anth$realtheo)
col_anth_split = col_anth_split[c(2,1)]

# calculate the hypervolume for each element of the splitted dataset
hv_list = mapply(function(x, y) 
    hypervolume_gaussian(x[, c(4,5,6)],
                         samples.per.point=100, name = y), 
    x = col_anth_split, 
    y = names(col_anth_split))

hv_list <- hypervolume_join(hv_list)

# calculate occupancy without groups
hv_occupancy <- hypervolume_n_occupancy(hv_list)

# calculate occupancy with groups
hv_occupancy_list_theo <- hypervolume_n_occupancy(hv_list, 
                                                 classification = c("Theoretical", "Observed"))

# plot
library(viridis)
jpeg("objects/hypervoume.jpeg", quality = 100, width = 7.5,
     height = 5, units = "in", res = 600)
plot(hv_occupancy_list_theo, 
     cex.random = .5, show.density = T, show.contour = F,
     num.points.max.random = 1000,
     cex.data=1,
     point.dark.factor=0.1,
     cex.axis=0.75,cex.names=2.0,
     cex.legend=1,
     point.alpha.min = .1,
     num.points.max.data = 272,
     colors = viridis::viridis(2))

dev.off()

# bootstrap the hypervolumes
hv_list_boot = hypervolume_n_resample(name = "example", hv_list)

# calculate occupancy on bootstrapped hypervolumes
hv_occupancy_boot_theo = hypervolume_n_occupancy_bootstrap(path = hv_list_boot,
                                                          name = "example_occ",
                                                          classification = c("Real", "Theoretical"))
# get the intersection
hv_output = get_occupancy_intersection_bootstrap(hv_occupancy_boot_theo)
hv_output

# make output table
observed_volume = hv_occupancy_list_theo@HVList[[2]]@Volume
theoretical_volume = hv_occupancy_list_theo@HVList[[1]]@Volume
overlap = hv_output[1,2]
overlap_2.5 = hv_output[1,5]
overlap_97.5 = hv_output[1,7]

percent_volume = observed_volume / theoretical_volume

percent_overlap = overlap / observed_volume
percent_overlap_2.5 = overlap_2.5 / observed_volume
percent_overlap_97.5 = overlap_97.5 / observed_volume

##### Calculate Hypervolumes and Intersection - Revised Submission #####


# split the dataset on theoretical vs. real
col_anth_split = split(colorspace, colorspace$realtheo)
col_anth_split = col_anth_split[c(2,1)]

# calculate the hypervolume for each element of the splitted dataset
hv_list = mapply(function(x, y) 
    hypervolume_gaussian(x[, c(4,5,6)],
                         samples.per.point=100, name = y), 
    x = col_anth_split, 
    y = names(col_anth_split))

hv_list <- hypervolume_join(hv_list)

# calculate occupancy without groups
hv_occupancy <- hypervolume_n_occupancy(hv_list)

# calculate occupancy with groups
hv_occupancy_list_theo <- hypervolume_n_occupancy(hv_list, 
                                                  classification = c("Theoretical", "Observed"))

# plot
library(viridis)
jpeg("objects/hypervolume.jpeg", quality = 100, width = 9,
     height = 6, units = "in", res = 600)
plot(hv_occupancy_list_theo, 
     cex.random = .5, show.density = T, show.contour = F,
     num.points.max.random = 1000,
     cex.data=1,
     point.dark.factor=0.1,
     cex.axis=0.75,cex.names=2.0,
     cex.legend=1,
     point.alpha.min = .1,
     num.points.max.data = 272,
     colors = rev(viridis::viridis(2)))

dev.off()

# bootstrap the hypervolumes
#hv_list_boot = hypervolume_n_resample(name = "example", hv_list)

# calculate occupancy on bootstrapped hypervolumes
#hv_occupancy_boot_theo = hypervolume_n_occupancy_bootstrap(path = hv_list_boot,
                                                           #name = "example_occ",
                                                           #classification = c("Real", "Theoretical"))
# get the intersection
#hv_output = get_occupancy_intersection_bootstrap(hv_occupancy_boot_theo)
#hv_output

hv_set <- hypervolume_set(hv_occupancy_list_theo@HVList[[2]], 
                          hv_occupancy_list_theo@HVList[[1]], 
                          check.memory = FALSE)

hv_stats <- hypervolume_overlap_statistics(hv_set)

hv_stats

# make output table
observed_volume = hv_occupancy_list_theo@HVList[[2]]@Volume
theoretical_volume = hv_occupancy_list_theo@HVList[[1]]@Volume
overlap = hv_output[1,2]
overlap_2.5 = hv_output[1,5]
overlap_97.5 = hv_output[1,7]

percent_volume = observed_volume / theoretical_volume
percent_volume

percent_overlap = (2 * overlap) / (observed_volume + theoretical_volume)
percent_overlap

percent_overlap_2.5 = overlap_2.5 / observed_volume
percent_overlap_2.5
percent_overlap_97.5 = overlap_97.5 / observed_volume
percent_overlap_97.5

## omit outliers and rerun analyses
colorspace_filt = colorspace |>
    dplyr::filter(if_all(c(Pelargonidin_standard, Cyanidin_standard,
                           Delphinidin_standard, Peonidin_standard,
                           Petunidin_standard, Malvidin_standard), ~ . <= 3))

# split the dataset on theoretical vs. real
col_anth_split2 = split(colorspace_filt, colorspace_filt$realtheo)
col_anth_split2 = col_anth_split2[c(2,1)]

# calculate the hypervolume for each element of the splitted dataset
hv_list2 = mapply(function(x, y) 
    hypervolume_gaussian(x[, c(4,5,6)],
                         samples.per.point=100, name = y), 
    x = col_anth_split2, 
    y = names(col_anth_split2))



hv_list2 <- hypervolume_join(hv_list2)

# calculate occupancy without groups
hv_occupancy2 <- hypervolume_n_occupancy(hv_list2)

# calculate occupancy with groups
hv_occupancy_list_theo2 <- hypervolume_n_occupancy(hv_list2, 
                                                  classification = c("Theoretical", "Observed"))

jpeg("objects/hypervolume2.jpeg", quality = 100, width = 9,
     height = 6, units = "in", res = 600)
plot(hv_occupancy_list_theo2, 
     cex.random = .5, show.density = T, show.contour = F,
     num.points.max.random = 1000,
     cex.data=1,
     point.dark.factor=0.1,
     cex.axis=0.75,cex.names=2.0,
     cex.legend=1,
     point.alpha.min = .1,
     num.points.max.data = 272,
     colors = rev(viridis::viridis(2)))

dev.off()

hv_set2 <- hypervolume_set(hv_occupancy_list_theo2@HVList[[2]], 
                          hv_occupancy_list_theo2@HVList[[1]], 
                          check.memory = FALSE)

hv_stats2 <- hypervolume_overlap_statistics(hv_set2)

hv_stats2


### calculate dispersion and edge occupancy statistics

# extract distance to centroid for random points within the hypervolume
library(hypervolume)
library(vegan)

observed_points = hv_list2@HVList[["real"]]@RandomPoints
theoretical_points = hv_list2@HVList[["theo"]]@RandomPoints

# points_obs, points_theo: matrices (n x p)
n1 <- nrow(observed_points)
n2 <- nrow(theoretical_points)
p  <- ncol(observed_points)

S1 <- cov(observed_points)
S2 <- cov(theoretical_points)

# pooled covariance
Sp <- ((n1-1)*S1 + (n2-1)*S2) / (n1 + n2 - 2)

# Centroids
centroid_obs = colMeans(observed_points)
centroid_theo = colMeans(theoretical_points)

# Mahalanobis using pooled covariance (distances on common scale)
dist_pooled_obs  <- mahalanobis(observed_points,  center = centroid_obs, cov = Sp)
dist_pooled_theo <- mahalanobis(theoretical_points, center = centroid_theo, cov = Sp)



# Mahalanobis distance to centroid
cov_obs = cov(observed_points)
cov_theo = cov(theoretical_points)

# Mahalanobis distances
distM_obs = mahalanobis(observed_points, center = centroid_obs, cov = cov_obs)
distM_theo = mahalanobis(theoretical_points, center = centroid_theo, cov = cov_theo)

# summarise
summary(distM_obs)
summary(distM_theo)

summary(dist_pooled_obs)
summary(dist_pooled_theo)

# KS & Anderson-Darling Test
k_test = ks.test(distM_obs, distM_theo)
k_test

k_testpooled = ks.test(dist_pooled_obs, dist_pooled_theo)
k_testpooled

library(kSamples)
ad_test = ad.test(distM_obs, distM_theo)

ad_test

ad_testpooled = ad.test(dist_pooled_obs, dist_pooled_theo)

ad_testpooled

# Combine into a single data frame
df <- data.frame(
    distM = c(distM_obs, distM_theo),
    group = factor(c(rep("Observed", length(distM_obs)),
                     rep("Theoretical", length(distM_theo))))
)
# Combine into a single data frame
df2 <- data.frame(
    distM = c(dist_pooled_obs, dist_pooled_theo),
    group = factor(c(rep("Observed", length(dist_pooled_obs)),
                     rep("Theoretical", length(dist_pooled_theo))))
)



library(ggplot2)
# Create CDF plot
ggplot(df, aes(x = distM, color = group)) +
    stat_ecdf(size = 1) +
    labs(
        x = "Mahalanobis distance to centroid",
        y = "Cumulative probability",
        title = "CDF of Mahalanobis distances: Observed vs Theoretical Hypervolumes"
    ) +
    theme_minimal() +
    theme(text = element_text(size = 14))

ggplot(df2, aes(x = distM, color = group)) +
    stat_ecdf(size = 1) +
    labs(
        x = "Mahalanobis distance to centroid",
        y = "Cumulative probability",
        title = "CDF of Mahalanobis distances: Observed vs Theoretical Hypervolumes"
    ) +
    theme_minimal() +
    theme(text = element_text(size = 14))


# Summary stats (mean, var, skew, kurtosis)
library(e1071) # for skewness, kurtosis
summary_obs = c(mean=mean(distM_obs), var=var(distM_obs), skew=skewness(distM_obs), kurt=kurtosis(distM_obs))
summary_theo = c(mean=mean(distM_theo), var=var(distM_theo), skew=skewness(distM_theo), kurt=kurtosis(distM_theo))

summary_obs
summary_theo

summary_obs2 = c(mean=mean(dist_pooled_obs), var=var(dist_pooled_obs), skew=skewness(dist_pooled_obs), kurt=kurtosis(dist_pooled_obs))
summary_theo2 = c(mean=mean(dist_pooled_theo), var=var(dist_pooled_theo), skew=skewness(dist_pooled_theo), kurt=kurtosis(dist_pooled_theo))

summary_obs2
summary_theo2

# plot raw distributions of mahalanobis distance
ggplot(df, aes(x = distM, color = group, fill = group)) +
    geom_density(aes(y = ..density..), alpha = 0.3, size = 1.2) +
    labs(
        x = "Mahalanobis Distance to Centroid",
        y = "Proportion",
    ) +
    scale_fill_manual(values = c("#440154FF", "#FDE725FF")) +
    scale_color_manual(values = c("#440154FF", "#FDE725FF")) +
    theme_classic() +
    theme(legend.position = "none",
        text = element_text(size = 15),
          axis.line = element_line(size = 1),
          axis.text = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 15),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 15))

# plot raw distributions of mahalanobis distance
ggplot(df2, aes(x = distM, color = group, fill = group)) +
    geom_density(aes(y = ..density..), alpha = 0.5, size = 1.2) +
    labs(
        x = "Mahalanobis Distance to Centroid",
        y = "Proportion",
    ) +
    scale_fill_manual(values = c("#440154FF", "#FDE725FF")) +
    scale_color_manual(values = c("#440154FF", "#FDE725FF")) +
    theme_classic() +
    theme(legend.position = "none",
          text = element_text(size = 15),
          axis.line = element_line(size = 1),
          axis.text = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 15),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 15))

ggsave(filename = "objects/mahalanobis_distance.jpg",
       device = "jpg", width = 6, height = 6, units = "in", dpi = 600)


##### Geographic Overlap, phylogenetic distance, and Color Distance - Original & Revised Submission #####


# pseudocode
# compile occurence dataset for all ruellia that we have color data for
# use overlap function from AM NAT
# write function to calculate distance in colorspace
# linear model predicting colorspace distance from Overlap and phylogenetic distance

# compile occurence dataset for all ruellia that we have color data for

# make species list
species_list = unique(antho_db$species)
# remove species and species nova
species_list = species_list[! species_list %in% c('species', 'species nova')]
species_list = paste("Ruellia", species_list)

# load in occurence data from Am Nat submission
gbif_man = readr::read_csv("occurence_data/GBIF_manual.csv") |>
    dplyr::select(species, decimalLatitude, decimalLongitude)
gbif_man$species = paste("Ruellia", gbif_man$species)

gbif = readr::read_csv("occurence_data/Ruellia_GBIF.22.csv") |>
    dplyr::select(species, decimalLatitude, decimalLongitude)
gbif = gbif[complete.cases(gbif),]

idig = readr::read_csv("occurence_data/Ruellia_IDIG.22.RAW.csv") |>
    dplyr::select(`dwc:specificEpithet`, `dwc:decimalLatitude`, `dwc:decimalLongitude`) |>
    dplyr::rename(species = `dwc:specificEpithet`, 
                  decimalLatitude = `dwc:decimalLatitude`, 
                  decimalLongitude = `dwc:decimalLongitude`)
idig = idig[complete.cases(idig),]
idig$species = paste("Ruellia", idig$species)

IMH = readr::read_csv("occurence_data/Ruellia_IMH.22.csv") |>
    dplyr::select(scientificName, decimalLatitude, decimalLongitude) |>
    dplyr::rename(species = scientificName)
IMH = IMH[complete.cases(IMH),]

TEX = readr::read_csv("occurence_data/Ruellia_TEX_22.csv") |>
    dplyr::select(scientificName, decimalLatitude, decimalLongitude) |>
    dplyr::rename(species = scientificName)
TEX = TEX[complete.cases(TEX),]

supp = readr::read_csv("occurence_data/supplemental_points_literature.csv") |>
    dplyr::select(species, decimalLatitude, decimalLongitude)
supp$species = paste("Ruellia", supp$species)
supp = supp[complete.cases(supp),]

tripp = readr::read_csv("occurence_data/Tripp_et_al_2021.csv") |>
    dplyr::select(species, lat_DD, long_DD) |>
    dplyr::rename(decimalLatitude = lat_DD,
                  decimalLongitude = long_DD)
tripp = tripp[complete.cases(tripp),]
tripp$species = paste("Ruellia", tripp$species)

ruellia_occ = rbind(supp, tripp, TEX, IMH, idig, gbif, gbif_man) |>
    dplyr::distinct() |>
    dplyr::rename(decimallongitude = decimalLongitude,
           decimallatitude = decimalLatitude)

# change names to accepted names
# load in names dataset
accepted_names = readr::read_csv("occurence_data/Ruellia_names_2022.csv")

matches = match(ruellia_occ$species, accepted_names$verbatum.name)
matches
acc_spec = c(accepted_names$accepted_species[matches])
acc_spec = paste("Ruellia", acc_spec)

ruellia_occ$species = ifelse(acc_spec == "Ruellia NA", ruellia_occ$species, acc_spec)

# remove species that we don't have flower color data for
ruellia_occ = ruellia_occ |>
    dplyr::filter(species %in% species_list)

# switch long and lat
ruellia_occ = ruellia_occ[, c(1,3,2)]
# filter incorrect longitude and latitude points
ruellia_occ = ruellia_occ |>
    dplyr::filter(decimallongitude > -360 & decimallongitude < 360 &
                    decimallatitude > -180 & decimallatitude < 180) |>
    dplyr::rename(decimalLatitude = decimallatitude,
                  decimalLongitude = decimallongitude)



# coordinate cleaner
library(CoordinateCleaner)
flags <- clean_coordinates(x = ruellia_occ, 
                          tests = c("capitals", 
                                    "centroids","equal",
                                    "institutions", 
                                    "outliers", "seas", 
                                    "zeros"))

#Exclude problematic records
ruellia_clean = ruellia_occ[flags$.summary,]


# use overlap function to compute overlap metric
#source("CalculateO.R")
#library(data.table)
#library(geosphere)
#rue_overlap = CalculateO(ruellia_clean)

# save
readr::write_csv(rue_overlap, "objects/rue_overlap.csv")

# calculate distance values in colorspace

# add species column and calculate means
species_simsamp = antho_db |>
    dplyr::select(species, hplc_id1) |>
    dplyr::rename(simsamp = hplc_id1)
species_simsamp$species = paste("Ruellia", species_simsamp$species)

real_col = colorspace |>
    dplyr::filter(realtheo == "real")

spec_col_anth = merge(real_col, species_simsamp, by = "simsamp") |>
    dplyr::group_by(species) |>
    dplyr::summarise(hue = mean(hue),
                     chroma = mean(chroma),
                     brightness = mean(brightness),
                     Pelargonidin = mean(Pelargonidin_standard),
                     Delphinidin = mean(Delphinidin_standard),
                     Cyanidin = mean(Cyanidin_standard),
                     Petunidin = mean(Petunidin_standard),
                     Peonidin = mean(Peonidin_standard),
                     Malvidin = mean(Malvidin_standard)) |>

# scale hue, chroma, and brightness to remove scale biases
    dplyr::mutate_if(is.numeric, ~(scale(.) %>% as.vector))

# calculate pairwise color distances
library(rdist)
color_dist = pdist(spec_col_anth[, c(2:4)])
color_dist

color_dists = c(color_dist)
color_dists

sp1 = rep(spec_col_anth$species, nrow(spec_col_anth))
sp2 = rep(spec_col_anth$species, each = nrow(spec_col_anth))

color_dist_df = data.frame(sp1, sp2, color_dists)


# merge color distances and overlap metric
rue_overlap = readr::read_csv("objects/rue_overlap.csv")

color_overlap = merge(color_dist_df, rue_overlap, by = c("sp1", "sp2"))

# calculate pairwise phylogenetic distances

# load in the phylogeny in the form of a .tre file.
ruellia_phylo <- read.tree('color data/color data/phylogeny/03-Ruellia_phylo.timed.tre')
plot(ruellia_phylo)
phylo_dist = c(cophenetic.phylo(ruellia_phylo))
phylo_dist

sp1 = rep(ruellia_phylo[["tip.label"]], length(ruellia_phylo[["tip.label"]]))
sp2 = rep(ruellia_phylo[["tip.label"]], each = length(ruellia_phylo[["tip.label"]]))

phylo_dist_df = data.frame(sp1, sp2, phylo_dist)

# change species labels to corresponding labels in antho_db
sp1_matches = match(phylo_dist_df$sp1, antho_db$phylo_ID)
sp2_matches = match(phylo_dist_df$sp2, antho_db$phylo_ID)

sp1_new = c(antho_db$species[sp1_matches])
sp1_new = paste("Ruellia", sp1_new)

phylo_dist_df$sp1 = ifelse(sp1_new == "Ruellia NA", NA, sp1_new)

sp2_new = c(antho_db$species[sp2_matches])
sp2_new = paste("Ruellia", sp2_new)

phylo_dist_df$sp2 = ifelse(sp2_new == "Ruellia NA", NA, sp2_new)

phylo_dist_df = phylo_dist_df[complete.cases(phylo_dist_df),]

# merge
color_ov_phylo = merge(color_overlap, phylo_dist_df, by = c("sp1", "sp2"))

# log O to make more normal
color_ov_phylo$log_O = log(color_ov_phylo$O + 0.001)

# Ensure unique pairs and keep other columns
color_ov_phylo = color_ov_phylo |>
    dplyr::group_by(sp1, sp2) |>
    dplyr::slice(1) |> # Keep the first occurrence
    dplyr::ungroup()
# run linear models
lm = lm(color_dists ~ log_O * phylo_dist, data = color_ov_phylo)
summary(lm)

lm = lm(color_dists ~ log_O + phylo_dist, data = color_ov_phylo)
summary(lm)

lm = lm(color_dists ~ phylo_dist, data = color_ov_phylo)
summary(lm)

lm = lm(color_dists ~ log_O, data = color_ov_phylo)
summary(lm)

color_phylo_plot = ggplot(data = color_ov_phylo,
                            mapping = aes(x = phylo_dist, y = color_dists)) +
    geom_point(alpha = .2) +
    xlab("Phylogenetic Distance") +
    ylab("Color Distance") +
    theme_classic() + 
    stat_smooth(method = "lm") +
    theme(axis.line = element_line(size = 1),
          axis.text = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 15),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 15))
color_phylo_plot

ggsave(filename = "objects/color_phylo.jpg", plot = color_phylo_plot,
       device = "jpg", width = 6, height = 6, units = "in", dpi = 600)

color_overlap_plot = ggplot(data = color_ov_phylo,
                            mapping = aes(x = log_O, y = color_dists)) +
    geom_point(alpha = .2) +
    xlab("log(Geographic Overlap)") +
    ylab("Color Distance") +
    theme_classic() +
    stat_smooth(method = "lm") +
    theme(axis.line = element_line(size = 1),
          axis.text = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 15),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 15))
color_overlap_plot

ggsave(filename = "objects/color_overlap.jpg", plot = color_overlap_plot,
       device = "jpg", width = 6, height = 6, units = "in", dpi = 600)

dists_ov_phylo = color_overlap_plot + 
    color_phylo_plot + 
    plot_annotation(tag_levels = "A"
)

dists_ov_phylo

ggsave(filename = "objects/color_ov_phylo.jpg", plot = dists_ov_phylo,
       device = "jpg", width = 7, height = 5, dpi = 600)

# calculate anthodist and compare to phylogenetic distance
antho_dist = pdist(spec_col_anth[, c(5:10)])
antho_dist

antho_dists = c(antho_dist)
antho_dists

sp1 = rep(spec_col_anth$species, nrow(spec_col_anth))
sp2 = rep(spec_col_anth$species, each = nrow(spec_col_anth))

antho_dist_df = data.frame(sp1, sp2, antho_dists)

antho_phylo = merge(antho_dist_df, phylo_dist_df, by = c("sp1", "sp2"))

antho_phylo_plot = ggplot(data = antho_phylo,
                          mapping = aes(x = phylo_dist, y = antho_dists)) +
    geom_point(alpha = .2) +
    xlab("Phylogenetic Distance") +
    ylab("Anthocyanin Distance") +
    theme_bw() + 
    stat_smooth(method = "lm")
antho_phylo_plot

ggsave(filename = "objects/antho_phylo.jpg", plot = antho_phylo_plot,
       device = "jpg", width = 7, height = 5, dpi = 600)


lm = lm(antho_dists ~ phylo_dist, data = antho_phylo)
summary(lm)
