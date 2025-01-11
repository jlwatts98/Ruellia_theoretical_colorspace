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
antho_db <- read_excel("color data/color data/antho data/antho_db.xlsx") %>%
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

spec_db <- read_excel("color data/color data/spec data/spec_db.xlsx")
load("color data/color data/spec data/03-mean.specs.RData")

# merge and filter data to only include species with anthos and corolla
antho_spec = merge(spec_db, antho_db, by = "spec_ID") 

# look at the structure of the dataframe
str(antho_spec)

# filter out the useless variables and rename them
antho_spec = antho_spec[, c(1:23, 30, 37, 44, 51, 58, 65)]

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
specs = antho_spec[, c(30:430)]

specs = stack(specs)
colnames(specs)[1] = "reflectance"
colnames(specs)[2] = "wavelength"
specs$wavelength = as.numeric(specs$wavelength) + 299

specs$sample = rep(antho_spec$spec_ID, 401)
specs$colorcode = rep(antho_spec$hex_color, 401)
specs$color = rep(antho_spec$color, 401)

# add anthocyanin concentrations to the specs data.frame
specs$Pelargonidin = rep(antho_spec$Pelargonidin, 401)
specs$Cyanidin = rep(antho_spec$Cyanidin, 401)
specs$Peonidin = rep(antho_spec$Peonidin, 401)
specs$Petunidin = rep(antho_spec$Petunidin, 401)
specs$Delphinidin = rep(antho_spec$Delphinidin, 401)
specs$Malvidin = rep(antho_spec$Malvidin, 401)

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

##### Model Predicting Reflectance Spectra from Anthocyanins #####

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

##### Theoretical Morphospace #####
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


##### Colorspace Comparisons and Trends #####

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
        dplyr::filter(wavelength >= 400) |>
        dplyr::summarise(brightness = sum(reflectance))
    
    spec_df$brightness = rep(brightness$brightness, each = 401)
    
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
    # calculate Y
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
# to avoid articial separation between 0 and 360)
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


##### ABP Branches #####
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
#
##### Hue and Brightness #####
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

##### Visualize Model Accuracy in Colorspace #####
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


##### Calculate Hypervolumes and Intersection #####
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

##### Geographic Overlap, phylogenetic distance, and Color Distance #####


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
                    decimallatitude > -180 & decimallatitude < 180)



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


spec_col_anth = merge(real_col_anth, species_simsamp, by = "simsamp") |>
    dplyr::group_by(species) |>
    dplyr::summarise(hue = mean(hue),
                     chroma = mean(chroma),
                     brightness = mean(brightness),
                     Pelargonidin = mean(Pelargonidin),
                     Delphinidin = mean(Delphinidin),
                     Cyanidin = mean(Cyanidin),
                     Petunidin = mean(Petunidin),
                     Peonidin = mean(Peonidin),
                     Malvidin = mean(Malvidin)) |>

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
    theme_bw() + 
    stat_smooth(method = "lm")
color_phylo_plot

color_overlap_plot = ggplot(data = color_ov_phylo,
                            mapping = aes(x = log_O, y = color_dists)) +
    geom_point(alpha = .2) +
    xlab("log(Geographic Overlap)") +
    ylab("Color Distance") +
    theme_bw() +
    stat_smooth(method = "lm")
color_overlap_plot

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
