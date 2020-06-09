library(tidyverse)
# devtools::install_github("ha0ye/rEDM")
library(rEDM)
# install.packages("plot3D")
library(plot3D)

# Analysis pipeline based on https://ha0ye.github.io/rEDM/articles/rEDM.html and Cenci, S., Sugihara, 
# G. and Saavedra, S., 2019. Regularized S‐map for inference and forecasting with noisy ecological 
# time series. Methods in Ecology and Evolution, 10(5), pp.650-660 for future projections

# Reading in the simulated data
temp <- as.data.frame(read_table2("data.txt",col_names=TRUE) %>% arrange(Time_since_Eemian))
temp
# First five rows
#Time_since_Eemian     Ne Area
#1                   0  262.5    1.50
#2                   1  312.5    1.50
#3                   2  362.5    1.50
#4                   3  412.5    1.00
#5                   4  462.5    0.00
str(temp)
#'data.frame':	121 obs. of  3 variables:
#  $ Time_since_Eemian: num  0 1 2 3 4 5 6 7 8 9 ...
#$ Ne               : num  262 312 362 412 462 ...
#$ Area          : num  1.5 1.5 1.5 1 0 0 0 0 0 0 ...
#- attr(*, "spec")=
#  .. cols(
#    ..   Time_since_Eemian = col_double(),
#    ..   Ne = col_double(),
#    ..   Area = col_double()
#    .. )

# Adding noise to the (overly smooth) simulated data to approximate a more realistic dataset
temp$Ne <- temp$Ne + rnorm(length(temp$Ne), sd = sd(temp$Ne) * 0.2)
temp$Area <- temp$Area + rnorm(length(temp$Area), sd = sd(temp$Area) * 0.2)

# Obtain the optimal values for E (embedding dimension) and tau 
# (the lag for time delay embedding)
opt_E_tau_output <- c("lib_prop_of_data","best_Ne_E","tau")

# Varying the amount of data used for “training” (lib), on which
# nearest neighbors can be identified
# versus “test” portion of the data (pred)
for (i in seq(0.1,0.9,0.1)) {
  lib <- c(1, floor(i*NROW(temp)))
  pred <- c(floor(i*NROW(temp))+1,NROW(temp))
  Ne_pred <- temp[,2]
  Ne_output <- simplex(Ne_pred, lib, pred, E = 1:10, tau = 1:30)
  opt_E_tau_output <- rbind(opt_E_tau_output,c(i,Ne_output$E[which((Ne_output$rho)==max(Ne_output$rho,na.rm=TRUE))][1],Ne_output$tau[which((Ne_output$rho)==max(Ne_output$rho,na.rm=TRUE))][1]))
}

as_tibble(opt_E_tau_output)
#V1               V2        V3   
#<chr>            <chr>     <chr>
#  1 lib_prop_of_data best_Ne_E tau  
#2 0.1              2         4    
#3 0.2              2         2    
#4 0.3              3         1    
#5 0.4              2         1    
#6 0.5              6         9    
#7 0.6              3         19   
#8 0.7              9         2    
#9 0.8              2         22   
#10 0.9              2         2 

# An E of 2 and tau of 1 gives the best Rho for the plurality of lib_prop
# During analyses of full datasets, all the returned combinations of E/tau
# and their influence on inferences should be evaluated

# Now we will evaluate whether our Ne data shows nonlinearity or not
Evalue <- 2
tauvalue <- 1
opt_theta_output <- c("lib_prop_of_data","E","theta","delta_MAE")
for (i in seq(0.1,0.9,0.1)) {
  lib <- c(1, floor(i*NROW(temp)))
  pred <- c(floor(i*NROW(temp))+1,NROW(temp))
  smap_output <- s_map(temp[,2], lib, pred, E = Evalue, tau=tauvalue)
    opt_theta_output <- rbind(opt_theta_output,c(i,Evalue,smap_output$theta[which(smap_output$rho==max(smap_output$rho))][1],(smap_output$rmse[which(smap_output$rho==max(smap_output$rho))][1]-smap_output$rmse[1])))
}

# Evidence for non-linearity (theta > 0), but dependent on lib proportion of data
as_tibble(opt_theta_output)
# A tibble: 10 x 4
#V1               V2    V3    V4               
#<chr>            <chr> <chr> <chr>            
#  1 lib_prop_of_data E     theta delta_MAE        
#2 0.1              2     0.75  7.36012741982967 
#3 0.2              2     0     0                
#4 0.3              2     0     0                
#5 0.4              2     0.5   -1.23210452838771
#6 0.5              2     0.75  -2.54115604484687
#7 0.6              2     1     -2.0629188376459 
#8 0.7              2     0.75  0.899547478183564
#9 0.8              2     1.5   -4.02583884126062
#10 0.9              2     2     -5.89040621819689

# Evaluating the prediction decay
for (i in seq(0.1,0.9,0.1)) {
  lib <- c(1, floor(i*NROW(temp)))
  pred <- c(floor(i*NROW(temp))+1,NROW(temp))
  simplex_output <- simplex(temp[,2], lib, pred, E = Evalue, tau = tauvalue,tp = 1:10)
  plot(simplex_output$tp, simplex_output$rho, type = "l", xlab = paste("Time to Prediction (tp) for",i), 
       ylab = "Forecast Skill (rho)")
}  
# Generally only have predictive power about two or so steps out

# Next step is to figure out the appropriate time lag 
# Given we are interested in predicting Ne, we will model up to a 
# max of 30 (1/4 of the total time series)
# We can now make a lagged block for Ne using these values.
maxtime_lag <- 30
temp2 <- make_block(temp,max_lag=maxtime_lag,restrict_to_lib = FALSE)
# Don't need to keep all columns
temp2 <- temp2[,-c(1,3:31,63:91)]

lib <- c(1, NROW(temp))
pred <- c(1,NROW(temp))
lag_rho <- NULL
for (i in 3:31) {
  block_lnlp_output <- block_lnlp(temp2, lib = lib, pred = pred, columns = c("Ne",names(temp2)[i],"Area"),target_column = "Ne", stats_only = FALSE, silent = TRUE, tp=1, theta=0, method="s-map",save_smap_coefficients = TRUE,first_column_time = TRUE)
  lag_rho <- rbind(lag_rho,c(names(temp2)[i],block_lnlp_output$rho))
}

Ne_rho <- as_tibble(lag_rho) %>% arrange(desc(V2))
Ne_rho$V1 <- factor(Ne_rho$V1,levels = unique(Ne_rho$V1))

ggplot() + geom_point(data=Ne_rho,aes(x=V1,y=as.numeric(V2))) + ylim(0.95,0.98)
Ne_rho

# Defining the appropriate lag
lag_cols <- "Ne_1"

# This suggests there is the highest correlation with the cols defined in lag_cols
# So we'll use this for setting up the lag. Let's see how good
# our prediction is at different amounts of training data
output <- NULL
for (i in seq(0.1,0.9,0.1)) {
  lib <- c(1, floor(i*NROW(temp)))
  pred <- c(floor(i*NROW(temp))+1,NROW(temp))
  block_lnlp_output <- block_lnlp(temp2, lib = lib, pred = pred, columns = c("Ne",lag_cols,"Area"),target_column = "Ne", stats_only = TRUE, silent = TRUE, tp=1, theta=0, first_column_time = TRUE)
  output <- rbind(output,cbind(i,block_lnlp_output))
}  

output

# Using as little as 20% of our data set to train gives a pretty good rho! (0.9663438)
i <- 0.2
lib <- c(1, floor(i*NROW(temp2)))
pred <- c(floor(i*NROW(temp2))+1,NROW(temp2))
block_lnlp_output <- block_lnlp(temp2, lib = lib, pred = pred, columns = c("Ne",lag_cols,"Area"),
                                target_column = "Ne", stats_only = FALSE, silent = TRUE, tp=1, theta=0, method="s-map",save_smap_coefficients = TRUE, first_column_time = TRUE)
list_of_model_predictions <- block_lnlp_output$model_output
first_data_frame_of_predictions <- list_of_model_predictions[[1]]
observed <- first_data_frame_of_predictions$obs
predicted <- first_data_frame_of_predictions$pred
plot_range <- range(c(observed, predicted), na.rm = TRUE)
plot(observed, predicted, xlim = plot_range, ylim = plot_range, xlab = "Observed", 
     ylab = "Predicted", asp = 1)
abline(a = 0, b = 1, lty = 2, col = "blue")

smap_coeffs <- block_lnlp_output$smap_coefficients[[1]]

predictions <- block_lnlp_output$model_output[[1]]
t <- predictions$time

plot(t, predictions$obs, type = "l", col = "black", ylab = "x", xlab = "")
lines((t), predictions$pred, lty = 2)
legend("topright", legend = c("observed", "predicted"), lty = c(1, 2), bty = "n")

plot(t, smap_coeffs[, 1], type = "l", col = "red", ylab = "effect of Ne", xlab = "")
plot(t, smap_coeffs[, 2], type = "l", col = "blue", ylab = paste("effect of",lag_cols), xlab = "")
plot(t, smap_coeffs[, 3], type = "l", col = "blue", ylab = "effect of area", xlab = "")

# Obtaining S-map coefficients: for our actual analyses use regularized s-maps
# e.g. Cenci, S., Sugihara, G. and Saavedra, S., 2019. Regularized S‐map for inference 
# and forecasting with noisy ecological time series. Methods in Ecology and Evolution, 
# 10(5), pp.650-660.
# https://zenodo.org/record/2535688#.Xt9dvWozbwc
# This will slightly increase rho and decrease RMSE instead of standard s-maps.
# Given we are satisfied that the data gives us good predictions, 
# going to generate s-map coefficients for full dataset that can be used in
# future forecasting
lib <- c(1, NROW(temp))
pred <- c(1,NROW(temp))
block_lnlp_output <- block_lnlp(temp2, lib = lib, pred = pred, columns = c("Ne",lag_cols,"Area"),
                                target_column = "Ne", method="s-map",stats_only = FALSE, silent = TRUE, tp=1,theta = 0,save_smap_coefficients = TRUE,first_column_time = TRUE)

smap_coeffs <- block_lnlp_output$smap_coefficients[[1]]

predictions <- block_lnlp_output$model_output[[1]]
t <- predictions$time

plot(t, predictions$obs, type = "l", col = "black", ylab = "x", xlab = "")
lines((t), predictions$pred, lty = 2)
legend("topright", legend = c("observed", "predicted"), lty = c(1, 2), bty = "n")

plot(t, smap_coeffs[, 1], type = "l", col = "red", ylab = "effect of Ne", xlab = "")
plot(t, smap_coeffs[, 2], type = "l", col = "blue", ylab = paste("effect of",lag_cols), xlab = "")
plot(t, smap_coeffs[, 3], type = "l", col = "blue", ylab = "effect of area", xlab = "")

# Adding this in to future projections of habitat suitability
temp_forecast <- as.data.frame(read_table2("forecasting.txt",col_names=TRUE) %>% arrange(Time_since_Eemian))

# Adding noise to the (overly smooth) simulated data to approximate a more realistic dataset
temp_forecast$Area <- temp_forecast$Area + rnorm(length(temp_forecast$Area), sd = sd(temp_forecast$Area) * 0.2)

# Building next Ne forecasts based on S-map coefficients
for (i in 1:dim(temp_forecast)[1]) {
  lag_index <- which(names(temp2)==lag_cols)
  Ne_t <- temp2$Ne[dim(temp2)[1]]
  Ne_lag_t <- temp2[dim(temp2)[1],lag_index]
  Area_t <- temp2$Area[dim(temp2)[1]]

  Ne_coeff <- smap_coeffs$c_1[(dim(temp2)[1]-1)]
  Ne_lag_coeff <- smap_coeffs$c_2[(dim(temp2)[1]-1)]
  Area_coeff <- smap_coeffs$c_3[(dim(temp2)[1]-1)]

  c0 <- smap_coeffs$c_0[(dim(temp2)[1]-1)]

  new_Ne <- c0 + Ne_t*Ne_coeff + Ne_lag_t*Ne_lag_coeff + Area_t*Area_coeff
  
  temp_forecast$Ne[i] <- new_Ne

  # Updating temp with the new prediction
  temp <- rbind(temp,temp_forecast[i,])

  # We can now make a lagged block for Ne using 
  # these values. 
  maxtime_lag <- 30
  temp2 <- make_block(temp,max_lag=maxtime_lag,restrict_to_lib = FALSE)
  # Don't need to keep all columns
  temp2 <- temp2[,-c(1,3:31,63:91)]
  
  # going to generate s-map coefficients for full dataset that 
  # can be used in future forecasting
  lib <- c(1, NROW(temp2))
  pred <- c(1,NROW(temp2))
  block_lnlp_output <- block_lnlp(temp2, lib = lib, pred = pred, columns = c("Ne",lag_cols,"Area"),
                                  target_column = "Ne", method="s-map",stats_only = FALSE, silent = TRUE, tp=1,theta = 0,save_smap_coefficients = TRUE,first_column_time = TRUE)
  
  smap_coeffs <- block_lnlp_output$smap_coefficients[[1]]
}

predictions <- block_lnlp_output$model_output[[1]]
t <- predictions$time

plot(t, predictions$obs, type = "l", col = "black", ylab = "x", xlab = "")
lines((t), predictions$pred, lty = 2)
legend("topright", legend = c("observed", "predicted"), lty = c(1, 2), bty = "n")

plot(t, smap_coeffs[, 1], type = "l", col = "red", ylab = "effect of Ne", xlab = "")
plot(t, smap_coeffs[, 2], type = "l", col = "blue", ylab = paste("effect of",lag_cols), xlab = "")
plot(t, smap_coeffs[, 3], type = "l", col = "blue", ylab = "effect of area", xlab = "")

# Looking at difference between incorporating area and not in terms of predictive power
block_lnlp_output <- block_lnlp(temp2, lib = lib, pred = pred, columns = c("Ne",lag_cols,"Area"),
                                target_column = "Ne", method="s-map",stats_only = FALSE, silent = TRUE, tp=1,theta = 0,save_smap_coefficients = TRUE,first_column_time = TRUE)

no_area_lnlp_output <- block_lnlp(temp2, lib = lib, pred = pred, columns = c("Ne",lag_cols),
                                target_column = "Ne", method="s-map",stats_only = FALSE, silent = TRUE, tp=1,theta = 0,save_smap_coefficients = TRUE,first_column_time = TRUE)

observed_no_area <- no_area_lnlp_output$model_output[[1]]$obs
predicted_no_area <- no_area_lnlp_output$model_output[[1]]$pred

observed_w_area <- block_lnlp_output$model_output[[1]]$obs
predicted_w_area <- block_lnlp_output$model_output[[1]]$pred

plot_range <- range(c(observed_no_area, predicted_no_area), na.rm = TRUE)
plot(observed_no_area, predicted_no_area, xlim = plot_range, ylim = plot_range, xlab = "Observed", 
     ylab = "Predicted")
abline(a = 0, b = 1, lty = 2, col = "darkgrey", lwd = 2)
abline(lm(predicted_no_area ~ observed_no_area), col = "black", lty = 3, lwd = 2)

points(observed_w_area, predicted_w_area, pch = 2, col = "red")
abline(lm(predicted_w_area ~ observed_w_area), col = "red", lty = 3, lwd = 2)

legend("bottom", legend = c(paste("(Ne only) rho =", round(no_area_lnlp_output$rho, 
                                                                 2)), paste("(Ne and area) rho =", round(block_lnlp_output$rho, 2))), lty = 3, 
       lwd = 2, col = c("black", "red"), box.col = NA, xpd = TRUE)

# More formally using CCM to evaluate the influence of Area on Ne
# lib_sizes will need to be adjusted based on total size of the dataset
# No time lag
Ne_xmap_Area <- ccm(temp, E = Evalue, lib_column = "Ne", 
                        target_column = "Area", lib_sizes = seq(10, 120, by = 10), num_samples = 100, 
                        random_libs = TRUE, replace = TRUE, silent = TRUE)
Area_xmap_Ne <- ccm(temp, E = Evalue, lib_column = "Area", target_column = "Ne", 
                        lib_sizes = seq(10, 120, by = 10), num_samples = 100, random_libs = TRUE, 
                        replace = TRUE, silent = TRUE)

Ne_xmap_Area_means <- ccm_means(Ne_xmap_Area)
Area_xmap_ne_means <- ccm_means(Area_xmap_Ne)

plot(Ne_xmap_Area_means$lib_size, pmax(0, Ne_xmap_Area_means$rho), type = "l", col = "red", 
     xlab = "Library Size", ylab = "Cross Map Skill (rho)", ylim = c(0, 1))
lines(Area_xmap_ne_means$lib_size, pmax(0, Area_xmap_ne_means$rho), col = "blue")
legend(x = "topleft", legend = c("Ne xmap Area", "Area xmap Ne"), col = c("red", 
                                                                                  "blue"), lwd = 1, bty = "n", inset = 0.02, cex = 0.8)

# This graph suggests that Area depends on Ne more strongly than Ne depends on Area 
# (potentially due to the use of simulated data, and/or the failure to account for time lags, the
# next exploratory analysis)
vars <- c("Ne","Area")

# generate all combinations of lib_column, target_column, tp
params <- expand.grid(lib_column = vars, target_column = vars, tp = -20:20)

# throw out cases where lib == target
params <- params[params$lib_column != params$target_column, ]

params$E <- Evalue

output <- do.call(rbind, lapply(seq_len(NROW(params)), function(i) {
  ccm(temp, E = params$E[i], lib_sizes = NROW(temp), 
      random_libs = FALSE, lib_column = params$lib_column[i], target_column = params$target_column[i], 
      tp = params$tp[i], silent = TRUE)
}))

output$direction <- paste(output$lib_column, "xmap to\n", output$target_column)

library(ggplot2)
time_delay_ccm_fig <- ggplot(output, aes(x = tp, y = rho, color = direction)) + 
  geom_line() + theme_bw()

print(time_delay_ccm_fig)

# This figure looks similar to Fig 2 of Ye, H., Deyle, E.R., Gilarranz, L.J. and Sugihara, G., 
# 2015. Distinguishing time-delayed causal interactions using convergent cross mapping. Scientific # reports, 5, p.14750., with x = Area, and y = Ne.
# The optimal cross-map lag from Ne to Area (Ne xmap to Area) is negative
output %>% filter(lib_column=="Ne") %>% arrange(desc(rho)) %>% slice(1)
# which suggests Area causes Ne, as estimating the historcal influence of Area from
# Ne is where the optimal peak occurs.
# Conversely, the optimal mapping of Area to Ne (Area xmap to Ne) is positive (11), as changes in Area # are not reflected in changes in Ne until some point in the future.
output %>% filter(lib_column=="Area") %>% arrange(desc(rho)) %>% slice(1)

# PLOTTING TO DEMONSTRATE S-MAP PREDICTIONS
# Smap coefficients + two distinct time points that can be used to demonstrate how s-map works
# following:
# Deyle, E.R., May, R.M., Munch, S.B. and Sugihara, G., 2016. Tracking and forecasting ecosystem interactions in real time
# Proceedings of the Royal Society B: Biological Sciences, 283(1822), p.20152258.
# temp2 is 0 indexed, while smap_coeff is not, hence the need to minus one
#line1 <- which(abs(smap_coeffs$c_1/mean(smap_coeffs$c_1,na.rm=TRUE)-1)+abs(smap_coeffs$c_2/mean(smap_coeffs$c_2,na.rm=TRUE)-1)+abs(smap_coeffs$c_3/mean(smap_coeffs$c_3,na.rm=TRUE)-1)==min(abs(smap_coeffs$c_1/mean(smap_coeffs$c_1,na.rm=TRUE)-1)+abs(smap_coeffs$c_2/mean(smap_coeffs$c_2,na.rm=TRUE)-1)+abs(smap_coeffs$c_3/mean(smap_coeffs$c_3,na.rm=TRUE)-1),na.rm=TRUE))
#line1 <- which(smap_coeffs$c_1==min(smap_coeffs$c_1[(lag_index+1):dim(smap_coeffs)[1]],na.rm=TRUE))-1
#line2 <- which(smap_coeffs$c_1==max(smap_coeffs$c_1,na.rm=TRUE))-1

# This is in dataframe position. Time to Eemian=line1-1
line1 <- which((abs(smap_coeffs$c_1/mean(smap_coeffs$c_1,na.rm=TRUE)-1)+abs(smap_coeffs$c_2/mean(smap_coeffs$c_2,na.rm=TRUE)-1)+abs(smap_coeffs$c_3/mean(smap_coeffs$c_3,na.rm=TRUE)-1))==min(abs(smap_coeffs$c_1/mean(smap_coeffs$c_1,na.rm=TRUE)-1)+abs(smap_coeffs$c_2/mean(smap_coeffs$c_2,na.rm=TRUE)-1)+abs(smap_coeffs$c_3/mean(smap_coeffs$c_3,na.rm=TRUE)-1),na.rm=TRUE))
line2 <- which((abs(smap_coeffs$c_1/mean(smap_coeffs$c_1,na.rm=TRUE)-1)+abs(smap_coeffs$c_2/mean(smap_coeffs$c_2,na.rm=TRUE)-1)+abs(smap_coeffs$c_3/mean(smap_coeffs$c_3,na.rm=TRUE)-1))==max(abs(smap_coeffs$c_1/mean(smap_coeffs$c_1,na.rm=TRUE)-1)+abs(smap_coeffs$c_2/mean(smap_coeffs$c_2,na.rm=TRUE)-1)+abs(smap_coeffs$c_3/mean(smap_coeffs$c_3,na.rm=TRUE)-1),na.rm=TRUE))

ggplot() + geom_vline(xintercept=line1-1,size=2, col="light grey") +
    geom_vline(xintercept=line2-1,size=2,col="light grey") +
    geom_line(data=smap_coeffs,aes(x=temp2$Time_since_Eemian,y=c_1/mean(smap_coeffs$c_1,na.rm=TRUE)-1), size=1,color="dark green") +
    geom_line(data=smap_coeffs,aes(x=temp2$Time_since_Eemian,y=c_2/mean(c_2,na.rm=TRUE)-1), size=1,color="red") + 
    geom_line(data=smap_coeffs,aes(x=temp2$Time_since_Eemian,y=c_3/mean(c_3,na.rm=TRUE)-1), size=1,color="blue") + 
  theme_bw(base_size=14) +
  xlab("Time since Eemian") +
  ylab("Relative interaction coefficient strength")

ggsave("relative_coefficients.pdf",device="pdf",width=12,height=4,units="in")
  
# Plotting original time series
ggplot() + geom_vline(xintercept=line1-1,size=2, col="light grey") +
  geom_vline(xintercept=line2-1,size=2,col="light grey") +
  geom_line(data=smap_coeffs,aes(x=temp2$Time_since_Eemian,y=temp2$Ne), size=1,color="dark green") +
  xlab("Time since Eemian") +
  ylab("Relative magnitude") + 
  theme_bw(base_size=14) +
  xlab("Time since Eemian") +
  ylab("Ne")

ggsave("Ne.pdf",device="pdf",width=12,height=2,units="in")

ggplot() + geom_vline(xintercept=line1-1,size=2, col="light grey") +
  geom_vline(xintercept=line2-1,size=2,col="light grey") +
  geom_line(data=smap_coeffs,aes(x=temp2$Time_since_Eemian,y=temp2$Area), size=1,color="blue") +
  xlab("Time since Eemian") +
  ylab("Relative magnitude") + 
  theme_bw(base_size=14) +
  xlab("Time since Eemian") +
  ylab("Relative area")

ggsave("Area.pdf",device="pdf",width=12,height=2,units="in")

# Will use line1 and line2 points as demonstration points based on contrasting smap coefficients
line1values <- c(temp2$Ne[line1],temp2[line1,lag_index],temp2$Area[line1],temp2$Ne[line1+1])
line2values <- c(temp2$Ne[line2],temp2[line2,lag_index],temp2$Area[line2],temp2$Ne[line2+1])

# 3D graph of time lag
# Values for axes increase in the direction of the arrows
# Time series coloured from past (blue) to present (red)
# Have to run the code twice to get the full extent of the line
# for some strange reason!
pdf("3D.pdf")

scatter3D(temp2$Ne,temp2[,lag_index],temp2$Area,type="l",
          colvar = as.integer(temp2$Time_since_Eemian),
          phi=10,theta=20,lwd=4,colkey=FALSE,
          xlab="Ne",ylab=gsub("_","-",lag_cols),zlab="Area")
scatter3D(temp2$Ne,temp2[,lag_index],temp2$Area,type="l",
          colvar = as.integer(temp2$Time_since_Eemian),
          phi=10,theta=20,lwd=4,colkey=FALSE,
          xlab="Ne",ylab=gsub("_","-",lag_cols),zlab="Area",add=TRUE)
scatter3D(line1values[1],line1values[2],line1values[3], pch=19,
          colvar = as.integer(temp2$Time_since_Eemian),
          phi=10,theta=20,lwd=4,colkey=FALSE, cex=2,
          xlab="Ne",ylab=gsub("_","-",lag_cols),zlab="Area",add=TRUE)
scatter3D(line2values[1],line2values[2],line2values[3], pch=19,
          colvar = as.integer(temp2$Time_since_Eemian),
          phi=10,theta=20,lwd=4,colkey=FALSE,cex=2,
          xlab="Ne",ylab=gsub("_","-",lag_cols),zlab="Area",add=TRUE)

dev.off()


# Likelihood surfaces
Ne_lag_range <- unique(temp2[,lag_index])
Area_t_range <- unique(temp2$Area)

# rows to remove because they have NA data
rows_to_keep <- which(!is.na(Ne_lag_range) & !is.na(Area_t_range))
Ne_lag_range <- Ne_lag_range[rows_to_keep]
Area_t_range <- Area_t_range[rows_to_keep]

# Obtaining the coordinates to plot the surface
surface_shape <- matrix(NA,nrow=length(Area_t_range),ncol=length(Ne_lag_range))
for (Ne_lag_t in 1:length(Ne_lag_range)) {
  for (Area_t in 1:length(Area_t_range)) {
    # For each of the combination of values in Ne_lag_range and Area_t_range
    # filtering temp2 to see if the combination exists in "actual" data
    matrixentry <- temp2 %>% filter(!!as.name(lag_cols)==Ne_lag_range[Ne_lag_t]) %>% filter(Area==Area_t_range[Area_t])
    # If it does
    if(dim(matrixentry)[1]>0) {
      # But it is NOT the last line in temp2!
      if(matrixentry$Time_since_Eemian!=temp2$Time_since_Eemian[(dim(temp2)[1])]) {
        # Then taking the observed value of Ne(t+1) corresponding to Ne-lag(t) and Area(t)
        surface_shape[Area_t,Ne_lag_t] <- temp2$Ne[which(temp2$Time_since_Eemian==(as.numeric(matrixentry$Time_since_Eemian+1)))]
      }
    } else {
      # Need to estimate the value of Ne(t+1) corresponding to the values given
      rowindex <- as.numeric(temp2 %>% filter(!!as.name(lag_cols)==Ne_lag_range[Ne_lag_t]) %>% select(Time_since_Eemian))+1
      Ne_coeff <- smap_coeffs$c_1[(rowindex-1)]
      # If Ne_coeff is not NA (i.e. not the first row etc)
      if (!is.na(Ne_coeff)) {
        Ne_lag_coeff <- smap_coeffs$c_2[(rowindex-1)]
        Area_coeff <- smap_coeffs$c_3[(rowindex-1)]
        c0 <- smap_coeffs$c_0[(rowindex-1)]
        Ne_t <- temp2$Ne[rowindex]
        surface_shape[Area_t,Ne_lag_t] <- c0 + Ne_t*Ne_coeff + Ne_lag_range[Ne_lag_t]*Ne_lag_coeff + Area_t_range[Area_t]*Area_coeff
      }
    }
  }  
}    

x <- matrix(Ne_lag_range,nrow=length(Area_t_range),ncol=length(Ne_lag_range),byrow = TRUE)
y <- matrix(Area_t_range,ncol=length(Area_t_range),nrow=length(Ne_lag_range))
z <- surface_shape

c(max(x,na.rm=TRUE),min(x,na.rm=TRUE))
max(x,na.rm=TRUE)-min(x,na.rm=TRUE)
c(max(y,na.rm=TRUE),min(y,na.rm=TRUE))
max(y,na.rm=TRUE)-min(y,na.rm=TRUE)
c(max(z,na.rm=TRUE),min(z,na.rm=TRUE))
max(z,na.rm=TRUE)-min(z,na.rm=TRUE)

# Plotting the full extent
pdf("full_correlation_surface.pdf")

scatter3D(x,y,z,col=NULL,cex=0.0001,theta=210,xlab=paste(gsub("_","-",lag_cols),"(t)",sep=""),ylab="Area(t)",zlab="Ne(t+1)",colkey=FALSE)
surf3D(x,y,z,border="black",add=TRUE,colkey=FALSE,col="grey")
scatter3D(line1values[2],line1values[3],line1values[4],col="white",cex=2,add=TRUE,pch=19)
scatter3D(line2values[2],line2values[3],line2values[4],col="white",cex=2,add=TRUE,pch=19)

dev.off()

points_to_plot <- rbind(line1values,line2values)

for (j in 1:2) {
  valueindex <- points_to_plot[j,]
  keep_rows <- NULL
  for (i in 1:dim(x)[1]) {
    if (any((x[i,] < (valueindex[2]+(max(temp2[,lag_index],na.rm = TRUE)/10))) & x[i,] > (valueindex[2]-(max(temp2[,lag_index],na.rm = TRUE)/10)))) {
      if(any((y[i,] < (valueindex[3]+(max(temp2$Area,na.rm = TRUE)/10))) & y[i,] > (valueindex[3]-(max(temp2$Area,na.rm = TRUE)/10)))) {
        keep_rows <- cbind(keep_rows,i)
      }
    }
  }

  keep_cols <- NULL
  for (i in 1:dim(x)[2]) {
    if (any((x[,i] < (valueindex[2]+(max(temp2[,lag_index],na.rm = TRUE)/10))) & x[,i] > (valueindex[2]-(max(temp2[,lag_index],na.rm = TRUE)/10)))) {
      if(any((y[,i] < (valueindex[3]+(max(temp2$Area,na.rm = TRUE)/10))) & y[,i] > (valueindex[3]-(max(temp2$Area,na.rm = TRUE)/10)))) {
        keep_cols <- cbind(keep_rows,i)
      }
    }
  }

  displayx <- x[keep_rows,keep_cols]
  displayy <- y[keep_rows,keep_cols]
  displayz <- z[keep_rows,keep_cols]

  c(max(displayx,na.rm=TRUE),min(displayx,na.rm=TRUE))
  max(displayx,na.rm=TRUE)-min(displayx,na.rm=TRUE)
  c(max(displayy,na.rm=TRUE),min(displayy,na.rm=TRUE))
  max(displayy,na.rm=TRUE)-min(displayy,na.rm=TRUE)
  c(max(displayz,na.rm=TRUE),min(displayz,na.rm=TRUE))
  max(displayz,na.rm=TRUE)-min(displayz,na.rm=TRUE)
  
  if (j==1) {
    circle_size <- 3
  } else {
    circle_size <- 3
  }
  
  pdf(paste("line",j,"values.pdf",sep=""))
  
  scatter3D(displayx,displayy,displayz,theta=180,cex=0.00001,xlab=paste(gsub("_","-",lag_cols),"(t)",sep=""),ylab="Area(t)",zlab="Ne(t+1)",colkey=FALSE)
  surf3D(displayx,displayy,displayz,add=TRUE,border="black",col="grey",colkey=FALSE)
  scatter3D(valueindex[2],valueindex[3],valueindex[4],col="black",cex=2,add=TRUE,pch=19)
  
  dev.off()
}
