ll <- function(x, op, om, pm,AP,AM) {
    # X is offset for male
    #      slope for male
    #      offset for female 
    #      slope for female

    # AP is the age of the father
    # AM                   mother
    pi <- x[1] + x[2]*AP
    mu <- x[3] + x[4]*AM
    l <- 0
    for (i in 1:length(op)) {
        J <- 0:(pm[i] - op[i] - om[i])
        pp <- dpois(op[i] + J, pi[i])
        mp <- dpois(pm[i] - (op[i] + J), mu[i])
        w <- pp*mp
        l <- l + log(sum(w))
    }
    return(l)
}

llf <- function(x, ...) {
    return(-ll(x, ...))
}
llpat <- function(x, ...) {
    # omitting the female slope
    return(-ll(c(x[1],x[2],x[3],0), ...))
}
llmat <- function(x, ...) {
    # omitting the male slope
    return(-ll(c(x[1], 0, x[2],x[3]), ...))
}

fit_model <- function(temp,starting_values,...){
    min_value <- Inf
    min_value_pat <- Inf
    min_value_mat <- Inf
    for (sv in 1:nrow(starting_values)){
        rfull_one_sv <- nlm(llf, 
                            as.numeric(starting_values[sv,]), 
                            print.level=0, 
                            hessian=T,
                            temp$father, 
                            temp$mother, 
                            temp$total,
                            temp$fathers_age,
                            temp$mothers_age,
                            ...)
        rfull_pat_one_sv <- nlm(llpat, 
                            as.numeric(starting_values[sv,]), 
                            print.level=0, 
                            hessian=T,
                            temp$father, 
                            temp$mother, 
                            temp$total,
                            temp$fathers_age,
                            temp$mothers_age,
                            ...)

        rfull_mat_one_sv <- nlm(llmat, 
                            as.numeric(starting_values[sv,]), 
                            print.level=0, 
                            hessian=T,
                            temp$father, 
                            temp$mother, 
                            temp$total,
                            temp$fathers_age,
                            temp$mothers_age,
                            ...)
        if (rfull_one_sv$minimum< min_value){
            rfull <- rfull_one_sv
            min_value <- rfull$minimum
        }
        if (rfull_pat_one_sv$minimum< min_value_pat){
            rfullpat <- rfull_pat_one_sv
            min_value_pat <- rfullpat$minimum
        }
        if (rfull_mat_one_sv$minimum< min_value_mat){
            rfullmat <- rfull_mat_one_sv
            min_value_mat <- rfullmat$minimum
        }
    }
    #
    maternal_slope_pvalue  <- pchisq(2*(min_value_pat - min_value),1,lower.tail=FALSE)
    paternal_slope_pvalue  <- pchisq(2*(min_value_mat - min_value),1,lower.tail=FALSE)
    covar_approx <- diag(solve(rfull$hessian))
    estimates <- data.frame(estimates = rfull$estimate,
                            conf_low = -1.96*sqrt(covar_approx)+rfull$estimate,
                            conf_high = +1.96*sqrt(covar_approx)+rfull$estimate,
                            pvalue= c(NA,paternal_slope_pvalue,NA,maternal_slope_pvalue))
    estimates$coefficent <- c('intercept_father','slope_father','intercept_mother','slope_mother')
    return(list(fitted_model=rfull,
                estimates=estimates,
                fitted_model_pat=rfullpat,
                fitted_model_mat=rfullmat))
}

#read data
dat <- read.table('./totalPerPn', header=TRUE)

#removing na if there are any
dat <- transform(dat,pat=ifelse(is.na(pat),0,pat))
dat <- transform(dat,mat=ifelse(is.na(mat),0,mat))

#Plotting the data to see where relationship with denom stops being linear
library(ggplot2)
ggplot(data=dat,aes(x=denom,y=total))+geom_point()+geom_smooth(col='blue')+geom_smooth(method='lm',col='red')+theme_bw()
sdat <- subset(dat,denom<4.5e5)
ggplot(data=sdat,aes(x=denom,y=total))+geom_point()+geom_smooth(col='blue')+geom_smooth(method='lm',col='red')+theme_bw()

#integrating out denom
correction_factor <- (lm(total~denom,data=sdat))
sdat$total_corrected <- as.integer(residuals(correction_factor)+abs(min(residuals(correction_factor),na.rm=T)))

#fit model
out <- fit_model(temp=data.frame(father=sdat$pat,
                                 mother=sdat$mat,
                                    total=sdat$total_corrected,
                                 fathers_age=sdat$AP,
                                 mothers_age=sdat$AM),
                 starting_values=data.frame(x=1,
                                            y=1.6,
                                            c=1,
                                            d=0.3))

#Extrapolate effect (have total of 1394292 microsatellites)
ratio <- mean(dat$denom)/1394292
gw_effect_paternal <- out$estimates[2,1]/ratio
gw_effect_maternal <- out$estimates[4,1]/ratio
