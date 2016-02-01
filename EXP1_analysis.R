# Thermocouple lab analysis
# ME 411
# 2016-1-18
setwd("~/Documents/Class Documents/Current Classes/Link to Instrumentation (ME 411)/EXP1")
library("scales")
# pdf("EXP1_plots.pdf")
png("EXP1_plot%03d.png")

#### read in the data ####
#all values are kept in SI
Tamb <- 22.5
T1 <- read.csv("Test12016-01-26 13-51.txt", skip=3, sep="\t")
T2 <- read.csv("Test22016-01-26 13-49.txt", skip=3, sep="\t")
T3 <- read.csv("Test32016-01-26 13-53.txt", skip=3, sep="\t")
T4 <- read.csv("Test42016-01-26 13-58.txt", skip=3, sep="\t")
T5 <- read.csv("Test52016-01-26 14-01.txt", skip=3, sep="\t")

#returns the indices of vec for which the value of vec is closest to targ
near <- function(vec, targ) sapply(targ, function(targ) which.min(abs(vec-targ)) )

#### initialize the main data frame ####
size <- c(rep("small",8), rep("large",8))
heat <- c(rep("water",4), rep("hotplate",4), rep("water",4), rep("hotplate",4))
thermo <- data.frame(
	size, heat, data=NA, tau=NA, adjR2=NA, model=NA, 
	rho=NA, heatCap=NA, R=NA, L=NA, h.tau.conversion=NA, h.trans= NA
	)
for (i in 1:length(thermo$size)){
	if ((thermo$size == "small")[i]){#average of copper and nickel
		thermo$R[i] <- 0.5e-3/2
		thermo$heatCap[i] <- (0.385e3+0.440e3)/2
		thermo$rho[i] <- (8960+8908)/2
		thermo$L[i] <- 1.37e-3
		thermo$h.tau.conversion[i] <- (
			thermo$rho[i]*thermo$heatCap[i]*thermo$R[i]*thermo$L[i]/
				(2*(thermo$R[i]+thermo$L[i]))
		)
	}
	else {#average of copper, nickel, and lead
		thermo$R[i] <- 3.34e-3/2
		thermo$heatCap[i] <- (0.385e3+0.440e3+0.129e3)/3
		thermo$rho[i] <- (8960+8908+11340)/3
		# thermo$L[i] <- 1.37e-3
		thermo$h.tau.conversion[i] <- (
			thermo$rho[i]*thermo$heatCap[i]*thermo$R[i]/3
		)
	}
}

#### organize each trial as an observation ####
# plot(T1$Time, log(T1$Ch-Tamb), type="l")# small on water
interval <-near(T1$Time,c(13.5, 14.5))
runholder <- list(T1[interval[1]:interval[2], ])
# thermo$data[1] <- list(runholder)

# plot(T2$Time, log(T2$Ch-Tamb), type="l")# small on water
plot(
	T2$Time, log(T2$Ch-Tamb), type="l", 
	xlim=c(15, 20), xlab="time (s)", ylab= "log(T-Tamb)"
	)
interval <- append(interval, near(T2$Time,c(15.3, 16.25, 39.8, 40.2, 68.6, 69.5)))
for (i in c(3,5,7))
	runholder <- append(runholder, list(T2[interval[i]:interval[i+1], ]))

# plot(T3$Time, log(T3$Ch-Tamb), type="l")# small on hotplate
interval <- near(T3$Time, c(13, 15, 37, 39, 59, 61, 82.5, 85))
for (i in c(1,3,5,7))
	runholder <- append(runholder, list(T3[interval[i]:interval[i+1], ]))

# plot(T4$Time, log(T4$Ch-Tamb), type="l")# large on water
interval <- near(T4$Time, c(28, 35, 85, 89, 115, 120, 150, 153))
for (i in c(1,3,5,7))
	runholder <- append(runholder, list(T4[interval[i]:interval[i+1], ]))

# plot(T5$Time, log(T5$Ch-Tamb), type="l")# large on hotplate
interval <- near(T5$Time, c(55, 70, 128, 140, 208, 220, 277.5, 290))
for (i in c(1,3,5,7))
	runholder <- append(runholder, list(T5[interval[i]:interval[i+1], ]))

for ( i in 1:16)#read the runs into the data frame
	thermo$data[i] <- runholder[i]

#### find the time constants ####
par.old <- par(mfrow= c(2,2))
for (i.trial in 1:16){
# i.trial <- 2
# {
	temps <- thermo$data[[i.trial]]$Ch
	times <- thermo$data[[i.trial]]$Time
	temps <- (temps-Tamb)/(max(temps)-Tamb)#just shift; normalizing only changes intercept
	dat <- data.frame(times, temps)
	x <- plot(
		times, log(temps),
		xlab= "time (s)", ylab= "log(T*)",
		bty="n", pch=19, col=alpha("black", alpha=0.1),
		main= paste(
			thermo$size[i.trial], 
			"thermocouple \nheated by", 
			thermo$heat[i.trial]
			)
		)
	model <- lm(log(temps)~times, dat)
	abline(model, col="red", lwd=3)
	legend("topright", 
	       legend=c(
		       	paste("tau= ", formatC(-1/model$coeff["times"], digits=3)),
		       	paste("adj. R sq. = ", formatC(summary(model)$adj, digits=3))
		       )
		)
	thermo$tau[i.trial] <- -1/model$coeff["times"]
	thermo$adjR2[i.trial] <- summary(model)$adj
	thermo$model[i.trial] <- list(model)
}
par(par.old)

#### post processing ####
for (i in 1:length(thermo$size))
	thermo$h.trans[i] <- thermo$h.tau.conversion[i]/thermo$tau[i]

plot(#box and whisker plot of time constants
	factor(paste(thermo$size, "\nwith", thermo$heat)), 
	thermo$tau,
	ylab= "time constant (s)",
	main= "distributions of time constants\nfor various thermocouple sizes\nand heating methods"
	)
plot(#box and whisker plot of transfer coefficients
	factor(paste(thermo$size, "\nwith", thermo$heat)), 
	thermo$h.trans,
	ylab= "convective transfer coefficient (W/m^2/K)",
	main= "distributions of transfer coefficients\nfor various thermocouple sizes\nand heating methods"
)
write.csv(x = data.frame(
		thermo$size, 
		thermo$heat, 
		tau= formatC(thermo$tau, digits=4), 
		adjR2= formatC(thermo$adjR2, digits=4),
		h.conv= formatC(thermo$h.trans, digits=4),
	       	rho= formatC(thermo$rho, digits=4),
		heatCap= formatC(thermo$heatCap, digits=4),
       		Radius= formatC(thermo$R, digits=4),
       		Length= formatC(thermo$L, digits=4),
       		h.to.tau.conversion= formatC(thermo$h.tau.conversion, digits=4)
		), file = "reduction.csv")
dev.off()

