rotate <- function(unit, theta, center=c(0,0)) {
	x <- (unit$x-center[1])*cos(theta)-(unit$y-center[2])*sin(theta)+center[1]
	y <- (unit$x-center[1])*sin(theta)+(unit$y-center[2])*cos(theta)+center[2]
	return(data.frame(x, y))
}

shift <- function(unit, r) data.frame(x=unit$x+r[1], y=unit$y+r[2])

hrev <- function(unit) return(data.frame(x=rev(unit$x), y=rev(unit$y)))

hilbert <- function(N){
	if(N==0) return(data.frame(x=0.5, y=0.5))
	else{
		unit <- hilbert(N-1)
		newUnit <- rbind(
			hrev(shift(rotate(unit, pi/2, c(0.5,0.5))/2, c(0,0.5))),#upper left
			unit/2,#lower left
			shift(unit/2, c(0.5,0)),#lower right
			hrev(shift(rotate(unit, -pi/2, c(0.5,0.5))/2, c(0.5,0.5)))#upper right
		)
		return(newUnit)
	}
}

colLines <- function(pts, ...){
	xh <- pts$x
	yh <- pts$y
	arrows(
		xh[-length(xh)], yh[-length(yh)], 
		xh[-1], yh[-1], 
		length=0, col=rainbow(n=length(xh)),
		...
	)
}

pdf("hilbertCurve.pdf")
par(mar=c(0,0,0,0), xpd=NA)

par(bg="black", mfrow=c(3,3))
plot(0.5,0.5, col="white", pch=20, bty="n", axes=F, xlab="", ylab="")
for (n in 1:8){
	plot.new()
	plot.window(xlim=c(0,1), ylim= c(0,1), asp=1)
	if (n<4) points(hilbert(n), col="white", pch=20)
	colLines(hilbert(n))
}

par(mfrow=c(2,2))
plot(0.5,0.5, col="white", pch=20, bty="n", axes=F, xlab="", ylab="")
i <- c(2,2,1)
j <- c(1,2,2)
for (n in 1:3){
	par(mfg=c(i[n],j[n]))
	plot.new()
	plot.window(xlim=c(0,1), ylim= c(0,1), asp=1)
	points(hilbert(n), col="white", pch=20)
	colLines(hilbert(n))
}

par(bg="white", mfrow=c(1,1))
plot.new()
plot.window(xlim=c(0,1), ylim= c(0,1), asp=1)
lines(hilbert(5), col=rgb(1,0,0), lwd=0.002*96)
border <- data.frame(x=c(0,1,1,0,0), y=c(0,0,1,1,0))
border <- shift(1.04*border, c(-0.02, -0.02))
lines(border, col=rgb(1,0,0), lwd=0.002*96)

dev.off()
