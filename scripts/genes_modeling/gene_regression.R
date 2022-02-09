



genedf <- rbind(genedf, data.frame(name=rownames(data), condition=cond, timepoint=timepoint,                                      slope=NA, monotonic=NA, treatment=line, time=time,
                                   mean=rowMeans(data, na.rm=TRUE), sd=rowSds(data, na.rm=TRUE)))
genelm <- genedf[genedf$name==gene & genedf$condition %in% c("UTT","UTTT","UTTTT"),]
linearMod <- lm(mean ~ time, data=genelm)
genedf[genedf$name==gene & genedf$condition %in% c("UTT","UTTT","UTTTT"),"slope"] <- linearMod$coefficients[2] 



?smooth.spline()
number_of_cycles = 2
max_y = 40

x = 1:500
a = number_of_cycles * 2*pi/length(x)

y = max_y * sin(x*a)
noise1 = max_y * 1/10 * sin(x*a*10)
y <- y + noise1

plot(x, y, type="l", ylim=range(-1.5*max_y,1.5*max_y,5), lwd = 5, col = "green")
lowpass.spline <- smooth.spline(x,y, spar = 0.6) ## Control spar for amount of smoothing
lowpass.loess <- loess(y ~ x, data = data.frame(x = x, y = y), span = 0.3) ## control span to define the amount of smoothing

lines(predict(lowpass.spline, x), col = "red", lwd = 2)
lines(predict(lowpass.loess, x), col = "blue", lwd = 2)
