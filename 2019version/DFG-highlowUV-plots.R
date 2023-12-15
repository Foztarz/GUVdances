#clean dataspace
rm(list = ls())
graphics.off()
#################################################################
#	Load & Install												#
#################################################################
chooseCRANmirror(graphics=FALSE, ind = which(getCRANmirrors()$Country =='Germany'))#'Sweden'))#
install.packages(c('rstan', 'circular', 'brms'))
library(rstan) # observe startup messages
library(circular)
library(brms)
options(mc.cores = parallel::detectCores()-1)

#################################################################
#	Useful Functions											#
#################################################################
# plot a circular histogram with stacked points, a mean vector and ±1 s.d. error bars
Cplot <- function(headings, sp, bt, ax, ...){
	#fit mean axis, fits mean direction unless otherwise specified
	if(missing(ax)){ax <- F}
	#spacing of stacked points, now automatically stacks towards centre unless otherwise specified
	if(missing(sp) & missing(bt)){sp <- 0.04}
	#bt specifies the stacking by a multipicative factor, 1 = stacked, 2 = 1 point's space between, 0.5 = half overlapping
	if( missing(sp) & !(missing(bt)) ){sp <- bt*.04}
	#	Get functions and packages
	if(sum(rownames(installed.packages()) %in% c('CircStats', 'circular'), na.rm = T)<2){
		install.packages(c('CircStats', 'circular'))
	}#if(sum(rownames(installed.packages()) %in% c('CircStats', 'circular'), na.rm = T)<2)
	if(!(	sum('CircCI'%in% ls())	)){
		CircCI <- function(mn, lci, uci, out, zro, drc, lng, ...){
			if(missing(lng)){lng<-10*360/5};	if(missing(drc)){drc<-'clock'}
			if(missing(zro)){zro <- pi/2};if(missing(out)){out <- 0.05}
			if(missing(uci)){uci <- lci}
			lwr <- mn - lci;	upr <- mn + uci
			circ.pos <- ( ((drc == 'clock')-1)*2 +1) * 
				-seq( pi*lwr/180, pi*upr/180, length.out = lng) + zro
			circ.x <- cos(circ.pos)*(1+out);	circ.y <- sin(circ.pos)*(1+out)
			lines(circ.x, circ.y, ...)
			lines.circular( as.circular(rep(lwr,2),units = 'degrees', 
				type = 'angles', modulo = '2pi', zero = zro, 
				rotation = drc, template = 'none'), 
				out*c(0.5, 1.5), modulo = '2pi', 
				zero = zro, rotation = drc, ...)
			lines.circular(as.circular(rep(upr,2),units = 'degrees', 
				type = 'angles', modulo = '2pi', zero = zro, 
				rotation = drc, template = 'none'),
			 	out*c(0.5, 1.5), modulo = '2pi', zero = zro, 
			 	rotation = drc, ...)
		}#CircCI <- function(mn, lci, uci, out, zro, drc, lng, ...)
	}#if(!(	sum('CircCI'%in% ls())	))
	if(!(	sum('mycirc'%in% ls())	)){
		mycirc <- function(angles, clock){
			if(missing(clock)){clock <- T}
			if(clock){
			return(		as.circular(angles,units='degrees',type='angles',
			modulo='2pi',zero=pi/2,rotation='clock',	template='none')	)
				}else{
				as.circular(angles,units='degrees',type='angles',
				 modulo='2pi',zero=pi/2,rotation='counter',template='none')
				}#if(clock)
			}#mycirc <- function(angles, clock)
	}#if(!(	sum('mycirc'%in% ls())	))
	#circular plot settings
	increments <- 5 #degrees
	zr <- pi/2 #start at top of screen (pi*	90	/180)
	bn <- 10*10*360/5 #bins 	
	degrad <- 180/pi #conversion from radians to degrees
	tcl <- rgb(1,1,1,0)#transparent colour
	pcl <- rgb(.3,.1,.1,.5)#point colour
	#plot characters
	lw <- 0.5 #line width
	pnt <- 2.5 #point size
	arw <- 10 #arrowhead angle
	arl <- 0.1 #arrowhead length
	#	set up input variables
	hd <- mycirc(headings)
	sm <- summary(hd)
	sv <- degrad*sd.circular(hd, na.rm=T)
	lbl <- 90*(1:4-1)
	plot(hd, col=tcl, main="", zero=zr, axes=F, shrink=1,tol=0.075)
	axis.circular(1, at = mycirc(lbl), labels = paste0(lbl, 'º'))
	par(new=T)
	plot.circular(hd, col=tcl,main="",zero=zr,axes=F,shrink=1.05,tol=0.075)
	points(hd,stack=T,bin=bn,sep=-sp,zero=zr,...)
	if(!(ax)){
		arrows.circular( mycirc(sm['Mean']),zero=zr,col='red4',lwd=3,
		 length=arl,angle=arw,shrink = sm['Rho'])
		 CircCI(sm['Mean'], sv, out = 0.15, zro=zr, drc='clock',col='red4',lwd=1)	}else{
		 sm2 <- summary(mycirc(hd*2))
		 sv2 <- degrad*sd.circular(hd*2, na.rm=T)/2
		 arrows.circular( mycirc(sm2['Mean']/2),zero=zr,col='red4',lwd=3,
		 length=arl,angle=arw,shrink = sm2['Rho'])
		 arrows.circular( mycirc(180+sm2['Mean']/2),zero=zr,col='red4',lwd=3,
		 length=arl,angle=arw,shrink = sm2['Rho'])
		 CircCI(sm2['Mean']/2, sv2, out = 0.15, zro=zr, drc='clock',col='red4',lwd=1)
		 CircCI(180+sm2['Mean']/2, sv2, out = 0.15, zro=zr, drc='clock',col='red4',lwd=1)
	 }#if(!(ax))
###################	END OF FUNCTION	###########################	
}

#################################################################
#	Read in the Data											#
#################################################################
cd <- read.table(paste0(Sys.getenv('HOME'),'/Dropbox/DFG2019-reapplication/colourdances', '.txt') , header = T, sep  = '\t')
cd$bee <- as.factor(cd$bee); cd$date <- as.factor(cd$date)
length(levels(cd$bee))#169 bees
head(cd)
#for the purposes of this analysis, select only the
# "Green low intensity" and "UV low intensity" conditions

MEmises_dat <- subset(cd, cd$light_type %in% c('ul', 'uh'))
#convert to radians
MEmises_dat$angles <- MEmises_dat$bearing * pi/180
		
#convert to (-pi, pi)
MEmises_dat1 <- MEmises_dat
#make it between -pi & pi
MEmises_dat1$angles[MEmises_dat1$angles > pi] <-
				-(2*pi - MEmises_dat1$angles[MEmises_dat1$angles > pi])
			
# a lot of animals didn't experience both conditions
#select only those that did
b_both <- unique(subset(MEmises_dat1$bee, MEmises_dat1$light_type == 'ul'))[unique(subset(MEmises_dat1$bee, MEmises_dat1$light_type == 'ul')) %in% unique(subset(MEmises_dat1$bee, MEmises_dat1$light_type == 'uh'))]

MEmises_dat1 <- subset(MEmises_dat1, MEmises_dat1$bee %in% b_both)
#factors should now only have levels found in the dataset
MEmises_dat1$bee <- factor(MEmises_dat1$bee)
MEmises_dat1$light_type <- factor(MEmises_dat1$light_type)

bees <- 1:length(levels(MEmises_dat1$bee))
cll <-  c('magenta', 'purple', 'green1', 'green4')
#################################################################
#	Load model													#
#################################################################
# form1 <- bf( angles ~ light_type + (1|bee), kappa ~ light_type + (1|bee))		
load(file = paste0(Sys.getenv('HOME'), '/Dropbox/DFG2019-reapplication/', 'br.data.uvMEmises-farbetanzen-UV.Rdata'))

#intercept should be something like
round(mean.circular(br.data$Y[!as.logical(br.data$X[,2])]),2)
round(A1inv(rho.circular(br.data$Y[!as.logical(br.data$X[,2])])),2)

#load fitted model
load(file = paste0(file = paste0(Sys.getenv('HOME'), '/Dropbox/DFG2019-reapplication/', '3rd_attempt.mod.Rdata')))#'longrunUV.mod.Rdata')))	
load(file = paste0(file = paste0(Sys.getenv('HOME'), '/Downloads/','longrunUV.mod-2.Rdata')))
#inspect
print(fit, pars = c('b_vec','b', 'mu_vec', 'b_Intercept', 'sd_1', 'b_kappa_Intercept', 'b_kappa', 'sd_2'))
#plot if needed
traceplot(fit, pars = c('b_vec','b', 'mu_vec', 'b_Intercept', 'sd_1', 'b_kappa_Intercept', 'b_kappa', 'sd_2'))
#chain 2 was problematic, it's possible that priors might avoid bad choices of mu intercept

#################################################################
#	Post-processing												#
#################################################################
fitted <- extract(fit)	#extract draws
#could consider excluding chain 3?
#circular median gives intercept angle
med.mu <- median(circular(fitted$b_Intercept))#robust estimate of mean angle
q.mu <- quantile(circular(fitted$b_Intercept), c(0.025, 0.975) )#CI of mean angle

#kappa converged well, let's focus on that
med.kappa <- median(fitted$b_kappa_Intercept)#robust estimate of kappa
q.k <- quantile(fitted$b_kappa_Intercept, c(0.025, 0.975) )#CI off kappa

#effect of change in condition on circular mean
med.b <- median(circular(fitted$b))
q.b <- quantile(circular(fitted$b), c(0.025, 0.975) )#CI of mean angle

#effect of change in condition on kappa
med.b.kappa <- median(fitted$b_kappa)#robust estimate of kappa
q.b.k <- quantile(fitted$b_kappa, c(0.025, 0.975) )#CI off kappa

#random effects arc for circular mean
med.sd <- median(circular(fitted$sd_1))#same as median(fitted$sd_1)
q.sd <- quantile(circular(fitted$sd_1), c(0.025, 0.975) )#CI of mean angle

#effect of change in condition on kappa
med.b.kappa.sd <- median(fitted$sd_2)#robust estimate of kappa
q.b.k.sd <- quantile(fitted$sd_2, c(0.025, 0.975) )#CI off kappa

levels(MEmises_dat1$light_type)
#intercept condition is 'uv high'
#contrast condition is 'uv low'
dev.new(width = 3); par(mai = c(0.9,0.9,0,0))
plot(1:2, A1(exp(med.kappa+c(0,med.b.kappa))), xlim = c(0,3)+c(1,-1)*c(0.7,0.7), ylim = c(0,1), col = c('magenta', 'purple'), lwd = 2, cex = 0.25, axes = F, xlab = 'Condition', ylab = 'Mean Vector Length')
polygon(c(1:2, rev(1:2)), A1(exp( c(c(q.k[1], q.k[1]+q.b.k[1]) -1.96* med.b.kappa.sd, rev(c(q.k[2], q.k[2]+q.b.k[2]))+1.96* med.b.kappa.sd) )), col = rgb(1,0,1,0.1), border = F)
for(b in bees){
	bb <- subset(MEmises_dat1, bee == levels(MEmises_dat1$bee)[b])
	k1 <- mle.vonmises(with(bb,subset(angles, light_type == 'uh')), bias=T)$kappa
	k2 <- mle.vonmises(with(bb,subset(angles, light_type == 'ul')), bias=T)$kappa
	points(1:2,A1(c(k1,k2)), col = rgb(0,0,0,0.1), pch = 19, cex = 0.5)
	lines(1:2, A1(c(k1,k2)), col = rgb(0,0,0,0.1), lwd = 0.25)
	# lines(1:2, A1(exp(med.kappa+c(0,med.b.kappa))+median(fitted$r_2_kappa_1[,b])), col = rgb(0,0,0,0.1), lwd = 0.25)
}#for(b in bees)
axis(2); axis(1, at = 1:2, labels = c('UV high', 'UV low'))
arrows(rep(1,2), A1(exp(q.k[1])), y1 = A1(exp(q.k[2])), length = 0.05, code = 3, angle = 90, col = 'magenta')
arrows(rep(2,2), A1(exp(q.k[1]+ q.b.k[1])), y1 = A1(exp(q.k[2]+ q.b.k[2])), length = 0.05, code = 3, angle = 90, col = 'purple')
lines(1:2, A1(exp(med.kappa+c(0,med.b.kappa))), col = c('purple4'), lwd = 2)
points(1:2, A1(exp(med.kappa+c(0,med.b.kappa))), col = c('magenta', 'purple'), pch = 21, bg = 'white')
legend('bottomleft', inset = 0.15, legend = c('Individual Bees', 'Prediction Interval', 'Population Mean', '95%CI'), pch = c(19,22,21,NA), lty = c(1,NA,NA,1), col = c(rgb(0,0,0,0.2), rgb(1,0,1,0.2), 'purple','purple'), lwd = c(1,3,1,1), cex = 0.5)
suppressWarnings(
dev.copy(pdf, paste0(Sys.getenv('HOME'),'/Dropbox/Plotmania/', 	
				'lowhighUV mean vector','.pdf'),
	width= par("din")[1], height= par("din")[2], useDingbats = F)	);
	dev.off();dev.set(dev.prev())
	
#drawing all values from the model
dev.new(width = 3); par(mai = c(0.9,0.9,0,0))
plot(1:2, A1(exp(med.kappa+c(0,med.b.kappa))), xlim = c(0,3)+c(1,-1)*c(0.7,0.7), ylim = c(0,1), col = c('magenta', 'purple'), lwd = 2, cex = 0.25, axes = F, xlab = 'Condition', ylab = 'Mean Vector Length')
polygon(c(1:2, rev(1:2)), A1(exp( c(c(q.k[1], q.k[1]+q.b.k[1]) -1.96* med.b.kappa.sd, rev(c(q.k[2], q.k[2]+q.b.k[2]))+1.96* med.b.kappa.sd) )), col = rgb(1,0,1,0.1), border = F)
for(b in bees){
	# bb <- subset(MEmises_dat1, bee == levels(MEmises_dat1$bee)[b])
	k1 <- med.kappa + median(fitted$r_2_kappa_1[,b])
	k2 <- med.kappa + med.b.kappa + median(fitted$r_2_kappa_1[,b])
	points(1:2,A1(exp(c(k1,k2))), col = rgb(0,0,0,0.1), pch = 19, cex = 0.5)
	lines(1:2, A1(exp(c(k1,k2))), col = rgb(0,0,0,0.1), lwd = 0.25)
	# lines(1:2, A1(exp(med.kappa+c(0,med.b.kappa))+median(fitted$r_2_kappa_1[,b])), col = rgb(0,0,0,0.1), lwd = 0.25)
}#for(b in bees)
axis(2); axis(1, at = 1:2, labels = c('UV high', 'UV low'))
arrows(rep(1,2), A1(exp(q.k[1])), y1 = A1(exp(q.k[2])), length = 0.05, code = 3, angle = 90, col = 'magenta')
arrows(rep(2,2), A1(exp(q.k[1]+ q.b.k[1])), y1 = A1(exp(q.k[2]+ q.b.k[2])), length = 0.05, code = 3, angle = 90, col = 'purple')
lines(1:2, A1(exp(med.kappa+c(0,med.b.kappa))), col = c('purple4'), lwd = 2)
points(1:2, A1(exp(med.kappa+c(0,med.b.kappa))), col = c('magenta', 'purple'), pch = 21, bg = 'white')
legend('bottomleft', inset = 0.15, legend = c('Individual Bees', 'Prediction Interval', 'Population Mean', '95%CI'), pch = c(19,22,21,NA), lty = c(1,NA,NA,1), col = c(rgb(0,0,0,0.2), rgb(1,0,1,0.2), 'purple','purple'), lwd = c(1,3,1,1), cex = 0.5)
suppressWarnings(
dev.copy(pdf, paste0(Sys.getenv('HOME'),'/Dropbox/Plotmania/', 	
				'lowhighUV mean vector modelled','.pdf'),
	width= par("din")[1], height= par("din")[2], useDingbats = F)	);
	dev.off();dev.set(dev.prev())
#the difference between these two implies that I should try an even more complicated mode :-/