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
			

# bees <- 1:length(levels(MEmises_dat1$bee))
# cll <-  c('magenta', 'purple', 'green1', 'green4')
# for(b in bees){
	# bb <- subset(MEmises_dat1, bee == levels(MEmises_dat1$bee)[b])
	# lvv <- levels(as.factor(paste0(subset(MEmises_dat1, bee == levels(MEmises_dat1$bee)[b])$light_type)))
	# dev.new(width = 7*length(lvv)/2, height = 5)
	# par(mfrow= c(1, length(lvv)), mai = rep(0,4))
	# for(ll in lvv){
		# cl <- ifelse(ll == 'ul', 'purple', ifelse(ll == 'uh', 'magenta', ifelse(ll == 'gl', 'green1', 'green4')))
		# Cplot(subset(bb, light_type == ll)$bearing, col = cl)
		# mtext(ll, line = -4)
		# }#for(ll in lvv)
		# mtext(levels(MEmises_dat1$bee)[b], outer = T, line = -4)
		# # save it
		# suppressWarnings(
		# dev.copy(pdf, paste0(Sys.getenv('HOME'),'/Dropbox/Plotmania/', 	
						# 'lowGreenUV bee',levels(cd$bee)[b],'.pdf'),
			# width= par("din")[1], height= par("din")[2], useDingbats = F)	);
			# dev.off();dev.set(dev.prev())
		# graphics.off()
# }#for(b in bees)

# a lot of animals didn't experience both conditions
#select only those that did
b_both <- unique(subset(MEmises_dat1$bee, MEmises_dat1$light_type == 'ul'))[unique(subset(MEmises_dat1$bee, MEmises_dat1$light_type == 'ul')) %in% unique(subset(MEmises_dat1$bee, MEmises_dat1$light_type == 'uh'))]

MEmises_dat1 <- subset(MEmises_dat1, MEmises_dat1$bee %in% b_both)
#factors should now only have levels found in the dataset
MEmises_dat1$bee <- factor(MEmises_dat1$bee)
MEmises_dat1$light_type <- factor(MEmises_dat1$light_type)

bees <- 1:length(levels(MEmises_dat1$bee))
cll <-  c('magenta', 'purple', 'green1', 'green4')
# for(b in bees){
	# bb <- subset(MEmises_dat1, bee == levels(MEmises_dat1$bee)[b])
	# lvv <- levels(as.factor(paste0(subset(MEmises_dat1, bee == levels(MEmises_dat1$bee)[b])$light_type)))
	# dev.new(width = 7*length(lvv)/2, height = 5)
	# par(mfrow= c(1, length(lvv)), mai = rep(0,4))
	# for(ll in lvv){
		# cl <- ifelse(ll == 'ul', 'purple', ifelse(ll == 'uh', 'magenta', ifelse(ll == 'gl', 'green1', 'green4')))
		# Cplot(subset(bb, light_type == ll)$bearing, col = cl)
		# mtext(ll, line = -4)
		# }#for(ll in lvv)
		# mtext(levels(MEmises_dat1$bee)[b], outer = T, line = -4)
		# # save it
		# suppressWarnings(
		# dev.copy(pdf, paste0(Sys.getenv('HOME'),'/Dropbox/Plotmania/', 	
						# 'lowhighUV bee',levels(cd$bee)[b],'.pdf'),
			# width= par("din")[1], height= par("din")[2], useDingbats = F)	);
			# dev.off();dev.set(dev.prev())
		# graphics.off()
# }#for(b in bees)


#################################################################
#	compile														#
#################################################################
form1 <- bf( angles ~ light_type + (light_type |bee), kappa ~ light_type + (light_type |bee))			
  
br.data <- make_standata(
							formula = form1,
							data = MEmises_dat1,
							family = von_mises
							)
save(br.data, file = paste0(Sys.getenv('HOME'), '/Dropbox/DFG2019-reapplication/', 'br.data.uvMEmises-farbetanzen-UV---.Rdata'))	
load(file = paste0(Sys.getenv('HOME'), '/Dropbox/DFG2019-reapplication/', 'br.data.uvMEmises-farbetanzen-UV---.Rdata'))

br.code <- make_stancode(
							formula = form1,
							data = MEmises_dat1,
							family = von_mises
							)	
#brms can't handle random effects in von_mises distributions (& others)
#save this stancode and rewrite it in the style of "handwritten.uvMEmises6-1.stan"
write.table(br.code, file = paste0(Sys.getenv('HOME'), '/Dropbox/DFG2019-reapplication/', 'br.uvMEmises-farbetanzen-UV---.stan'), col.names = F, row.names = F, quote = F)					

#"handwritten.uvMEmises6-1" is based on brms code for the formula:
# mu ~ 1 + (1|animal), kappa ~ 1 + (1|animal)
# i.e. each animal may have a different mean and precision
# mu is now estimated as a 2-component unit_vector
# and its transform to an angle: 
# temp_angle = atan2(mu_vec[2], mu_vec[1]);
# instead of random effects standard deviation,
# a ±arc containing animal means is estimated
# vector<lower=0, upper=pi()>[M_1] sd_1;
# target += uniform_lpdf(sd_1 | 0, pi())  //no bias towards smaller SD.
# - 1 * uniform_lccdf(0 | 0, pi()); // N.B. sd_1 is already bounded to (0,pi())
# this is scaled for each animal by z_1,
# vector<lower=-1, upper=1>[N_1] z_1[M_1];
# target += uniform_lpdf(z_1[1] | -1, 1) - normal_lpdf(z_1[1]| 0, 0.25);
# N.B. the prior on mu_vec is now set according to pop. estimate of mu
round(mean.circular(br.data$Y),2)
round(mean.circular(br.data$Y[!as.logical(br.data$X[,2])]),2)
# could set kappa this way, but safer to have a very weak prior 0.1
round(A1inv(rho.circular(br.data$Y)),2)
round(A1inv(rho.circular(br.data$Y[!as.logical(br.data$X[,2])])),2)
# target += von_mises_unitvector_lpdf(mu_vec | mean_circular(Y,N), 5);
# let's see if that helps anything
# 

# # a few possible priors for Z
# plot(seq(-1,1,length.out = 10^3)*180, dunif(seq(-1,1,length.out = 10^3), -1,1, log = T)-dnorm(seq(-1,1,length.out = 10^3), 0,0.25, log = T), type = 'l', ylim = c(-3,5), xlab = 'Maximum Offset from Population Mean (°)', ylab = 'Probability Density')
# lines(seq(-1,1,length.out = 10^3)*180, dvonmises(seq(-1,1,length.out = 10^3)*pi, pi,1, log = T), col = 'blue')
# lines(seq(-1,1,length.out = 10^3)*180, dvonmises(seq(-1,1,length.out = 10^3)*2*pi, pi,1, log = T), col = 'red')
# legend('top', inset = 0.1, legend = c('InverseNormal(0,0.25)', 'vonMises(pi,1)', 'Doubled vonMises(2*pi/2,1)'), col = c('black','blue','red'), lty = 1, cex = 0.75)

# print(true.mu)
# print(round(true.offset))
# print(true.kappa)

#can I get this to run 4 chains?
#yes, use brms default von_mises_real_lpdf
fit <- stan(file =
				paste0(Sys.getenv('HOME'), '/Dropbox/DFG2019-reapplication/',
					'br.uvMEmises-farbetanzen-UV-1---.stan'),#I wrote this myself
			data = br.data,#data object brms made
			# iter = 400, warmup = 200,
			iter = 4000, warmup = 2000,
			init = "random", cores = getOption("mc.cores", 1L),
			chains = 4,#+1,
			control = list(adapt_delta = 0.95))

save(fit, file = paste0(file = paste0(Sys.getenv('HOME'), '/Dropbox/DFG2019-reapplication/', 'longrunUV-1---.mod.Rdata')))
load(file = paste0(file = paste0(Sys.getenv('HOME'), '/Dropbox/DFG2019-reapplication/', 'longrunUV-1---.mod.Rdata')))
			
#why are there divergent transitions?
#N.B. Some initialisation problems, using von Mises lpdf on mu_vec
# numeric overflow in cyl_bessel_i<double>(double,double)
# at least weeds out the bad chains, but would be nice to avoid
# actually good for convergence, perhaps run more chains with a strong prior	
#update, can be avoided by using normal_lpdf for kappa>100
print(fit, pars = c('b_vec','b', 'mu_vec', 'b_Intercept', 'sd_1', 'b_kappa_Intercept', 'b_kappa', 'sd_2'))
traceplot(fit, pars = c('b_vec','b', 'mu_vec', 'b_Intercept', 'sd_1', 'b_kappa_Intercept', 'b_kappa', 'sd_2'))

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
dev.new(width = 3, height = 5);par(mai = c(0.5,0.9,0,0))
plot(1:2, A1(exp(med.kappa+c(0,med.b.kappa))), xlim = c(0,3)+c(1,-1)*0.7, ylim = c(0,1), col = c('magenta', 'purple'), lwd = 2, cex = 0.25, axes = F, xlab = '', ylab = 'Mean Vector Length'); axis(1, at = 1:2, labels = c('UV high', 'UV low'));axis(2)

polygon(c(1,2,2,1), A1(exp( c(q.k[1]-1.96* med.b.kappa.sd, q.k[1]+ q.b.k[1]-1.96* med.b.kappa.sd, rev(c(q.k[2]+1.96* med.b.kappa.sd, q.k[2]+ q.b.k[2]+1.96* med.b.kappa.sd))) )), col = rgb(1,0,1,0.1), border = F)

bees <- 1:length(levels(MEmises_dat1$bee))
cll <-  c('magenta', 'purple', 'green1', 'green4')
for(b in bees){
	bb <- subset(MEmises_dat1, bee == levels(MEmises_dat1$bee)[b])
	k1 <- mle.vonmises(with(bb,subset(angles, light_type == 'uh')), bias=T)$kappa
	k2 <- mle.vonmises(with(bb,subset(angles, light_type == 'ul')), bias=T)$kappa
	points(1:2, A1(c(k1,k2)), col = rgb(0,0,0,0.1), pch = 19, cex = 0.5)
	lines(1:2, A1(c(k1,k2)), col = rgb(0,0,0,0.1), lwd = 0.5)
	# lines(1:2, A1(exp(med.kappa+c(0,med.b.kappa))+median(fitted$r_2_kappa_1[,b])), col = rgb(0,0,0,0.1), lwd = 0.25)
}#for(b in bees)
lines(1:2, A1(exp(med.kappa+c(0,med.b.kappa))), col = 'purple4', lwd = 2)
arrows(rep(1,2), A1(exp(q.k[1])), y1 = A1(exp(q.k[2])), length = 0.05, code = 3, angle = 90, col = 'magenta', lwd = 2)
arrows(rep(2,2), A1(exp(q.k[1]+ q.b.k[1])), y1 = A1(exp(q.k[2]+ q.b.k[2])), length = 0.05, code = 3, angle = 90, col = 'purple', lwd = 2)
points(1:2, A1(exp(med.kappa+c(0,med.b.kappa))), col = c('magenta', 'purple'), pch = 21, bg = 'white')
legend('bottomleft', inset = 0.05,
		legend = c('Individual bees',
					'Prediction Interval',
					'Model Estimate',
					'95% CI'),
		pch = c(19, 15, 21, NA),
		lty = c(1,NA,NA,1),
		lwd = c(1,3,1,1),
		col = c(rgb(0,0,0,0.2), rgb(1,0,1,0.2), 'purple', 'purple'),
		cex = 0.65)
suppressWarnings(
dev.copy(pdf, paste0(Sys.getenv('HOME'),'/Dropbox/Plotmania/', 	
				'lowhighUV bee','longrun---model','.pdf'),
	width= par("din")[1], height= par("din")[2], useDingbats = F)	);
	dev.off();dev.set(dev.prev())
		
		
dev.new(width = 3, height = 5)
par(mfrow = c(2,1), mai = c(0,0,0,0))
plot.circular(circular(pi/2), axes = F, col = 'magenta', stack = T, bins = 360)
for(b in bees){
	bb <- with(MEmises_dat1, subset(bearing, bee == levels(MEmises_dat1$bee)[b] & light_type == 'uh'))
	arrows.circular(mean(circular(bb,  units = 'degrees', zero = pi/2, rotation = 'clock')), length = 0, shrink = rho.circular(circular(bb,  units = 'degrees', zero = pi/2, rotation = 'clock')), col = rgb(1,0,1,0.2), lwd = 3, lend = 'butt')
}#for(b in bees)
plot.circular(circular(pi/2), axes = F, col = 'purple', stack = T, bins = 360)
for(b in bees){
	bb <- with(MEmises_dat1, subset(bearing, bee == levels(MEmises_dat1$bee)[b] & light_type == 'ul'))
	arrows.circular(mean(circular(bb,  units = 'degrees', zero = pi/2, rotation = 'clock')), length = 0, shrink = rho.circular(circular(bb,  units = 'degrees', zero = pi/2, rotation = 'clock')), col = rgb(0.5,0,0.5,0.2), lwd = 3, lend = 'butt')
}#for(b in bees)
suppressWarnings(
dev.copy(pdf, paste0(Sys.getenv('HOME'),'/Dropbox/Plotmania/', 	
				'lowhighUV bee','circle mean vec','.pdf'),
	width= par("din")[1], height= par("din")[2], useDingbats = F)	);
	dev.off();dev.set(dev.prev())


