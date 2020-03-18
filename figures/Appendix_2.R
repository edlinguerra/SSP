##Example 3

library(SSP)
library(tidyr)
library(ggplot2)
library(dplyr)

data(pilot)
sectors <- levels(pilot$Sector)

#Defining arguments for simulation and sampling
N = 30
sites = 30
cases = 20
n = 20
m = 20
k = 10

#Lists to store results
sum.l <- opt.l <- qua.l <- vector(mode = "list", length = 3)

#Loop SSP at each sector
for ( i in 1:length(sectors)){
  dat <- pilot[pilot$Sector==sectors[i],2:length(pilot)]

  #parameters for simulation
  par <- assempar(data = dat, type = "cover", Sest.method = "chao")

  # Simulation of data
  sim <- simdata(Par = par, cases = cases, N = N, sites = sites)

  # Quality of simulated data
  qua <- datquality(data = dat, dat.sim = sim, Par = par, transformation = "fourth root", method = "bray")
  qua$sector <- rep(sectors[i], nrow(qua))
  qua.l[[i]] <- qua

  # Sampling and estimation of multse for each data set
  samp <- sampsd(dat.sim = sim, Par = par, transformation = "fourth root",
                 method = "bray", n = n, m = m, k = k)

  # average of multse for each potential sampling design
  sum <- summary_ssp(samp, multi.site = TRUE)

  #Optimal sample sizes
  opt <- ioptimum(sum)
  opt <- as.data.frame(opt)
  opt$sv <- c("sites", "samples")
  opt <- pivot_longer(opt, cols = c("c1", "c2", "c3"), names_to = "cut", values_to = "effort")
  opt$sector <- rep(sectors[i], nrow(opt))
  opt.l[[i]] <- opt

  #arrangement to plot
  sum$sector <- rep(sectors[i], nrow(sum))
  sum.l[[i]] <- sum

  }

#combine summary into a data frame
sum.df <- do.call(rbind.data.frame, sum.l)
sum.df$sector <- factor(sum.df$sector, levels = c("E", "M", "I"))

#combine optimal sample sizes into a data frame
opt.df <- do.call(rbind.data.frame, opt.l)
opt.df$sector <- factor(opt.df$sector, levels =  c("E", "M", "I"))

#Combine quality features into a data frame
qua.df <- do.call(rbind.data.frame, qua.l)

# Generation of plot
my_breaks <- function(x) {
  y <- seq(min(x), max(x), 1)
  y <- round(y,0)
}

#Definition of values for shade areas
shade.opt <- opt.df %>%
             group_by(sector, sv) %>%
             filter(cut != "c1") %>%
             summarise(xmin = min(effort), xmax = max(effort))

shade.sub <- opt.df %>%
             group_by(sector, sv) %>%
             filter(cut != "c3") %>%
             summarise(xmin = min(effort), xmax = max(effort))

#plot
fig.3 <- ggplot(sum.df, aes(x=samples, y=mean))+
         geom_point()+
         geom_errorbar(aes(ymin=lower, ymax=upper), width=.1)+
         facet_grid(sector~sv, scales = "free_x")+
         theme_bw(base_size=16) +
         ylab ("Multivariate pseudo SE")+
         xlab("Sampling effort")+
         scale_y_continuous(breaks=seq(0.0, max(sum.df$upper), 0.025))+
         scale_x_continuous(breaks = my_breaks)+
         theme(axis.text.x = element_text(colour="black", size=rel(0.7)),
         axis.text.y = element_text(colour="black", size=rel(0.7)),
         axis.title.x = element_text(colour="black", size=rel(0.9)),
         axis.title.y = element_text(colour="black", size=rel(0.9)),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.border = element_rect(size=0.4),
         axis.ticks= element_line(size=0.2))+
         geom_rect(data = shade.opt, aes_(x = NULL,y = NULL,
                                   xmin=~xmin, xmax=~xmax, ymin=min(sum.df$lower), ymax=max(sum.df$upper)), alpha=0.5, fill="grey10")+
         geom_rect(data = shade.sub, aes_(x = NULL,y = NULL,
                                   xmin=~xmin, xmax=~xmax, ymin=min(sum.df$lower), ymax=max(sum.df$upper)), alpha=0.5, fill="grey50")+
         geom_text(aes_(x = ~samples, y = ~upper + 0.01, label = ~cum), na.rm = TRUE)
fig.3



