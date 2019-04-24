# Created by MRU on March 31, 2019
# To regress type of CCF relationship of pixels with climate data

library("ncdf4") 
library(zoo)
library(raster)
library(rasterVis)
library(gpclib)
library(maptools)
library(maps)
library(forecast)
library(astsa)
library(TSA)
library(ggplot2)
library(ggpubr)  #for several ggplots together
library("MultNonParam")  #for Kruskal Wallis
library("FSA")  #for the Dunn Test
library(rcompanion) #for organizing with the non-parametric test results
library(multcompView)  #for organizing with the non-parametric test results

# 1. Read climate data (ppt and rad for now) and relationship data --------
####  and calculate annual means for each climate data for each pixel
####  and add the calculated stuff to a data frame for everrrrything
# relationship:
ccfprecrad <- read.table('/Users/uribem/Google Drive/PhD/Time series/ccf/ccf_precrad_everything_sif.csv', sep=",", header=TRUE)
#prec:
evi.resampled <- list.files(path="/Volumes/Maria/maiac_vi/maiacevi_resampledwithrad/", pattern = "*.gri$", full.names = TRUE)
evi.resampled.ras <- stack(evi.resampled)

#evi.arrays <- as.array(evi.resampled.ras[[2:179]]) 
#Create the stack using the resampled rasters:
#trmmgri <- raster("/Volumes/MRUDATA/Time series work/precipitation_trmm3b43/resampling/trmm_resampled/3B43.20000101.7A.HDF_resampled.gri")
trmmres <- list.files(path="/Volumes/depot/trmm_tmpa_download/prec_trmmresampled_withrad", pattern = "*_resampledlpjguesswithrad.gri$", full.names = TRUE)
trmmres_st <- stack(trmmres[3:180]) #don't include the first file which is for Jan 2000 as there's no EVI data
trmmres_st.trops <- crop(trmmres_st, evi.resampled.ras)
#origprec <- raster("/Volumes/Maria/trmm_tmpa_download/downloads/3B43.20000101.7A.HDF.nc4")  # take a look at the original nc file
trmmlayers2 <- split(trmmres, ceiling(seq_along(trmmres)/12))  # to divide the file list (months) by year
trmmlayers3 <- lapply(trmmlayers2, function(x) stack(x))  # stack each group of months (i.e. each year)
trmmlayers5 <- lapply(trmmlayers3, function(x) stack(x[[1]]*24*31, x[[2]]*24*28, x[[3]]*24*31, x[[4]]*24*30,  # calculate total month prec
                                                     x[[5]]*24*31, x[[6]]*24*30, x[[7]]*24*31, x[[8]]*24*31,
                                                     x[[9]]*24*30, x[[10]]*24*31, x[[11]]*24*30, x[[12]]*24*31))
trmm.year.mean <- lapply(trmmlayers5, function(x) calc(x, mean))  # mean monthly prec per year # NOT USED!
trmm.year.total <- lapply(trmmlayers5, function(x) calc(x, sum))  # total year prec (now I have one raster for each year annual prec)
trmm.year.total.stack <- stack(trmm.year.total)  # stack all years
trmm.year.total.mean <- mean(trmm.year.total.stack)  # all years mean (i.e., MAP 2000-2014)
plt <- levelplot(trmm.year.total.mean, margin = F, colorkey=TRUE) # plot
plt + layer(sp.lines(world.outlines.sp, col = "black", lwd = 0.5))
trmm.year.total.mean <- crop(trmm.year.total.mean, evi.resampled.ras)
trmm.year.total.mean <- mask(trmm.year.total.mean, continents)
trmm.year.total.mean <- mask(trmm.year.total.mean, meanevi.ras)  #mask for low EVI values too
map.df <- as.data.frame(trmm.year.total.mean, xy=TRUE)
colnames(map.df) <- c("lon","lat","MAP")

ccf_precrad_map <- merge(ccfprecrad, map.df, by=c("lon","lat"))
ccf_precrad_map$relationship <- as.factor(ccf_precrad_map$relationship)

#rad:
radnc <- "/Volumes/depot/radiation/CERES_EBAF-Surface_Ed4.0_Subset_200003-201501.nc" # clm gpp netcdf file name
radin <- brick(radnc)  # read as a brick of rasters because it has a third dimension time
radin <- raster::rotate(radin[[1:178]])
radin.trops <- crop(radin, evi.resampled.ras)
allmonths <- c(1:178)
allmonths.div <- split(allmonths, ceiling(seq_along(allmonths)/12))  # to divide the file list (months) by year
allmonths.div2 <- list(c(1:10), c(11:22), c(23:34), c(35:46), c(47:58),
                       c(59:70), c(71:82), c(83:94), c(95:106), c(107:118),
                       c(119:130), c(131:142), c(143:154), c(155:166), c(167:178))
radlayers3 <- lapply(allmonths.div2, function(x) subset(radin.trops, x[1]:x[length(x)]))  #subset by group of months of each year
radlayers.years.mean <- lapply(radlayers3, function(x) calc(x, mean))  #mean of each year
rad.year.mean.stack <- stack(radlayers.years.mean)  # stack all years
rad.year.mean.mean <- mean(rad.year.mean.stack)  # all years mean (i.e., MAP 2000-2014)
rad.year.mean.mean <- mask(rad.year.mean.mean, continents)  #mask with continents
rad.year.mean.mean <- mask(rad.year.mean.mean, meanevi.ras)  #mask for low EVI values too
radmean.df <- as.data.frame(rad.year.mean.mean, xy=TRUE)  
colnames(radmean.df) <- c("lon","lat","radmean")
ccf_precrad_map <- merge(ccf_precrad_map, radmean.df, by=c("lon","lat"))

# Dry season length:
dslfile <- "/Volumes/depot/RADS/duration.dry.season.TRMM.1998-2016.nc" # clm gpp netcdf file name
dslin <- brick(dslfile)  # read as a brick of rasters because it has a third dimension time
dslin20002014 <- dslin[[3:16]]
dslin20002014 <- raster::rotate(dslin20002014)
dslin20002014 <- crop(dslin20002014, evi.resampled.ras)
dslin200020141deg <- resample(dslin20002014, evi.resampled.ras)
plt <- levelplot(dslin20002014[[1:6]], margin = F, colorkey=TRUE) # plot
plt + layer(sp.lines(world.outlines.sp, col = "black", lwd = 0.5))
plot(is.na(dslin20002014[[3]]))
# check the nc file:
ncin <- nc_open(dslfile)    #open the file
print(ncin) 
testlon <- ncvar_get(ncin, "lon")  #get each dimension variable
testlat <- ncvar_get(ncin, "lat")
testlev <- ncvar_get(ncin, "lev")
testt <- ncvar_get(ncin, "time")
LonIdx <- which( ncin$dim$lon$vals == 124.875)  #to validate with some random coordinates
LatIdx <- which( ncin$dim$lat$vals == 0.625)
MyVariable <- ncvar_get(ncin, 'durdry')[ LonIdx, LatIdx, ]  #extract for coordinates and validate the results
# Ok, everything looks good, so keep going:
dslin200020141deg[dslin200020141deg > 365] <- 365
dsl.mean <- mean(dslin200020141deg)  # all years mean (i.e., Mean DSL 2000-2014)
dsl.mean <- mask(dsl.mean, continents)  #mask with continents
dsl.mean <- mask(dsl.mean, meanevi.ras)  #mask for low EVI values too
plt <- levelplot(dsl.mean.months, margin = F, colorkey=TRUE) # plot
plt + layer(sp.lines(world.outlines.sp, col = "black", lwd = 0.5))
dsl.mean.months <- dsl.mean/30
dsl.df <- as.data.frame(dsl.mean, xy=TRUE)
colnames(dsl.df) <- c("lon","lat","dsl")
ccf_precrad_map <- merge(ccf_precrad_map, dsl.df, by=c("lon","lat"))

#prec-rad correlation:
precradcorr.month <-  corLocal(trmmres_st.trops,radin.trops)   # correlation between raster stacks of prec and rad
#continents <- readOGR("/Volumes/Maria/LandCover/TM_WORLD_BORDERS_SIMPL-0.3/TM_WORLD_BORDERS_SIMPL-0.3.shp")  
precradcorr.month <- mask(precradcorr.month, continents)  # mask continental area
precradcorr.month <- mask(precradcorr.month, meanevi.ras)  #mask for low EVI values too
precradcorr.df <- as.data.frame(precradcorr.month, xy=TRUE)  # 
colnames(precradcorr.df) <- c("lon","lat","precradcorr")
ccf_precrad_map <- merge(ccf_precrad_map, precradcorr.df, by=c("lon","lat"))

# Temperature:
tempfile19912000 <- "cru_ts4.02.1991.2000.tmp.dat.nc"
tempfile20012010 <- "cru_ts4.02.2001.2010.tmp.dat.nc"
tempfile20112017 <- "cru_ts4.02.2011.2017.tmp.dat.nc"
tempras19912000 <- stack(paste0("/Volumes/depot/temperature/",tempfile19912000))
tempras20012010 <- stack(paste0("/Volumes/depot/temperature/",tempfile20012010))
tempras20112017 <- stack(paste0("/Volumes/depot/temperature/",tempfile20112017))
# check one of the nc files:
ncin <- nc_open(paste0("/Volumes/Maria/temperature/",tempfile19912000))    #open the file
print(ncin) 
testlon <- ncvar_get(ncin, "lon")  #get each dimension variable
testlat <- ncvar_get(ncin, "lat")
testt <- ncvar_get(ncin, "time")
LonIdx <- which( ncin$dim$lon$vals == -54.25)  #to validate with some random coordinates
LatIdx <- which( ncin$dim$lat$vals == 2.75)
MyVariable <- ncvar_get(ncin, 'tmp')[ LonIdx, LatIdx, ] 
# Everything looks ok, so let's keep going...
tempras20002014 <- stack(tempras19912000[[109:120]], tempras20012010, tempras20112017[[1:48]])
tempras20002014.trops <- crop(tempras20002014, evi.resampled.ras)
tempras20002014.trops.1deg <- resample(tempras20002014.trops, evi.resampled.ras)
temp.mean <- mean(tempras20002014.trops.1deg)
temp.mean <- mask(temp.mean, meanevi.ras)  #mask for low EVI values too
plot(temp.mean)
temp.mean.df <- as.data.frame(temp.mean, xy=TRUE)  # 
colnames(temp.mean.df) <- c("lon","lat","meantemp")
ccf_precrad_map <- merge(ccf_precrad_map, temp.mean.df, by=c("lon","lat"))

write.csv(ccf_precrad_map, "/Users/uribem/Google Drive/PhD/Time series/ccf/ccf_precrad_everything-allregs_sif.csv", row.names = FALSE)

# Prec seasonality:
precsi <- raster("/Volumes/Maria/seasonality/prec/prec_si_2000-2014.tif")
precsi.m <- mask(precsi, continents)
precsi.m <- mask(precsi.m, meanevi.ras)  #mask for low EVI values too
precseas.df <- as.data.frame(precsi.m, xy=TRUE)
colnames(precseas.df) <- c("lon","lat","prec_si")
ccf_precrad_map <- merge(ccf_precrad_map, precseas.df, by=c("lon","lat"))

write.csv(ccf_precrad_map, "/Users/uribem/Google Drive/PhD/Time series/ccf/ccf_precrad_everything-allregs_sif2.csv", row.names = FALSE)


# 2. First set of plots: ---------------------------------------------------------------
# Boxplots:
boxplot(ccf_precrad_map$meantemp ~ ccf_precrad_map$relationship, ylab="Temperature", xlab="Type of relationship")
# Violin plots:
myviolpal <- c("#4c0000","#033159","#c1a100","#281d0b","#a3a3a3","#a50098","#310059","#002104","#a34300") #color palette for rel type
templab <- "Mean Temperature (°C)"
mapviolin <- ggplot(ccf_precrad_map, aes(x=relationship, y=MAP, fill=relationship)) +  # MAP
  geom_violin(trim=TRUE, lwd=0.2, color="black") + scale_x_discrete(limits=c(1:9)) +  
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.03) + 
  geom_boxplot(width=0.05, outlier.size=0.05, fill="#ffffff", lwd=0.3, color="black") + 
  labs(y=expression(paste("MAP (mm ", yr^-1, ")", sep = ""))) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"), 
        axis.title.x=element_blank(), axis.text.x=element_blank(), legend.position="none")
#mapviolin + scale_fill_manual(values=mypal)  
meanradviolin <- ggplot(ccf_precrad_map, aes(x=relationship, y=radmean, fill=relationship)) +   # MEAN RAD
  geom_violin(trim=TRUE, lwd=0.2, color="black") + scale_x_discrete(limits=c(1:9)) +  
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.03) + 
  geom_boxplot(width=0.05, outlier.size=0.05, fill="#ffffff", lwd=0.3, color="black") + 
  labs(y=expression(paste("Mean Daily Radiation (W ", m^-2, ")", sep = ""))) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"), 
        axis.title.x=element_blank(), axis.text.x=element_blank(), legend.position="none") 
#meanradviolin + scale_fill_manual(values=mypal)  
meandslviolin <- ggplot(ccf_precrad_map, aes(x=relationship, y=dsl, fill=relationship)) +   # MEAN RAD
  geom_violin(trim=TRUE, lwd=0.2, color="black") + scale_x_discrete(limits=c(1:9)) +  
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.03) + 
  geom_boxplot(width=0.05, outlier.size=0.05, fill="#ffffff", lwd=0.3, color="black") + 
  labs(y="Dry season length (days)") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"), 
        axis.title.x=element_blank(), axis.text.x=element_blank(), legend.position="none")  #axis.text.x=element_blank(), 
#meanradviolin + scale_fill_manual(values=mypal)  
precradviolin <- ggplot(ccf_precrad_map, aes(x=relationship, y=precradcorr, fill=relationship)) +  # PREC-RAD CORR
  geom_violin(trim=TRUE, lwd=0.2, color="black") + scale_x_discrete(limits=c(1:9)) +  
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.03) + 
  geom_boxplot(width=0.05, outlier.size=0.05, fill="#ffffff", lwd=0.3, color="black") + 
  theme(panel.background = element_rect(fill = "white", colour = "grey50"), 
        axis.title.x=element_blank(), axis.text.x=element_blank(), legend.position="none") 
#precradviolin + scale_fill_manual(values=mypal)  
tempviolin <- ggplot(ccf_precrad_map, aes(x=relationship, y=meantemp, fill=relationship)) +  # MEAN TEMP
  geom_violin(trim=TRUE, lwd=0.2, color="black") + scale_x_discrete(limits=c(1:9)) +  
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.03) + 
  geom_boxplot(width=0.05, outlier.size=0.05, fill="#ffffff", lwd=0.3, color="black") + 
  labs(y=templab) +  #x="Type of relationship", 
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.position="none",  axis.title.x=element_blank()) 
#tempviolin + scale_fill_manual(values=mypal)  
ccf_precrad_map$latabs <- abs(ccf_precrad_map$lat)
latviolin <- ggplot(ccf_precrad_map, aes(x=relationship, y=latabs, fill=relationship)) +  # MEAN TEMP
  geom_violin(trim=TRUE, lwd=0.2, color="black") + scale_x_discrete(limits=c(1:9)) +  
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.03) + 
  geom_boxplot(width=0.05, outlier.size=0.05, fill="#ffffff", lwd=0.3, color="black") + 
  labs(y=templab) +  #x="Type of relationship", 
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.position="none",  axis.title.x=element_blank()) 
latviolin + scale_fill_manual(values=myviolpal)  
precsiviolin <- ggplot(ccf_precrad_map, aes(x=relationship, y=prec_si, fill=relationship)) +  # PREC SEAS
  geom_violin(trim=TRUE, lwd=0.2, color="black") + scale_x_discrete(limits=c(1:9)) +  
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.03) + 
  geom_boxplot(width=0.05, outlier.size=0.05, fill="#ffffff", lwd=0.3, color="black") + 
  labs(y=templab) +  #x="Type of relationship", 
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.position="none",  axis.title.x=element_blank()) 
precsiviolin + scale_fill_manual(values=myviolpal)  

ggarrange(mapviolin + scale_fill_manual(values=mypal),
          meanradviolin+ scale_fill_manual(values=mypal),
          meandslviolin + scale_fill_manual(values=mypal),
          tempviolin+ scale_fill_manual(values=mypal),
          labels = c("", "", ""),
          ncol = 1, nrow = 4, align="v")
mapviolin + scale_fill_manual(values=mypal)  
meanradviolin + scale_fill_manual(values=mypal)  
tempviolin + scale_fill_manual(values=mypal)  

# other violin stuff:
mapviolin <- ggplot(subset(ccf_precrad_map, ccf_precrad_map$MAP > 500), aes(x=relationship, y=MAP)) + 
  geom_violin() + scale_x_discrete(limits=c(1:9)) +  
  geom_jitter(shape=16, position=position_jitter(0.2))  #points in different locations

library("ggpubr")
ggline(ccf_precrad_map, x = "relationship", y = "MAP", 
       add = c("mean_se", "jitter"))
ccf_precrad_map2 <- subset(ccf_precrad_map, ccf_precrad_map$MAP > 1000)
ggline(ccf_precrad_map2, x = "relationship", y = "MAP", 
       add = c("mean_se", "jitter"))

# 3. ANOVAS ---------------------------------------------------------------
## Ojo!!! I can't use this cause they assume normality of the distribution for the variable of interest 
## and assume same sample sizes for the different groups.
#Compute the analysis of variance #Summary of the analysis #TukeyHSD for comparisons between pairs of groups
map.aov <- aov(MAP ~ relationship, data = ccf_precrad_map) ; summary(map.aov) ; TukeyHSD(map.aov)
rad.aov <- aov(radmean ~ relationship, data = ccf_precrad_map) ; summary(rad.aov) ; TukeyHSD(rad.aov)
tmp.aov <- aov(meantemp ~ relationship, data = ccf_precrad_map) ; summary(tmp.aov) ; TukeyHSD(tmp.aov)
dsl.aov <- aov(dsl ~ relationship, data = ccf_precrad_map) ; summary(dsl.aov) ; TukeyHSD(dsl.aov)
lat.aov <- aov(lat ~ relationship, data = ccf_precrad_map) ; summary(lat.aov) ; TukeyHSD(lat.aov)
prec_si.aov <- aov(prec_si ~ relationship, data = ccf_precrad_map) ; summary(prec_si.aov) ; TukeyHSD(prec_si.aov)
# Non-parametric with different sample sizes (but post hoc for ANOVA):
pt <- pairwise.wilcox.test(x=ccf_precrad_map$MAP, g=ccf_precrad_map$relationship, p.adjust.method = "bonferroni")
pt <- pt$p.value
pt1 <- fullPTable(pt)
multcompLetters(pt1, compare="<", threshold=0.05, Letters=letters, reversed = FALSE)
relationship_count <- table(ccf_precrad_map$relationship)  #number of pixels per relationship type
barplot(relationship_count)

# 4. Kruskal Wallis (used instead of ANOVA for non-normal a --------
ccf_precrad_map_three <- ccf_precrad_map %>% filter(ccf_precrad_map$relationship == c(1,2,3))
KW.map <- kruskal.test(MAP ~ relationship, data = ccf_precrad_map)   #formula=response~group
DT.map <- dunnTest(MAP ~ relationship,data=ccf_precrad_map,    # Dunn test, see ?dunnTest for options
              method="bonferroni", kw=TRUE) #Adjusts p-values for multiple comparisons with Bonferroni
PT.map <- DT.map$res    # Letter display
DT.map2 <- cldList(P.adj ~ Comparison, data = PT.map, threshold = 0.05)
DT.map2[,2]  # contains the values for all groups
print(DT.map,dunn.test.results=TRUE)
# for mean radiation:
KW.rad <- kruskal.test(radmean ~ relationship, data = ccf_precrad_map) 
DT.rad <- dunnTest(radmean ~ relationship,data=ccf_precrad_map, method="bonferroni", kw=TRUE) 
PT.rad <- DT.rad$res    # Letter display
DT.rad2 <- cldList(P.adj ~ Comparison, data = PT.rad, threshold = 0.05)
print(DT.rad,dunn.test.results=TRUE)
# for mean temperature:
KW.meantemp <- kruskal.test(meantemp ~ relationship, data = ccf_precrad_map)
DT.meantemp <- dunnTest(meantemp ~ relationship,data=ccf_precrad_map, method="bonferroni", kw=TRUE) 
PT.meantemp <- DT.meantemp$res    # Letter display
DT.meantemp2 <- cldList(P.adj ~ Comparison, data = PT.meantemp, threshold = 0.05)
print(DT.rad,dunn.test.results=TRUE)
# for dry season length:
KW.dsl <- kruskal.test(dsl ~ relationship, data = ccf_precrad_map)
DT.dsl <- dunnTest(dsl ~ relationship,data=ccf_precrad_map, method="bonferroni", kw=TRUE) 
PT.dsl <- DT.dsl$res    # Letter display
DT.dsl2 <- cldList(P.adj ~ Comparison, data = PT.dsl, threshold = 0.05)
print(DT.dsl,dunn.test.results=TRUE)
# for latitude:
KW.lat <- kruskal.test(latabs ~ relationship, data = ccf_precrad_map)
DT.lat <- dunnTest(latabs ~ relationship,data=ccf_precrad_map, method="bonferroni", kw=TRUE) 
PT.lat <- DT.lat$res    # Letter display
DT.lat2 <- cldList(P.adj ~ Comparison, data = PT.lat, threshold = 0.05)
print(DT.lat,dunn.test.results=TRUE)
# for precipitation seasonality:
KW.prec_si <- kruskal.test(prec_si ~ relationship, data = ccf_precrad_map)
DT.prec_si <- dunnTest(prec_si ~ relationship,data=ccf_precrad_map, method="bonferroni", kw=TRUE) 
PT.prec_si <- DT.prec_si$res    # Letter display
DT.prec_si2 <- cldList(P.adj ~ Comparison, data = PT.prec_si, threshold = 0.05)
print(DT.prec_si,dunn.test.results=TRUE)
# for prec-rad correlation:
KW.precradcorr <- kruskal.test(precradcorr ~ relationship, data = ccf_precrad_map)
DT.precradcorr <- dunnTest(precradcorr ~ relationship,data=ccf_precrad_map, method="bonferroni", kw=TRUE) 
PT.precradcorr <- DT.precradcorr$res    # Letter display
DT.precradcorr2 <- cldList(P.adj ~ Comparison, data = PT.precradcorr, threshold = 0.05)
print(DT.precradcorr,dunn.test.results=TRUE)

# 5. New plots W/ SIGNIFICANCE --------------------------------------------
# Violin plots with significance:
maprel.wilcox <- compare_means(MAP ~ relationship, data=ccf_precrad_map, 
                               method = "kruskal.test", p.adjust.method="bonferroni")
meanradrel.wilcox <- compare_means(radmean ~ relationship, data=ccf_precrad_map, method = "kruskal.test")
tmprel.wilcox <- compare_means(meantemp ~ relationship, data=ccf_precrad_map, method = "kruskal.test")
dslrel.wilcox <- compare_means(dsl ~ relationship, data=ccf_precrad_map, method = "kruskal.test")
mycomparisons <- list(c("1","3"),c("1","5"),c("2","3"))
mapviolin <- ggplot(ccf_precrad_map, aes(x=relationship, y=MAP, fill=relationship)) +  # MAP
  geom_violin(trim=TRUE, lwd=0.2, color="black") + scale_x_discrete(limits=c(1:9)) +  
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.03) + 
  geom_boxplot(width=0.05, outlier.size=0.05, fill="#ffffff", lwd=0.3, color="black") + 
  labs(y=expression(paste("MAP (mm ", yr^-1, ")", sep = ""))) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"), 
        axis.title.x=element_blank(), axis.text.x=element_blank(), legend.position="none") +
  stat_compare_means(comparisons = mycomparisons, label.y = c(6250, 6760), label = "p.signif") # Add pairwise comparisons p-value
mapviolin <- ggplot(ccf_precrad_map, aes(x=relationship, y=MAP, fill=relationship)) +  # MAP
  geom_violin(trim=TRUE, lwd=0.2, color="black") + scale_x_discrete(limits=c(1:9)) +  
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.03) + 
  geom_boxplot(width=0.05, outlier.size=0.05, fill="#ffffff", lwd=0.3, color="black") + 
  ylim(0,8500) +
  labs(y=expression(paste("MAP (mm ", yr^-1, ")", sep = ""))) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"), 
        axis.title.x=element_blank(), legend.position="none") +  #,axis.text.x=element_blank(),
  annotate(geom="text", size=2.8, x=c(0.9,1.1,0.9,1.1,1), y=c(7700,7700,7000,7000,6300), label=c("3","5","6","7","9")) +
  annotate(geom="text", size=2.8,x=c(1.9,2.1,1.9,2.1,2), y=c(5600,5600,4900,4900,4200), label=c("3","5","6","7","9")) +
  annotate(geom="text", size=2.8,x=c(2.9,3.1,2.9,3.1,2.9,3.1,3), y=c(7900,7900,7200,7200,6500,6500,5800), label=c("1","2","4","5","6","7","8")) +
  annotate(geom="text", size=2.8,x=c(3.9,4.1,4), y=c(4300,4300,3700), label=c("3","5","9")) +
  annotate(geom="text", size=2.8,x=c(4.9,5.1,4.9,5.1,4.9,5.1,5), y=c(6650,6650,5950,5950,5250,5250,4450), label=c("1","2","3","4","6","7","8")) +
  annotate(geom="text", size=2.8,x=c(5.9,6.1,5.9,6.1), y=c(8200,8200,7500,7500), label=c("1","2","3","9")) +
  annotate(geom="text", size=2.8,x=c(6.9,7.1,6.9,7.1,6.9,7.1), y=c(8800,8800,8100,8100,7400,7400), label=c("1","2","3","5","8","9")) +
  annotate(geom="text", size=2.8,x=c(7.9,8.1,7.9,8.1), y=c(5600,5600,4900,4900), label=c("3","5","7","9")) +
  annotate(geom="text", size=2.8,x=c(8.9,9.1,8.9,9.1,8.9,9.1), y=c(6500,6500,5800,5800,5200,5200), label=c("1","2","4","6","7","8"))
mapviolin <- ggplot(ccf_precrad_map, aes(x=relationship, y=MAP, fill=relationship)) +  # MAP
  geom_violin(trim=TRUE, lwd=0.2, color="black") + scale_x_discrete(limits=c(1:9)) +  
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.03) + 
  geom_boxplot(width=0.05, outlier.size=0.05, fill="#ffffff", lwd=0.3, color="black") + 
  #ylim(0,8100) +
  labs(y=expression(paste("MAP (mm ", yr^-1, ")", sep = ""))) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"), 
        axis.title.x=element_blank(), legend.position="none") +  #,axis.text.x=element_blank(),
  annotate(geom="text", size=2.7, x=c(0.9,1.1,0.9,1.1,1), y=c(6700,6700,6400,6400,6100), label=c("3","5","6","7","9")) +
  annotate(geom="text", size=2.7, x=c(1.9,2.1,1.9,2.1,2), y=c(4500,4500,4200,4200,3900), label=c("3","5","6","7","9")) +
  annotate(geom="text", size=2.7, x=c(2.9,3.1,2.9,3.1,2.9,3.1,3), y=c(6400,6400,6100,6100,5800,5800,5500), label=c("1","2","4","5","6","7","8")) +
  annotate(geom="text", size=2.7, x=c(3.9,4.1,4), y=c(3700,3700,3400), label=c("3","5","9")) +
  annotate(geom="text", size=2.7, x=c(4.9,5.1,4.9,5.1,4.9,5.1,5), y=c(5050,5050,4750,4750,4450,4450,4150), label=c("1","2","3","4","6","7","8")) +
  annotate(geom="text", size=2.7, x=c(5.9,6.1,5.9,6.1), y=c(7500,7500,7200,7200), label=c("1","2","3","9")) +
  annotate(geom="text", size=2.7, x=c(6.9,7.1,6.9,7.1,6.9,7.1), y=c(7700,7700,7400,7400,7100,7100), label=c("1","2","3","5","8","9")) +
  annotate(geom="text", size=2.7, x=c(7.9,8.1,7.9,8.1), y=c(4900,4900,4600,4600), label=c("3","5","7","9")) +
  annotate(geom="text", size=2.7, x=c(8.9,9.1,8.9,9.1,8.9,9.1), y=c(5500,5500,5200,5200,4900,4900), label=c("1","2","4","6","7","8"))
meanradviolin <- ggplot(ccf_precrad_map, aes(x=relationship, y=radmean, fill=relationship)) +   # MEAN RAD
  geom_violin(trim=TRUE, lwd=0.2, color="black") + scale_x_discrete(limits=c(1:9)) +  
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.03) + 
  geom_boxplot(width=0.05, outlier.size=0.05, fill="#ffffff", lwd=0.3, color="black") + 
  labs(y=expression(paste("Mean Daily Radiation (W ", m^-2, ")", sep = ""))) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"), 
        axis.title.x=element_blank(), legend.position="none") +
  #ylim(150,340) +
  annotate(geom="text", size=2.7, x=c(0.9,1.1,0.9,1.1,1), y=c(299,299,293,293,287), label=c("3","5","6","7","9")) +
  annotate(geom="text", size=2.7, x=c(1.9,2.1,1.9,2.1,2), y=c(296,296,290,290,284), label=c("3","5","6","7","9")) +
  annotate(geom="text", size=2.7, x=c(2.9,3.1,2.9,3.1,2.9,3.1,3), y=c(303,303,297,297,291,291,285), label=c("1","2","4","5","6","7","8")) +
  annotate(geom="text", size=2.7, x=c(3.9,4.1,3.9,4.1), y=c(289,289,283,283), label=c("3","5","7","9")) +
  annotate(geom="text", size=2.7, x=c(4.9,5.1,4.9,5.1,4.9,5.1), y=c(291,291,285,285,279,279), label=c("1","2","3","4","6","8")) +
  annotate(geom="text", size=2.7, x=c(5.9,6.1,5.9,6.1,6), y=c(294,294,288,288,282), label=c("1","2","3","5","7")) +
  annotate(geom="text", size=2.7, x=c(6.9,7.1,6.9,7.1,6.9,7.1), y=c(295,295,289,289,283,283), label=c("1","2","3","4","6","8")) +
  annotate(geom="text", size=2.7, x=c(7.9,8.1,7.9,8.1), y=c(281,281,275,275), label=c("3","5","7","9")) +
  annotate(geom="text", size=2.7, x=c(8.9,9.1,8.9,9.1), y=c(292,292,286,286), label=c("1","2","4","8"))
tempviolin <- ggplot(ccf_precrad_map, aes(x=relationship, y=meantemp, fill=relationship)) +  # MEAN TEMP
  geom_violin(trim=TRUE, lwd=0.2, color="black") + scale_x_discrete(limits=c(1:9)) +  
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.03) + 
  geom_boxplot(width=0.05, outlier.size=0.05, fill="#ffffff", lwd=0.3, color="black") + 
  labs(y=templab) +  #x="Type of relationship", 
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.position="none",  axis.title.x=element_blank()) +
  annotate(geom="text", size=2.7, x=c(0.9,1.1), y=c(31.3), label=c("2","3")) +
  annotate(geom="text", size=2.7, x=c(1.9,2.1), y=c(30.9), label=c("1","8")) +
  annotate(geom="text", size=2.7, x=c(2.9,3.1,2.9,3.1), y=c(31.6,31.6,30.6,30.6), label=c("1","4","6","8")) +
  annotate(geom="text", size=2.7, x=c(4), y=c(29.7), label=c("3")) +
  annotate(geom="text", size=2.7, x=c(6), y=c(30.7), label=c("3")) +
  annotate(geom="text", size=2.7, x=c(7.9,8.1), y=c(30.8), label=c("2","3"))
meandslviolin <- ggplot(ccf_precrad_map, aes(x=relationship, y=dsl, fill=relationship)) +   # MEAN RAD
  geom_violin(trim=TRUE, lwd=0.2, color="black") + scale_x_discrete(limits=c(1:9)) +  
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.03) + 
  geom_boxplot(width=0.05, outlier.size=0.05, fill="#ffffff", lwd=0.3, color="black") + 
  labs(y="Dry season length (days)") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"), 
        axis.title.x=element_blank(), legend.position="none") + #axis.text.x=element_blank(), 
  annotate(geom="text", size=2.7, x=c(0.9,1.1,0.9,1.1,0.9,1.1,1), y=c(344,344,336,336,328,328,320), label=c("2","3","5","6","7","8","9")) +
  annotate(geom="text", size=2.7, x=c(1.9,2.1,1.9,2.1), y=c(323,323,315,315), label=c("1","3","6","9")) +
  annotate(geom="text", size=2.7, x=c(2.9,3.1,2.9,3.1,2.9,3.1), y=c(327,327,319,319,311,311), label=c("1","2","4","5","6","8")) +
  annotate(geom="text", size=2.7, x=c(3.9,4.1), y=c(291), label=c("3","9")) +
  annotate(geom="text", size=2.7, x=c(4.9,5.1,5), y=c(283,283,275), label=c("1","3","9")) +
  annotate(geom="text", size=2.7, x=c(5.9,6.1,5.9,6.1,5.9,6.1), y=c(317,317,309,309,301,301), label=c("1","2","3","7","8","9")) +
  annotate(geom="text", size=2.7, x=c(6.9,7.1), y=c(306), label=c("1","6")) +
  annotate(geom="text", size=2.7, x=c(7.9,8.1,7.9,8.1), y=c(310,310,302,302), label=c("1","3","6","9")) +
  annotate(geom="text", size=2.7, x=c(8.9,9.1,8.9,9.1,8.9,9.1), y=c(253,253,245,245,237,237), label=c("1","2","4","5","6","8"))
ggarrange(mapviolin + scale_fill_manual(values=myviolpal),
          meanradviolin+ scale_fill_manual(values=myviolpal),
          meandslviolin + scale_fill_manual(values=myviolpal),
          tempviolin+ scale_fill_manual(values=myviolpal),
          labels = c("", "", ""),
          ncol = 2, nrow = 2, align="v")


# #5.1. with a different coloring idea: -----------------------------------
mapviolin <- ggplot(ccf_precrad_map, aes(x=relationship, y=MAP, fill=relationship, color=relationship)) +  # MAP
  geom_violin(trim=TRUE, lwd=0.2, color="black", alpha = .7) + scale_x_discrete(limits=c(1:9)) +  
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.03) + 
  geom_boxplot(width=0.05, outlier.size=0.05, fill="#ffffff", lwd=0.3, color="black") + 
  #ylim(0,8100) +
  labs(y=expression(paste("MAP (mm ", yr^-1, ")", sep = ""))) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"), 
        axis.title.x=element_blank(), legend.position="none") +  #,axis.text.x=element_blank(),
  annotate(geom="text", size=2.7, x=c(0.9,1.1,0.9,1.1,1), y=c(6700,6700,6400,6400,6100), label=c("3","5","6","7","9")) +
  annotate(geom="text", size=2.7, x=c(1.9,2.1,1.9,2.1,2), y=c(4500,4500,4200,4200,3900), label=c("3","5","6","7","9")) +
  annotate(geom="text", size=2.7, x=c(2.9,3.1,2.9,3.1,2.9,3.1,3), y=c(6400,6400,6100,6100,5800,5800,5500), label=c("1","2","4","5","6","7","8")) +
  annotate(geom="text", size=2.7, x=c(3.9,4.1,4), y=c(3700,3700,3400), label=c("3","5","9")) +
  annotate(geom="text", size=2.7, x=c(4.9,5.1,4.9,5.1,4.9,5.1,5), y=c(5050,5050,4750,4750,4450,4450,4150), label=c("1","2","3","4","6","7","8")) +
  annotate(geom="text", size=2.7, x=c(5.9,6.1,5.9,6.1), y=c(7500,7500,7200,7200), label=c("1","2","3","9")) +
  annotate(geom="text", size=2.7, x=c(6.9,7.1,6.9,7.1,6.9,7.1), y=c(7700,7700,7400,7400,7100,7100), label=c("1","2","3","5","8","9")) +
  annotate(geom="text", size=2.7, x=c(7.9,8.1,7.9,8.1), y=c(4900,4900,4600,4600), label=c("3","5","7","9")) +
  annotate(geom="text", size=2.7, x=c(8.9,9.1,8.9,9.1,8.9,9.1), y=c(5500,5500,5200,5200,4900,4900), label=c("1","2","4","6","7","8"))
meanradviolin <- ggplot(ccf_precrad_map, aes(x=relationship, y=radmean, fill=relationship, color=relationship)) +   # MEAN RAD
  geom_violin(trim=TRUE, lwd=0.2, color="black", alpha = .7) + scale_x_discrete(limits=c(1:9)) +  
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.03) + 
  geom_boxplot(width=0.05, outlier.size=0.05, fill="#ffffff", lwd=0.3, color="black") + 
  labs(y=expression(paste("Mean Daily Radiation (W ", m^-2, ")", sep = ""))) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"), 
        axis.title.x=element_blank(), legend.position="none") +
  #ylim(150,340) +
  annotate(geom="text", size=2.7, x=c(0.9,1.1,0.9,1.1,1), y=c(299,299,293,293,287), label=c("3","5","6","7","9")) +
  annotate(geom="text", size=2.7, x=c(1.9,2.1,1.9,2.1,2), y=c(296,296,290,290,284), label=c("3","5","6","7","9")) +
  annotate(geom="text", size=2.7, x=c(2.9,3.1,2.9,3.1,2.9,3.1,3), y=c(303,303,297,297,291,291,285), label=c("1","2","4","5","6","7","8")) +
  annotate(geom="text", size=2.7, x=c(3.9,4.1,3.9,4.1), y=c(289,289,283,283), label=c("3","5","7","9")) +
  annotate(geom="text", size=2.7, x=c(4.9,5.1,4.9,5.1,4.9,5.1), y=c(291,291,285,285,279,279), label=c("1","2","3","4","6","8")) +
  annotate(geom="text", size=2.7, x=c(5.9,6.1,5.9,6.1,6), y=c(294,294,288,288,282), label=c("1","2","3","5","7")) +
  annotate(geom="text", size=2.7, x=c(6.9,7.1,6.9,7.1,6.9,7.1), y=c(295,295,289,289,283,283), label=c("1","2","3","4","6","8")) +
  annotate(geom="text", size=2.7, x=c(7.9,8.1,7.9,8.1), y=c(281,281,275,275), label=c("3","5","7","9")) +
  annotate(geom="text", size=2.7, x=c(8.9,9.1,8.9,9.1), y=c(292,292,286,286), label=c("1","2","4","8"))
tempviolin <- ggplot(ccf_precrad_map, aes(x=relationship, y=meantemp, fill=relationship, color=relationship)) +  # MEAN TEMP
  geom_violin(trim=TRUE, lwd=0.2, color="black", alpha = .7) + scale_x_discrete(limits=c(1:9)) +  
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.03) + 
  geom_boxplot(width=0.05, outlier.size=0.05, fill="#ffffff", lwd=0.3, color="black") + 
  labs(y=templab) +  #x="Type of relationship", 
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.position="none",  axis.title.x=element_blank()) +
  annotate(geom="text", size=2.7, x=c(0.9,1.1), y=c(31.3), label=c("2","3")) +
  annotate(geom="text", size=2.7, x=c(1.9,2.1), y=c(30.9), label=c("1","8")) +
  annotate(geom="text", size=2.7, x=c(2.9,3.1,2.9,3.1), y=c(31.6,31.6,30.6,30.6), label=c("1","4","6","8")) +
  annotate(geom="text", size=2.7, x=c(4), y=c(29.7), label=c("3")) +
  annotate(geom="text", size=2.7, x=c(6), y=c(30.7), label=c("3")) +
  annotate(geom="text", size=2.7, x=c(7.9,8.1), y=c(30.8), label=c("2","3"))
meandslviolin <- ggplot(ccf_precrad_map, aes(x=relationship, y=dsl, fill=relationship, color=relationship)) +   # MEAN RAD
  geom_violin(trim=TRUE, lwd=0.2, color="black", alpha = .7) + scale_x_discrete(limits=c(1:9)) +  
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.03) + 
  geom_boxplot(width=0.05, outlier.size=0.05, fill="#ffffff", lwd=0.3, color="black") + 
  labs(y="Dry season length (days)") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"), 
        axis.title.x=element_blank(), legend.position="none") + #axis.text.x=element_blank(), 
  annotate(geom="text", size=2.7, x=c(0.9,1.1,0.9,1.1,0.9,1.1,1), y=c(344,344,336,336,328,328,320), label=c("2","3","5","6","7","8","9")) +
  annotate(geom="text", size=2.7, x=c(1.9,2.1,1.9,2.1), y=c(323,323,315,315), label=c("1","3","6","9")) +
  annotate(geom="text", size=2.7, x=c(2.9,3.1,2.9,3.1,2.9,3.1), y=c(327,327,319,319,311,311), label=c("1","2","4","5","6","8")) +
  annotate(geom="text", size=2.7, x=c(3.9,4.1), y=c(291), label=c("3","9")) +
  annotate(geom="text", size=2.7, x=c(4.9,5.1,5), y=c(283,283,275), label=c("1","3","9")) +
  annotate(geom="text", size=2.7, x=c(5.9,6.1,5.9,6.1,5.9,6.1), y=c(317,317,309,309,301,301), label=c("1","2","3","7","8","9")) +
  annotate(geom="text", size=2.7, x=c(6.9,7.1), y=c(306), label=c("1","6")) +
  annotate(geom="text", size=2.7, x=c(7.9,8.1,7.9,8.1), y=c(310,310,302,302), label=c("1","3","6","9")) +
  annotate(geom="text", size=2.7, x=c(8.9,9.1,8.9,9.1,8.9,9.1), y=c(253,253,245,245,237,237), label=c("1","2","4","5","6","8"))
ggarrange(mapviolin + scale_fill_manual(values=myviolpal) + scale_color_manual(values=myviolpal),
          meanradviolin+ scale_fill_manual(values=myviolpal) + scale_color_manual(values=myviolpal),
          meandslviolin + scale_fill_manual(values=myviolpal) + scale_color_manual(values=myviolpal),
          tempviolin+ scale_fill_manual(values=myviolpal) + scale_color_manual(values=myviolpal),
          labels = c("", "", ""),
          ncol = 2, nrow = 2, align="v")


# #5.2. with the new significance data: -----------------------------------
samplesizes <- table(ccf_precrad_map$relationship)
samplesizes.lab <- paste0()
mapviolin <- ggplot(ccf_precrad_map, aes(x=relationship, y=MAP, fill=relationship, color=relationship)) +  # MAP
  geom_violin(trim=TRUE, lwd=0.2, color="black", alpha = .8) + scale_x_discrete(limits=c(1:9)) +  
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.03) + 
  geom_boxplot(width=0.1, outlier.size=0.05, fill="#ffffff", lwd=0.3, color="black") + 
  #ylim(0,8100) +
  labs(y=expression(paste("MAP (mm ", yr^-1, ")", sep = ""))) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"), 
        axis.title.x=element_blank(), axis.title.y=element_text(size=8), axis.text=element_text(size=6), legend.position="none") +  #,axis.text.x=element_blank(),
  annotate(geom="text", size=2.5, x=c(1:9), y=c(7270), label=as.character(DT.map2[,2])) +
  annotate(geom="text", size=2.5, x=c(1:9), y=c(7900), label=c(as.vector(samplesizes)))
meanradviolin <- ggplot(ccf_precrad_map, aes(x=relationship, y=radmean, fill=relationship, color=relationship)) +   # MEAN RAD
  geom_violin(trim=TRUE, lwd=0.2, color="black", alpha = .8) + scale_x_discrete(limits=c(1:9)) +  
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.03) + 
  geom_boxplot(width=0.1, outlier.size=0.05, fill="#ffffff", lwd=0.3, color="black") + 
  labs(y=expression(paste("Mean Daily Radiation (W ", m^-2, ")", sep = ""))) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"), 
        axis.title.x=element_blank(), axis.title.y=element_text(size=8), axis.text=element_text(size=6), legend.position="none") +
  #ylim(150,340) +
  annotate(geom="text", size=2.5, x=c(1:9), y=c(289), label=as.character(DT.rad2[,2]))
tempviolin <- ggplot(ccf_precrad_map, aes(x=relationship, y=meantemp, fill=relationship, color=relationship)) +  # MEAN TEMP
  geom_violin(trim=TRUE, lwd=0.2, color="black", alpha = .8) + scale_x_discrete(limits=c(1:9)) +  
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.03) + 
  geom_boxplot(width=0.1, outlier.size=0.05, fill="#ffffff", lwd=0.3, color="black") + 
  labs(y=templab) +  #x="Type of relationship", 
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.position="none",  axis.title.x=element_blank(), axis.title.y=element_text(size=8), axis.text=element_text(size=6)) +
  annotate(geom="text", size=2.5, x=c(1:9), y=c(31.8), label=as.character(DT.meantemp2[,2]))
meandslviolin <- ggplot(ccf_precrad_map, aes(x=relationship, y=dsl, fill=relationship, color=relationship)) +   # MEAN RAD
  geom_violin(trim=TRUE, lwd=0.2, color="black", alpha = .8) + scale_x_discrete(limits=c(1:9)) +  
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.03) + 
  geom_boxplot(width=0.1, outlier.size=0.05, fill="#ffffff", lwd=0.3, color="black") + 
  labs(y="Dry season length (days)") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"), 
        axis.title.x=element_blank(), axis.title.y=element_text(size=8), axis.text=element_text(size=6), legend.position="none") + #axis.text.x=element_blank(), 
  annotate(geom="text", size=2.5, x=c(1:9), y=c(322), label=as.character(DT.dsl2[,2]))
latabsviolin <- ggplot(ccf_precrad_map, aes(x=relationship, y=latabs, fill=relationship, color=relationship)) +   # MEAN RAD
  geom_violin(trim=TRUE, lwd=0.2, color="black", alpha = .8) + scale_x_discrete(limits=c(1:9)) +  
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.03) + 
  geom_boxplot(width=0.1, outlier.size=0.05, fill="#ffffff", lwd=0.3, color="black") + 
  labs(y="Distance from the Equator (degrees)") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"), 
        axis.title.x=element_blank(), axis.title.y=element_text(size=8), axis.text=element_text(size=6), legend.position="none") + #axis.text.x=element_blank(), 
  annotate(geom="text", size=2.5, x=c(1:9), y=c(20.5), label=as.character(DT.lat2[,2]))
precsiviolin <- ggplot(ccf_precrad_map, aes(x=relationship, y=prec_si, fill=relationship, color=relationship)) +   # MEAN RAD
  geom_violin(trim=TRUE, lwd=0.2, color="black", alpha = .8) + scale_x_discrete(limits=c(1:9)) +  
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.03) + 
  geom_boxplot(width=0.1, outlier.size=0.05, fill="#ffffff", lwd=0.3, color="black") + 
  labs(y="Precipitation Seasonality Index") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"), axis.title.x=element_blank(), 
        axis.title.y=element_text(size=8), axis.text=element_text(size=6), legend.position="none") + #axis.text.x=element_blank(), 
  annotate(geom="text", size=2.5, x=c(1:9), y=c(1.5), label=as.character(DT.prec_si2[,2]))
precradcorrviolin <- ggplot(ccf_precrad_map, aes(x=relationship, y=precradcorr, fill=relationship, color=relationship)) +   # MEAN RAD
  geom_violin(trim=TRUE, lwd=0.2, color="black", alpha = .8) + scale_x_discrete(limits=c(1:9)) +  
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.03) + 
  geom_boxplot(width=0.1, outlier.size=0.05, fill="#ffffff", lwd=0.3, color="black") + 
  labs(y="Precipitation-Radiation Correlation Coefficient") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"), axis.title.x=element_blank(), 
        axis.title.y=element_text(size=7.5), axis.text=element_text(size=6), legend.position="none") + #axis.text.x=element_blank(), 
  annotate(geom="text", size=2.5, x=c(1:9), y=c(0.8), label=as.character(DT.precradcorr2[,2]))
precradcorrviolin + scale_fill_manual(values=myviolpal) + scale_color_manual(values=myviolpal)

all <- ggarrange(mapviolin + scale_fill_manual(values=myviolpal) + scale_color_manual(values=myviolpal),  #myviolpal
          meanradviolin+ scale_fill_manual(values=myviolpal) + scale_color_manual(values=myviolpal),
          meandslviolin + scale_fill_manual(values=myviolpal) + scale_color_manual(values=myviolpal),
          tempviolin+ scale_fill_manual(values=myviolpal) + scale_color_manual(values=myviolpal),
          latabsviolin + scale_fill_manual(values=myviolpal) + scale_color_manual(values=myviolpal),
          precradcorrviolin + scale_fill_manual(values=myviolpal) + scale_color_manual(values=myviolpal),
          labels = c("", "", ""),
          ncol = 3, nrow = 2, align="v")
ggexport(all, plotlist = NULL, filename = "violinplot_all_rel_sif_sig8.png", ncol = 2,
         nrow = 3, width = 2500, height = 1200, pointsize = 1, res = 300,
         verbose = TRUE)
one <- ggarrange(precradcorrviolin + scale_fill_manual(values=myviolpal) + scale_color_manual(values=myviolpal),  #myviolpal,
                 labels = c("", "", ""),
                 ncol = 1, nrow = 1, align="v")
ggexport(one, plotlist = NULL, filename = "/Users/uribem/violinplot_precradcorr_sif_sig.png", ncol = 1,
         nrow = 1, width = 800, height = 400, pointsize = 1, res = 300,
         verbose = TRUE)


 logistic <- multinom(relationship2 ~ MAP + radmean + meantemp + dsl, data=ccf_precrad_map)
z <- summary(logistic)$coefficients/summary(logistic)$standard.errors
p <- (1 - pnorm(abs(z), 0, 1))*2
ccf_precrad_map$relationship2 <- relevel(ccf_precrad_map$relationship, ref = "2")

pairasbio <- pairw.anova(ccf_precrad_map$MAP, ccf_precrad_map$relationship)
plot(pairw.anova(ccf_precrad_map$MAP, ccf_precrad_map$relationship), type=2)


# 6. Prec-Rad corr map ----------------------------------------------------
## read and organize data:
ccfresults <- read.table('/Users/uribem/Google Drive/PhD/Time series/ccf/results/ccf_prec_rad_1deg_1360_4.csv', sep=",", header=TRUE)
corratpeak <- data.frame(x=lons, y=lats, corrcoeff=ccfresults[,3])
corratpeak$corrcoeff[is.na(corratpeak$corrcoeff)] <- -99.99
## create raster from data frame:
corratpeakrast <- rasterFromXYZ(corratpeak, res=c(1,1),
                                crs="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0",
                                digits=5)
## crop and mask
corr_cropped <- crop(corratpeakrast, extent(tropicsaspoly))
corr_masked <- mask(corr_cropped, tropicsaspoly)
corr_masked2 <- mask(corr_masked, meanevi.ras)   # low evi mask
# plot and export:
breakpointsforcorrs <- c(-100,-34,-1.0,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,1.0)  #layers to differentiate colors in the lag map
myKeycorrs <- list(rectangles=list(col = c("gray88","gray60",brewer.pal(8,"PuOr"))), 
                        text=list(c('NA','NS','-1','-0.6','-0.4','-0.2','0','0.2','0.4','0.6'),
                                  at=breakpointsforcorrs, cex=0.6,font=0.4),columns=1)
setEPS()
postscript(filename, width = 1000)
plt <- levelplot(corr_masked2, margin = F, at=breakpointsforcorrs, 
                 col.regions=c("gray88","gray60",brewer.pal(8,"PuOr")),
                 colorkey=FALSE, key=myKeycorrs) #, key=myKeycorrs
print(plt + layer(sp.lines(world.outlines.sp, col = "black", lwd = 0.5)))
dev.off()

display.brewer.all()
