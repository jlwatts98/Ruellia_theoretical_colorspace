`CalculateO` <-
function(data, data.names=NULL)
{

#data should be a dataframe containing all UNIQUE locality data points in at least three columns
#with the first column being species, second longitude, third latitude.
#This fuction first, subsets by records with more than 1 record for each taxon.
#This function calculates the point-proximity metric (O) introduced
#by Cardillo & Warren Global Ecol. Biogeogr. 2016 25:951–963.
#It calculates great circle distance in meters (rather than Euclidian distance) between points
#as Euclidian distance is not appropriate for points on a sphere.
#The function returns a dataframe of the species combination and O-value

	coos <- setDT(data)[, if (.N > 1) .SD, by = species]
	coos$species <- as.character(coos$species)

##get all combinations of taxa
	combos_O <- combn(unique(coos$species),2, simplify=T)

#meat of the function
	O_vals <- list(); sp_i<-list(); sp_j<-list()

for(k in seq_along(combos_O[1,]))
	{
	i <- combos_O[1,k]; sp_i[k]<- i
  j <- combos_O[2,k]; sp_j[k]<- j

#set up dataframes for each species in species pair "k"	
	sp1 <- coos[which(coos$species==i), ]
	sp1 <- sp1[,c(2,3)]
	sp2 <- coos[which(coos$species==j), ]
	sp2 <- sp2[,c(2,3)]

#calculate nearest distance between each point of species 2 and 1 using great circle distance
	heterosp_sp1_sp2 <- geosphere::distm(sp1, sp2, fun = distHaversine)
#calculate nearest distance between each point of species 2 and 2 using great circle distance
	consp_sp1_sp2 <- geosphere::distm(sp1, sp1, fun = distHaversine)
#remove identical point comparisons
	consp_sp1_sp2[consp_sp1_sp2 == 0] <- NA

#find minimum distance for each unique point comparisons
#first for between species, then for within species
	between_sp1_sp2 <- apply(heterosp_sp1_sp2, 1, min, na.rm=T)
	within_sp1_sp2  <- apply(consp_sp1_sp2, 2, min, na.rm=T)

#calculate O for sp2 
	O_stat_sp1sp2 <- (within_sp1_sp2 - between_sp1_sp2)
	O_stat_sp1sp2 <- sum(O_stat_sp1sp2 > 0)/length(O_stat_sp1sp2)

#do the same for sp1
	heterosp_sp2_sp1 <- geosphere::distm(sp2, sp1, fun = distHaversine)
	consp_sp2_sp1    <- geosphere::distm(sp2, sp2, fun = distHaversine)
	consp_sp2_sp1[consp_sp2_sp1 == 0] <- NA
	between_sp2_sp1 <- apply(heterosp_sp2_sp1, 1, min, na.rm=T)
	within_sp2_sp1  <- apply(consp_sp2_sp1, 2, min, na.rm=T)

	O_stat_sp2sp1 <- (within_sp2_sp1 - between_sp2_sp1)
	O_stat_sp2sp1 <- sum(O_stat_sp2sp1 > 0)/length(O_stat_sp2sp1)

#calcualte O-statistic for the pair and add to output list
	O_XY <- (O_stat_sp1sp2 + O_stat_sp2sp1)/2

	O_vals[k]<-O_XY
	
}
	O_vals<-data.frame(matrix(unlist(O_vals)))
	sp_i<-data.frame(matrix(unlist(sp_i)))
	sp_j<-data.frame(matrix(unlist(sp_j)))
	O_df<-cbind.data.frame(sp_i,sp_j,O_vals);colnames(O_df)=c("sp1","sp2","O")
	O_df2<-O_df[,c(2,1,3)]; colnames(O_df2)=colnames(O_df)
	O_df<-rbind.data.frame(O_df,O_df2)
	return(O_df)
}