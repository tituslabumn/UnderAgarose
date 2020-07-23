rm(list=ls())
library(svDialogs)
library(ggplot2)
library(plyr)
library(stringr)
library(xlsx)
############################################################
#######################START CODE HERE######################
############################################################

setwd(choose.dir(default = "S:/LAB-titus004/AnnikaSchroder/Analysis/no cells touching/Ax2/", caption = "Select folder"))
#make raw and normalized folders
dir.create("Output")
#set frame rate 
framerate <- 15 # dlg_input("Enter frame interval (seconds): ", Sys.info()[""])$res

# manually set experiment number/ID and frame rate
concentration <- dlg_input("Enter concentration of chemoattractant loaded in middle well: ", Sys.info()[""])$res
G120
celltype <- dlg_input("Enter type of cells/genotype: ", Sys.info()[""])$res

# LOAD DATA
file <- list.files(pattern="Tracker_Results*",recursive=FALSE)

df <- read.csv(file, header=T) # open data

#######################
####assign cell IDs####
#######################
#split data to isolate first frame only, and assign cell IDs to each cell
totalframes <- max(df$Slice)
n.slice1 <- count(df$Slice==1)[2,2]
df$ID <- ""
df$distance <- ""

for (i in 1:n.slice1) {          # assign cell IDs to cells in first slice
  df[i,"ID"] <- paste("cell",i)
}

listcellIDs <- unique(df$ID)[1:n.slice1]

dfcells <- ""

for (k in 1:totalframes) {  # totalframe is the same as slice
  if (k == 1) {
    prev <- subset(df, Slice == 1)
  } else {
    prev <- subset(dfcells, Slice == (k-1))
  }
  current <- subset(df, Slice == k)
  
  # IF YOU EXPORTED CENTROID:
   current$X <- as.numeric(current$X)
   current$Y <- as.numeric(current$Y)
   prev$X <- as.numeric(prev$X)
   prev$Y <- as.numeric(prev$Y)
  # IF YOU EXPORTED CENTER OF MASS
  # current$XM <- as.numeric(current$XM)
  # current$YM <- as.numeric(current$YM)
  # prev$XM <- as.numeric(prev$XM)
  # prev$YM <- as.numeric(prev$YM)
   
  #now loop through current slice and find nearest neighbor of the previous slice to assign cell ID
  for(i in 1:nrow(current)) {
    
    #make new distance column
    prev$distance <- NA
    
    #find distance between first point and all cells from previous slice
    for (j in 1:nrow(prev)) {   
      # IF YOU EXPORTED CENTROID:
       prev[j,"distance"] <- sqrt((current[i,"X"] - prev[j,"X"])^2 + (current[i,"Y"] - prev[j,"Y"])^2)
      # IF YOU EXPORTED CENTER OF MASS:
      # prev[j,"distance"] <- sqrt((current[i,"XM"] - prev[j,"XM"])^2 + (current[i,"YM"] - prev[j,"YM"])^2)
    }
    
    #select smallest distance for each point
    cell <- which.min(prev[,"distance"]) #finds row number for smallest value
    if (prev[cell,"distance"] < 41.7) { #if nearest neighbor is less than 15um assign cell ID
      #if (prev[cell,"distance"] < 15) { #if exported as um
      current[i,"ID"] <- prev[cell,"ID"]
      current[i,"distance"] <- prev[cell,"distance"]
    } else { # if nearest neighbor is more than 5um assign new cell name
      lastcell <- unlist(strsplit(tail(listcellIDs,1),split = " "))[2] #identify last cell in the list
      newcell <- as.numeric(lastcell) + 1
      newcellID <- paste("cell",newcell)
      listcellIDs <- c(listcellIDs,newcellID)
      current[i,"ID"] <- newcellID
      # current[i,"ID"] <- "cell x"
    }
  }
  dfcells <- rbind(dfcells,current)
}
rm(current,prev,i,j,k)
dfcells <- dfcells[2:nrow(dfcells),] #remove the first row which is empty

#make sure all numeric columns are numeric so we can make nice plots
dfcells$`Time (min)` <- (as.numeric(dfcells$Slice) - 1)* as.numeric(framerate) / 60 # in minutes
# dfcells$XM <- as.numeric(as.character(dfcells$XM))
# dfcells$YM <- as.numeric(as.character(dfcells$YM))
dfcells$X <- as.numeric(as.character(dfcells$X))
dfcells$Y <- as.numeric(as.character(dfcells$Y))

#plot the change in xy distance for each cell
ggplot(dfcells, aes(x=X, y=Y, color = ID)) + geom_point()
                     
ggsave(paste0("Output/",concentration,"_",celltype,"_raw_cell_xyscatter.png"))

####verify that the number of each cell is not more than the number of frames

n_cells <- aggregate(data.frame(count = dfcells), list(value = dfcells$ID), length)
n_cells <- n_cells[c("value","count.ID")]
n_cells$FLAG <- ""
for (i in 1:nrow(n_cells)) {
  if (n_cells[i,2] <= totalframes) {
    n_cells[i,3] = "false"
  } else {
    n_cells[i,3] = "TRUE"
  }
}

n_cells #check for any "true" flags

# exclude all cells that have fewer than 8 time points (2min)
keepcells <- subset(n_cells,n_cells$count.ID > 7)
list.keepcells <- unlist(keepcells[1])

dfcellstrimmed <- dfcells[dfcells$ID %in% list.keepcells,]

#############################################
######Calculate chemotaxis parameters########
#############################################

# first adjust xy to origin for spider plots
data <- ""
for (j in 1:length(list.keepcells)) { # subset data to loop thru all tracks
  current <- subset(dfcellstrimmed, ID == list.keepcells[j])
  X.ori <- current[1,3]
  Y.ori <- current[1,4]
  current$X.origin <- current$X - X.ori # calculate x and y adjusted to origin
  current$Y.origin <- current$Y - Y.ori
  data <- rbind(data,current) # bind all cell tracks together
}
rm(current,X.ori,Y.ori,j)
data <- data[!apply(data == "", 1, all),] # remove the first row which is empty
colkeep <- c("ID","Time (min)","X.origin","Y.origin", "Circ.")
data <- data[,colkeep]

data$X.origin <- as.numeric(as.character(data$X.origin))
data$Y.origin <- as.numeric(as.character(data$Y.origin))

#plot the change in xy distance for each cell
ggplot(data, aes(x=X.origin, y=Y.origin, color = ID)) + geom_point()

ggsave(paste0("Output/",concentration,"_",celltype,"_trimmed_cell_xyscatter.png"))

# add columns for other metrics we want
data$`[cAMP] (nM)` <- concentration
data$Cells <- celltype
data$`Net Displacement (um)` <- ""
data$`Speed (um/min)` <- ""
data$Directionality <- ""
data$`Chemotactic Index` <- ""
data$`MSD (um^2)` <- ""
data2 <- ""

# calculate distance and speed per cell by subsetting and binding
for (cell in list.keepcells) {
  temp <- subset(data, ID == cell)
  # set first time point as zeroes
  temp[1,8:12] <- 0
  temp[,2:4] <- as.numeric(unlist(temp[,2:4]))
  
  # reset timer to 0 for all cells
  temp[1,2] <- 0
  for (i in 2:nrow(temp)) {
    temp[i,2] <- temp[i-1,2] + 0.25
    rm(i)
  }
  
  # now loop through all other time points + calculate
  for (i in 2:nrow(temp)) {
    # distance
    #                  current x       prev x      current y    prev y 
    temp[i,8] <- sqrt((temp[i,3] - temp[i-1,3])^2 + (temp[i,4] - temp[i-1,4])^2)
    # speed
                  #  net displacement    change in time from prev step
    temp[i,9] <- as.numeric(temp[i,8]) / (temp[i,2] - temp[i-1,2])
  }
  
  for (i in 2:nrow(temp)) {
    # directionality
    totaldistance <- sum(as.numeric(temp[2:i,8]))
    #           final euclidean distance traveled / sum of path length
    temp[i,10] <- sqrt(temp[i,3]^2 + temp[i,4]^2) / totaldistance
    
    # chemotactic index
    #net movement in x / euclidean distance
    #sqrt(temp[i,3]^2 + temp[i,4]^2)
    temp[i,11] <- temp[i,3] / sqrt(temp[i,3]^2 + temp[i,4]^2)
    
    # instantenous mean square displacement
    temp[i,12] <- totaldistance / sqrt(i-1)
    
  }
  
  data2 <- rbind(data2,temp) # bind all cell tracks together
}
rm(temp,cell,i)
data2 <- data2[2:nrow(data2),] #remove the first row which is empty

#save final data frame
write.csv(data2,paste0("./Output/",concentration,"_",celltype,"_processing_results.csv"))

