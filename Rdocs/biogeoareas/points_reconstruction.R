#install.packages("openxlsx")
#install.packages("rgplates")
library(openxlsx)
library(rgplates)

#set working directory
wd = "/Users/moba692/Documenting_R/GPlates_R_for_Biologists"
setwd(wd)

fn = "world_area_points_v1.xlsx"
xls = openxlsx::read.xlsx(fn, sheet="Dupin2017fig5")
xls

ages = seq(from=1, to=1000, by=1)
reconstructed_points_table = xls

cat("\nReconstructing paleolongitudes (x) & paleolatitudes (y) for ", length(ages), " ages, ending at ", max(ages), "...", sep="")
i=1
for (i in 1:length(ages))
{
  # Nice to print each step so you can see timing
  cat(ages[i], ",", sep="")
  
  tmpresult = rgplates::reconstruct(x=xls[,c("long","lat")], age=ages[i], model="MERDITH2021")
  
  # Change the names of the columns to specify age
  colnames(tmpresult) = paste0(c("long_","lat_"), ages[i])
  tmpresult
  
  # Add tmpresult to table
  reconstructed_points_table = cbind(reconstructed_points_table, tmpresult)
  
}

reconstructed_points_table

head(reconstructed_points_table[,1:7])
tail(reconstructed_points_table[,1:7])


# Create a new workbook
wb <- createWorkbook()

# Add a sheet and write the data
addWorksheet(wb, "Reconstructed_Points_1000mya")
writeData(wb, "Reconstructed_Points_1000mya", reconstructed_points_table)

# Save the workbook
saveWorkbook(wb, "reconstructed_points.xlsx", overwrite = TRUE)
