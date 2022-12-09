# Randomizing the population map
# Sutherland Bioinformatics
# 2022-12-09

pop.map <-  read.delim(file = "~/Documents/pyes/stacks_workflow/01-info_files/population_map.txt"
                       , header = F
                       , sep = "\t"
                       )

head(pop.map)

# Create randomized vector
pop.map$rand.order <- sample(pop.map$V1)
pop.map

# Create random pop map
pop.map <- pop.map[,c("rand.order", "V2")]
head(pop.map)

# Write out
write.table(x = pop.map, file = "~/Documents/pyes/stacks_workflow_randomized_samples/01-info_files/population_map_rand.txt"
            , sep = "\t"
            , quote = F
            , row.names = F
            , col.names = F
            )
