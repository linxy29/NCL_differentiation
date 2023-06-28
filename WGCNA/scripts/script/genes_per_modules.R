# Packages
library(rjson)
library(stringr)
library(ggplot2)

# Functions
source("https://gitlab.univ-nantes.fr/E114424Z/veneR/raw/master/loadFun.R?inline=false")
source("https://gitlab.univ-nantes.fr/E198672Y/embryo-functions/-/raw/master/functions.R?inline=false")

res_all <- data.frame()
for (file in snakemake@input){
	print(paste0("file : ", file))
	tmp <- readTable(file, rnames=FALSE)
	n_mods <- table(tmp$Module)
	res <- listToDf(n_mods)
	iter_name <- str_match(file, "SF-([0-9]+)")[2]
	res$iters <- rep(iter_name, nrow(res))
	res_all <- rbind(res_all, res)
}

pdf(snakemake@output)
print("output : ", snakemake@output)
print(ggplot(res_all, aes(x=iters, y=genes)) + geom_point())
dev.off()
