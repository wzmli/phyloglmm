
################# Data ------------
data("dune", package = "vegan")
data("dune.env", package = "vegan")
data("dune.taxon", package = "vegan")
dune.taxon$abr = rownames(dune.taxon)
# to get full species names in order to build the phylogeny
dune.taxon = left_join(dune.taxon, 
                       data.frame(sp = c("Achillea millefolium", "Agrostis stolonifera", 
                                         "Aira praecox", "Alopecurus geniculatus", "Anthoxanthum odoratum",
                                         "Bellis perennis", "Bromus hordaceus", "Chenopodium album", 
                                         "Cirsium arvense", "Eleocharis palustris", "Elymus repens", 
                                         "Empetrum nigrum", "Hypochaeris radicata", "Juncus articulatus",
                                         "Juncus bufonius", "Scorzoneroides autumnalis", "Lolium perenne", 
                                         "Plantago lanceolata", "Poa pratensis", "Poa trivialis",
                                         "Comarum palustris", "Ranunculus flammula", "Rumex acetosa", 
                                         "Sagina procumbens", "Salix repens", "Trifolium pratense", 
                                         "Trifolium repens", "Vicia lathyroides", 
                                         "Brachythecium rutabulum", "Calliergonella cuspidata")) %>% 
                         mutate(abr = str_replace(sp, "^([A-Za-z]{3,4})[a-z]* (.{4}).*$", "\\1\\2"),
                                abr3 = str_replace(sp, "^([A-Za-z]{3})[a-z]* (.{3}).*$", "\\1_\\2")),
                       by = "abr")
splist_dune = dune.taxon$sp
splist_dune_phylocom = phylomatic_names(taxa = str_replace(splist_dune, "_", " "), format = "isubmit")
# phylomatic(taxa = splist_dune_phylocom, taxnames = F, get = "GET", storedtree = "R20120829")
tree_dune_phylocom = phylomatic(taxa = splist_dune_phylocom, taxnames = F, get = "GET", storedtree = "zanne2014")
# missing mosses Calliergonella_cuspidata. But we do not have traits for mosses neither.
# so it is fine as we will remove mosses from analyses later.

plot(tree_dune_phylocom,  show.node.label = F)
# nodelabels()
add.scale.bar(length = 100)
edgelabels(text = round(tree_dune_phylocom$edge.length, 2), frame = "n")

simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
} # to make the first letter upper case
tree_dune_phylocom = ladderize(tree_dune_phylocom)
tree_dune_phylocom$tip.label = unname(sapply(tree_dune_phylocom$tip.label, simpleCap))

# read traits
dune.traits = read.table(file = "hacked_code/data_clean/dune_traits_Z.txt", header = T)
names(dune.traits)[1] = "abr3"
# update species names
dune.traits$abr3[dune.traits$abr3 == "Pot_pal"] = "Com_pal" # Potentilla palustris -- Compalu pal...
dune.traits$abr3[dune.traits$abr3 == "Leo_aut"] = "Sco_aut" # Leontodom autumnalis -- Scorzoneroides aut...
# merge speciese full names and remove mosses (with NA data)
dune.traits2 = left_join(dune.taxon, dune.traits, by = "abr3") %>% 
  select(sp, SLA:Lifespan) %>% 
  mutate(sp = str_replace(sp, " ", "_")) %>% 
  na.omit()# no traits for mosses.

hist(dune.traits2$SLA)
hist(log(dune.traits2$SLA))
hist(dune.traits2$Height)
hist(log(dune.traits2$Height))
hist(dune.traits2$LDMC)
hist(log(dune.traits2$LDMC))
hist(dune.traits2$Seedmass)
hist(log(dune.traits2$Seedmass))
dune.traits2$SLA = log(dune.traits2$SLA)
dune.traits2$Height = log(dune.traits2$Height)
dune.traits2$LDMC = log(dune.traits2$LDMC)
dune.traits2$Seedmass = log(dune.traits2$Seedmass)
dune.traits2 = rename(dune.traits2, annual = Lifespan, log.sla = SLA,
                      log.height = Height, log.ldmc = LDMC, log.seed.mass = Seedmass)
# annual = 1, perennial = 0
dune.traits2$annual = (as.numeric(as.factor(dune.traits2$annual)) - 2) * -1
# scale traits
dune.traits2[, 2:5] = scale(dune.traits2[, 2:5])

# veg data: wide to long
dune.veg2 = mutate(dune, site = 1:20) %>%
  gather("abr", "freq", -site) %>%
  left_join(select(dune.taxon, sp, abr), by = "abr") %>%
  mutate(sp = str_replace(sp, " ", "_")) %>%
  select(site, sp, freq) %>% tbl_df()
colSums(dune>0)

# envi data
dune.env2 = mutate(dune.env, site = 1:20)
hist(dune.env2$A1)
hist(log(dune.env2$A1))
dune.env2$A1 = log(dune.env2$A1)
names(dune.env2)[1] = "log.A1"
dune.env2[,1] = scale(dune.env2[,1])

# all data together
dune.all = left_join(dune.veg2, dune.traits2, by = "sp") %>%
  left_join(dune.env2, by = "site") %>%
  na.omit() # remove mosses, who do not have traits data
str(dune.all)
dune.all$Moisture = as.numeric(as.character(dune.all$Moisture)) # ordered factor, treat as numeric
dune.all$Use = as.numeric(dune.all$Use)
dune.all$Manure = as.numeric(as.character(dune.all$Manure))
select(dune.all, site, log.A1:Manure) %>% unique
names(dune.all)
# scale numeric traits and envi variables
dune.all[, c(4:7, 9)] = scale(dune.all[, c(4:7, 9)])

# update the phylogeny
dune.phylo2 = drop.tip(tree_dune_phylocom, 
                       tip = tree_dune_phylocom$tip.label[
                         !tree_dune_phylocom$tip.label %in% unique(dune.all$sp)])
plot(dune.phylo2, use.edge.length = T, show.node.label = F)
edgelabels(text = round(dune.phylo2$edge.length,2), frame = "n")

# plot the data along with the phylogeny
dune.veg2.wide = dcast(dune.veg2, site~sp, value.var = "freq")
row.names(dune.veg2.wide) = dune.veg2.wide$site; dune.veg2.wide$site = NULL

#pdf(file = "dune_phylo.pdf", width = 10, height = 12)
par(mar = c(0.1,0.1,1.1,0.0))
plot(dune.phylo2,
     show.tip.label=T, cex=1, x.lim=1200,
     label.offset=4, edge.width=1, direction = "rightwards")
for(i in 1:20){
  tiplabels(tip = which(dune.phylo2$tip.label %in% 
                          names(dune.veg2.wide)[
                            dune.veg2.wide[i, ] > 0]), pch = 20, 
            cex=dune.veg2.wide[i, ][dune.veg2.wide[i, ]>0]/3, 
            adj = c(950-i*30, 0.5),
            col = "blue")
}
# mtext(text = "1958", side = 2, line = 1, at = 1000, col = "blue")
mtext(text = "Sites", side = 3, line = 0, at = 600, col = "black")
dev.off()
