#' ---
#' title: "Diversity analysis using rdiversity and vegan packages"
#' author: "Rory Thomson"
#' date: '`r format(Sys.Date(), "%B %d %Y")`'
#' output: html_document
#' ---
#'
#' ## Load the packages needed

library(Rory646BCI)
library(rdiversity)
library(magrittr)
library(vegan)

#' ## Setting up dataframes and metacommunity for rdiversity analysis

# Copy of taxonomy dataframe

taxonomy <- as.data.frame(bci_taxa)

# Order the species correctly

rownames(taxonomy) <- taxonomy$GenusSpecies

taxonomy <- taxonomy[row.names(bci_2010),]

# Creating metapopulation

# Copy of 2010 dataframe

meta <- metacommunity(bci_2010)

#' ## Gamma Diversity
#'
## Diversity of whole metacommunity

results <- meta_gamma(meta, qs = seq(from = 0, to = 5))

plot(diversity ~ q, type = "l", data = results)

## Taking into account similarity between species

taxSim <- taxonomy %>%
  tax2dist(c(GenusSpecies = 0, Genus = 1, Family = 2, Other = 3)) %>%
  dist2sim("linear")

metatax <- metacommunity(bci_2010, taxSim)

results4tax <- meta_gamma(metatax, qs = seq(from = 0, to = 5))

# Add line to previous graph

lines(diversity ~ q, col = 2, data = results4tax)

## Taking into account similarity between genera

taxSim <- taxonomy %>%
  tax2dist(c(GenusSpecies=0, Genus=0, Family= 1, Other= 2)) %>%
  dist2sim("linear")

metatax <- metacommunity(bci_2010, taxSim)

results4tax <- meta_gamma(metatax, qs = 0:5)

lines(diversity ~ q, col = 3, data = results4tax)

## Taking into account similarity between families

taxSim <- taxonomy %>%
  tax2dist(c(GenusSpecies=0, Genus=0, Family=0, Other=1)) %>%
  dist2sim("linear")

metatax <- metacommunity(bci_2010, taxSim)

results4tax <- meta_gamma(metatax, qs = 0:5)

lines(diversity ~ q, col = 4, data = results4tax)

#' ### Conclusions
#'
#' As you can see in each analysis as q increases, and therefore weighting of abundance in analysis, the
#' diversity decreases. You can also see that diversity decreases when calculating the analysis of the
#' metacommunity (black), taking account similarity of species (red), similarity between genera (green)
#' and similarity between families (blue).
#'

#' ## Alpha diversity
#'
## Diversity of whole metacommunity

results <- norm_meta_alpha(meta, qs = seq(from = 0, to = 5))

plot(diversity ~ q, type = "l", data = results)

## Taking into account similarity between species

taxSim <- taxonomy %>%
  tax2dist(c(GenusSpecies = 0, Genus = 1, Family = 2, Other = 3)) %>%
  dist2sim("linear")

metatax <- metacommunity(bci_2010, taxSim)

results4tax <- norm_meta_alpha(metatax, qs = 0:5)

# Add line to previous graph

lines(diversity ~ q, col =2, data = results4tax)

## Taking into account similarity between genera

taxSim <- taxonomy %>%
  tax2dist(c(GenusSpecies=0, Genus=0, Family= 1, Other= 2)) %>%
  dist2sim("linear")

metatax <- metacommunity(bci_2010, taxSim)

results4tax <- norm_meta_alpha(metatax, qs = 0:5)

lines(diversity ~ q, col = 3, data = results4tax)

## Taking into account similarity between families

taxSim <- taxonomy %>%
  tax2dist(c(GenusSpecies=0, Genus=0, Family=0, Other=1)) %>%
  dist2sim("linear")

metatax <- metacommunity(bci_2010, taxSim)

results4tax <- norm_meta_alpha(metatax, qs = 0:5)

lines(diversity ~ q, col = 4, data = results4tax)

#' ## Conclusions
#'
#' As you can see in each analysis as q increases, and therefore weighting of abundance in analysis, the
#' diversity decreases. You can also see that diversity decreases when calculating the analysis of the
#' metacommunity (black), taking account similarity of species (green), similarity between genera (red) and
#' similarity between families (blue).
#'
#' When compared to the gamma diversity graph, you can see that the diversity calculated is lower is this
#' due to alpha diversity being a measure of the average diversity across subcommunities whereas the gamma
#' diversity is a measure of the whole metacommunity.
#'

#' ## Beta Diversity
#'
#' Diversity measure for metapopulation
#'
beta.meta <- raw_meta_beta(meta, qs = 0:5)

plot(diversity ~ q, data = beta.meta, type = "line")

#' ## Conclusion
#'
#' As you increase the weighting on abundance the beta diversity, a measure of the difference between
#' subcommunities, increases.

#' Difference in diversity between quadrats

beta.div <- raw_sub_beta(meta, qs = 1)

hist(beta.div$diversity)

#' The histogram has the majority of subcommunities in the same interval in the histogram, showing that
#' there is little difference between the diversity of the majority of subcommunities.
#'
#'
#' ## Analysis using Vegan Package
#'
#'  Preparing dataset for vegan package

# Swap rows and columns, for format suitable to use with vegan

abundances <- t(bci_2010)

# Calculate species richness (gamma diversity)

specnumber(abundances, groups = 1)

#' Jaccard
beta.jac <- betadiver(abundances, measure = "j")

jac.data <- as.matrix(beta.jac$a)

jac.means <- rowMeans(jac.data)

plot(jac.means)


#' Sorensen

beta.sor <- betadiver(abundances, measure = "sor")

sor.data <- as.matrix(beta.sor$a)

sor.means <- rowMeans(sor.data)

plot(sor.means)


#' ## Conclusions
#'
#' As you can see both diversity measures follow the same general trend and the lowest value for both is
#' found when the index = 100.
