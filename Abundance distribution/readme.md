# Estimating Species Relative Abundance distributions
most eBirders also specify how many individuals of each species were observed. So, in this chapter, we’ll take advantage of these counts to model a relative measure of species abundance.

# Multiple steps involved:
1. extract habitat covariayes
2.make a mosaic of the tiles
3. calculate landcover metrics
4. make the prediction surface raster
5. test-train the model
6. assess the models
7. predict the relative abundances

#Most of the scipts are adapted from Ebird Best Practices eBook: https://cornelllabofornithology.github.io/ebird-best-practices/
