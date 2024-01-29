# PyramidScheme.jl

PyramidScheme.jl is a package to easily and efficiently compute pyramids of a given datacube which might be larger than RAM.
PyramidScheme.jl provides the Pyramid struct to handle the layers of a pyramid. 
A pyramid is a collection of subarrays of a larger array so that every layer is half the size in every dimension so that pixels are combined to get to the next level of the pyramid. 
These different layers allow to lazily give an overview of the dataset without having to load the whole array into memory. 

PyramidScheme.jl allows to efficiently compute these Pyramids and to interactively plot them with the Makie.jl ecosystem. 
PyramidScheme.jl is based on the DiskArrays.jl for the handling of larger than memory arrays and on DiskArrayEngine.jl for the computations on these arrays. 
In the future, we plan to be able to lazily define computations on DiskArrays so that these computations are only applied to the pixels that are currently used in the interactive plot.

In this talk we will show how to use PyramidScheme.jl to easily and efficiently compute the Pyramids from a given DiskArray and how to then interactively plot these Pyramids using Makie.jl.
