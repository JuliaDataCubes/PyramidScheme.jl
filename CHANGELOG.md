
File filter
# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Links are updated with Changelog.jl using the command:

```julia
Changelog.generate(
    Changelog.CommonMark(),             # output type
    "CHANGELOG.md";                     # input and output file
    repo = "JuliaDataCubes/PyramidScheme.jl", # default repository for links
)
```

## [0.1.2]

### Added 

- Implemented proper MakieRecipies this should make Observable(::Pyramid) work
- Added `cat` of Pyramids to concatenate on all levels

### Fixed

- Improved computation of the dimensions of the aggregated levels
