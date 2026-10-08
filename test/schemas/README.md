# Convention schemas

Copies of the JSON schemas of the zarr-conventions used by the multiscale writer,
taken unchanged from the `v0.1` tag of each repository:

- `multiscales-v0.1.json`: https://raw.githubusercontent.com/zarr-conventions/multiscales/refs/tags/v0.1/schema.json
- `proj-v0.1.json`: https://raw.githubusercontent.com/zarr-conventions/proj/refs/tags/v0.1/schema.json
- `spatial-v0.1.json`: https://raw.githubusercontent.com/zarr-conventions/spatial/refs/tags/v0.1/schema.json

When moving to a new tag, add the new files and update the `Convention` constants in `src/conventions.jl`.

`examples/geospatial-pyramid.json` and `examples/dem-multiresolution.json` are the examples of the
same name from the `v0.1` tag of https://github.com/zarr-conventions/multiscales. The generated layout
is checked against the first, the second has an upsampled level and is used to test reading layouts.
