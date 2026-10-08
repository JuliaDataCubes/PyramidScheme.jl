# Convention schemas

Copies of the JSON schemas of the zarr-conventions used by the multiscale writer,
taken unchanged from the `v0.1` tag of each repository:

- `multiscales-v0.1.json`: https://raw.githubusercontent.com/zarr-conventions/multiscales/refs/tags/v0.1/schema.json
- `proj-v0.1.json`: https://raw.githubusercontent.com/zarr-conventions/proj/refs/tags/v0.1/schema.json
- `spatial-v0.1.json`: https://raw.githubusercontent.com/zarr-conventions/spatial/refs/tags/v0.1/schema.json

When moving to a new tag, add the new files and update the `Convention` constants in `src/conventions.jl`.

`examples/geospatial-pyramid.json` is the example of the same name from the `v0.1` tag of
https://github.com/zarr-conventions/multiscales, used to check that the generated layout matches it.
