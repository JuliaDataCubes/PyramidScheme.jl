# topozarr reference metadata

Root group attributes written by [topozarr](https://github.com/carbonplan/topozarr) 0.2.0
(commit 8635079) for the cases in `generate.py`, and a small store written by it (`store.zarr`).
The tests compare them against the metadata written by PyramidScheme, which differs only in
the identity `transform` of level 0 and the name of the resampling method (`"average"` for `mean`).

Regenerate with:

```sh
uv run --no-project --python 3.12 --with topozarr==0.2.0 python generate.py . store.zarr
```

`store-0.1.9.zarr` was written by topozarr 0.1.9 (commit d661b5a), which measured the
`translation` between levels between cell centers. It checks that such stores can still be read.
