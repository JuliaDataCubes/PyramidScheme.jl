import json, sys
import numpy as np, xarray as xr, xproj  # noqa: F401
from topozarr import create_pyramid

out = sys.argv[1]

def dataset(nx, ny, epsg, time=None, res=1.0):
    x = np.arange(nx) * res
    y = (np.arange(ny) * res)[::-1]  # north up
    dims, shape, coords = ("y", "x"), (ny, nx), {"x": x, "y": y}
    if time:
        dims, shape = ("time",) + dims, (time,) + shape
        coords["time"] = np.arange(time)
    ds = xr.Dataset({"v": (dims, np.zeros(shape, "float32"))}, coords=coords)
    return ds.proj.assign_crs(spatial_ref=epsg)

cases = {
    "square": (dataset(32, 32, "EPSG:4326"), dict(levels=3)),
    "non_square": (dataset(48, 32, "EPSG:4326"), dict(levels=3)),
    "extra_dim": (dataset(32, 32, "EPSG:4326", time=3), dict(levels=3)),
    "projected": (dataset(32, 32, "EPSG:3857", time=2, res=100.0), dict(levels=2)),
    "sparse_factors": (dataset(64, 48, "EPSG:4326"), dict(factors=[1, 2, 8])),
    "sum": (dataset(32, 32, "EPSG:4326"), dict(levels=2, method="sum")),
    "odd_size": (dataset(50, 30, "EPSG:4326"), dict(levels=3)),
}
for name, (ds, kw) in cases.items():
    p = create_pyramid(ds, x_dim="x", y_dim="y", **kw)
    with open(f"{out}/{name}.json", "w") as f:
        json.dump(p.attrs, f, indent=1)
    print(name, list(p.attrs["multiscales"]["layout"][1]["transform"].items()))

# A small store written by topozarr, read back by the tests
if len(sys.argv) > 2:
    ds = dataset(32, 32, "EPSG:4326")
    ds["v"].values[:] = np.arange(32 * 32, dtype="float32").reshape(32, 32)
    create_pyramid(ds, levels=2, x_dim="x", y_dim="y").write(sys.argv[2])
