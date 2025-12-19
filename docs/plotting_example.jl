using DimensionalData
using PyramidScheme: PyramidScheme as PS
using GLMakie
using YAXArrays

data = rand(20000, 20000)
dd = DimArray(data, (X(1:size(data, 1)), Y(1:size(data, 2))))
dd = YAXArray((X(1:size(data, 1)), Y(1:size(data, 2))), data)
savecube(dd, "test/data/test.zarr")
PS.buildpyramids("test/data/test.zarr")
pyramid = PS.Pyramid("test/data/test.zarr")
# pyramid = PS.Pyramid(dd)

fig = Figure()
sli = Slider(fig[2, 1], range = LinRange(extrema(data)..., 100), startvalue = (maximum(data) + minimum(data))/2, tellwidth = false, tellheight = true)
pyr_obs = lift(sli.value) do v
    pyramid .< v
end

ax, plt = heatmap(fig[1, 1], pyr_obs)
fig