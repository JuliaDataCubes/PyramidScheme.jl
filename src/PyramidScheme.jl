module PyramidScheme
using DiskArrayEngine
using DiskArrays
using Zarr

d = zopen("test/data/E048N018T3_rqatrend_VH_A117__thresh_4_8_8_2023-03-22T16:00:28.438.zarr/")

a = d.arrays["layer"]
testdata = a[1000:1063,1000:1063]
n_level = ceil(Int,log2(15000/256))
pixel_base_size = 2^n_level
a.metadata

output_sizes = [ceil(Int,64 / 2^i) for i in 1:5]
output_arrays = [zeros(i,i) for i in output_sizes]



using Statistics
function aggregate_by_factor(xout,x)
    fac = ceil(Int,size(x,1)/size(xout,1))
    for j in 1:size(xout,2)
        for i in 1:size(xout,1)
            xout[i,j] = mean(skipmissing(view(x,((i-1)*fac+1):min(size(x,1),(i*fac)),((j-1)*fac+1):min(size(x,2),(j*fac)))))
        end
    end
end

xout = zeros(10,10)
x = reshape(1:361,19,19)
aggregate_by_factor(xout,x)



xout
end

function all_pyramids!(xout,x)
    xnow = x
    for i in 1:length(xout)
        aggregate_by_factor(xout[i],xnow)
        xnow = xout[i]
    end
end

all_pyramids!(output_arrays,testdata)

output_arrays[1]

heatmap(output_arrays[3])

function generate_pyramid()
# Write your package code here.


end
