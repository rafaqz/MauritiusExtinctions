

using CUDA
CUDA.allowscalar(false)
includet("map_file_list.jl")
includet("learn.jl")

files = get_map_files()
# imgname = files.mus.surveyor_general_1872_from_gleadow.filename
imgname = files.mus[2].filename
side = 64 
datasets = project_to_training(imgname; side, min_threshold=300, batchsize=256, distcheck=20);
eachindex(datasets.images.train)
i = rand(eachindex(datasets.images.train)); Plots.plot(datasets.images.train[i]; title=datasets.labels.train[i]); Plots.scatter!([(side ÷ 2, side ÷ 2)]) 
# model = Flux.Chain(
#     InstanceNorm(3),
#     resnet.layers[1:end-1]
#     ,
#     Flux.Chain(
#         AdaptiveMeanPool((1, 1)) ,
#         Flux.flatten,
#         Dense(512, datasets.nclasses),
#     )
# )

nbands = 3
model = Flux.Chain(
    # Convolution layers
    # 1st convolution: on nbands * 64 * 64 layers
	# InstanceNorm(nbands),
    Conv((7, 7), nbands => 32, selu; pad=(1, 1)),
	# Conv((7, 7), 32 => 32, selu; pad=SamePad()),
    # MaxPool((2, 2)),
    # # 2nd convolution
    # Conv((7, 7), 32 => 32, selu; pad=(1, 1)),
	# Conv((7, 7), 32 => 32, selu; pad=SamePad()),
    # MaxPool((2, 2)),
    # # 3rd convolution
    # Conv((7, 7), 32 => 64, selu; pad=(1, 1)),
	# Conv((7, 7), 64 => 64, selu; pad=SamePad()),
    # MaxPool((2, 2)),
    # Dropout(0.25),
	# reshape the array
    # Flux.flatten,
    # Dense(64 * 4 * 4 => 500, selu),
    # Dropout(0.25),
	# Finish with `nclasses` on the axis for each observation
    # Dense(500 => datasets.nclasses),
    # Nicer properties for crossentropy loss
    # softmax,
)

x = gpu(model)
typeof(Flux.params(x))

train!(gpu(model), datasets, "train_save"; n_epochs=10)
