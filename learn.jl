# using Flux
using StatsBase
using Random
using Flux
using Metalhead
using InvertedIndices
using Images
import Pipe as P
using Pipe: @pipe
using JuliennedArrays
using CUDA
using ProgressLogging
using BSON
using Flux: Losses.logitbinarycrossentropy
using Plots
includet("map_file_list.jl")

function project_to_training(imgname::String; 
    frac_test = 0.3,
    frac_validate = 0.0,
    frac_train = 1 - (frac_test + frac_validate),
    cuda=true,
    min_threshold=200,
    batchsize=64,
    side=224,
    step=8,
    distcheck=4,
)
    projectname::String=string(splitext(imgname)[1], ".json")
    original_image = load_image(imgname)
    project = open_output(MapRasterization.MapSelection, projectname)
    category_image = reshape(project.output, size(original_image))
    keep = BitArray(project.settings.category_keep)
    labels = (1:length(keep))[keep]
    half = side ÷ 2
    centers = CartesianIndices((half+distcheck+1:step:size(original_image, 1)-half-distcheck, half+1+distcheck:step:size(original_image, 2) - half-distcheck))
    images = map(collect(centers)) do I
        i, j = Tuple(I)
        label = category_image[I]
        label == 0 && return missing
        nh = (I + CartesianIndex(0, distcheck), I + CartesianIndex(distcheck, 0), I + CartesianIndex(-distcheck, 0), I + CartesianIndex(0, -distcheck))
        for h in nh
            label == category_image[h] || return missing
        end

        img = view(original_image, i-half:i+half-1, j-half:j+half-1)
        label => img
    end |> skipmissing |> collect
    lab, img = rand(images) 
    # Plots.heatmap(img; title=lab) |> display
    counts = [l => count(x -> first(x) == l, images) for l in labels]
    large_counts = filter(x -> x[2] > min_threshold, counts) # Ignore small classes?
    used_labels = map(first, large_counts)
    used_counts = [l => count(x -> first(x) == l, images) for l in used_labels]
    images_per_label = minimum(map(last, used_counts))
    shuffled = shuffle(images)
    equal_groups = map(used_labels) do l
        group = filter(x -> first(x) == l, shuffled)[1:images_per_label]
        map(group) do (label, image)
            image .= rand(eltype(image))
            (label => image)
        end
    end |> Iterators.flatten |> collect |> shuffle
    n = length(equal_groups)
    ntrain = floor(Int, n * frac_train)
    ntest = floor(Int, n * frac_test)
    nvalidate = floor(Int, n * frac_validate)
    inds = eachindex(equal_groups)
    test_indices = sample(inds, ntest; replace=false)
    validate_indices = sample(inds[Not(test_indices)], nvalidate; replace=false)
    train_indices = inds[Not(vcat(test_indices, validate_indices))]
    indices = (test=test_indices, train = train_indices, validate=validate_indices)
    image_sets = map(indices) do I
        map(last, equal_groups[I])
    end
    label_sets = map(indices) do I
        vcat(map(x -> first(x), equal_groups[I]))
    end
    loaders = map(image_sets, label_sets) do i, l
        if length(i) > 1
            data = permutedims(Float32.(Align(map(channelview, i), 1, 2, 3)), [2, 3, 1, 4]) 
            labels = Flux.onehotbatch(l, used_labels)
            Flux.DataLoader((data, labels); batchsize, shuffle=true)
        else
            nothing
        end
    end
    len = map(length, indices)
    return (; images=image_sets, labels=label_sets, loaders, length=len, nclasses=length(used_labels))
end

"Stores the history through all the epochs of key training/validation performance metrics."
mutable struct TrainingMetrics
    val_acc::Vector{AbstractFloat}
    val_loss::Vector{AbstractFloat}
    TrainingMetrics(n_epochs::Integer) = new(zeros(n_epochs), zeros(n_epochs))
end

"Trains given model for a given number of epochs and saves the model that performs best on the validation set."
function train!(model, data, filename::String; n_epochs::Integer=10)
    model = model |> gpu
    optimizer = ADAM()
    params = Flux.params(model) # transfer learning, so only training last layers

    metrics = TrainingMetrics(n_epochs)

    # zero init performance measures for epoch
    epoch_acc = 0.0
    epoch_loss = 0.0

    # so we can automatically save the model with best val accuracy
    best_acc = 0.0

    # X and y are already in the right shape and on the gpu
    # if they weren't, Zygote.jl would throw a fit because it needs to be able to differentiate this function
    loss(X, y) = logitbinarycrossentropy(model(X), y)

    @info "Beginning training loop..."
    for epoch_idx ∈ 1:n_epochs
        @info "Training epoch $(epoch_idx)..."
        # train 1 epoch, record performance
        @withprogress for (batch_idx, (imgs, labels)) ∈ enumerate(data.loaders.train)
            X = @pipe imgs |> gpu |> float32.(_)
            y = @pipe labels |> gpu |> float32.(_)

            gradients = gradient(() -> loss(X, y), params)
            Flux.Optimise.update!(optimizer, params, gradients)

            @logprogress batch_idx / length(enumerate(data.loaders.train))
        end

        # reset variables
        epoch_acc = 0.0
        epoch_loss = 0.0

        @info "Validating epoch $(epoch_idx)..."
        # val 1 epoch, record performance
        @withprogress for (batch_idx, (imgs, labels)) ∈ enumerate(data.loaders.test)
            X = @pipe imgs |> gpu
            y = @pipe labels |> gpu

            # feed through the model to create prediction
            ŷ = model(X)

            # calculate the loss and accuracy for this batch, add to accumulator for epoch results
            batch_acc = @pipe ((((σ.(ŷ) .> 0.5) .* 1.0) .== y) .* 1.0) |> cpu |> reduce(+, _)
            epoch_acc += batch_acc
            batch_loss = logitbinarycrossentropy(ŷ, y)
            @show batch_loss
            epoch_loss += (batch_loss |> cpu)

            @logprogress batch_idx / length(enumerate(data.loaders.test))
        end
        # add acc and loss to lists
        metrics.val_acc[epoch_idx] = epoch_acc / data.length.test
        metrics.val_loss[epoch_idx] = epoch_loss / data.length.test

        # automatically save the model every time it improves in val accuracy
        if metrics.val_acc[epoch_idx] >= best_acc
            @info "New best accuracy: $(metrics.val_acc[epoch_idx])! Saving model out to $(filename).bson"
            BSON.@save joinpath(@__DIR__, "$(filename).bson")
            best_acc = metrics.val_acc[epoch_idx]
        end
    end

    return model, metrics
end

