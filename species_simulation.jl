include("species_common.jl")

uncleared = gpu_cleanup(Rasters.modify(BitArray, Raster("uncleared.nc")))

# Build auxiliary rasters
lc_predictions_paths = (
    mus="$outputdir/lc_predictions_mus.nc",
    reu="$outputdir/lc_predictions_reu.nc",
    rod="$outputdir/lc_predictions_rod.nc",
)

# Set up and run simulations
k = :reu
k = :rod
k = :mus
pred_keys = (:cat, :black_rat, :norway_rat, :mouse, :pig)

include("species_rules.jl")
pred_pops_aux = isdefined(Main, :pred_pops_aux) ? pred_pops_aux : map(_ -> nothing, dems)
(; ruleset, rules, pred_ruleset, endemic_ruleset, islands, pred_response) = def_syms(
    pred_df, introductions_df, island_endemic_tables, auxs, aggfactor; 
    replicates=nothing, pred_pops_aux, pred_keys
);
(; output, endemic_output, pred_output, init, output_kw) = islands[k]
# @time sim!(output, ruleset; proc=SingleCPU(), printframe=true);
# Optimisation
pred_pops_aux = map(islands) do island
    (; pred_output, init) = island
    @time sim!(pred_output, pred_ruleset; proc=SingleCPU(), printframe=true);
    A = cat(pred_output...; dims=3)
    DimArray(A, (dims(init.pred_pop)..., dims(pred_output)...))
end
pred_pops_aux

# Store so we don't have to run the above
jldsave("sym_setup.jld2";
    auxs, pred_pops_aux, pred_response
);

(; ruleset, rules, pred_ruleset, endemic_ruleset, islands) = def_syms(
    pred_df, introductions_df, island_endemic_tables, auxs, aggfactor; 
    replicates=nothing, pred_pops_aux
);
(; output, endemic_output, pred_output, init, output_kw) = islands[k]
@time sim!(endemic_output, endemic_ruleset; proc=SingleCPU(), printframe=true);

mkoutput = mk(init, ruleset; landcover=lc_all[k], output_kw..., ncolumns=5)
mkoutput = mk_pred(init, pred_ruleset; landcover=lc_all[k], output_kw...)
sim!(mkoutput, ruleset; proc=SingleCPU(), printframe=true);
display(mkoutput)

# k = :mus
# k = :rod
# k = :reu
# (; output, endemic_output, pred_output, init, output_kw) = islands[k]

# p = Rasters.rplot(lc_all.mus[Ti=At(1700:2018)]; colorrange=(1, 6))
# save("images/landcover_simulation.png", p)
# mkoutput = mk_endemic(init, endemic_ruleset; ncolumns=5, maxpops, pred_pop=pred_pops_aux[k], landcover=lc_all[k], output_kw...)
# display(mkoutput)

@time sim!(endemic_output, endemic_ruleset; proc=SingleCPU(), printframe=true);
cu_endemic_output = Adapt.adapt(CuArray, endemic_output)
CUDA.@profile sim!(cu_endemic_output, endemic_ruleset; proc=CuGPU(), printframe=true)
@time preds = predict_extinctions(endemic_ruleset, islands)
# @time sim!(output, endemic_ruleset; proc=CPUGPU(), printframe=true);
sum(map(xs -> xs .> 10, preds.reu)) .+ first(DynamicGrids.tspan(output))
