include("species_common.jl")

# 600 feral pigs in mauritius 2009 !! Lubisi et al

# Set up and run simulations
k = :reu
k = :rod
k = :mus
# pred_keys = (:cat, :black_rat, :norway_rat, :mouse, :pig, :macaque)
# pred_keys = (:cat, :black_rat, :norway_rat, :pig, :mouse)
pred_keys = (:cat, :black_rat, :norway_rat)

include("makie.jl")
include("species_rules.jl")

pred_funcs = (;
    cat =        p -> 2.0f0p.black_rat + 0.5f0p.norway_rat + 10f0p.urban + 2f0p.cleared,
    black_rat  = p -> -0.2f0p.cat - 0.1f0p.norway_rat + 0.5f0p.native + 0.3f0p.abandoned + 0.3f0p.forestry + 1p.urban,
    norway_rat = p -> -0.1f0p.cat - 0.1f0p.black_rat + 1.5f0p.urban - 0.2f0p.native,
    mouse =      p -> -0.3f0p.cat - 0.2f0p.black_rat - 0.2f0p.norway_rat + 0.8f0p.cleared + 1.5f0p.urban,
)[pred_keys]
# pred_pops_aux = isdefined(Main, :pred_pops_aux) ? pred_pops_aux : map(_ -> nothing, dems)
(; ruleset, rules, pred_ruleset, endemic_ruleset, islands, pred_response) = def_syms(
    pred_df, introductions_df, island_endemic_tables, auxs, aggfactor; 
    replicates=nothing, pred_keys, first_year, last_year, extant_extension,
    pred_pops_aux = map(_ -> nothing, dems),
    pred_funcs,
);
(; output, endemic_output, pred_output, init, output_kw) = islands[k];
# @time sim!(output, ruleset; proc=SingleCPU(), printframe=true);
# Optimisation
pred_pops_aux = map(islands) do island
    (; pred_output, init) = island
    @time sim!(pred_output, pred_ruleset; proc=SingleCPU(), printframe=true);
    A = cat(pred_output...; dims=3)
    DimArray(A, (dims(init.pred_pop)..., dims(pred_output)...))
end
sum(getproperty.(pred_pops_aux.mus[Ti=At(2009)], :pig))
sum(getproperty.(pred_pops_aux.mus[Ti=At(2009)], :macaque))
sum(getproperty.(pred_pops_aux.rod[Ti=At(2009)], :cat))

Makie.plot(getproperty.(pred_pops_aux.mus[Ti=At(2009)], :pig))
Makie.plot(getproperty.(pred_pops_aux.mus[Ti=At(2009)], :cat))
Makie.plot(getproperty.(pred_pops_aux.mus[Ti=At(2009)], :norway_rat))
Makie.plot(getproperty.(pred_pops_aux.mus[Ti=At(2009)], :black_rat))
Makie.plot(getproperty.(pred_pops_aux.mus[Ti=At(2009)], :mouse))
Makie.plot(getproperty.(pred_pops_aux.rod[Ti=At(2009)], :macaque))
p = Makie.plot(masks.rod)
Makie.scatter!(p.axis, [(63.43,-19.69)])

# Store so we don't have to run the above
jldsave("sym_setup2_$aggfactor.jld2";
    auxs, pred_pops_aux, pred_response
);


(; ruleset, rules, pred_ruleset, endemic_ruleset, islands) = def_syms(
    pred_df, introductions_df, island_endemic_tables, auxs, aggfactor; 
    replicates=nothing, first_year, last_year, extant_extension,
    pred_pops_aux=map(_ -> nothing, dems),
    pred_keys,
);
k = :rod
k = :mus
(; output, endemic_output, pred_output, init, output_kw) = islands[k];
# @time sim!(endemic_output, endemic_ruleset; proc=SingleCPU(), printframe=true);

@time sim!(output, ruleset; proc=SingleCPU(), printframe=true);
mkoutput = mk_pred(init, pred_ruleset; landcover=lc_all[k], output_kw...)
mkoutput = mk(init, ruleset; landcover=lc_all[k], output_kw..., ncolumns=3)
display(mkoutput)
maximum(output[end].causes)

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
