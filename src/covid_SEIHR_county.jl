module covid_SEIHR_county
using Turing
using AxisArrays
using MCMCChains
using Optim
using LineSearches
using Random
"""
get_vector_param_chain
"""
function get_vector_param_chain(gq_chain, param_name)
  n_samples, n_chains = size(gq_chain)
  l = length(gq_chain[1][param_name]) - 1
  Chains(permutedims(reshape(hcat(getindex.(gq_chain, param_name)...), l + 1, n_samples, n_chains), [2, 1, 3]),
  Symbol.(string.(string(param_name,"["), string.(0:l), "]")))
end

"""
get_gq_chains
"""
function get_gq_chains(model, sample_chains)
  gq = generated_quantities(model, sample_chains)
  scalar_gq_names = collect(filter(key -> gq[1][key] isa Float64, keys(gq[1])))
  vector_gq_names = collect(filter(key -> gq[1][key] isa Vector{Float64}, keys(gq[1])))
  
  gq_subchains = vcat(
    Chains(permutedims(reinterpret(reshape, Float64, NamedTuple{Tuple(scalar_gq_names)}.(gq)), [2, 1, 3]), scalar_gq_names),
    get_vector_param_chain.([gq], vector_gq_names))

  cat(gq_subchains..., dims = 2)
end
export get_gq_chains

# function ChainsCustomIndex(c::Chains, indices_to_keep::BitMatrix)
#   min_length = minimum(mapslices(sum, indices_to_keep, dims = 1))
#   Chains(cat([c.value[indices_to_keep[:, i], :, i][1:min_length, :] for i in 1:size(c.value, 3)]..., dims = 3), c.logevidence, c.name_map, c.info)
# end
function ChainsCustomIndex(c::Chains, indices_to_keep::BitMatrix)
    min_length = minimum(mapslices(sum, indices_to_keep, dims = 1))
  v = c.value
  new_v = copy(v.data)
  new_v_filtered = cat([new_v[indices_to_keep[:, i], :, i][1:min_length, :] for i in 1:size(v, 3)]..., dims = 3)
  aa = AxisArray(new_v_filtered; iter = v.axes[1].val[1:min_length], var = v.axes[2].val, chain = v.axes[3].val)

  Chains(aa, c.logevidence, c.name_map, c.info)
end
export ChainsCustomIndex

# https://gist.github.com/sethaxen/b18e4101a7fd0b2cfe986dd240e99213
# utilities for working with Turing model parameter names using only the DynamicPPL API

"""
    flattened_varnames_list(model::DynamicPPL.Model) -> Vector{Symbol}

Get a vector of varnames as `Symbol`s with one-to-one correspondence to the
flattened parameter vector. 

```julia
julia> @model function demo()
           s ~ Dirac(1)
           x = Matrix{Float64}(undef, 2, 4)
           x[1, 1] ~ Dirac(2)
           x[2, 1] ~ Dirac(3)
           x[3] ~ Dirac(4)
           y ~ Dirac(5)
           x[4] ~ Dirac(6)
           x[:, 3] ~ arraydist([Dirac(7), Dirac(8)])
           x[[2, 1], 4] ~ arraydist([Dirac(9), Dirac(10)])
           return s, x, y
       end
demo (generic function with 2 methods)

julia> flattened_varnames_list(demo())
10-element Vector{Symbol}:
 :s
 Symbol("x[1,1]")
 Symbol("x[2,1]")
 Symbol("x[3]")
 Symbol("x[4]")
 Symbol("x[:,3][1]")
 Symbol("x[:,3][2]")
 Symbol("x[[2, 1],4][1]")
 Symbol("x[[2, 1],4][2]")
 :y
```
"""
function flattened_varnames_list(model::DynamicPPL.Model)
    varnames_ranges = varnames_to_ranges(model)
    nsyms = maximum(maximum, values(varnames_ranges))
    syms = Vector{Symbol}(undef, nsyms)
    for (var_name, range) in varnames_to_ranges(model)
        sym = Symbol(var_name)
        if length(range) == 1
            syms[range[begin]] = sym
            continue
        end
        for i in eachindex(range)
            syms[range[i]] = Symbol("$sym[$i]")
        end
    end
    return syms
end
export flattened_varnames_list

# code snippet shared by @torfjelde
"""
    varnames_to_ranges(model::DynamicPPL.Model)
    varnames_to_ranges(model::DynamicPPL.VarInfo)
    varnames_to_ranges(model::DynamicPPL.Metadata)

Get `Dict` mapping variable names in model to their ranges in a corresponding parameter vector.

# Examples

```julia
julia> @model function demo()
           s ~ Dirac(1)
           x = Matrix{Float64}(undef, 2, 4)
           x[1, 1] ~ Dirac(2)
           x[2, 1] ~ Dirac(3)
           x[3] ~ Dirac(4)
           y ~ Dirac(5)
           x[4] ~ Dirac(6)
           x[:, 3] ~ arraydist([Dirac(7), Dirac(8)])
           x[[2, 1], 4] ~ arraydist([Dirac(9), Dirac(10)])
           return s, x, y
       end
demo (generic function with 2 methods)

julia> demo()()
(1, Any[2.0 4.0 7 10; 3.0 6.0 8 9], 5)

julia> varnames_to_ranges(demo())
Dict{AbstractPPL.VarName, UnitRange{Int64}} with 8 entries:
  s           => 1:1
  x[4]        => 5:5
  x[:,3]      => 6:7
  x[1,1]      => 2:2
  x[2,1]      => 3:3
  x[[2, 1],4] => 8:9
  x[3]        => 4:4
  y           => 10:10
```
"""
function varnames_to_ranges end

varnames_to_ranges(model::DynamicPPL.Model) = varnames_to_ranges(DynamicPPL.VarInfo(model))
varnames_to_ranges(varinfo::DynamicPPL.UntypedVarInfo) = varnames_to_ranges(varinfo.metadata)
function varnames_to_ranges(varinfo::DynamicPPL.TypedVarInfo)
    offset = 0
    dicts = map(varinfo.metadata) do md
        vns2ranges = varnames_to_ranges(md)
        vals = collect(values(vns2ranges))
        vals_offset = map(r -> offset .+ r, vals)
        offset += reduce((curr, r) -> max(curr, r[end]), vals; init=0)
        Dict(zip(keys(vns2ranges), vals_offset))
    end

    return reduce(merge, dicts)
end
export varnames_to_ranges
function varnames_to_ranges(metadata::DynamicPPL.Metadata)
    idcs = map(Base.Fix1(getindex, metadata.idcs), metadata.vns)
    ranges = metadata.ranges[idcs]
    return Dict(zip(metadata.vns, ranges))
end
export varnames_to_ranges

function ChainsCustomOrder(c::Chains, sorted_var_names)
  v = c[sorted_var_names].value
  x, y, z = size(v)
  unsorted = v.axes[2].val
  sorted = collect(zip(indexin(sorted_var_names, unsorted), sorted_var_names))

  new_axes = (v.axes[1], Axis{:var}([n for (_, n) in sorted]), v.axes[3])
  new_v = copy(v.data)
  for i in eachindex(sorted)
      new_v[:, i, :] = v[:, sorted[i][1], :]
  end

  aa = AxisArray(new_v, new_axes...)

  # Sort the name map too:
  namemap = deepcopy(c.name_map)
  namemap = (parameters = sorted_var_names,)
  # for names in namemap
  #     sort!(names, by=string, lt=lt)
  # end

  return Chains(aa, c.logevidence, namemap, c.info)
end
export ChainsCustomOrder

function augment_chains_with_forecast_samples(original_chains::Chains, model, model_forecast, augment_type)::Chains
  n_samples = size(original_chains)[1]
  n_chains = size(original_chains)[3]
  forecast_params_in_order = flattened_varnames_list(model_forecast)  
  new_forecast_params = setdiff(forecast_params_in_order, flattened_varnames_list(model))
  n_new_forecast_params = length(new_forecast_params)
  
  if augment_type == "zeros"
    new_forecast_params_chain = setrange(Chains(zeros(n_samples, n_new_forecast_params, n_chains), new_forecast_params), original_chains.value.axes[1].val)
  elseif augment_type == "randn"
    new_forecast_params_chain = setrange(Chains(randn(n_samples, n_new_forecast_params, n_chains), new_forecast_params), original_chains.value.axes[1].val)
  else
    error("Invalid augment_type")
  end
  
  ChainsCustomOrder(hcat(original_chains, new_forecast_params_chain), forecast_params_in_order) 
end
export augment_chains_with_forecast_samples

"""
Try n_reps different initializations to get MAP estimate.
"""
function optimize_many_MAP(model, n_reps = 100, top_n = 1, verbose = true)
  lp_res = repeat([-Inf], n_reps)
  for i in eachindex(lp_res)
      if verbose
          println(i)
      end
      Random.seed!(i)
      try
          lp_res[i] = optimize(model, MAP(), LBFGS(linesearch = LineSearches.BackTracking())).lp
      catch
      end
  end
  eligible_indices = findall(.!isnan.(lp_res) .& isfinite.(lp_res))
  best_n_seeds =  eligible_indices[sortperm(lp_res[eligible_indices], rev = true)][1:top_n]
  
  map(best_n_seeds) do seed
    Random.seed!(seed)
    optimize(model, MAP(), LBFGS(linesearch = LineSearches.BackTracking())).values.array
  end
end
export optimize_many_MAP

end # module
