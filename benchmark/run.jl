# julia --project=benchmark benchmark/run.jl [TrussTopology ...]

import Pkg
Pkg.activate(@__DIR__)
# As a package, not `include`d, so that Sequoia's Clarabel extension loads.
Pkg.develop(path=joinpath(@__DIR__, ".."); io=devnull)

using Clarabel, Dates, JuMP, Printf, Sequoia

# Chordal decomposition crashes on these sparse PSD blocks (Clarabel indexes an
# OrderedSet). The blocks are small, so the dense factorization is the reliable path.
const OPTIMIZER = optimizer_with_attributes(Clarabel.Optimizer,
                                            "verbose" => false,
                                            "chordal_decomposition_enable" => false)
const INSTANCES = normpath(joinpath(@__DIR__, "..", "ref", "instances"))

struct BenchRow
    instance::String
    variables::Union{Int,Nothing}
    integers::Union{Int,Nothing}
    status::String
    objective::Float64
    bound::Float64
    gap::Float64
    nodes::Union{Int,Nothing}
    failed::Union{Int,Nothing}
    seconds::Float64
end

"Directory names under `ref/instances`, in alphabetical order."
function instance_groups()
    names = [name for name in readdir(INSTANCES) if isdir(joinpath(INSTANCES, name))]
    return sort!(names)
end

"Sparse SDPA files in one group, in filename order."
function instance_files(group::AbstractString)
    dir = joinpath(INSTANCES, group)
    isdir(dir) || throw(ArgumentError("no instance group named $group"))
    files = String[]
    for name in readdir(dir)
        lower = lowercase(name)
        if endswith(lower, ".dat-s") || endswith(lower, ".dat-s.gz")
            push!(files, joinpath(dir, name))
        end
    end
    return sort!(files, by=basename)
end

"Filename with the `.dat-s` or `.dat-s.gz` suffix removed."
function instance_label(path::AbstractString)
    name = basename(path)
    lower = lowercase(name)
    for suffix in (".dat-s.gz", ".dat-s")
        endswith(lower, suffix) && return name[1:end-length(suffix)]
    end
    return name
end

dash(x::Nothing) = "-"
dash(x::Integer) = string(x)
dash(x::Real) = isfinite(x) ? @sprintf("%.12e", x) : "-"

format_gap(g::Real) = isfinite(g) ? @sprintf("%.4f%%", 100 * g) : "-"
format_seconds(t::Real) = isfinite(t) ? @sprintf("%.2f", t) : "-"

function row_from_result(label, problem, result)
    return BenchRow(label, length(problem.c), length(problem.integers),
                    string(result.status), result.objective, result.bound, gap(result),
                    result.nodes, result.unknown, result.time)
end

error_row(label, seconds) = BenchRow(label, nothing, nothing, "error", NaN, NaN, NaN, nothing, nothing, seconds)

"Load one file, solve it, and return a table row. A thrown instance is recorded as `error`."
function solve_file(path::AbstractString; time_limit::Real)
    label = instance_label(path)
    started = time()
    try
        problem = read_dat_s(path)
        model = jump_model(problem)
        result = solve(model; optimizer=OPTIMIZER,
                       options=Options(time_limit=Float64(time_limit), log_every=Inf))
        return row_from_result(label, problem, result)
    catch err
        err isa InterruptException && rethrow()
        println("  error: ", sprint(showerror, err))
        return error_row(label, time() - started)
    end
end

function write_log(group::AbstractString, rows; time_limit::Real)
    path = joinpath(@__DIR__, group, "README.md")
    mkpath(dirname(path))
    limit = isinf(time_limit) ? "none" : "$(time_limit) s"
    date = Dates.format(Dates.today(), "yyyy-mm-dd")
    open(path, "w") do io
        println(io, "# ", group)
        println(io)
        println(io, "Sequoia with Clarabel, chordal decomposition off. Time limit ", limit, ". ", date, ".")
        println(io)
        println(io, "| Instance | Variables | Integers | Status | Objective | Bound | Gap | Nodes | Failed | Time (s) |")
        println(io, "| --- | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |")
        for row in rows
            cells = (row.instance, dash(row.variables), dash(row.integers), row.status,
                     dash(row.objective), dash(row.bound), format_gap(row.gap),
                     dash(row.nodes), dash(row.failed), format_seconds(row.seconds))
            println(io, "| ", join(cells, " | "), " |")
        end
    end
    return path
end

"Solve every instance in one group, rewriting its README after each one."
function benchmark_group(group::AbstractString; time_limit::Real=60)
    rows = BenchRow[]
    for path in instance_files(group)
        println(group, "/", instance_label(path))
        row = solve_file(path; time_limit)
        push!(rows, row)
        write_log(group, rows; time_limit)
        println("  ", row.status, "  ", format_seconds(row.seconds), "s  nodes=", dash(row.nodes))
    end
    return rows
end

"""
    benchmark(groups=nothing; time_limit=60)

Solve the sparse SDPA groups under `ref/instances` and write `benchmark/<group>/README.md`.
`groups` is a name or a list of names; the default is every group directory. Each instance
stops after `time_limit` seconds (`Inf` runs to completion).
"""
function benchmark(groups=nothing; time_limit::Real=60)
    names = if groups === nothing
        instance_groups()
    elseif groups isa AbstractString
        [String(groups)]
    else
        String[String(g) for g in groups]
    end
    isempty(names) && throw(ArgumentError("no instance groups under $INSTANCES"))
    for name in names
        benchmark_group(name; time_limit)
    end
    return names
end

if !isempty(PROGRAM_FILE) && abspath(PROGRAM_FILE) == @__FILE__
    benchmark(isempty(ARGS) ? nothing : ARGS)
end
