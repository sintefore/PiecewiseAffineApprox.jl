using Documenter, JuMP, PiecewiseAffineApprox

# TODO: Consider generating figures as part of build
docdir = @__DIR__
cp(
    joinpath(docdir, "approx.svg"),
    joinpath(docdir, "src/assets/approx.svg"),
    force = true,
)
cp(
    joinpath(docdir, "approx_3D.png"),
    joinpath(docdir, "src/assets/approx_3D.png"),
    force = true,
)
cp(
    joinpath(docdir, "approxanim.mp4"),
    joinpath(docdir, "src/assets/approxanim.mp4"),
    force = true,
)
cp(
    joinpath(docdir, "rotation.mp4"),
    joinpath(docdir, "src/assets/rotation.mp4"),
    force = true,
)

pages = [
    "Introduction" => "index.md",
    "How to" => [
        "Approximate a function" => "howto/approximate_function.md",
        "Approximate points" => "howto/approximate_points.md",
        "Add constraints to a model" => "howto/add_to_model.md",
        "Plot approximations" => "howto/plot_approximation.md",
    ],
    "Background" => "background/background.md",
    "API reference" => "reference/api.md",
]

# TODO: Enable modules
# TODO: Enable doctests
Documenter.makedocs(
    sitename = "PiecewiseAffineApprox",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        edit_link = "main",
        assets = String[],
    ),
    doctest = true,
    # modules = [PiecewiseAffineApprox],
    pages = pages,
    remotes = nothing,
)

Documenter.deploydocs(;
    repo = "github.com/sintefore/PiecewiseAffineApprox.jl.git",
)
