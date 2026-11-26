using Documenter, JuMP, PiecewiseAffineApprox

# The figures of the background pages are static images in src/assets/figures,
# the plots of the manual are made when the documentation is built.

pages = [
    "Introduction" => "index.md",
    "Manual" => [
        "Approximate a function" => "howto/approximate_function.md",
        "Approximate data" => "howto/approximate_points.md",
        "Control the error" => "howto/control_error.md",
        "Choose the number of planes" => "howto/choose_planes.md",
        "Non-convex functions" => "howto/nonconvex.md",
        "Add to a JuMP model" => "howto/add_to_model.md",
        "Plot approximations" => "howto/plot_approximation.md",
        "Save and load" => "howto/save_load.md",
    ],
    "Methods" => [
        "Overview" => "background/overview.md",
        "MILP" => "background/milp.md",
        "Cluster" => "background/cluster.md",
        "Interpol" => "background/interpol.md",
        "Progressive and FullOrder" => "background/progressive.md",
        "Difference of convex" => "background/dc.md",
        "Formulations" => "background/formulations.md",
    ],
    "Extending" => [
        "Add an algorithm" => "howto/add_algorithm.md",
        "Add a formulation" => "howto/add_formulation.md",
    ],
    "API reference" => "reference/api.md",
]

Documenter.makedocs(
    sitename = "PiecewiseAffineApprox",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        edit_link = "main",
        assets = ["assets/figures.css"],
    ),
    doctest = true,
    modules = [PiecewiseAffineApprox],
    checkdocs = :exports,
    pages = pages,
    repo = Remotes.GitHub("sintefore", "PiecewiseAffineApprox.jl"),
)

Documenter.deploydocs(;
    repo = "github.com/sintefore/PiecewiseAffineApprox.jl.git",
)
