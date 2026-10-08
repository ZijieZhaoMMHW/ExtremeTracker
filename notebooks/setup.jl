using Pkg
cd(@__DIR__) do
    Pkg.activate(".")
    Pkg.develop(path="..")
    Pkg.instantiate()
end
println("Notebook environment ready. Select a Julia kernel in Jupyter.")
