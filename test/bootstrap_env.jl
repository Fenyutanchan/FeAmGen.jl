import Pkg

root_directory      =   dirname(@__FILE__)
path_project_toml   =   joinpath(root_directory, "Project.toml")
path_manifest_toml  =   joinpath(root_directory, "Manifest.toml")

Pkg.activate(root_directory)

needed_packages_registered  =   String[
    "SymEngine"
]
needed_packages_without_registered  =   String[
    "AmpTools",
    "FeAmGen"
]

if !isfile(path_manifest_toml)
    try rm(path_project_toml) catch end
else
    Pkg.instantiate()
    Pkg.update()
end

if !isfile(path_project_toml)
    for package_name ∈ needed_packages_without_registered
        Pkg.add(url="https://github.com/zhaoli-IHEP/$package_name.jl.git")
    end
    
    for package_name ∈ needed_packages_registered
        Pkg.add(package_name)
    end
else
    installed_packages  =   Pkg.project().dependencies
    for package_name ∈ setdiff(needed_packages_without_registered, keys(installed_packages))
        Pkg.add(url="https://github.com/zhaoli-IHEP/$package_name.jl.git")
    end
    
    for package_name ∈ setdiff(needed_packages_registered, keys(installed_packages))
        Pkg.add(package_name)
    end
end
