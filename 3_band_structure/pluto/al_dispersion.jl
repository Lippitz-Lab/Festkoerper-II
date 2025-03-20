### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 377c3c6c-b464-11ed-1b81-839fde44a868
begin
	using LinearAlgebra

	using Unitful 
	import Unitful: Å, eV
	
	import RangeHelpers:  range, around
	using PGFPlotsX

end

# ╔═╡ f2e2610a-6718-489a-a086-0dbf1424a7e2
md"""
Dies folgt, leicht modifiziert,  **Computational Physics for the Masses. Part 3: Solid Stuff** by Tomas Polakovic (2022) 
 [https://tpolakovic.github.io/posts/cmpm3/](https://tpolakovic.github.io/posts/cmpm3/) unter CC BY NC 4.0
"""

# ╔═╡ 2d96c428-1a40-449f-9a88-fad9c478db49
# ╠═╡ disabled = true
#=╠═╡
	pgfsave("../Alu_free_electron.tikz.tex",myaxis;  include_preamble= false)

  ╠═╡ =#

# ╔═╡ ad76edc3-9970-498a-9b3f-ced94133496d
# ╠═╡ disabled = true
#=╠═╡
	pgfsave("../Alu_empty_lattice.tikz.tex",myaxis2; include_preamble= false)

  ╠═╡ =#

# ╔═╡ 8cb70cd9-bd96-4d3c-8d4d-5bb1ad0d8025
md"""
## Here the real work starts
"""

# ╔═╡ 6864d7b4-7806-4643-9492-360a19b85fe3
begin
	#a₀ = 1.889726125Å
	a₀ = 0.5291772109Å
	Ha = 27.2eV;
end;

# ╔═╡ 6b1193d6-441a-4792-916d-25a4abd258bb
begin 

struct Lattice{T<:Real, N}
    R::Union{T, Array{T}}
    G::Union{T, Array{T}}
    V::T

    function Lattice(R::Union{<:Real,Matrix{<:Real}})
        G = 2π * inv(R)
        R,G = promote(R,G)
        V = det(R)
        dim = isempty(size(R)) ? 1 : first(size(R))
        new{eltype(R), dim}(R, G, V)
    end
end
	
function Lattice(R::Union{<:Quantity, Matrix{<:Quantity}})
    Unitful.NoUnits.(R ./ a₀) |> Lattice
end

function Lattice(a::T, b::T, γ::Real) where T <: Union{Real, Quantity}
    γ = deg2rad(γ)
    R = [a*sin(γ) zero(T);
         a*cos(γ) b]
    
    Lattice(R)
end

function Lattice(a::T, b::T, c::T, α::Real, β::Real, γ::Real) where T <: Union{Real, Quantity}
    α, β, γ = deg2rad.((α, β, γ))
    γ = (cos(α) * cos(β) - cos(γ)) / (sin(α) * sin(β))
    γ = clamp(γ, -1, 1) |> acos # floating point arithmetic could theoretically push the value above ±1

    R = [a*sin(β) -b*sin(α)*cos(γ) zero(T);
         zero(T)  b*sin(α)*sin(γ)  zero(T);
         a*cos(β) b*cos(α)         c]

    Lattice(R)
end;

end

# ╔═╡ fa04d65a-283d-44a2-97cd-d5472786af62
begin
	struct UnitCell{T, N}
	    positions::Vector{T}
	    species::Vector{Symbol}
	
	    function UnitCell(species::Vector{Symbol}, rs::T...) where T
	        positions = collect(rs)
	            if length(species) == length(rs)
	                new{T, length(first(rs))}(positions, species)
	            else
	                throw(DimensionMismatch("Number of species and positions does not match"))
	            end
	    end
	
	end
	
	function UnitCell(species::Symbol, rs...)
	    species = fill(species, length(rs))
	    UnitCell(species, rs...)
	end
	
	function UnitCell(rs...)
	    UnitCell(:none, rs...)
	end;
end

# ╔═╡ ffca9658-4a51-4eab-aa92-82add5e55276
struct Crystal{N}
    lattice::Lattice
    cell::UnitCell

    function Crystal(l::Lattice{T,N}, c::UnitCell) where {T,N}
        new{N}(l, c)
    end
end

# ╔═╡ de4d2eb9-208d-4754-94d8-661b1111d524
Al = Crystal(
    Lattice(2.856Å, 2.856Å, 2.856Å, 60,60,60),
    UnitCell(:Al, [0.0,0.0,0.0])
);

# ╔═╡ fa6955b8-e19f-4a9b-b16c-59f1a3b00e8c
Base.length(c::UnitCell) = c.positions |> length

# ╔═╡ 5a95170b-cf51-4bdf-ab58-1f31ebde9351
function elH(k, n, crystal::Crystal)
    G = crystal.lattice.G
	si = sort(Iterators.product(fill(-n:n,length(k))...) |> collect |> vec, by=norm)
    map(si) do Gn
        q = G' * collect(k .- Gn)
        1/2*norm(q)^2
    end
end;

# ╔═╡ 74dd28e2-da9d-44cc-bad7-c1e86797fedd
function elH2(k, n, crystal::Crystal)
    G = crystal.lattice.G
	si = sort(Iterators.product(fill(-n:n,length(k))...) |> collect |> vec, by=norm)
    map(si) do Gn
        q = G' * collect(k .- Gn)
        1/2*norm(q)^2
    end
end;

# ╔═╡ 8d77a2dc-729d-4fb7-9b8d-fc788fdd9e65

function nfH(k, n::Integer, V::Function, crystal::Crystal)
	G = crystal.lattice.G
	e = sort(Iterators.product(fill(-n:n,length(k))...) |> collect |> vec, by=norm) |> collect
	k = G' * k
	Gs = (G' * collect(g) for g ∈ e)
	H = [V(j - i) for i ∈ Gs, j ∈ Gs]
	H .+= diagm((1/2*norm(k .+ g)^2 for g ∈ Gs) |> collect)
end


# ╔═╡ a937b581-9972-49dc-8de3-8931413fbba9
function kpath(kpoints, dk)
    vertices = reverse([point.second for point ∈ kpoints])
    labels = [point.first for point ∈ kpoints]
    path = [last(vertices)]
    plength = zeros(typeof(dk),1)
    idxs = [1]

    while length(vertices) >= 2
        pop!(path)
        v1 = pop!(vertices)
        v2 = last(vertices)
        dir = v2 .- v1
        dirm = norm(dir)
        segment = [v1 .+ dir .* dd 
            for dd ∈ range(start = 0, stop = 1, step = around(dk/dirm))]
        path = append!(path,segment)
        idxs = push!(idxs, last(idxs) + length(segment) - 1)
    end
    
    plength = append!(plength,
        Iterators.accumulate(
            +,[norm(v2 .- v1) for (v1,v2) ∈ zip(path[1:end-1], path[2:end])]
            ))
    points = [lab => plength[i] for (lab,i) ∈ zip(labels, idxs)]
    (path=path, plength=plength, ppoints=points)
end;

# ╔═╡ e397dc0c-feb0-4d4c-adf0-c9036f3caece
Alks = kpath([
        :Γ => [0,0,0],
        :X => [1/2,0,1/2],
        :W => [1/2,1/4,3/4],
	  	:L => [1/2,1/2,1/2],
		:Γ => [0,0,0],
        :K => [3/8,3/8,3/4],
        :X => [1/2,0,1/2]
        ], 0.01);

# ╔═╡ 0797e31e-5a0c-4016-84b6-d3444338d22b
begin
	alH_el(k) = elH(k, 1, Al)
	es_al_el = sort.(alH_el.(Alks.path));
	
	myaxis = @pgf PGFPlotsX.Axis(
	    {
	        ymin = 0, ymax = 22,
	        xmin = 0, xmax = Alks.plength[end],
			ylabel = raw"E (eV)",
			xtick = [p.second for p ∈ Alks.ppoints], 
			xticklabels = [string(p.first) for p ∈ Alks.ppoints],
			xmajorgrids,
			width= "55mm", height= "55mm", 
	    },	
		);
	
	@pgf for n ∈ 1:6 #length(es_al_el[1])
	    p = Plot(
	        {
	            no_marks,
	        },
	        Coordinates(Alks.plength, [e[n] * 27.2 for e ∈ es_al_el])
	    )
		push!(myaxis, p)
		end

	myaxis
end

# ╔═╡ 74f2cde2-543d-4d66-8970-cf3a77687d2c
length(es_al_el[1])

# ╔═╡ f4f952d6-dd79-4481-a308-e3f9b8569273
begin
	#V(g, Q, q) = ifelse(norm(g) ≈ 0, 0, 4π * Q/(norm(g)^2 .+ q^2));
	V(g, Q, q) = ifelse( (norm(g) ≈ 0) | (norm(g) > 2), 0, 3e-2);
	
	AlV(k) = V(k, 3, 3)
	
	AlH(k) = nfH(k, 1, AlV, Al);
	
	Ales = AlH.(Alks.path) .|> eigvals;
	E0 = Ales[1][1];
	
	
	myaxis2 = @pgf PGFPlotsX.Axis(
	    {
	        ymin = 0, ymax = 16,
	        xmin = 0, xmax = Alks.plength[end],
			ylabel = "E (eV)",
			xtick = [p.second for p ∈ Alks.ppoints], 
			xticklabels = [string(p.first) for p ∈ Alks.ppoints],
			xmajorgrids,
		width= "105mm", height= "55mm", 
	    },	
		);

	@pgf for n ∈ 1:6 #length(es_al_el[1])
	    p = Plot(
	        {
	            no_marks, dashed, gray
	        },
	        Coordinates(Alks.plength, [e[n] * 27.2 for e ∈ es_al_el])
	    )
		push!(myaxis2, p)
		end
	
	@pgf for n ∈ 1:6 #length(Ales[1])
	    p = Plot(
	        {
	            no_marks,
	        },
	        Coordinates(Alks.plength, [(e[n] - E0) * 27.2 for e ∈ Ales])
	    )
		push!(myaxis2, p)
		end

	myaxis2
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PGFPlotsX = "8314cec4-20b6-5062-9cdb-752b83310925"
RangeHelpers = "3a07dd3d-1c52-4395-8858-40c6328157db"
Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[compat]
PGFPlotsX = "~1.6.2"
RangeHelpers = "~0.1.14"
Unitful = "~1.22.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.4"
manifest_format = "2.0"
project_hash = "9c2a00ec62e702dcd86efdf2c13d2ecbda29deac"

[[deps.ArgCheck]]
git-tree-sha1 = "680b3b8759bd4c54052ada14e52355ab69e07876"
uuid = "dce04be8-c92d-5529-be00-80e4d2c0e197"
version = "2.4.0"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.DefaultApplication]]
deps = ["InteractiveUtils"]
git-tree-sha1 = "c0dfa5a35710a193d83f03124356eef3386688fc"
uuid = "3f0dd361-4fe0-5fc6-8523-80b14ec94d85"
version = "1.1.0"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.7.2+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.MacroTools]]
git-tree-sha1 = "72aebe0b5051e5143a079a4685a46da330a40472"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.15"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.6+0"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.OrderedCollections]]
git-tree-sha1 = "cc4054e898b852042d7b503313f7ad03de99c3dd"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.0"

[[deps.PGFPlotsX]]
deps = ["ArgCheck", "Dates", "DefaultApplication", "DocStringExtensions", "MacroTools", "OrderedCollections", "Parameters", "Requires", "Tables"]
git-tree-sha1 = "e5df51ffc01f8771d94c8db2d164be1f6927849c"
uuid = "8314cec4-20b6-5062-9cdb-752b83310925"
version = "1.6.2"

    [deps.PGFPlotsX.extensions]
    ColorsExt = "Colors"
    ContourExt = "Contour"
    MeasurementsExt = "Measurements"
    StatsBaseExt = "StatsBase"

    [deps.PGFPlotsX.weakdeps]
    Colors = "5ae59095-9a9b-59fe-a467-6f913c188581"
    Contour = "d38c429a-6771-53c6-b99e-75d170b6e991"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.RangeHelpers]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "c4d05f8331a339b144a0df913782eb76ffd4f40b"
uuid = "3a07dd3d-1c52-4395-8858-40c6328157db"
version = "0.1.14"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "598cd7c1f68d1e205689b1c2fe65a9f85846f297"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.12.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "c0667a8e676c53d390a09dc6870b3d8d6650e2bf"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.22.0"

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    InverseFunctionsUnitfulExt = "InverseFunctions"

    [deps.Unitful.weakdeps]
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"
"""

# ╔═╡ Cell order:
# ╟─f2e2610a-6718-489a-a086-0dbf1424a7e2
# ╠═de4d2eb9-208d-4754-94d8-661b1111d524
# ╠═e397dc0c-feb0-4d4c-adf0-c9036f3caece
# ╠═5a95170b-cf51-4bdf-ab58-1f31ebde9351
# ╠═0797e31e-5a0c-4016-84b6-d3444338d22b
# ╠═74f2cde2-543d-4d66-8970-cf3a77687d2c
# ╠═2d96c428-1a40-449f-9a88-fad9c478db49
# ╠═74dd28e2-da9d-44cc-bad7-c1e86797fedd
# ╠═8d77a2dc-729d-4fb7-9b8d-fc788fdd9e65
# ╠═f4f952d6-dd79-4481-a308-e3f9b8569273
# ╠═ad76edc3-9970-498a-9b3f-ced94133496d
# ╟─8cb70cd9-bd96-4d3c-8d4d-5bb1ad0d8025
# ╠═377c3c6c-b464-11ed-1b81-839fde44a868
# ╠═6864d7b4-7806-4643-9492-360a19b85fe3
# ╠═6b1193d6-441a-4792-916d-25a4abd258bb
# ╠═fa04d65a-283d-44a2-97cd-d5472786af62
# ╠═ffca9658-4a51-4eab-aa92-82add5e55276
# ╠═fa6955b8-e19f-4a9b-b16c-59f1a3b00e8c
# ╠═a937b581-9972-49dc-8de3-8931413fbba9
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
