### A Pluto.jl notebook ###
# v0.19.5

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ d5444f6e-d82f-11ec-27a2-c14197ad75ac
begin
	import Pkg

	Pkg.activate(Base.current_project())
end

# ╔═╡ 72279e98-4e1f-4af9-a844-92069f7823e5
using Revise

# ╔═╡ 40a7ffee-f2b1-4425-bcf0-ffeabab92392
begin
	using OrdinaryDiffEq, ModelingToolkit
end

# ╔═╡ 4d4db583-7cea-457b-9906-bffc00ffb9e8
begin
	using CardiovascularModels 
	using Plots, StatsPlots
	using Turing
	using CSV, DataFrames, DataFramesMeta
end

# ╔═╡ 671f72f8-a4f5-4ddc-9b8e-e6f5c2d8adeb
using PlutoUI, HypertextLiteral

# ╔═╡ 047fa4a0-88f1-4fa6-9636-1c49687fe229
begin
	abp = CSV.File("../parameter_estimation/data/abp.csv") |> DataFrame 
	qrs = CSV.File("../parameter_estimation/data/qrs.csv") |> DataFrame
	# Convert pressure to Pascal
	mmHg2pascal(mmHg) = mmHg * 133.322
	@transform!(abp, :ABP_pascal = mmHg2pascal.(:ABP))
end


# ╔═╡ 1365e025-358b-426a-9567-4d836f44e260
begin
	abp_plot = plot(abp.seconds, abp.ABP, labels = "ABP");
	vline!(abp_plot, qrs.seconds, labels = "QRS")
end


# ╔═╡ 1c5ba960-86ea-4c97-8bca-e5e70f28f29c
begin
@variables t
D = Differential(t)

# Ventricle
@named ventricle = Ventricle(Ees = 200e6,
    Vd = 0., V0 = 0., λ = 33e3, P0 = 10.)

# This cardiac driver function uses a vector of qrs complexes 
# to indicate when a contraction starts.
@named card_driver = qrsDriver(qrs_times = qrs.seconds)

@named aortic_valve = Valve(R = 1e6)

# Aorta
@named aorta1 = ElasticVessel(Ees = 100e6, Vd = 5e-6)
@named aorta2 = ElasticVessel(Ees = 100e6, Vd = 100e-6)

@named systemic_resistance = Resistor(R = 50e6)

eqs_driver = [
    ventricle.drv ~ card_driver.contraction
]

@named vein = ElasticVessel(Ees = 5e6, Vd = 500e-6)
@named mitral_valve = Valve(R = 10e6)
end

# ╔═╡ 697cf823-1650-49a6-99cb-647c3d1c563d


# ╔═╡ d48bca66-835e-43ce-b1dc-7b9c11b378b7
md"### Vessel with Inertia"

# ╔═╡ 380ad5b0-c01d-4895-9ba8-218985e780c4
@named aortaInertia = InertiaCompartment(Area = 4e-4, V = 0.1e-3)

# ╔═╡ a96db9c2-635a-471c-b447-4f10ee9c5981
@named aortaResistance = Resistor(R=1e6)

# ╔═╡ b41261bc-31dd-44ef-9cef-20be2da3149c
begin
	# Connect systems in series.
	# vein.out -> mitral_valve.in
	# mitral_valve.out -> ventricle.in
	# etc.
	eqs_con = serial_connect(vein, 
	    mitral_valve, 
	    ventricle,
	    aortic_valve,
	    aorta1,
		aortaInertia,
		aorta2,
	    systemic_resistance,
	    vein)
	
	eqs_comb = [eqs_driver; eqs_con]
end

# ╔═╡ a8a67950-deb8-41a2-9369-7955a973faec
begin
	# ODE parameters
	# We put 2.25 liters of blood in circulation.
	volume_start = [ventricle.V => 1e-4, 
	    aorta1.V => 0.5e-4,
	    aorta2.V => 1.5e-4,
	    vein.V => 2e-3]
	time_span = (0.0, 10.0)
	
end

# ╔═╡ c6d0eedb-2ed1-42d9-b323-f218052f7335
@named hemo_sys = ODESystem(eqs_comb, t, [], []; 
		systems = [
			ventricle, card_driver, aortic_valve, 
			aorta1, 
			aortaInertia, 
			aorta2, systemic_resistance, vein, mitral_valve]);

# ╔═╡ 3db57c48-9451-42d6-882e-61f8268c972f
problem = ODEProblem(structural_simplify(hemo_sys), volume_start, time_span, [])

# ╔═╡ 3d7129ad-c697-460b-82b0-9c30fe4b1dc7
sol_def = solve(problem, Tsit5(); dtmax = 0.01, reltol = 1e-6);

# ╔═╡ 1222e0ba-d300-4925-8678-6d87b1f3b354
begin
	pascal2mmhg(t, pascal) = t, pascal * 0.00750062
	m32ml(t, m3) = t, m3 * 1e6
end

# ╔═╡ 7dd851ea-6382-4226-9be2-1bb7924fdfc9
md"## Update parameters"

# ╔═╡ c509ac08-0fd3-4b99-bad9-1508d013bd44
begin
	# Extract default parameter values for the ModelingToolkit model.
	p_def = ModelingToolkit.defaults(hemo_sys)

	# Get the names of the parameters in correct order. Used to pass a vector of 
	# parameters correctly.
	p_order = parameters(hemo_sys)
end;

# ╔═╡ 872d51c1-9ed8-4fe8-b305-1cfde2d417eb
md"## Fit model with new parameters"

# ╔═╡ b44ab3e5-5b3c-4fb0-a3a9-da2b6a62b31d


# ╔═╡ 44110d26-b774-4b2a-8861-b6ab4cb794c7
function make_param_sliders(params_ordered, default_values)

	function itr(d_val)
		if (d_val == 0)
			return -1:0.01:1
		else
			return (d_val/2.):(d_val/100.):(d_val*2.)
		end
	end
	
	PlutoUI.combine() do Child
		@htl("""
		<h6>Model parameters</h6>
		<ul>
		$([
			@htl("<li>$(name): $(Child(Slider(
				itr(default_values[name]),
				default = default_values[name],
				show_value = true
				)))</li>")
			for name in params_ordered
		])
		</ul>
		""")
		
	end
end

# ╔═╡ c6f9f403-4214-44d2-8082-7a997abb90b2
@bind reset_params Button("Reset")

# ╔═╡ 68447bce-2fc9-48e8-b75f-64d5379c8be3
reset_params; @bind new_params confirm(make_param_sliders(p_order, p_def))


# ╔═╡ 872e270c-de6f-49f8-bd77-fb4725602adf
sol = solve(problem, Tsit5(); 
                      p = new_params, 
                      dtmax = 0.01, reltol = 1e-4);

# ╔═╡ b238205d-2df9-48f4-9f46-79c4daf9339e
md"## Observed vs fitted"

# ╔═╡ b929bbf0-ac75-4ea2-82b4-03329b635cf2
let
	p1 = plot(sol, vars=[(pascal2mmhg, 0,aorta2.P)], tspan =  (1,6),
                ylabel = "mmHg", label = "Model")
	plot!(p1, abp.seconds, abp.ABP, label = "Observed", 
	legend=:outertop)
end


# ╔═╡ 82f482d0-5d2b-4f4a-a1dc-4a5d8cc16a32
md"## Pressure plots [mmHg]"

# ╔═╡ 9a4724b1-50a0-4fae-8cae-fa4bf82f5be5
let
	plot_p1 = plot(sol, vars=[(pascal2mmhg, 0,ventricle.P), (pascal2mmhg, 0,aorta1.P)], tspan =  (1,6));
	plot_p2 = plot(sol, vars=[(pascal2mmhg, 0, vein.P)], tspan =  (1,6));
	plot(plot_p1, plot_p2, layout = (2,1), ylabel = "mmHg", legend = :outertop)
	
end

# ╔═╡ 85ed3cef-30cb-4f28-bacb-1e605244e5ba
md"## Volume plots [ml]"

# ╔═╡ 7d684352-d2d7-4c1e-90ae-dfdcce7fe96c
let
	plot_v1 = plot(sol, vars=[(m32ml, 0, ventricle.V), (m32ml, 0, aorta1.V)], tspan =  (1,6));
	plot_v2 = plot(sol, vars=[(m32ml, 0, vein.V)], tspan =  (1,6));
	plot(plot_v1, plot_v2, layout = (2,1), ylabel = "ml", legend = :outertop)
end

# ╔═╡ 2f20d7a7-4f08-4ac7-a8d8-6739ec7a4dc2
md"Unused code"

# ╔═╡ bf79667a-29ab-41b4-ad11-61c9bb8dfc1a
md"""
Model parameteres:

Ventricle.Ees: $(@bind ventricle_Ees Slider(10e6:10e6:1000e6, default = 200e6, show_value = true))


Aorta.Ees: $(@bind aorta_Ees Slider(10e6:10e6:1000e6, default = 200e6, show_value = true))

systemic\_resistance.R: $(@bind systemic_resistance_R Slider(10e6:1e6:100e6, default = 50e6, show_value = true))

C (cardiac driver): $(@bind C Slider(0:0.01:1, default = 0.2, show_value = true))

"""

# ╔═╡ ba8b5f1d-dc0d-48fa-b7ac-d81250393525
# Update parameters.
    # varmap_to_vars converts a dictionary of parameter values to a correctly ordered vector.
pnew = ModelingToolkit.varmap_to_vars([
	aorta1.Ees => aorta1_Ees,
 	ventricle.Ees => ventricle_Ees,
    systemic_resistance.R => systemic_resistance_R,       
    card_driver.C => C],
		p_order, 
        defaults = p_def)

# ╔═╡ Cell order:
# ╠═d5444f6e-d82f-11ec-27a2-c14197ad75ac
# ╠═40a7ffee-f2b1-4425-bcf0-ffeabab92392
# ╠═4d4db583-7cea-457b-9906-bffc00ffb9e8
# ╠═671f72f8-a4f5-4ddc-9b8e-e6f5c2d8adeb
# ╠═72279e98-4e1f-4af9-a844-92069f7823e5
# ╠═047fa4a0-88f1-4fa6-9636-1c49687fe229
# ╠═1365e025-358b-426a-9567-4d836f44e260
# ╠═1c5ba960-86ea-4c97-8bca-e5e70f28f29c
# ╠═697cf823-1650-49a6-99cb-647c3d1c563d
# ╠═d48bca66-835e-43ce-b1dc-7b9c11b378b7
# ╠═380ad5b0-c01d-4895-9ba8-218985e780c4
# ╠═a96db9c2-635a-471c-b447-4f10ee9c5981
# ╠═b41261bc-31dd-44ef-9cef-20be2da3149c
# ╠═a8a67950-deb8-41a2-9369-7955a973faec
# ╠═c6d0eedb-2ed1-42d9-b323-f218052f7335
# ╠═3db57c48-9451-42d6-882e-61f8268c972f
# ╠═3d7129ad-c697-460b-82b0-9c30fe4b1dc7
# ╠═1222e0ba-d300-4925-8678-6d87b1f3b354
# ╟─7dd851ea-6382-4226-9be2-1bb7924fdfc9
# ╠═c509ac08-0fd3-4b99-bad9-1508d013bd44
# ╟─872d51c1-9ed8-4fe8-b305-1cfde2d417eb
# ╠═872e270c-de6f-49f8-bd77-fb4725602adf
# ╠═b44ab3e5-5b3c-4fb0-a3a9-da2b6a62b31d
# ╠═44110d26-b774-4b2a-8861-b6ab4cb794c7
# ╟─68447bce-2fc9-48e8-b75f-64d5379c8be3
# ╟─c6f9f403-4214-44d2-8082-7a997abb90b2
# ╟─b238205d-2df9-48f4-9f46-79c4daf9339e
# ╠═b929bbf0-ac75-4ea2-82b4-03329b635cf2
# ╟─82f482d0-5d2b-4f4a-a1dc-4a5d8cc16a32
# ╠═9a4724b1-50a0-4fae-8cae-fa4bf82f5be5
# ╟─85ed3cef-30cb-4f28-bacb-1e605244e5ba
# ╠═7d684352-d2d7-4c1e-90ae-dfdcce7fe96c
# ╠═2f20d7a7-4f08-4ac7-a8d8-6739ec7a4dc2
# ╟─bf79667a-29ab-41b4-ad11-61c9bb8dfc1a
# ╠═ba8b5f1d-dc0d-48fa-b7ac-d81250393525
