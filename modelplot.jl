using SpectralFitting, XSPECModels, Relxill, Warmabs, Plots,  LaTeXStrings

begin
column = 4.0
abs_xi = 1.0
disc_xi = 1.0
disc_dens = 17
Γ = 2.5
index = 1.0
a = 0.998
r_abs = 3.0*ISCO(a)
θ = 70.0
A_Fe = 1.0
K_in = 0.05
K_out = 0.009
K_pl = 15000.0
K_abs = -1000000
line_width = 0.1
end

begin
# ISCO function for manual range setting
function ISCO(a::Float64)
    Z_1 = 1.0+((1-a^2)^(1/3))*((1+a)^(1/3)+(1-a)^(1/3))
    Z_2 = (3*a^2+Z_1^2)^(1/2)
    if a >= 0 
 r = 3+Z_2-((3-Z_1)*(3+Z_1+2*Z_2))^(1/2)
    else
 r = 3+Z_2+((3-Z_1)*(3+Z_1+2*Z_2))^(1/2)
    end
end

function SpectralFitting._invoke_guard!(output, domain, model::XS_Relconv{<:Number})
    for i in eachindex(output)
        if output[i] <= 0
            output[i] = eps(Float64)
        end
        #throw("BAD")
    end
    SpectralFitting.invoke!(output, domain, model)
end


inner_disk_abs = XS_Relconv()(GaussianLine()+XillverD5())
outer_disk = XS_Relconv()(XillverD5())
comp_model_abs  = PhotoelectricAbsorption()*(inner_disk_abs+outer_disk+PowerLaw())

function patcher!(p)
    p.c1.a = clamp(p.c1.a, 0, 0.998)
    p.c2.a = clamp(p.c2.a, 0, 0.998)
    p.c1.inner_r = ISCO(p.c1.a)
    p.c1.outer_r = p.c1.inner_r > p.c1.outer_r ? p.c1.inner_r*1.1 : p.c1.outer_r
    p.c2.inner_r = p.c1.outer_r
    p.c1.r_break = (p.c1.inner_r+p.c1.outer_r)/2
    p.c2.r_break = (p.c2.inner_r+p.c2.outer_r)/2
end

patched_comp_model_abs = ParameterPatch(comp_model_abs; patch = patcher!)

patched_comp_model_abs.m1.ηH = 0.168

patched_comp_model_abs.c1.index1 = index
patched_comp_model_abs.c1.index2 = index
patched_comp_model_abs.c1.r_break = (ISCO(a)+r_abs)/2
patched_comp_model_abs.c1.a = a
patched_comp_model_abs.c1.θ_obs = θ
patched_comp_model_abs.c1.inner_r = ISCO(a)
patched_comp_model_abs.c1.outer_r = r_abs
6.9
patched_comp_model_abs.a1.K = K_abs
patched_comp_model_abs.a1.μ.frozen = true
 patched_comp_model_abs.a1.σ = line_width
 

patched_comp_model_abs.a2.K = K_in
patched_comp_model_abs.a2.Γ = Γ
patched_comp_model_abs.a2.A_Fe = A_Fe
patched_comp_model_abs.a2.logXi = disc_xi
patched_comp_model_abs.a2.density = disc_dens
patched_comp_model_abs.a2.inclination = θ

patched_comp_model_abs.c2.index1 = index
patched_comp_model_abs.c2.index2 = index
patched_comp_model_abs.c2.r_break = (400+r_abs)/2
patched_comp_model_abs.c2.a = a
patched_comp_model_abs.c2.θ_obs = θ
patched_comp_model_abs.c2.inner_r = r_abs
patched_comp_model_abs.c2.outer_r = 400

patched_comp_model_abs.a3.K = K_out
patched_comp_model_abs.a3.Γ = Γ
patched_comp_model_abs.a3.A_Fe = A_Fe
patched_comp_model_abs.a3.logXi = disc_xi
patched_comp_model_abs.a3.density = disc_dens
patched_comp_model_abs.a3.inclination = θ

patched_comp_model_abs.a4.K = K_pl
patched_comp_model_abs.a4.a = Γ

patched_comp_model_abs

energy = collect(range(1.0, 12.0, 150))

full_model_abs = invokemodel(energy,patched_comp_model_abs)


inner_disk = XS_Relconv()(XillverD5())

comp_model  = PhotoelectricAbsorption()*(inner_disk+outer_disk+PowerLaw())

patched_comp_model = ParameterPatch(comp_model; patch = patcher!)


patched_comp_model.m1.ηH = 0.168

patched_comp_model.c1.index1 = index
patched_comp_model.c1.index2 = index
patched_comp_model.c1.r_break = (ISCO(a)+r_abs)/2
patched_comp_model.c1.a = a
patched_comp_model.c1.θ_obs = θ
patched_comp_model.c1.inner_r = ISCO(a)
patched_comp_model.c1.outer_r = r_abs

patched_comp_model.a1.K = K_in
patched_comp_model.a1.Γ = Γ
patched_comp_model.a1.A_Fe = A_Fe
patched_comp_model.a1.logXi = disc_xi
patched_comp_model.a1.density = disc_dens
patched_comp_model.a1.inclination = θ

patched_comp_model.c2.index1 = index
patched_comp_model.c2.index2 = index
patched_comp_model.c2.r_break = (400+r_abs)/2
patched_comp_model.c2.a = a
patched_comp_model.c2.θ_obs = θ
patched_comp_model.c2.inner_r = r_abs
patched_comp_model.c2.outer_r = 400

patched_comp_model.a2.K = K_out
patched_comp_model.a2.Γ = Γ
patched_comp_model.a2.A_Fe = A_Fe
patched_comp_model.a2.logXi = disc_xi
patched_comp_model.a2.density = disc_dens
patched_comp_model.a2.inclination = θ

patched_comp_model.a3.K = K_pl
patched_comp_model.a3.a = Γ

patched_comp_model

full_model = invokemodel(energy,patched_comp_model)
end
begin
    plot(energy[1:end-1], energy[1:end-1].*energy[1:end-1].*full_model./diff(energy), linewidth=5,color="#cd1e69", label="Model without absorber")
    plot!(energy[1:end-1], energy[1:end-1].*energy[1:end-1].*full_model_abs./diff(energy), linewidth=5 ,formatter=(_...) -> "", ylabel=L"Counts (s$^{-1}$keV$^{-1}$)",color="#4dfa00",xscale=:log10, yscale=:log10, xlim=(2, 12.0), xlabel="Energy (keV)", xticks = ([2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], ["2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"]), label="Model with absorber", legend = :bottomleft,legendfontcolor="white",tickfontcolor="white",guidefontcolor="white",background_color= RGBA(1, 1, 1, 0),background_color_outside = RGBA(1, 1, 1, 0))
    patched_comp_model_abs.a1.K = 0.0
	patched_comp_model_abs.a2.K = 0.0
	patched_comp_model_abs.a3.K = 0.0
    patched_comp_model_abs.a4.K = K_pl
	inner_model_abs = invokemodel(energy, patched_comp_model_abs)
	plot!(energy[1:end-1], energy[1:end-1].*energy[1:end-1].*inner_model_abs./diff(energy), linewidth=5,color="#6500ff",label = "Power law")
	plot!(x_foreground_color_axis=:white, y_foreground_color_axis=:white)
    plot!(x_foreground_color_border=:white, y_foreground_color_border=:white)
end