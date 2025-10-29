using SpectralFitting, XSPECModels, Relxill, Warmabs, Plots, LaTeXStrings

function SpectralFitting._invoke_guard!(output, domain, model::XS_Relconv{<:Number})
    for i in eachindex(output)
        if output[i] <= 0
            output[i] = eps(Float64)
        end
        #throw("BAD")
    end
    SpectralFitting.invoke!(output, domain, model)
end



comp_model_abs = AutoCache(XS_WarmAbsorber(),abstol=1e-9)

#= using Colors
c1 = colorant"red"
c2 = colorant"green"
colors = range(c1, stop=c2, length=21) =#

begin
    comp_model_abs  = PhotoelectricAbsorption()*(XS_Relconv()(AutoCache(XS_WarmAbsorber(),abstol=1e-9)*XillverD5())+PowerLaw())
    energy = collect(range(1.0, 12.0, 150))
    A_Fe = 2
    Γ = 3
    θ = 60
    K_PL = 1e5
    K_Disc = 1

#PhotoelectricAbsorption
comp_model_abs.m1.ηH = 0.168

#Relconv
comp_model_abs.c1.θ_obs = θ 
comp_model_abs.c1.r_break = 4 
comp_model_abs.c1.outer_r = 5 

#WarmAbsorber
comp_model_abs.m2.model.column = 1.7
comp_model_abs.m2.model.vturb = 200
comp_model_abs.m2.model.redshift = 0.0658
comp_model_abs.m2.model.Feabund = A_Fe 
comp_model_abs.m2.model.rlogxi = 3.5

#XillverD5
comp_model_abs.a1.K = K_Disc
comp_model_abs.a1.Γ = Γ
comp_model_abs.a1.A_Fe = A_Fe
comp_model_abs.a1.inclination = 60

#powerlaw
comp_model_abs.a2.K = K_PL
comp_model_abs.a2.a = Γ

end

#begin
function DarkPlot(;
        xscale=:log10,
        yscale=:log10,
        xlabel="Energy (keV)",
        ylabel=L"Counts (s$^{-1}$keV$^{-1}$)",
        xlim=(2, 12.0),
        xticks=([2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],["2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"]),
        x_foreground_color_axis=:white, 
        y_foreground_color_axis=:white, 
        x_foreground_color_border=:white, 
        y_foreground_color_border=:white, 
        legend = :bottomright,
        legendfontcolor="white",
        tickfontcolor="white",
        guidefontcolor="white",
        background_color=RGBA(1, 1, 1, 0),
        background_color_outside=RGBA(1, 1, 1, 0))   
    plot(
        xscale=xscale,
        yscale=yscale,
        xlabel=xlabel,
        ylabel=ylabel,
        xlim=xlim,
        xticks=xticks,
        x_foreground_color_axis=x_foreground_color_axis, 
        y_foreground_color_axis=y_foreground_color_axis, 
        x_foreground_color_border=x_foreground_color_border, 
        y_foreground_color_border=y_foreground_color_border, 
        legend = legend,
        legendfontcolor=legendfontcolor,
        tickfontcolor=tickfontcolor,
        guidefontcolor=guidefontcolor,
        background_color=background_color,
        background_color_outside=background_color_outside)
end

DarkPlot()
        
full_model_abs = invokemodel(energy,comp_model_abs) 
plot!(energy[1:end-1], energy[1:end-1].*energy[1:end-1].*full_model_abs./diff(energy), linewidth=1, color="green", label="full model")

comp_model_abs.a2.K = 0
disk_model_abs = invokemodel(energy,comp_model_abs) 
plot!(energy[1:end-1], energy[1:end-1].*energy[1:end-1].*disk_model_abs./diff(energy), linewidth=1, color="blue", label="disc model")

comp_model_abs.a2.K = K_PL
comp_model_abs.a1.K = 0
PL_model_abs = invokemodel(energy,comp_model_abs) 
plot!(energy[1:end-1], energy[1:end-1].*energy[1:end-1].*PL_model_abs./diff(energy), linewidth=1, color="red", label="power law")

#end

#= k=21
begin
plot(x_foreground_color_axis=:white,x_foreground_color_border=:white, y_foreground_color_border=:white, y_foreground_color_axis=:white,ylabel=L"Counts (s$^{-1}$keV$^{-1}$)",xscale=:log10, yscale=:log10, xlim=(2, 12.0), xlabel="Energy (keV)", xticks = ([2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], ["2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"]), legend = :bottomright, legendfontcolor="white",tickfontcolor="white",guidefontcolor="white",background_color= RGBA(1, 1, 1, 0),background_color_outside = RGBA(1, 1, 1, 0))
k=1
for i in collect(range(0,5,21))   
comp_model_abs.model.rlogxi = i
full_model_abs = invokemodel(energy,comp_model_abs) 
plot!(energy[1:end-1], energy[1:end-1].*energy[1:end-1].*full_model_abs./diff(energy), linewidth=1, color=colors[k], label=i)
k+=1
end
plot!()
end
 =#
