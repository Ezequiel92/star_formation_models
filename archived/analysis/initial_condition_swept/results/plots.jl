using Glob
using DelimitedFiles
using Plots
pyplot()

filelist_basic = glob("basic_model_*.dat", "data/")
filelist_simp = glob("simplified_model_*.dat", "data/")
filelist_full_const = glob("full_model_constant_tau_no_IO_*.dat", "data/")
filelist_full_var = glob("full_model_variable_tau_no_IO_*.dat", "data/")
filelist_full_kennicutt = glob("full_model_kennicutt_no_IO_*.dat", "data/")

function plotColumn(data, N, name)
	labels = ["ionized gas fraction", 
			"atomic gas fraction", 
			"molecular gas fraction", 
			"metallicity", 
			"SFR", 
			"star fraction"]
	models = ["basic model" "simplified model" "full model constant τ - no IO" "full model variable τ - no IO" "full model Kennicutt - no IO"]
	
	x = data[1][:,1]
	if N == 7
		y = [data[2][2:end,N], 
			(data[3][2:end,N] ./ data[3][2:end,end]), 
			(data[4][2:end,N] ./ data[4][2:end,end]), 
			(data[5][2:end,N] ./ data[5][2:end,end])]
		x = x[2:end]
		models = models[:,2:end]
	else
		y = [data[1][:,N], 
			data[2][:,N], 
			(data[3][:,N] ./ data[3][:,end]), 
			(data[4][:,N] ./ data[4][:,end]), 
			(data[5][:,N] ./ data[5][:,end])]
	end
	
	plot(x, y, 
		xaxis=(:log, (1e-4, Inf)), 
		yaxis=:log,
		title=head,
		label=models,
		legend=:outerbottom,
		xlabel="t [Gyr]",
		ylabel=labels[N - 1],
		size=(800, 800),
		lw=3,
		xtickfontsize=10,
		ytickfontsize=10,
		titlefontsize=15,
		xguidefontsize=12,
		yguidefontsize=12,
		legendfontsize=12
	)
	savefig(replace(replace(replace(replace(name, 
					"basic_model" => labels[N - 1]), 
					" " => "_"), 
					".dat" => ".png"), 
					"data" => "plots"))
end


for (fl_basic, fl_simp, fl_fullC, fl_fullV, fl_fullK) in zip(filelist_basic, filelist_simp, filelist_full_const, filelist_full_var, filelist_full_kennicutt)
	multi_model_data = [readdlm(fl_basic, '\t', Float64, '\n'; skipstart=1),
						readdlm(fl_simp, '\t', Float64, '\n'; skipstart=1),
						readdlm(fl_fullC, '\t', Float64, '\n'; skipstart=1),
						readdlm(fl_fullV, '\t', Float64, '\n'; skipstart=1),
						readdlm(fl_fullK, '\t', Float64, '\n'; skipstart=1)
						]
	open(fl_basic) do fl_head
        global head = readline(fl_head)
    end
	
	for param in 2:7
		plotColumn(multi_model_data, param, fl_basic)
	end

end 







