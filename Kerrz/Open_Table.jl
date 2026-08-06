r = FitParam(10.,lower_limit = 1.5, upper_limit = 10., frozen = false)
h = FitParam(5.,lower_limit = 1.5, upper_limit = 50., frozen = false)

rs = collect(range(1.5,10.0,10))
hs = collect(logrange(1.5,50.0,10))

run(`$kerrz emissivity --velocity corotate --photon-index $(model.Γ) --nthreads $(Threads.nthreads()) --ring-like h:$(model.h),x:$(model.r) --output $cur_dir/$emisivity_out_file`)

run(`$kerrz lineprof  --nradii 100 --nangles 200 --spin $(model.a) --incl $(model.θ) --ng $domain_size --rin $R_In --rout $R_Out --emissivity-profile  $cur_dir/$emisivity_out_file --output $cur_dir/$lineprof_out_file`)