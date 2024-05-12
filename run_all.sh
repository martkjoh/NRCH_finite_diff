julia -t 4 --project=. article/run/run_1.jl&
julia -t 4 --project=. article/run/run_2.jl&
julia -t 4 --project=. article/run/run_3.jl&
julia -t 3 --project=. article/run/run_4.jl&
julia -t 1 --project=. article/run/run_5.jl&

wait
