using MPI

function runmpi(file; ntasks = 1)
    # Some MPI runtimes complain if more resources are requested
    # than available.
    if MPI.MPI_LIBRARY == MPI.OpenMPI
        oversubscribe = `--oversubscribe`
    else
        oversubscribe = ``
    end
    MPI.mpiexec() do cmd
        Base.run(
            `$cmd $oversubscribe -n $ntasks $(Base.julia_cmd()) --startup-file=no --project=$(Base.active_project()) $file`;
            wait = true,
        )
        true
    end
end

if !Sys.iswindows()
    runmpi(joinpath(@__DIR__, "distributed", "dlimiter.jl"), ntasks = 3)
end
