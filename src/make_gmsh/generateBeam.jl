import BenchmarkExample: TimoshenkoBeam

TimoshenkoBeam.generateMsh("msh/beam.msh", transfinite = 101, order = 1)

for n in 2:25
    TimoshenkoBeam.generateMsh("msh/beam_seg2_$n.msh", transfinite = n + 1, order = 1)
    TimoshenkoBeam.generateMsh("msh/beam_seg3_$n.msh", transfinite = n + 1, order = 2)
end
