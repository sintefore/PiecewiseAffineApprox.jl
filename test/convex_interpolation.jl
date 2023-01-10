# Points to approximate
x = [i for i in -2:0.1:2]
z = x .^ 2

# Test interpolation routine
pwl = PWL.interpolatepw(x, z, optimizer; planes = 5);

@test length(pwl.x) == 6
@test length(pwl.z) == 6

@test pwl.z[2] ≈ 1.44
@test pwl.z[3] ≈ 0.16

# Test convexification
cpwl = PWL.convexify(pwl, optimizer)

@test cpwl.c[2] ≈ -1.6
@test cpwl.c[3] ≈ 0.0
