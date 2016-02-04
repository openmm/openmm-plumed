__kernel void addForces(__global const real* restrict forces, __global real4* restrict forceBuffers, __global int* restrict atomIndex) {
    for (int atom = get_global_id(0); atom < NUM_ATOMS; atom += get_global_size(0)) {
        int index = atomIndex[atom];
        real4 f = forceBuffers[atom];
        f.xyz += (real3) (forces[3*index], forces[3*index+1], forces[3*index+2]);
        forceBuffers[atom] = f;
    }
}

