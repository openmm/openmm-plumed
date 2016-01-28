real4 delta = pos2-pos1;
real r = SQRT(delta.x*delta.x + delta.y*delta.y + delta.z*delta.z);
float2 bondParams = PARAMS[index];
real deltaIdeal = r-bondParams.x;
real deltaIdeal2 = deltaIdeal*deltaIdeal;
energy += bondParams.y*deltaIdeal2*deltaIdeal2;
real dEdR = 4*bondParams.y*deltaIdeal2*deltaIdeal;
dEdR = (r > 0.0f) ? (dEdR/r) : 0.0f;
delta.xyz *= dEdR;
real4 force1 = delta;
real4 force2 = -delta;

