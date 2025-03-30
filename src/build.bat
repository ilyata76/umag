Rem Build script for MICROMAG5

ifort /O3  magmod.f90 -c magmod.obj
ifort /O3 micromag.f90 -c micromag.obj
ifort /O3 neighbour.f90 -c neighbour.obj
ifort /O3 llg_step.f90 -c llg_step.obj
ifort /O3 norm_vector.f90 -c norm_vector.obj
ifort /O3 anisotropy.f90 -c anisotropy.obj
ifort /O3 ash.f90 -c ash.obj
ifort /O3 create_mesh.f90 -c create_mesh.obj
ifort /O3 create_regions.f90 -c create_regions.obj
ifort /O3 demag.f90 -c demag.obj
ifort /O3 demag_tensor.f90 -c demag_tensor.obj
ifort /O3 exchange.f90 -c exchange.obj
ifort /O3 f.f90 -c f.obj
ifort /O3 f1.f90 -c f1.obj
ifort /O3 f2.f90 -c f2.obj
ifort /O3 fun.f90 -c fun.obj
ifort /O3 g.f90 -c g.obj
ifort /O3 g1.f90 -c g1.obj
ifort /O3 g2.f90 -c g2.obj
ifort /O3 gun.f90 -c gun.obj
ifort /O3 initial_magnetization.f90 -c initial_magnetization.obj
ifort /O3 nxx.f90 -c nxx.obj
ifort /O3 nxy.f90 -c nxy.obj
ifort /O3 output_magfield.f90 -c output_magfield.obj
ifort /O3 vector_cross_vector.f90 -c vector_cross_vector.obj
ifort /O3 vector_dot_scalar.f90 -c vector_dot_scalar.obj
ifort /O3 vector_dot_vector.f90 -c vector_dot_vector.obj
ifort /O3 vector_equal.f90 -c vector_equal.obj
ifort /O3 vector_plus_vector.f90 -c vector_plus_vector.obj
ifort /O3 zeeman.f90 -c zeeman.obj
ifort /O3 impulse.f90 -c impulse.obj

link micromag magmod neighbour llg_step norm_vector anisotropy ash create_mesh create_regions demag demag_tensor exchange f f1 f2 fun g g1 g2 gun initial_magnetization nxx nxy output_magfield vector_cross_vector vector_dot_scalar vector_dot_vector vector_equal vector_plus_vector zeeman impulse