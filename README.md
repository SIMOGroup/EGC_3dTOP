# Introduction
We provide here a Matlab-based 3D topology optimization package (EGC-3dTOP) based extruded geometric components (EGCs).
Each EGC is constructed by extruding a convex/non-convex polygon along the axis of the EGC and rounding the ends of the EGC. 
Using an adaptive mapping technique which allows mapping each ECG onto a support domain, the EGCs are mapped onto an analytical 
grid to obtain an effective density field for material interpolation. Moreover, 2D-plane calculations can be utilized to replace 
3D-space calculations to enhance computing efficiency. The positions and the cross-sectional areas of the ECGs are simultaneously 
optimized through the determination of an optimum set of geometry parameters. Some structural benchmark problems were investigated to 
verify the applicability of the proposed approach. Compared with the solid isotropic material penalization (SIMP) approach, 
the underlying approach does not require any filtering or projection techniques. Hence, it can produce a stiffer optimum design with 
an explicit boundary description whilst the number of design variables dramatically reduces. 

1. Structure of EGC-3dTOP software package: 
- main_EGC_3dTOP.m: the main function for implementing EGC-3dTOP to solve 3D topology optimization problems.
- FEM.m: the FEM function.
- density.m: Density.
- other functions
2. How to solve 3D topology optimization problems: 
- Define parameter of problems as given in main_EGC_3dTOP.m
- Run main_EGC_3dTOP.m and wait until the optimization process done.
- Get output

# Contributors
- Van-Nam Hoang
- Hung Nguyen-Xuan

# References:
 Van-Nam Hoang, H. Nguyen-Xuan, Extruded-geometric-component-based 3D topology optimization, Computer Methods in Applied Mechanics and Engineering, in press, 2020 (Q1)
