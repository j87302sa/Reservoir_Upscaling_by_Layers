The contents uploaded for this project are concerned with upscaling of reservoir from fine-scaled to coarse-scaled. The softwares used for this projects are MATLAB, ECLIPSE100, PETREL.
These 3 softwares were dependent on each other to produce a desire result which is to simulate the hydrocarbon production data with optimised RF/UR by producing relevant data such as GOR, Water Cut,
Saturation Profile, Production Plateau etc.

The contents uploaded here are mainly focused on upscaling of the reservoir using MATLAB. Upscaling is executed for both porosity and permeability data. Porosity is upscaled using simple "Arithmetic Averaging" method since it is static data. Upscaling permeability can be done by both "Arithmetic Averaging" method and "Flow Based Upscaling" method. The latter was used to executed upscaling of permeability with the aim of having better accuracy and narrow error margin between fine-scaled and coarse-scaled (upscaled) data. The aims and objectives are described as follow:

Aims 
* To shrink the dimension of 120:60:10 into 60:30:10 (Z-direction layers are input as unchanged due to the importance of reservoir thickness for simulation)
* To generate a MATLAB code that can execuete Flow-based or Analytical upscaling method by generating coarse-scale reservoir layers
* To transform the resutant datum into GRDECL format, use them in ECLIPSE100 software to execute field production simulation, develop the model further with
  modification of wells-placement and a choice of secondary reservoir depletion method
* To generate field data in PETREL (GOR, watercut, production plateu, saturation profile)
* To narrow down the error margin of fine-scaled simulation data Vs coarse-scaled simulation data

Objectives
* To tackle the complication of flow-based upscaling coding in MATLAB
* To challenge the error margin to be as narrow as possible
* To examine geological uncertainties
* To analyse comparisons between different secondary reservoir depletion methods

What is Arithematic Averaging/Upscaling Method?

It is simple averaging methods, harmonic and arithmetic averaging. For 
permeability(k) values in series, harmonic averaging method is used whereas arithmetic averaging is used for k values in parallel. The concept is similar to averaging of circuits in electricity currents.  
For x-direction flow, it was assumed that there was a barrier restricting vertical flow and only horizontal flow is possible. Therefore, harmonic averaging is applied for horizontal flow between 2 k values in series. Then, arithmetic averaging is applied to 2 averaged values (harmonic) and final k value is obtained. 
However, for y-direction flow, a barrier restricting horizontal flow is assumed to present. Therefore, arithmetic averaging is applied to 2 parallel k values and then harmonic averaging is applied to 2 averaged k values (arithmetic) to obtain final k value.

What is Flow Based Upscaling Method?

Flow based upscaling or pressure solver method concerns with solving Darcy’s equation twice, first for solving flow (Q) and second for solving permeability (K), by solving the pressure fields with different boundary conditions for x and y in the porous medium. 

The equation used for solving pressure field is 𝑃𝑅𝐸𝑆 = 𝐴\𝐵 where A is the coefficient matrix and B is the right hand side matrix which we need to write the discretized equation to figure out the right-hand side value.  

Coefficient matrix (A) has dimension of (𝑁𝑥 × 𝑁𝑦)×(𝑁𝑥 ×𝑁𝑦) = (2×2)×(2×2) = 4×4. 
Right-hand-side matrix (B) has dimension of 4 × 1 and is equal to 𝑘𝑙𝑒𝑓𝑡 × 𝑖𝑛𝑙𝑒𝑡 𝑝𝑟𝑒𝑠𝑠𝑢𝑟𝑒 for x direction flow and equal to 𝑘𝑡𝑜𝑝 × 𝑖𝑛𝑙𝑒𝑡 𝑝𝑟𝑒𝑠𝑠𝑢𝑟𝑒 for y-direction flow. 

For boundary condition at x-direction flow, the pressure values at top and bottom are set zero and only arbitrarily assumed pressure values at left and right are present and vice versa for boundary condition at y-direction flow. By writing discretization codes to determine permeability values of 2 x 2 matrix (ktop left, ktop right, kbottom right, kbottom left), pressure fields was determined. 

From that determined pressure values and permeability values, flow Q was calculated which later was used to figure out final upscaled permeability values for x and y direction flows
