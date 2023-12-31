# Direct_sensitivity_method_demonstration_on_benchmark_problems

Implementation of the direct sensitivity method to evaluate the gradient information and solved with fmincon solver in Matlab with the "Simple Reactor Problem". 

The implementation steps are based on the journal paper entitled "Time Scaling Transformation Avoiding Sensitivity Discontinuity for Nonlinear Optimal Control" by Hainan Wang, Xinhong Liu, and Dr. Edward P. Gatzke. 

The "Simple Reactor Problem" is based on the 6.1 Case 1 Simple Batch Reactor problem from Section 6 titled "Case Studies" from the above mentioned paper, inside. This "Simple Reactor Problem" optimization problem is addressed in three distinct ways as detailed in the paper, each termed as a separate 'approach'. Correspondingly, 

The folder named "N = 3 Approach 1" corresponds to the implementation steps of the Approach 1 in the above mentioned paper. 

The folder named "N = 3 Approach 2" corresponds to the implementation steps of the Approach 2 in the above mentioned paper. 

The folder named "N = 3 Approach 3" corresponds to the implementation steps of the Approach 3 in the above mentioned paper. 

To execute the MATLAB code in the aforementioned three folders, in addition to installing MATLAB main body, it is also necessary to install the Symbolic Math Toolbox and the Global Optimization Toolbox. The author used MATLAB version 2023a for developing these MATLAB codes.




The paper can be downloaded at https://pubs.acs.org/doi/epdf/10.1021/acs.iecr.2c02238. Should you find this code segment beneficial, kindly consider citing it as follows:

Wang, H.; Liu, X.; Gatzke, E. P. Time Scaling Transformation Avoiding Sensitivity Discontinuity for Nonlinear Optimal Control. Ind. Eng. Chem. Res. 2023, 62 (36), 14407–14426.


