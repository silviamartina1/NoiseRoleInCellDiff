# NoiseRoleInCellDiffentiation 
This repository collects the code scripts written to perform stochastic simulations, vector fields stability analysis on the toggle-switch model studied in the paper : "The Construction Role of the Noise in Cell Fate Induction". The codes are two scripts in Python and other five scripts in Julia programming languagues. The script `Bifurcation_diagram_Fig2a.py` contains the code to generate Fig.2a. in the main text. 
The script `VectorFields` is able to generate the Fig.2b in the main text and as well supplementary figures present in the SM (Supplementary Material Text) related to the stability analysis by setting the two parameters K<sub>X</sub> and K<sub>Y</sub>. The files .jl written in Julia programming language are referred to the stochastic simulations. Specifically the file `Toggle_Switch_parameters_SYM.jl` and `Toggle_Switch_parameters_ASYM.jl` are those where the parameters are set up, respectively for the symmetric case (K<sub>X</sub>=K<sub>Y</sub>) and for the asymmetric case (K<sub>X</sub> not equal to K<sub>Y</sub>). 
The main script is in the file `Toggle_Switch_main` where a function to perform stochastic simulation is implemented. This main script has been used to perform dynamics simulations, like those used to create the Fig.2c.
This function calls other functions at each time step that are called from the script `Toggle_Switch_functions`. Those are functions performing specific single operation each time step along the total simulation time such as the update of the molecular concentrations and the computation of the mean and the variance of cell distribution along X and Y directions. 
The file `Toggle_Switch_main_variance_SYM` describes the analysis of the variance on the stable manifold. 


## License
This project is licensed under the Apache License 2.0. For more details, see the [LICENSE](./LICENSE) file.
