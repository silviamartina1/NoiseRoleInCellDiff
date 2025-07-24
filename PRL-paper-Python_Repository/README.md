This repository contains two scripts in Python written to perform the complementary study on the stability analysis regarding the work related to the manuscript: `The Constructive Role of Noise in Cell Fate Induction`. 
The package dependencies are manage with Poetry. Installation steps are given below: 

```
 poetry lock 
 poetry install
```
you can run the scripts with all the steps and see how the whole pipeline can be executed with the following command 

 ```poetry run python VectorFields.py```

According to the parameter settings, the script will generate an output as in Fig.2b in the main text where the vector field is for K<sub>X</sub>=K<sub>Y</sub>=20. This code allows to visualize the vector field of the dynamics when the initial state of the dynamics is on the each point on the grid.

The second script is executed by the following command: 
 
 ```poetry run python Bifurcation_diagram_Fig.2a.py```

It allows to create the Fig.2a in the main text. It allows to study the nature of the steady states by checking the eigenvalues of the system dynamics. 

