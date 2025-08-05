# Multi-Phase-Field-for-FRC-using-Puck-theory

The codes presented here are used in the paper "Revisiting Multi-Phase field for FRCs using Puck failure theory". The codes can be used with AT2, AT1 models. We have tested the files for version Abaqus2020, Abaqus2022, Abaqus2024. 
All the material properties required are presented in the pfall.f file, and the input file corresponding to each example in the paper are given. 

Please note that if your initial mesh is bigger than 70000, third line in the module Kvisual "UserVar(70000,16,4)" has to be changed to the upper limit of the number of the mesh. 

If you are using this code for the academic research or industrial purpose, please cite our papers: 

1. P.K.A.V. Kumar, R. Fleischhaker, A. Dean, H.E. Pettermann "Revisiting multi-phase field model for FRCs using Puck failure theory", Accepted for publication in Composite structures, 2025.
2. A. Dean, P.K.A.V. Kumar, J. Reinoso, C. Gerendt, E. Mahdi, M. Paggi, R. Rolfes, "A multi-phase field fracture model for long fibre reinforced composites based on the Puck theory of failure", Composite structures, Volume 251, 2020, 112446.
3. P.K.A.V. Kumar, A. Dean, J. Reinoso, M. Paggi, "A multi phase-field-cohesive zone model for laminated composites: Application to delamination migration", Composite structures, Volume 276, 2021, 114471.

Please report any problems and suggestions to pavan.kumar@ilsb.tuwien.ac.at

Authors:

Pavan Kumar Asur: asurpavankumar@gmail.com , pavan.kumar@tuwien.ac.at

Aamir Dean, "aamir-dean", a.dean@isd.uni-hannover.de
