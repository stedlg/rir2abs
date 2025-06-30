# rir2abs

**Paper:** Geometry informed estimation of surface absorption profiles from room impulse responses, EUSIPCO 2022 (https://hal.archives-ouvertes.fr/hal-03636502/document)

**Authors:** Stéphane Dilungana, Antoine Deleforge, Cédric Foy, Sylvain Faisan 

main.m allows to perform the geometry informed estimation of surface absorption profiles from room impulse responses.
In this approach, absorption coefficients are estimated through nonlinear constrained optimization on magnitude spectrograms.

Please consider the following steps before running main.m :

- Install the ROOMSIM shoebox room acoustics simulator : https://sourceforge.net/projects/roomsim/
- Create results/ and data/ folders inside the rir2abs/ folder

In main.m, set path_package to rir2abs/ and path_roomsim to the ROOMSIM folder.
In main.m, the "simulate data" section allows to choose parameters for the ROOMSIM simulator.

In main.m, the "main parameters" allows to choose the main model parameters. Default parameters correspond to those used in the paper. 

**Affiliations:** Inria Nancy Grand Est (MULTISPEECH), Cerema (UMRAE), ICube (IMAGeS)




