# A3 Lab
This project models different kinds of thin-film optical filters using a transfer matrix model. The versatility of the model allows multiple kinds of optical filters to be modelled including anti-reflection coatings, neutral density filters, highly reflective DBR stacks and single wavelength transmission VCSELs.
The data used for refractive index is in the repository also, with the files named after the materials.
## Files and their functions
Each file in this repository is named after the task in the A3 lab script that they are deisgned to provide answers for.
* transfer_matrix function (first used in file task7_together) takes as input:
  * incident angle
  * incident light wavelength
  * light polarisation
  * materials in stack (including air and substrate)
  * depths of layers in stack (including air and substrate for array size agreement, although these depths are unused)
  * tt returns the reflection coefficient r and the transmission coefficient t
* stack function (first used in file task 12) takes as input:
  * incident wangle
  * incident wavelength
  * light polarisation (in task 13)
  * number of periods N
  * for the MgF2-Ta2O5 DBR, it returns the list of materials and their optimised depths in task 12
  * In task 13, it returns reflectance R
 * stack_fixed_d (first used task 14)
   * same as stack, but this time can specify different materials and depths
   * output now reflectance R
