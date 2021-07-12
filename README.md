# The Automated Preprocessing Pipe-Line for the Estimation of Scale-wise Entropy from EEG Data (APPLESEED)


## About
APPLESEED (Puglia, Slobin, & Williams, 2021) preprocesses task-based or resting-state EEG data and then computes scale-wise entropy on the data. Scale-wise entropy is an adaptation of multiscale entropy (Costa et al., 2002) in which sample entropy (Richman & Moorman, 2000) is calculated across discontinuous data epochs (Grandy et al., 2016) on the residuals (i.e., after subtracting the within-person average response across trials/epochs) of the preprocessed EEG data, with the entropy similarity criterion parameter, r, recalculated at each scale. 

APPLESEED is implemented as a set of MATALB scripts and functions, including APPLESEED_setup(), which prepares raw EEG data for use in APPLESEED, APPLESEED(), which executes the preprocessing and entropy estimation pipeline on the prepared dataset, and APPLESEED_batch, which demonstrates how APPLESEED can be run as a batch across multiple subjects’ data.

See Puglia, Slobin, & Williams (2021) for more information on the development and optimization of APPLESEED.

* Costa, M., Goldberger, A.L., Peng, C.-K., 2002. Multiscale Entropy Analysis of Complex Physiologic Time Series. Phys. Rev. Lett. 89, 068102. 
* Grandy, T.H., Garrett, D.D., Schmiedek, F., Werkle-Bergner, M., 2016. On the estimation of brain signal entropy from sparse neuroimaging data. Sci. Rep. 6, 23073. 
* Puglia, M.H., Slobin, J.S., Williams, C.L., 2021. The Automated Preprocessing Pipe-Line for the Estimation of Scale-wise Entropy from EEG Data (APPLESEED): Development and validation for use in pediatric populations. bioRxiv. https://doi.org/10.1101/2021.07.10.450198.
* Richman, J.S., Moorman, J.R., 2000. Physiological time-series analysis using approximate entropy and sample entropy. Am. J. Physiol. Circ. Physiol. 278, H2039–H2049. 


## Accompanying dataset
An example dataset is available for download from https://openneuro.org/datasets/ds003710. The user should place all APPLESEED scripts within the MATLAB path, download the example dataset to an "APPLESEED_Example_Dataset" directory, and change into this directory in MATLAB. The example code in APPLESEED_batch will then run APPLESEED_setup() and APPLESEED() from the MATLAB command line on the example dataset which includes: EEG data from 48 recording sessions from 13 infants and a channel location file and a bin file (located in the "APPLESEED_Example_Dataset > code" directory). 

See the following publications for more information on the tutorial dataset:

* Puglia, M.H., Krol, K.M., Missana, M., Williams, C.L., Lillard, T.S., Morris, J.P., Connelly, J.J., & Grossmann, T. (2020). Epigenetic tuning of brain signal entropy in emergent human social behavior. BMC medicine, 18(1), 1-24. https://doi.org/10.1186/s12916-020-01683-x
* Puglia, M.H., Slobin, J.S., Williams, C.L., 2021. The Automated Preprocessing Pipe-Line for the Estimation of Scale-wise Entropy from EEG Data (APPLESEED): Development and validation for use in pediatric populations. bioRxiv. https://doi.org/10.1101/2021.07.10.450198.


## License
Copyright 2021 Meghan H. Puglia, University of Virginia (meghan.puglia@virignia.edu)

APPLESEED is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

APPLESEED is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with APPLESEED.  If not, see <https://www.gnu.org/licenses/>.