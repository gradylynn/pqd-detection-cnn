# Power Quality Disturbance Prediction with Convolutional Neural Networks
## _Research by Grady Lynn, advised by Dr. Talayeh Razzaghi from the University of Oklahoma_

This repository is for the research that I began in the Spring 2020 semester with Dr. Razzaghi at the University of Oklahoma.
The goal of our research is to use convolutional nerual networks to detect and classify Power Quality Disturbances (PQDs) in power grids.

Some of the starting places for our research are papers listed below:
- Shouxiang Wang, Haiwen Chen, "A novel deep learning method for the classification of power quality disturbances using deep convolutional neural network"
- R.Igual, C.Medrano, F.J.Arcega, G.Mantescu, "Integral mathematical model of power quality disturbances"

The pqmodel.m was found [here](https://data.mendeley.com/datasets/6kmkk9bjdx/1) and is used to simulate different classes of power quality disturbances used for training our model.


## Data Generation
We use simulated data to train our CNN. We follow the parameters for data gereration outlined by R. Igual in "Integral mathematical model of power quality disturbances". This paper outlines parameters for generating samples of normal sine waves along with abnormal sine waves and provides a matlab script under the GNU General Public License that we use to generate wave-sample data.

## CNN Design
We implement a CNN structure using tensorflow.keras in python. 
