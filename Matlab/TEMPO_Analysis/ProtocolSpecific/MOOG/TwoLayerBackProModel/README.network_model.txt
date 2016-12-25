2 layer Backpropagation Network Model
-------------------------------------

This folder contains an implementation of a two layer backpropagation model
to learn a coordinate transformation (rotation) of heading direction based
where heading direction is encoded by input units distributed evenly about
the sphere in spherical coordinates.  The model simulates the three protocols
of visual only inputs (rotation), vestibular only inputs (linear transformation),
and combined (both inputs).

The input data and target data are generated using the 5 parameter cosine
tuned model used for fitting tuning curves in MST neurons.  The curve fitting
code must be included in the path (folder containing Curvefit_* modules) for
the network model to run.  The network is built using the Matlab Neural Network
Toolbox (user must also have this installed).

Major Modules Description
-------------------------

run_network_model
   top level script, execute this to train and simulate the model
create_sample_points
   creates the array and matrices of spherical coordinates distributed
   evenly about the sphere that determine the preferred direction of
   input and output units
init_params
   creates parameters for the tuning function models that represent
   the input and output units
make_training_set
   most important module besides run_network_model.  creates the input
   layer data and output layer target data.  
simulate_inputs
   simulates the output of the trained network for each of the training
   inputs.  this is what creates the hidden layer data.
run_network_model_trials
   iterates run_network_model for some number of iterations with the
   option to use different numbers of hidden layer units and numbers
   of input and output units distributed about the sphere.  saves results
   for each iteration.
make_chris_data
   creates a clean matlab file with the output of the hidden layer units
   for a single network simulation.  
   a sample is included in network_data.mat

Parameters in run_network_model
-------------------------------

remake_training_set
   1 = rebuild the training data in make_training_set (must be done once)
retrain_network
   1 = train a new network starting with random weights and biases
simulate_network
   1 = simulate the trained network for each of the training inputs.
   creates the hidden layer outputs.
plot_training_results (must have simulate_network=1)
   1 = plots the output of the trained network against the target data.
   used to verify visually that network is learning properly.
plot_hidden_layer
   1 = plots the outputs of the hidden layer in the trained network

These parameters control the first 3 phases of the network model:
   - generate input and target data
   - initialize network and train using Neural Network Toolbox
   - simulate the network outputs with the training data set

The last phase is analyzing and plotting results.
The code currently generates a slew of plots, some of them redundant.
The code could definitely use a clean up...

Questions?
pvw1@cec.wustl.edu
