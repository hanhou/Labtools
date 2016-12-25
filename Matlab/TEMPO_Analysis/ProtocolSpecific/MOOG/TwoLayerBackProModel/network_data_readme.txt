Hi Chris,

Attached is a v6 matlab file containing outputs of hidden units for a
network trained at gaze angles -20,0,20.  Here is a description of the
variables in the mat file

hidden_unit_outputs - the raw output of the hidden layer units in [-1 1]
  dim 1 - outputs for different stim types
  dim 2 - outputs for different gaze angles
  dim 3 - outputs for different spherical points
  dim 4 - outputs for different hidden units

hidden_unit_amplitude - the max value when hidden_unit_outputs are
normalized so that min response is 0.
  dim 1 - outputs for different stim types
  dim 2 - outputs for different gaze angles
  dim 3 - outputs for different hidden units

hidden_unit_pref_dir_az - preferred direction in azimuth for each hidden
unit, dims are same as above

hidden_unit_pref_dir_el - preferred direction in azimuth for each hidden
unit, dims are same as above

num_stim_types - size of stim type dimension (3)

hidden_unit_num_gaze_angles - size of gaze angle dimension (3)

num_unique_points - size of the points sampled about the sphere
dimension (26)

num_hidden_units - size of the hidden unit dimension (150)

hidden_unit_gaze_angles_r - parallel to the gaze angle dimension, giving
the gaze angle for each offest [-20 0 20], in radians

unique_point_azimuth_r - parallel to the points about the sphere
dimension, giving the azimuth angle for each offset, in radians

unique_point_elevation - same as above for elevation angle

-----------------------

Defines for the stim type directions are in the MOOG directory in
Curvefit_defines.  I usually put this macro at the top of all my scripts,
to import all defines relevant to curve fitting or the models.
If you don't want this, the defines are
CURVEFIT_OCCULAR_STIM = 2;
CURVEFIT_VESTIBULAR_STIM = 1;
CURVEFIT_OCCULAR_AND_VESTIBULAR_STIM = 3;

As an example
hidden_unit_outputs(2,3,10,23)
is the output of the 23rd hidden unit, at the point (0,0) in azimuth and
elevation, for the occular only stim type, and at gaze angle 20 degrees.

The corresponding angles would be in
hidden_unit_gaze_angles_r(3)
unique_point_azimuth_r(10)
unique_point_elevation_r(10)

let me know if you have any questions, or you think the data is bunk.
paul
