# FlowAnalysis

Install FlowCal with: _pip install FlowCal_

To run the script: _python mainFlow.py X_

X is the only accepted argument. It's the gating percentage (values between 20 and 90).
If no argument is given the default percentage is 50% 

From [FlowCal wiki](https://taborlab.github.io/FlowCal/python_tutorial/gate.html#density-gate): _Density gating automatically identifies the region with the highest density of events in a two-dimensional diagram (FSC vs SSC), and calculates how big it should be to capture a certain percentage of the total event count._

The script outputs one csv file per inducer. e.g. data_ahl.txt. You can change this in the _processData_ function
