%***File paths below are for Windows 10 64 bit***

%Different Ways to Start STK or attached to running STK
%------------------------------------------------------

%% COM method of attached to running STK

uiApplication = actxGetRunningServer('STK12.application');
root = uiApplication.Personality2;
scen = root.CurrentScenario;

clear uiApplication root scen
%% COM method of launching STK
uiApplication = actxserver('STK12.application');    %Change STK11 to start different STK version
uiApplication.Visible = 1;
root=uiApplication.Personality2;

%CREATE A NEW SCENARIO
%using Connect commands
root.ExecuteCommand('New / Scenario MatlabConnect');
pause(5);
root.CloseScenario;

%using the Object Model API
root.NewScenario('MatlabOM');
pause(5);
root.CloseScenario;

uiApplication.Quit;
uiApplication.release;
clear uiApplication root 
