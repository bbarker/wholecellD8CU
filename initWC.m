%%%%%%             Configuration                %%%%%%
%%%%%%  please customize for your installation  %%%%%%

% where the whole cell model is installed:
WCDIR = '/home/brandon/DREAM8/WholeCell';
% where our custom scripts are located:
D8DIR = '/home/brandon/DREAM8/wholecellD8CU';

userData = struct;
userData.email = 'beb82@cornell.edu';
userData.firstName = 'Brandon';
userData.lastName = 'Barker';
userData.affiliation = 'Cornell';
userData.userName = 'beb82';
userData.hostName = 'gulab.cornell.edu';
userData.ip = '128.0.0.1';

%%%%%%            End of Configuration          %%%%%%


cd(WCDIR);
setWarnings();
setPath();

cd(D8DIR);
addpath(D8DIR);
