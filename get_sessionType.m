%
% 2015-Apr-13
%

function session_type = get_sessionType(thisRID, thisSID)

thisSession = ['rat' thisRID, '-' thisSID];

% setting
standard4scene = {'rat232-04','rat232-05','rat232-06','rat232-07', ...
    'rat219-08','rat219-09','rat219-10','rat219-11','rat219-12','rat219-13','rat214-05','rat214-06','rat214-07','rat214-08', ...
    'rat221-01','rat221-02','rat221-03','rat221-04','rat234-01','rat234-02','rat234-03','rat234-04', ...
    'rat285-01','rat285-02','rat285-03','rat285-04','rat295-01','rat295-02','rat295-03','rat295-04', ...
    'rat415-10', 'rat415-11', 'rat415-12', 'rat415-13','rat427-01', 'rat427-02', 'rat427-03', ...
    'rat561-01','rat561-02','rat561-03','rat561-04','rat561-05'...
    'rat027-05', 'rat037-05','rat052-08','rat057-05','rat057-14','rat057-15'};
ambiguity_zebra_pebbles = {'rat214-11','rat219-14','rat232-08','rat234-05','rat285-05','rat295-05', 'rat415-14', 'rat427-04','rat561-06'};
ambiguity_bamboo_mountain = {'rat219-15','rat232-09','rat234-06','rat285-06','rat295-06', 'rat415-16'};
trace_zebra_pebbles = {'rat214-12','rat219-16','rat219-17','rat232-10','rat232-11','rat234-07', 'rat285-07', 'rat295-07', 'rat415-21','rat561-07'};
blocked_trace = {'rat232-12','rat234-08','rat285-08','rat295-09', 'rat415-22', 'rat415-23', 'rat427-16', 'rat427-17','rat561-08'};

new_scene_learning = {'rat232-13','rat234-09','rat285-09','rat295-10', 'rat415-17', 'rat415-18', 'rat427-11', 'rat427-12','rat561-09','rat561-10'};
standard6scene = {'rat232-14', 'rat232-15','rat234-10','rat234-11','rat295-11','rat295-12', 'rat415-19', 'rat415-20', 'rat427-13', 'rat427-14', 'rat427-15','rat561-11'};
%

% determine session type
sessionTot = {standard4scene, ambiguity_zebra_pebbles, ambiguity_bamboo_mountain, trace_zebra_pebbles, blocked_trace, new_scene_learning, standard6scene};
session_type = 6;
for iter = 1:length(sessionTot)
    if sum(strcmp(thisSession, sessionTot{iter})) == 1
        session_type = iter;
        return;
    end
end
%

end