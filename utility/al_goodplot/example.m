set(gcf,'color',[0.98 0.98 1]);
subplot(1,4,1)
load hospital
tbl = table(hospital.Age,hospital.Sex,hospital.BloodPressure(:,2),'VariableNames',{'Age','Sex','BloodPressure'});
al_goodplot(tbl.BloodPressure(tbl.Sex=='Male'));
al_goodplot(tbl.BloodPressure(tbl.Sex=='Female'),2);
% Bilateral plots for two different groups.
% Default options for color, width, laterality and all are set if not provided.
xticks([1 2])
xticklabels({'Male', 'Female'})
title('Blood pressure')

%%
subplot(1,4,[2 3])
x=normrnd(2,1,100,3);
y=[x; normrnd(1,0.5,100,3)];
z=[x; normrnd(3.5,0.5,100,3)];
al_goodplot(y,[],0.5,[0 0 0.75; 0 0 1; 0.25 0.25 1],'left',[],std(x(:))/1000);
al_goodplot(z,[],0.5,[0 0.5 0],'right',[],std(x(:))/1000);
% Unilateral plots for 2 timepoints (left: before, right: after), 3 groups.
% One can produce multiple plots at once using a NxP input, P plots (1 per column).
% One can use different options for each column.
% If options are given only for 1 column, it is replicated for the others.
xlim([0.4 3.6])
xticks([1 2 3])
xticklabels({'grp A', 'grp B', 'grp C'})
title('before | after')

%%
subplot(1,4,4)
al_goodplot([],0.5,[],[0.2 0.2 0.2]);
% One can display the manual by using the type 'man' or not providing any data.
axis off