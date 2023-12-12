function Obj = estimate_Obj(p)
%%   实现变参数星座覆盖计算
% 对于Walker星座 可优化参数有轨道参数，以及星座内部参数


% 打开STK

uiap = actxGetRunningServer('STK11.application');
root = uiap.Personality2;
% root.NewScenario('walker');
sc = root.CurrentScenario;

% 设定场景时间（可以不用设置直接默认）

% strBegTime = '01 Dec 2023 09:12:00.000';
% strEndTime = '02 Dec 2023 09:12:00.000';
% scenario.SetTimePeriod(strBegTime,strEndTime);
% scenario.Epoch = strBegTime;
% root.ExecuteCommand('Animate * Reset');   % 运行时在解析SetTimePeriod时出现错误

% 建立种子卫星

seedSat = sc.Children.New(18,'mysat');
seedSat.SetPropagatorType('ePropagatorJ4Perturbation');  % 轨道模型可以修改

% 设置卫星根数

a = 6371 + p(:,1); e = 0.0; i = p(:,2); w = 0; Raan = p(:,3); M = 0;
propagator = seedSat.Propagator;

% propagator.StarTime = sc.StartTime;
% propagator.StopTiem = sc.StopTime;
% propagator.Step = 60.0;

propagator.InitialState.Representation.AssignClassical('eCoordinateSystemJ2000',a,e,i,w,Raan,M);


% 创建卫星传感器

seedSen = seedSat.Children.New('eSensor','Sen');
% 设置传感器参数
seedSen.CommonTasks.SetPatternSimpleConic(44.85, 0.1);
seedSen.VO.ProjectionType = 'eProjectionEarthIntersections';


% 创建walker星座

nPlan = 4; % 轨道平面数
nPerPlan = 3; % 每个轨道平面卫星数量
% nRANNSpreed = 1; % 相位因子

walkercommand = ['Walker */Satellite/mysat Type Delta NumPlanes ' num2str(nPlan) ' NumSatsPerPlane ' num2str(nPerPlan) ' InterPlanePhaseIncrement ' num2str(p(:,4)) ' ColorByPlane Yes'];
root.ExecuteCommand(walkercommand);



% 创建覆盖区域对象 
% 针对的覆盖区域可以修改
% 美军夏威夷：
% 珍珠港海军基地：北纬21°21′，西经157°58′
% 斯科菲尔德兵营：北纬21°29′，西经158°03′
% 怀厄勒空军基地：北纬21°19′，西经157°56′
% 美国本土： 北纬24°50′到北纬49°   西经66°50′到西经125°
% 巴拿马运河：北纬8°到北纬9°，西经79°到80°
% 苏伊士运河：北纬30°到北纬31°，东经32°到33°
% 马六甲海峡：北纬2°到7°，东经32°到33°
% 北冰洋航道：北纬65°到北纬85°，东经12°到13°
%

covdef = sc.Children.New('eCoverageDefinition','mycov');

covdef.Grid.BoundsType = 'eBoundsLatLonRegion';
covdef.Grid.Bounds.MinLongitude = 12;   % 最小经度
covdef.Grid.Bounds.MaxLongitude = 13;    % 最大经度
covdef.Grid.Bounds.MinLatitude = 2;     % 最小纬度
covdef.Grid.Bounds.MaxLatitude = 7;      % 最大纬度

covdef.Grid.Resolution.LatLon = 0.25; % 网点间隔为0.25°
covdef.Graphics.Static.IsPointsVisible = 0; % 将网格设置为不可见
covdef.Advanced.AutoRecompute = 0; % 取消覆盖区域对象的自动更新


% 选择航天器/传感器对象，并激活

covdef.AssetList.AvailableAssets;
for i = 4 : 3 + nPlan * nPerPlan
covdef.AssetList.Add(covdef.AssetList.AvailableAssets{i});
end

% 建立覆盖品质参数（最大重访时间）
figmerit1 = covdef.Children.New('eFigureofmerit','revisittime');
figmerit1.SetDefinitionType('eFmrevisittime');
covdef.ComputeAccesses();

% 直接生成报告
root.ExecuteCommand('ReportCreate */CoverageDefinition/mycov/FigureOfMerit/revisittime Type Save Style "Value By Grid Point" File "D:\revis.txt"');

% 利用DataProviders提取数据到工作区

retimeDP = figmerit1.DataProviders.Item('Value by Longitude').Exec;
retimedataMax = cell2mat(retimeDP.DataSets.GetDataSetByName('Average').GetValues);



Obj = max(retimedataMax);



covdef.Unload
seedSat.Unload

satItems = root.ExecuteCommand('ShowNames * Class Satellite');
satPaths = strsplit(strtrim(satItems.Item(0)),' ');

for i = 1 : nPlan * nPerPlan
        
        sat = root.GetObjectFromPath(char(satPaths(i)));
        sat.Unload
end



end