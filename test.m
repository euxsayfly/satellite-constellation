% test 
% 测试matlab与stk的互联


%% 启动stk

uiapplication = actxserver('STK11.application'); % stk未运行
% uiapplication = actxGetRunningServer('STK11.application'); % stk已经运行

% 建立对象

root = uiapplication.Personality2;

% 新建场景

root.NewScenario('scnew');

% 关闭当前场景

root.CloseScenario;

%% 场景设置

% 拾取当前场景，sc是以后当前场景的handle

sc = root.CurrentScenario;

% 设置仿真时长

sc.SetTimePeriod('Today','+24hr') % 仿真时间设置格式1（不推荐）

sc.StartTime = '23 Nov 2023 04:00:00.000';
sc.StopTime = '24 Nov 2023 04:00:00.000';
sc.SetTimePeriod(sc.StartTime,sc.StopTime); % 仿真时间设置格式2（推荐使用，后续建立对象会用到相关参数）

% 设置场景中的单位
root.UnitPreferences.SetCurrentUnit('Distance','km');
root.UnitPreferences.SetCurrentUnit('Latitude','deg');
root.UnitPreferences.SetCurrnetUnit('Longitude','deg');

%% 创建对象

sat = sc.Children.New('eSatellite','mysat'); % 新建卫星对象。控制句柄为sat

sat.Propagator.Propagate; % 执行Propagate命令后，卫星将运行，即在STK中可以看到卫星运行轨迹

sat.get % 即可查看对象的所有属性

% 几个常见的卫星属性
sat.PropagatorType % 可以获取卫星对象的轨道动力学模型
sat.Propagator     % 卫星轨道生成器
sat.InstanceName   % 获取sat对应在STK场景中卫星对象的名字

% 查看场景中有多少个卫星
satCollection = sc.Children.GetElements('eSatellite');
satCollection.Count
sat1 = satCollection.Item('sat_1'); % 获取卫星对象控制句柄


%% 覆盖特性
uiap = actxserver('STK11.application');
root = uiap.Personality2;
root.NewScenario('Coverage');
sc = root.CurrentScenario;
sat = sc.Children.New(18,'mysat');

kep = sat.Propagator.InitialState.Representation.ConvertTo('eOrbitStateClassical');
kep.SizeShapeType = 'eSizeShapeAltitude';
kep.SizeShape.ApogeeAltitude = 500;
kep.SizeShape.PerigeeAltitude = 500;
kep.Orientation.Inclination = 50;
sat.Propagator.InitialState.Representation.Assign(kep);
sat.Propagator.Propagate;

% 新建传感器
sen = sat.Children.New('eSensor','mysen');
sen.CommonTasks.SetPatternSimpleConic(77,1);
sen.VO.ProjectionType = 'eProjectionEarthIntersections';

% 新建CoverageDefinition对象
covdef = sc.Children.New('eCoverageDefinition','mycov');

covdef.Grid.BoundsType = 'eBoundsLatLonRegion'; 
covdef.Grid.Bounds.MinLongitude = -120;
covdef.Grid.Bounds.MaxLongitude = 120;
covdef.Grid.Bounds.MinLatitude = -30;
covdef.Grid.Bounds.MaxLatitude = 30;

covdef.Grid.Resolution.LatLon = 1; % 网点间隔为1°

covdef.Graphics.Static.IsPointsVisible = 0; % 将网格设置为不可见
covdef.Advanced.AutoRecompute = 0;    % 取消Coveragedefinition对象的自动更新

% 选择航天器/传感器对象，并激活
covdef.AssetList.AvailableAssets
covdef.AssetList.Add(covdef.AssetList.AvailableAssets{2});

% Figureofmerit对象包含覆盖性的参数
figmerit1 = covdef.Children.New('eFigureofmerit','revisittime');
figmerit1.SetDefinitionType('eFmrevisittime');
covdef.ComputeAccesses();

% 直接生成报告
% root.ExecuteCommand('ReportCreate */CoverageDefinition/mycov/FigureOfMerit/revisittime Type Save Style "Value By Grid Point" File "D:\revis.txt"');

% 利用DataProviders提取数据到matlab工作区

retimeDP = figmerit1.DataProviders.Item('Value by Longitude').Exec;
retimedataMax = cell2mat(retimeDP.DataSets.GetDataSetByName('Maximum').GetValues);
retimedataAve = cell2mat(retimeDP.DataSets.GetDataSetByName('Average').GetValues);
retimedataMin = cell2mat(retimeDP.DataSets.GetDataSetByName('Minimum').GetValues);

%% sensor对象操作设置
sen = sat.Children.New('eSensor','mysen');

sen.PatternType % 查看sensor的形状及参数获取，获得的参数为：eSnSimpleConic

sen.Pattern.get % 设置sensor形状及角度

sen.CommonTasks.SetPatternSimpleConic(60,1); % 角度设置，半张角60°，角分辨率1°


%% 数据获取 DataProviders的使用方法

% 读取轨道六根数

uiap = actxserver('STK11.application');
root = uiap.Personality2;
root.NewScenario('elements');
sc = root.CurrentScenario;

for j = 1:20
satname = ['str_',num2str(j)];
sat = sc.Children.New(18,satname);
%设置卫星轨道动力学类型，这里选用高精度HPOP
sat.SetPropagatorType('ePropagatorHPOP');

%设置卫星轨道参数
kep = sat.Propagator.InitialState. Representation.ConvertTo('eOrbitStateClassic');
kep.SizeShapeType = 'eSizeShapeSemimajorAxis';
kep.LocationType = 'eLocationTrueAnomaly';
kep.Orientation.AscNodeType = 'eAscNodeLAN';

kep.SizeShape.SemimajorAxis = 6378.14 + 1000 + 5000 * rand(1);
kep.SizeShape.Eccentricity = 0.2 * rand(1);
kep.Orientation.Inclination = 50 * rand(1);
kep.Orientation.ArgOfPerigee = 360 * rand(1);
kep.Orientation.AscNode.Value = 360 * rand(1);
kep.Location.Value = 360 * rand(1);
sat.Propagator.InitialState.Representation.Assign(kep);
sat.Propagator.Propagate;
end

outputData = cell(21,7);
outputData{1,1} = {'卫星名称'};
outputData{1,2} = {'半长轴'};
outputData{1,3} = {'偏心率'};
outputData{1,4} = {'倾角'};
outputData{1,5} = {'升交点赤经'};
outputData{1,6} = {'近地点辐角'};
outputData{1,7} = {'真近点角'};

%获取所有卫星的路径
satpathcollection = root.ExecuteCommand('ShowNames * Class Satellite');
satpathcollection.Item(0);
%Item(0)中包含了所有卫星路径，但是是用一个字符串包含所有路径，
%中间用空格隔开，下面的语句是先将字符串分割，然后去除空格元素，
%则satPaths的每个元素都只包含一个卫星路径
satPaths = regexp(satpathcollection.Item(0),' ','split');
satPaths(cellfun(@isempty,satPaths)) = [];

for i = 1:length(satPaths)
sat = root.GetObjectFromPath(satPaths{i});%注意是大括号
%将时间改为历元形式
root.UnitPreferences.Item('DateFormat').SetCurrentUnit('EpSec');
%计算多个历元语句为times = {0;1000;2000};
times = {0;1000;2000};
elems = {'Semi-major Axis';'Eccentricity';'Inclination';'RAAN';'Arg of Perigee';'True Anomaly'};

satDP = sat.DataProviders.Item('Classical Elements').Group.Item('ICRF').ExecSingleElementsArray(times,elems);

outputData{i+1,1} = sat.InstanceName;
outputData{i+1,2} = cell2mat(satDP.GetArray(int32(0)));
outputData{i+1,3} =cell2mat(satDP.GetArray(int32(1)));
outputData{i+1,4} = cell2mat(satDP.GetArray(int32(2)));
outputData{i+1,5} = cell2mat(satDP.GetArray(int32(3)));
outputData{i+1,6} = cell2mat(satDP.GetArray(int32(4)));
outputData{i+1,7} = cell2mat(satDP.GetArray(int32(5)));
end




