function [param] = funConfig2(filename)

param.config(1).a = 460;
param.config(1).b = 650;
param.config(1).Mratio = param.config(1).b/param.config(1).a;
param.config(1).d = 278;
param.config(1).c = param.config(1).d/param.config(1).Mratio;
param.config(1).Dd = 9.148262329;
param.config(1).pixelL = 0.20;

param.config(2).a = 460;
param.config(2).b = 650;
param.config(2).Mratio = param.config(2).b/param.config(2).a;
param.config(2).d = 278;
param.config(2).c = param.config(2).d/param.config(2).Mratio;
param.config(2).Dd = 9.784784777;
param.config(2).pixelL = 0.152;

param.config(3).a = 460;
param.config(3).b = 650;
param.config(3).Mratio = param.config(3).b/param.config(3).a;
param.config(3).d = 278;
param.config(3).c = param.config(3).d/param.config(3).Mratio;
param.config(3).Dd = 9.042558996;
param.config(3).pixelL = 0.1;

param.roi(2).tanLMTF = 5;
param.roi(2).normLMTF = 4;
param.roi(2).tanLCNI = 10;
param.roi(2).normLCNI = 7;

param.roi(1).tanLMTF = 5;
param.roi(1).normLMTF = 4;
param.roi(1).tanLCNI = 10;
param.roi(1).normLCNI = 6;

param.roi(3).tanLMTF = 5;
param.roi(3).normLMTF = 4;
param.roi(3).tanLCNI = 10;
param.roi(3).normLCNI = 5;

param.roi(4).tanLMTF = 5;
param.roi(4).normLMTF = 4;
param.roi(4).tanLCNI = 10;
param.roi(4).normLCNI = 5;

xml_write(filename,param);

end