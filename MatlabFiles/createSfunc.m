clear
clc
close all
clear mex

%% Import structures to Matlab
Simulink.importExternalCTypes('../inc/SolverStruct.h');

%% Legacy Code Tool
def = legacy_code('initialize');
def.SFunctionName = 'magicSolver';
def.SourceFiles = {'SolverAlgorithm.c', 'MathFunction.c', 'SystemDynamic.c', 'PrintResult.c', 'BoundaryLine.c', 'AdaptiveGrid.c'};
def.HeaderFiles = {'SolverAlgorithm.h', 'MathFunction.h', 'SystemDynamic.h', 'SolverStruct.h', 'PrintResult.h', 'BoundaryLine.h', 'AdaptiveGrid.h'};
def.IncPaths = {'../inc/'};
def.SrcPaths = {'../src'};
def.OutputFcnSpec = 'void MagicBox(SolverInput u1[1], DynParameter u2[1], EnvFactor u3[1], SolverOutput y1[1], double u4, double u5)';
legacy_code('GENERATE_FOR_SIM', def)

%% Input settings
%inputSolver
