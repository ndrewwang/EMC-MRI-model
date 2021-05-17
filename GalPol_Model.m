function out = GalPol_Model(L_param,cavg,ipulse,current_profile_array,t_eval_string,pos_eval_pts,Ve_fn,V0_fn,D_fn,tp0_fn,chi_fn,lambda_fn,kappa_fn,partfrac_fn);
%
% GalPol_Model.m
%
% Model exported on Nov 5 2020, 15:07 by COMSOL 5.4.0.346.
L_param = string(L_param)+' [m]';
cavg = string(cavg) + ' [mol/m^3]';
pos_eval_pts = string(pos_eval_pts);

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.modelPath('/COMSOL Simulations/MATLAB Link');

model.label('GalPol_Model.mph');

model.comments(['Galpol\n\n']);

model.param.set('L', L_param, 'electrode separation');
model.param.set('F', '96487[C/mol]', 'Faraday''s constant');
model.param.set('R', '8.3145[J/mol/K]', 'gas constant');
model.param.set('T', '298.15[K]', 'ambient temperature');
model.param.set('nup', '1', 'cation stoich');
model.param.set('num', '1', 'anion stoich');
model.param.set('zp', '1', 'cation charge');
model.param.set('zm', '-zp*nup/num', 'anion charge');
model.param.set('sp', '-1', 'cation rxn stoich');
model.param.set('sm', '0', 'anion rxn stoich');
model.param.set('s0', '0', 'solvent rxn stoich');
model.param.set('cavg', cavg, 'avg salt concentration');
model.param.set('ipulse', ipulse, 'pulse current');
model.param.set('pos_eval_pts', pos_eval_pts, 'positions for export');

model.component.create('comp1', false);

model.component('comp1').geom.create('geom1', 1);

model.component('comp1').defineLocalCoord(false);

model.result.table.create('tbl1', 'Table');

model.func.create('pw1', 'Piecewise');
model.func('pw1').set('smooth', 'contd1');
model.func('pw1').set('smoothzone', '0.001');
model.func('pw1').set('smoothends', true);
model.func('pw1').set('pieces', current_profile_array);

model.component('comp1').mesh.create('mesh1');

model.component('comp1').geom('geom1').repairTolType('relative');
model.component('comp1').geom('geom1').create('i1', 'Interval');
model.component('comp1').geom('geom1').feature('i1').set('coord', {'0' 'L'});
model.component('comp1').geom('geom1').feature('fin').set('repairtoltype', 'relative');
model.component('comp1').geom('geom1').run;

model.component('comp1').variable.create('var1');

model.component('comp1').variable('var1').set('Ve', Ve_fn, 'solute partial molar volume');
model.component('comp1').variable('var1').set('Vp', 'Ve*(1-tp0)', 'cation partial molar volume');
model.component('comp1').variable('var1').set('Vm', 'Ve*tp0', 'anion partial molar volume');
model.component('comp1').variable('var1').set('V0', V0_fn, 'solvent partial molar volume');
model.component('comp1').variable('var1').set('D', D_fn, 'thermo diffusion coeff');
model.component('comp1').variable('var1').set('tp0', tp0_fn, 'cation transference number');
model.component('comp1').variable('var1').set('kappa', kappa_fn, 'conductivity');
model.component('comp1').variable('var1').set('chi', chi_fn, 'darken factor');
model.component('comp1').variable('var1').set('Nm', 'i/F/zm-zp*Np/zm');
model.component('comp1').variable('var1').set('N0', '(vbox/V0)-Np*(Vp/V0)-Nm*(Vm/V0)');
model.component('comp1').variable('var1').set('cT', '(1/V0)+(nup+num-(nup*Vp+num*Vm)/V0)*c');
model.component('comp1').variable('var1').set('lambda', lambda_fn);
model.component('comp1').variable('var1').set('y', partfrac_fn);
model.component('comp1').variable('var1').set('c0', 'cT-2*c');

model.component('comp1').physics.create('g', 'GeneralFormPDE', 'geom1');
model.component('comp1').physics('g').identifier('cation_balance');
model.component('comp1').physics('g').field('dimensionless').field('cation_flux');
model.component('comp1').physics('g').field('dimensionless').component({'Np'});
model.component('comp1').physics('g').prop('Units').set('DependentVariableQuantity', 'molarflux');
model.component('comp1').physics('g').create('dir1', 'DirichletBoundary', 0);
model.component('comp1').physics('g').feature('dir1').selection.set([1]);
model.component('comp1').physics('g').create('cons1', 'Constraint', 0);
model.component('comp1').physics('g').feature('cons1').selection.set([2]);
model.component('comp1').physics.create('g2', 'GeneralFormPDE', 'geom1');
model.component('comp1').physics('g2').identifier('charge_continuity');
model.component('comp1').physics('g2').field('dimensionless').field('current_density');
model.component('comp1').physics('g2').field('dimensionless').component({'i'});
model.component('comp1').physics('g2').prop('Units').set('DependentVariableQuantity', 'currentdensity');
model.component('comp1').physics('g2').create('dir1', 'DirichletBoundary', 0);
model.component('comp1').physics('g2').feature('dir1').selection.set([1]);
model.component('comp1').physics('g2').create('cons1', 'Constraint', 0);
model.component('comp1').physics('g2').feature('cons1').selection.set([2]);
model.component('comp1').physics.create('g3', 'GeneralFormPDE', 'geom1');
model.component('comp1').physics('g3').identifier('volume_continuity');
model.component('comp1').physics('g3').field('dimensionless').field('vol_avg_vel');
model.component('comp1').physics('g3').field('dimensionless').component({'vbox'});
model.component('comp1').physics('g3').prop('Units').set('DependentVariableQuantity', 'velocity');
model.component('comp1').physics('g3').create('cons1', 'Constraint', 0);
model.component('comp1').physics('g3').feature('cons1').selection.set([1]);
model.component('comp1').physics('g3').create('cons2', 'Constraint', 0);
model.component('comp1').physics('g3').feature('cons2').selection.set([2]);
model.component('comp1').physics.create('g4', 'GeneralFormPDE', 'geom1');
model.component('comp1').physics('g4').identifier('macinness_equation');
model.component('comp1').physics('g4').field('dimensionless').field('cat_echem_pot');
model.component('comp1').physics('g4').field('dimensionless').component({'mup'});
model.component('comp1').physics('g4').prop('Units').set('DependentVariableQuantity', 'activationenergy');
model.component('comp1').physics('g4').create('cons1', 'Constraint', 0);
model.component('comp1').physics('g4').feature('cons1').selection.set([2]);
model.component('comp1').physics('g4').create('cons2', 'Constraint', 0);
model.component('comp1').physics('g4').feature('cons2').selection.set([1]);
model.component('comp1').physics.create('g5', 'GeneralFormPDE', 'geom1');
model.component('comp1').physics('g5').identifier('ficks_law');
model.component('comp1').physics('g5').field('dimensionless').field('conc');
model.component('comp1').physics('g5').field('dimensionless').component({'c'});
model.component('comp1').physics('g5').prop('Units').set('DependentVariableQuantity', 'concentration');
model.component('comp1').physics('g5').create('cons1', 'Constraint', 0);
model.component('comp1').physics('g5').feature('cons1').selection.set([2]);
model.component('comp1').physics('g5').create('cons2', 'Constraint', 0);
model.component('comp1').physics('g5').feature('cons2').selection.set([1]);

model.component('comp1').mesh('mesh1').create('edg1', 'Edge');

model.capeopen.label('Thermodynamics Package');

model.component('comp1').view('view1').axis.set('xmin', -60.000244140625);
model.component('comp1').view('view1').axis.set('xmax', 6060);

model.component('comp1').physics('g').prop('Units').set('SourceTermQuantity', 'reactionrate');
model.component('comp1').physics('g').feature('gfeq1').set('f', 'nup*d(c,t)+d(Np,x)');
model.component('comp1').physics('g').feature('gfeq1').set('Ga', 0);
model.component('comp1').physics('g').feature('gfeq1').set('da', 0);
model.component('comp1').physics('g').feature('dir1').set('r', '-sp*i/F');
model.component('comp1').physics('g2').prop('Units').set('SourceTermQuantity', 'currentsource');
model.component('comp1').physics('g2').feature('gfeq1').set('f', 'd(i,x)');
model.component('comp1').physics('g2').feature('gfeq1').set('Ga', 0);
model.component('comp1').physics('g2').feature('gfeq1').set('da', 0);
model.component('comp1').physics('g2').feature('dir1').set('r', 'ipulse*pw1(root.t[1/s])');
model.component('comp1').physics('g3').prop('Units').set('CustomSourceTermUnit', '1/s');
model.component('comp1').physics('g3').feature('gfeq1').set('f', 'd(tp0,x)*Ve*i/F + d(c,x)*d(Ve,x)*D/(1-Ve*c) + d(vbox,x)');
model.component('comp1').physics('g3').feature('gfeq1').set('Ga', 0);
model.component('comp1').physics('g3').feature('gfeq1').set('da', 0);
model.component('comp1').physics('g3').feature('cons1').set('R', 'vbox + Vp*sp*i/F');
model.component('comp1').physics('g4').prop('Units').set('CustomSourceTermUnit', 'J/C/m');
model.component('comp1').physics('g4').feature('gfeq1').set('f', 'd(mup,x)/F/zp+i/kappa-(nup+num)*R*T*(1-tp0)*chi/cT/V0/F/zp/nup/c*d(c,x)');
model.component('comp1').physics('g4').feature('gfeq1').set('Ga', 0);
model.component('comp1').physics('g4').feature('gfeq1').set('da', 0);
model.component('comp1').physics('g4').feature('cons1').set('R', 'mup');
model.component('comp1').physics('g5').prop('Units').set('SourceTermQuantity', 'molarflux');
model.component('comp1').physics('g5').feature('gfeq1').set('f', 'Np+nup*D*d(c,x)-tp0*i/F/zp-nup*c*vbox');
model.component('comp1').physics('g5').feature('gfeq1').set('Ga', 0);
model.component('comp1').physics('g5').feature('gfeq1').set('da', 0);
model.component('comp1').physics('g5').feature('init1').set('c', 'cavg');
model.component('comp1').physics('g5').feature('cons1').set('R', 'Np+sp*i/F');

model.component('comp1').mesh('mesh1').feature('size').set('hauto', 1);
model.component('comp1').mesh('mesh1').feature('size').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('size').set('hmin', 2.0E-7);
model.component('comp1').mesh('mesh1').run;

model.study.create('std1');
model.study('std1').create('time', 'Transient');
model.study('std1').feature('time').set('activate', {'g' 'on' 'g2' 'on' 'g3' 'on' 'g4' 'on' 'g5' 'on'});

model.sol.create('sol1');
model.sol('sol1').study('std1');
model.sol('sol1').attach('std1');
model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').create('t1', 'Time');
model.sol('sol1').feature('t1').create('fc1', 'FullyCoupled');
model.sol('sol1').feature('t1').feature.remove('fcDef');

model.result.dataset.create('cpt1', 'CutPoint1D');
model.result.numerical.create('pev1', 'EvalPoint');
model.result.numerical('pev1').set('probetag', 'none');
model.result.create('pg6', 'PlotGroup1D');
model.result.create('pg1', 'PlotGroup1D');
model.result.create('pg2', 'PlotGroup1D');
model.result.create('pg3', 'PlotGroup1D');
model.result.create('pg4', 'PlotGroup1D');
model.result.create('pg5', 'PlotGroup1D');
model.result('pg6').create('ptgr1', 'PointGraph');
model.result('pg6').feature('ptgr1').selection.set([1]);
model.result('pg1').create('lngr1', 'LineGraph');
model.result('pg1').feature('lngr1').set('xdata', 'expr');
model.result('pg1').feature('lngr1').selection.all;
model.result('pg2').create('lngr1', 'LineGraph');
model.result('pg2').feature('lngr1').set('xdata', 'expr');
model.result('pg2').feature('lngr1').selection.all;
model.result('pg3').create('lngr1', 'LineGraph');
model.result('pg3').feature('lngr1').set('xdata', 'expr');
model.result('pg3').feature('lngr1').selection.all;
model.result('pg4').create('lngr1', 'LineGraph');
model.result('pg4').feature('lngr1').set('xdata', 'expr');
model.result('pg4').feature('lngr1').selection.all;
model.result('pg5').create('lngr1', 'LineGraph');
model.result('pg5').feature('lngr1').set('data', 'dset1');
model.result('pg5').feature('lngr1').set('xdata', 'expr');
model.result('pg5').feature('lngr1').selection.all;
model.result.export.create('tbl1', 'Table');

model.study('std1').feature('time').set('tlist', t_eval_string);
model.study('std1').feature('time').set('usertol', true);
model.study('std1').feature('time').set('rtol', 0.0001);
model.study('std1').feature('time').set('plot', true);
model.study('std1').feature('time').set('discretization', {'g' 'physics' 'g2' 'physics' 'g3' 'physics' 'g4' 'physics' 'g5' 'physics'});

model.sol('sol1').attach('std1');
model.sol('sol1').label('Solver 1');
model.sol('sol1').feature('v1').set('resscalemethod', 'auto');
model.sol('sol1').feature('v1').set('clist', {t_eval_string '6.0[s]'});
model.sol('sol1').feature('t1').set('tlist', t_eval_string);


model.sol('sol1').feature('t1').set('rtol', 0.0001);
model.sol('sol1').feature('t1').set('atolglobalvaluemethod', 'manual');
model.sol('sol1').feature('t1').set('atolglobal', '0.0010');
model.sol('sol1').feature('t1').set('fieldselection', 'comp1_mup');
model.sol('sol1').feature('t1').set('atolmethod', {'comp1_mup' 'global' 'comp1_c' 'global' 'comp1_i' 'global' 'comp1_Np' 'global' 'comp1_vbox' 'global'});
model.sol('sol1').feature('t1').set('atolvaluemethod', {'comp1_mup' 'manual' 'comp1_c' 'manual' 'comp1_i' 'manual' 'comp1_Np' 'manual' 'comp1_vbox' 'manual'});
model.sol('sol1').feature('t1').set('atolfactor', {'comp1_mup' '0.1' 'comp1_c' '0.1' 'comp1_i' '0.1' 'comp1_Np' '0.1' 'comp1_vbox' '0.1'});
model.sol('sol1').feature('t1').set('atol', {'comp1_mup' '1e-3' 'comp1_c' '1e-3' 'comp1_i' '1e-3' 'comp1_Np' '1e-3' 'comp1_vbox' '1e-3'});
model.sol('sol1').feature('t1').set('atoludot', {'comp1_mup' '1e-3' 'comp1_c' '1e-3' 'comp1_i' '1e-3' 'comp1_Np' '1e-3' 'comp1_vbox' '1e-3'});
model.sol('sol1').feature('t1').set('initialstepbdf', '0.0010');
model.sol('sol1').feature('t1').set('bwinitstepfrac', '0.0010');
model.sol('sol1').feature('t1').set('plot', true);
model.sol('sol1').feature('t1').feature('dDef').set('ooc', false);
model.sol('sol1').feature('t1').feature('dDef').set('rhob', 400);
model.sol('sol1').feature('t1').feature('fc1').set('dtech', 'auto');
model.sol('sol1').runAll;

model.result.dataset('dset1').label('Solution 1');
model.result.dataset('cpt1').set('pointx', 'range(0,L/pos_eval_pts,L)');
model.result.numerical('pev1').set('data', 'cpt1');
model.result.numerical('pev1').set('table', 'tbl1');
model.result.numerical('pev1').set('expr', {'Np' 'i' 'vbox' 'mup' 'c' 'c0'});
model.result.numerical('pev1').setResult;
model.result('pg6').label('Voltage (mup/F)');
model.result('pg6').set('xlabel', 'Time (s)');
model.result('pg6').set('ylabel', '-mup/F (V)');
model.result('pg6').set('xlabelactive', false);
model.result('pg6').set('ylabelactive', false);
model.result('pg6').feature('ptgr1').set('expr', '-mup/F');
model.result('pg6').feature('ptgr1').set('unit', 'V');
model.result('pg6').feature('ptgr1').set('descr', '-mup/F');
model.result('pg1').set('xlabel', 'x-coordinate (m)');
model.result('pg1').set('ylabel', 'Dependent variable Np (mol/(m<sup>2</sup>*s))');
model.result('pg1').set('xlabelactive', false);
model.result('pg1').set('ylabelactive', false);
model.result('pg1').feature('lngr1').set('descr', 'Dependent variable Np');
model.result('pg1').feature('lngr1').set('xdataexpr', 'x');
model.result('pg1').feature('lngr1').set('xdataunit', 'm');
model.result('pg1').feature('lngr1').set('xdatadescr', 'x-coordinate');
model.result('pg1').feature('lngr1').set('smooth', 'internal');
model.result('pg1').feature('lngr1').set('resolution', 'normal');
model.result('pg2').set('xlabel', 'x-coordinate (mm)');
model.result('pg2').set('ylabel', 'Dependent variable i (A/m<sup>2</sup>)');
model.result('pg2').set('xlabelactive', false);
model.result('pg2').set('ylabelactive', false);
model.result('pg2').feature('lngr1').set('expr', 'i');
model.result('pg2').feature('lngr1').set('unit', 'A/m^2');
model.result('pg2').feature('lngr1').set('descr', 'Dependent variable i');
model.result('pg2').feature('lngr1').set('xdataexpr', 'x');
model.result('pg2').feature('lngr1').set('xdataunit', 'mm');
model.result('pg2').feature('lngr1').set('xdatadescr', 'x-coordinate');
model.result('pg2').feature('lngr1').set('smooth', 'internal');
model.result('pg2').feature('lngr1').set('resolution', 'normal');
model.result('pg3').set('xlabel', 'x-coordinate (m)');
model.result('pg3').set('ylabel', 'Dependent variable vbox (m/s)');
model.result('pg3').set('xlabelactive', false);
model.result('pg3').set('ylabelactive', false);
model.result('pg3').feature('lngr1').set('expr', 'vbox');
model.result('pg3').feature('lngr1').set('unit', 'm/s');
model.result('pg3').feature('lngr1').set('descr', 'Dependent variable vbox');
model.result('pg3').feature('lngr1').set('xdataexpr', 'x');
model.result('pg3').feature('lngr1').set('xdataunit', 'm');
model.result('pg3').feature('lngr1').set('xdatadescr', 'x-coordinate');
model.result('pg3').feature('lngr1').set('smooth', 'internal');
model.result('pg3').feature('lngr1').set('resolution', 'normal');
model.result('pg4').set('xlabel', 'x-coordinate (m)');
model.result('pg4').set('ylabel', 'Dependent variable mup (J/mol)');
model.result('pg4').set('xlabelactive', false);
model.result('pg4').set('ylabelactive', false);
model.result('pg4').feature('lngr1').set('expr', 'mup');
model.result('pg4').feature('lngr1').set('unit', 'J/mol');
model.result('pg4').feature('lngr1').set('descr', 'Dependent variable mup');
model.result('pg4').feature('lngr1').set('xdataexpr', 'x');
model.result('pg4').feature('lngr1').set('xdataunit', 'm');
model.result('pg4').feature('lngr1').set('xdatadescr', 'x-coordinate');
model.result('pg4').feature('lngr1').set('smooth', 'internal');
model.result('pg4').feature('lngr1').set('resolution', 'normal');
model.result('pg5').set('xlabel', 'x-coordinate (m)');
model.result('pg5').set('ylabel', 'Dependent variable c (mol/m<sup>3</sup>)');
model.result('pg5').set('axislimits', true);
model.result('pg5').set('xmin', -7.766152266412974E-6);
model.result('pg5').set('xmax', 0.001007766152266413);
model.result('pg5').set('ymin', 0);
model.result('pg5').set('ymax', 4000);
model.result('pg5').set('xlabelactive', false);
model.result('pg5').set('ylabelactive', false);
model.result('pg5').feature('lngr1').set('expr', 'c');
model.result('pg5').feature('lngr1').set('unit', 'mol/m^3');
model.result('pg5').feature('lngr1').set('descr', 'Dependent variable c');
model.result('pg5').feature('lngr1').set('xdataexpr', 'x');
model.result('pg5').feature('lngr1').set('xdataunit', 'm');
model.result('pg5').feature('lngr1').set('xdatadescr', 'x-coordinate');
model.result('pg5').feature('lngr1').set('smooth', 'internal');
model.result('pg5').feature('lngr1').set('resolution', 'normal');
model.result.export('tbl1').set('filename', '/COMSOL Simulations/MATLAB Link/comsol_output.csv');
model.result.export('tbl1').run

out = model;
