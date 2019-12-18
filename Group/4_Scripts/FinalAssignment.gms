Sets
i                Units           /Coal, Nuclear, CCGT, Hydro, Geothermal, Diesel/
j                Renewable Units /Wind_Onshore, Wind_Offshore, PV/
k                Hydrogen Units  /Electrolyzer, Storage, FuelCell/
t                Timeperiods     /t1*t8760/
;

Parameters

inv(i)                   Annualized investment costs ($_MW)
invR(j)                  Annualized investment costs of renewables ($_MW)
invH(k)                  Annualized investment costs of hydroge ($_MW)
size(i)                  Capacity size (MW)
sizeR(j)                 Renewables capacity size (MW)
sizeH(k)                 Hydrogen capacity size (MW)
ra(i)                    Resource available (MW)
raR(j)                   Renewables resource available (MW)
p_min(i)                 Minimum production level (MW)
p_minH(k)                Minimum charge_discharge_storage level (MW)
c_e(i)                   Emission costs ($_kgCO2)
c_vom(i)                 Variable O&M costs of conventional plants ($_MWh)
c_fom(k)                 Fixed O&M costs of hydrogen ($_MW)
q_e(i)                   Efficiency adjusted emission quantity (kgCO2_MWh)
c_fn(i)                  Efficiency adjusted fuel cost ($_MWh)
n(k)                     Hydrogen system efficiency per component (Percentage)

d(t)                     Demand (MW)
wind_Onshore(t)          Wind Onshore profile
wind_Offshore(t)         Wind Offshore profile
pv(t)                    PV profile
hydro(t)                 Maximum hourly hydro generation based on monthly data (MWh)

fuelPrice(i)             Fuel price change (Percentage of total)
emisPrice(i)             Emission price change (Percentage of total)
subsidy                  Wind Offshore subsidy ($_MWh)
l0(k)                    Storage level at t=0 (MWh)
HInv(k)                  Hydrogen investment consts change (Percentage of total)

$onEcho > TechnologyData.txt
par=inv          rng=inv!A2:B7        rdim=1
par=invR         rng=invR!A2:B4       rdim=1
par=size         rng=size!A2:B7       rdim=1
par=sizeR        rng=sizeR!A2:B4      rdim=1
par=ra           rng=ra!A2:B7         rdim=1
par=raR          rng=raR!A2:B4        rdim=1
par=p_min        rng=p_min!A2:B7      rdim=1
par=c_e          rng=c_e!A2:B7        rdim=1
par=c_vom        rng=c_vom!A2:B7      rdim=1
par=q_e          rng=q_e!A2:B7        rdim=1
par=c_fn         rng=c_fn!A2:B7       rdim=1
par=fuelPrice    rng=fuelPrice!A2:B7  rdim=1
par=emisPrice    rng=emisPrice!A2:B7  rdim=1
$offEcho
$call gdxxrw.exe TechnologyData.xlsx maxDupeErrors=8760 @TechnologyData.txt
$gdxin TechnologyData.gdx
$load inv invR size sizeR ra raR p_min c_e c_vom q_e c_fn fuelPrice emisPrice
$gdxin

$onEcho > Demand_Profiles.txt
par=d                    rng=d!A2:B8761                  rdim=1
par=wind_Onshore         rng=wind_Onshore!A2:B8761       rdim=1
par=wind_Offshore        rng=wind_Offshore!A2:B8761      rdim=1
par=pv                   rng=pv!A2:B8761                 rdim=1
par=hydro                rng=hydro!A2:B8761              rdim=1
$offEcho
$call gdxxrw.exe Demand_Profiles.xlsx maxDupeErrors=8760 @Demand_Profiles.txt
$gdxin Demand_Profiles.gdx
$load d wind_Onshore wind_Offshore pv hydro
$gdxin
;

$onEcho > TechnologyData_Hydrogen.txt
par=invH         rng=invH!A2:B4       rdim=1
par=sizeH        rng=sizeH!A2:B4      rdim=1
par=p_minH       rng=p_minH!A2:B4     rdim=1
par=c_fom        rng=c_fom!A2:B4      rdim=1
par=n            rng=n!A2:B4          rdim=1
par=HInv         rng=HInv!A2:B4       rdim=1
$offEcho
$call gdxxrw.exe TechnologyData_Hydrogen.xlsx maxDupeErrors=8760 @TechnologyData_Hydrogen.txt
$gdxin TechnologyData_Hydrogen.gdx
$load invH sizeH p_minH c_fom n HInv
$gdxin

subsidy=0;
l0('Storage')=0;


Free Variables
Z                Objective function value
;

Integer Variables
p_max(i)         Units built
pR_max(j)        Renewables units built
pH_max(k)        Hydrogen units built

;

Positive Variables
p(i,t)           Production level
l(k,t)           Hydrogen level
;

Equations
obj              Objective function
min(i,t)         Minimum production level constraint
max(i,t)         Maximum production level constraint
max_C(i)         Maximum capacity
min_RG           Minimum renewables generation
max_hydro(i,t)   Maximum hydro generation
dem(t)           Demand constraint

minH(k,t)        Minimum level constraint for hydrogen
maxH(k,t)        Maximum level constraint for hydrogen
storageH0(k,t)   Hydrogen storage level constraint at t=0
storageH(k,t)    Hydrogen storage level constraint
discharge(k,t)   Fuel cell level constraint
charge(k,t)      Electrolyzer level constraint
;

*****************************************************************************************************
********************************************Normal Model*********************************************
*****************************************************************************************************
*obj ..                   Z =e= sum((i,t),(c_fn(i)*fuelPrice(i)+(emisPrice(i)*c_e(i)*q_e(i))+c_vom(i))*p(i,t)) + sum(i,inv(i)*p_max(i)*size(i)) + sum(j,invR(j)*pR_max(j)*sizeR(j)) - sum(t,pR_max('Wind_Offshore')*sizeR('Wind_Offshore')*wind_Offshore(t)*subsidy);
*dem(t) ..                sum(i,p(i,t)) =e= d(t) - pR_max('Wind_Onshore')*sizeR('Wind_Onshore')*wind_Onshore(t) - pR_max('Wind_Offshore')*sizeR('Wind_Offshore')*wind_Offshore(t) - pR_max('PV')*sizeR('PV')*pv(t);

min(i,t) ..              p(i,t) =g= p_min(i);
max(i,t) ..              p(i,t) =l= p_max(i)*size(i);
max_C(i) ..              p_max(i) =l= (ra(i))/size(i);
min_RG ..                sum(t, pR_max('Wind_Onshore')*sizeR('Wind_Onshore')*wind_Onshore(t) + pR_max('Wind_Offshore')*sizeR('Wind_Offshore')*wind_Offshore(t) + pR_max('PV')*sizeR('PV')*pv(t)) =g= 27663851.875;
max_hydro('Hydro',t) ..  p('Hydro',t) =l= hydro(t);

*****************************************************************************************************
****************************************Hydrogen extension*******************************************
*****************************************************************************************************

obj ..                   Z =e= sum((i,t),(c_fn(i)*fuelPrice(i)+(emisPrice(i)*c_e(i)*q_e(i))+c_vom(i))*p(i,t)) + sum(i,inv(i)*p_max(i)*size(i)) + sum(j,invR(j)*pR_max(j)*sizeR(j)) + sum(k, (invH(k)*HInv(k)+c_fom(k))*pH_max(k)*sizeH(k)) - sum(t,pR_max('Wind_Offshore')*sizeR('Wind_Offshore')*wind_Offshore(t)*subsidy);

dem(t) ..                sum(i,p(i,t)) + l('FuelCell',t) - l('Electrolyzer',t) =e= d(t) - pR_max('Wind_Onshore')*sizeR('Wind_Onshore')*wind_Onshore(t) - pR_max('Wind_Offshore')*sizeR('Wind_Offshore')*wind_Offshore(t) - pR_max('PV')*sizeR('PV')*pv(t);

minH(k,t) ..             l(k,t) =g= p_minH(k);
maxH(k,t) ..             l(k,t) =l= pH_max(k)*sizeH(k);

*storageH0('Storage',t)$(ord(t)=1) ..     l('Storage',t) =e= l0('Storage')+ l('Electrolyzer',t) - l('FuelCell',t);
*storageH('Storage',t)$(ord(t)>1) ..      l('Storage',t) =e= l('Storage',t-1)+ l('Electrolyzer',t) - l('FuelCell',t);
*discharge('FuelCell',t) ..               l('FuelCell',t) =l= l('Storage',t);
*charge('Electrolyzer',t) ..              l('Electrolyzer',t) =l= pR_max('Wind_Onshore')*sizeR('Wind_Onshore')*wind_Onshore(t) + pR_max('Wind_Offshore')*sizeR('Wind_Offshore')*wind_Offshore(t) + pR_max('PV')*sizeR('PV')*pv(t);


storageH0('Storage',t)$(ord(t)=1) ..     l('Storage',t) =e= l0('Storage')+ l('Electrolyzer',t)*n('Electrolyzer') - l('FuelCell',t)/n('FuelCell');
storageH('Storage',t)$(ord(t)>1) ..      l('Storage',t) =e= l('Storage',t-1)+ l('Electrolyzer',t)*n('Electrolyzer') - l('FuelCell',t)/n('FuelCell');
discharge('FuelCell',t) ..               l('FuelCell',t) =l= l('Storage',t)*n('FuelCell');
charge('Electrolyzer',t) ..              l('Electrolyzer',t) =l= pR_max('Wind_Onshore')*sizeR('Wind_Onshore')*wind_Onshore(t) + pR_max('Wind_Offshore')*sizeR('Wind_Offshore')*wind_Offshore(t) + pR_max('PV')*sizeR('PV')*pv(t);

option intVarUp = 0;

Model Carolina /all/ ;
Carolina.OptCR=0.001;
*$onecho > cplex.opt
*epmrk = 0.6
*$offecho
*Carolina.OptFile=1;
Solve Carolina using MIP minimizing Z ;

file results / results.txt /;
*Specifiying that we want a comma separated file:
results.pc = 5;
put results;
put "Time period", "Coal", "Nuclear", "CCGT", "Hydro", "Geothermal", "Diesel"/;
loop(t,
put t.tl,
         loop(i,
                 put p.l(i,t);
         );
put /;
);
putclose;

file resultsH / resultsH.txt /;
*Specifiying that we want a comma separated file:
resultsH.pc = 5;
put resultsH;
put "Time period", "Electrolyzer", "Storage", "FuelCell"/;
loop(t,
put t.tl,
         loop(k,
                 put l.l(k,t);
         );
put /;
);
putclose;

execute_unload "results.gdx" z.L p_max.L pR_max.L pH_max.L dem.M

execute 'gdxxrw.exe results.gdx o=results.xlsx var=z.L rng=System!'
execute 'gdxxrw.exe results.gdx o=results.xlsx var=p_max.L rng=Units! rDim=1'
execute 'gdxxrw.exe results.gdx o=results.xlsx var=pR_max.L rng=RenewableUnits! rDim=1'
execute 'gdxxrw.exe results.gdx o=results.xlsx var=pH_max.L rng=HydrogenUnits! rDim=1'
execute 'gdxxrw.exe results.gdx o=results.xlsx equ=dem.M rng=Price! rDim=1'
