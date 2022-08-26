within FreeFluids.LMedia;
package Tests
        partial model FluidTestingA
          replaceable package Medium = Modelica.Media.Interfaces.PartialTwoPhaseMedium;
          parameter Medium.AbsolutePressure p = 1.0e5;
          parameter Medium.Temperature initialT = 273.19;
          parameter Medium.Temperature finalT = 372.15;
          parameter Real fract = 0.1 "fraction of the density of stateP to be used in state Dlow";
          Medium.Temperature T(start = initialT) "We will ramp the temperature";
          Medium.MolarMass MM;
          //creation of state from p and T and its properties
          Medium.ThermodynamicState StateP "original state from p and T";
          Real H;
          Real D;
          Real S;
          Real U "internal energy";
          Real Cp;
          Real Cv;
          Real Gamma;
          Real SS "speed of sound";
          Real Beta "isobaric expansion coefficient";
          Real Kappa "isothermal compressibility";
          //Real DerD_p_T;
          //Real DerD_T_p;
          Real DerD_p_h;
          Real DerD_h_p;
          Real Mu;
          Real Th;
          //reconstruction of the state from other variables
          Medium.ThermodynamicState StateD "state reproduced from d,T";
          Medium.ThermodynamicState StateH "state reproduced from p,h";
          Medium.ThermodynamicState StateS "state reproduced from p,s";
          //BaseProperties model
          Medium.BaseProperties BaseProp;
        algorithm
        //Construction of StateP and calculation of properties
          for i in 1:100 loop
            StateP := Medium.setState_pTX(p, T);
            MM := Medium.molarMass(StateP);
            H := Medium.specificEnthalpy(StateP);
            D := Medium.density(StateP);
            S := Medium.specificEntropy(StateP);
            U := Medium.specificInternalEnergy(StateP);
            Cp := Medium.specificHeatCapacityCp(StateP);
            Cv := Medium.specificHeatCapacityCv(StateP);
            Gamma := Cp / Cv;
            SS := Medium.velocityOfSound(StateP);
            Beta := Medium.isobaricExpansionCoefficient(StateP);
            Kappa := Medium.isothermalCompressibility(StateP);
            DerD_p_h := Medium.density_derp_h(StateP);
            DerD_h_p := Medium.density_derh_p(StateP);
            //DerD_p_T := Medium.density_derp_T(StateP);
            //DerD_T_p := Medium.density_derT_p(StateP);
            Mu := Medium.dynamicViscosity(StateP);
            Th := Medium.thermalConductivity(StateP);
            
            StateD := Medium.setState_dTX(D, T);
            StateH := Medium.setState_phX(p, H);
            StateS := Medium.setState_psX(p, S);
          end for;
        equation
        //Construction of BaseProperties
          BaseProp.p = p;
          BaseProp.T = T;
          der(T) = finalT - initialT;
          annotation(
            Documentation(info = "<html>
            <body>
            <p>This test maintains a constant pressure, and ramps temperature between the selected values. At each temperature, a thermodynamic state is created from p and T and, from it, different properties are calculated. The state is recovered also from d_T, p_H, and p_S, in order to check if the state is correctly reconstructed. </p>
            </body>
            </html>"));
        end FluidTestingA;

  model TestA1A
    extends FluidTestingA(redeclare replaceable package Medium = FreeFluids.LMedia.Fluids.Water(refState = "User", reference_T = 273.15, highPressure = true, inputChoice = "pT"), p = 50e5, initialT = 0.1 + 273.15, fract = 0.1, finalT = 260 + 273.15);
  end TestA1A;

        model TestA1B
          extends TestA1A(redeclare package Medium = FreeFluids.TMedia.Fluids.Water(refState = "User", reference_T = 273.15, highPressure = true, inputChoice = "pT"));
        end TestA1B;

        model TestA1C
          extends TestA1A(redeclare package Medium = Modelica.Media.Water.WaterIF97_pT);
        end TestA1C;

    model TestA2A
      extends FluidTestingA(redeclare replaceable package Medium = LMedia.Fluids.R134A(highPressure = true, refState = "IIR", reference_T = 100.0, inputChoice = "pT"), p = 20.0e5, initialT = (-50.0) + 273.15, finalT = 65.0 + 273.15, fract = 0.5);
    end TestA2A;

        model TestA2B
          extends TestA2A(redeclare package Medium = FreeFluids.TMedia.Fluids.R134A(highPressure = true, refState = "IIR", reference_T = 100.0, inputChoice = "pT"));
        end TestA2B;

  model TestModel
    extends Modelica.Media.Examples.Utilities.PartialTestModel(redeclare package Medium = FreeFluids.LMedia.Fluids.Water(inputChoice = "pT"));
  equation
  
  annotation(
      Documentation(info = "<html><head></head><body>Run with the old frontend. this is done by checking the box at \"Simulation Setup/Translation flags/Enable old frontend for code generation\"</body></html>"));end TestModel;

    model TestModel2
      extends Modelica.Media.Examples.Utilities.PartialTestModel2(redeclare package Medium = FreeFluids.LMedia.Fluids.Water(inputChoice = "pT"));
    equation

    end TestModel2;

    model ThreeTanks
      extends Modelica.Fluid.Examples.Tanks.ThreeTanks(redeclare package Medium = FreeFluids.LMedia.Fluids.Water(inputChoice = "ph"));
    equation
    
      annotation(
        experiment(StartTime = 0, StopTime = 140, Tolerance = 1e-06, Interval = 0.4));
    end ThreeTanks;

  model HeatingSystem
    extends Modelica.Fluid.Examples.HeatingSystem(redeclare package Medium = FreeFluids.LMedia.Fluids.Water(inputChoice = "ph"));
  equation
  
    annotation(
      experiment(StartTime = 0, StopTime = 6000, Tolerance = 1e-6, Interval = 12),
      Documentation(info = "<html><head></head><body>Run with the old frontend. this is done by checking the box at \"Simulation Setup/Translation flags/Enable old frontend for code generation\"</body></html>"));
  end HeatingSystem;
    
    model TestModelS
        extends Modelica.Media.Examples.Utilities.PartialTestModel(redeclare package Medium = Modelica.Media.Water.StandardWater);
    equation
    
    end TestModelS;
    
    model TestModel2S
      extends Modelica.Media.Examples.Utilities.PartialTestModel2(redeclare package Medium = Modelica.Media.Water.StandardWater);
    equation
    
    end TestModel2S;
    
  model ThreeTanksS
    extends Modelica.Fluid.Examples.Tanks.ThreeTanks(redeclare package Medium = Modelica.Media.Water.StandardWater);
  equation
  
    annotation(
      experiment(StartTime = 0, StopTime = 140, Tolerance = 1e-06, Interval = 0.4));
  end ThreeTanksS;
    
    model HeatingSystemS
      extends Modelica.Fluid.Examples.HeatingSystem(redeclare package Medium = Modelica.Media.Water.StandardWater);
    equation
    
      annotation(
        experiment(StartTime = 0, StopTime = 6000, Tolerance = 1e-6, Interval = 12));
    end HeatingSystemS;
  end Tests;