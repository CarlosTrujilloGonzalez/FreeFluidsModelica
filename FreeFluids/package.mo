within ;
package FreeFluids
  import SI = Modelica.SIunits;
  import pi = Modelica.Constants.pi;
  import g_n = Modelica.Constants.g_n;
  import R = Modelica.Constants.R "8.3144 Pa·m3/(K·mol)";

  type Fraction = Real(min = 0.0, max = 1.0, nominal = 0.5);
  type MassFlux = Real(quantity = "MassFlux", final unit = "kg/(s.m2)");
  type Possitive = Real(min = 0.0);
  annotation(
    uses(Modelica(version = "3.2.3")));
end FreeFluids;
